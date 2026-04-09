"""Pooled multi-sample MIL training.

Loads embeddings from all 12 patient samples and trains a single shared
ProportionGuidedMIL model. This addresses data starvation (9-44 spots per
sample) by pooling to ~1,300 training spots, preventing cell-type collapse.

Usage:
  python train_pooled_mil.py [--epochs 300] [--lr 1e-4] [--device cuda]
"""

import argparse
import json
import logging
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
MODEL_DIR = REPO_ROOT / "CITEgeist" / "model"
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(MODEL_DIR))

OUTPUT_ROOT = REPO_ROOT / "output" / "morphology_assignment"

SAMPLES = [
    "HCC22-088-P1-S1", "HCC22-088-P1-S2",
    "HCC22-088-P2-S1", "HCC22-088-P2-S2",
    "HCC22-088-P3-S1_A", "HCC22-088-P3-S2",
    "HCC22-088-P4-S1", "HCC22-088-P4-S2_1i_rep",
    "HCC22-088-P5-S1", "HCC22-088-P5-S2_F_rep",
    "HCC22-088-P6-S1", "HCC22-088-P6-S2_D",
]

CELL_TYPES = [
    "Endothelial", "Fibroblasts", "B_Cells", "Macrophages", "Monocytes",
    "CD8_T_Cells", "CD4_T_Cells", "Cancer_Luminal", "Cancer_Basal",
    "Dendritic_Cells",
]

CELL_TYPES_FOR_PROPS = [
    "Endothelial", "Fibroblasts", "B_Cells", "Macrophages", "Monocytes",
    "CD8_T_Cells", "CD4_T_Cells", "Cancer_Luminal", "Cancer_Basal",
    "Dendritic_Cells",
]

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
logger = logging.getLogger(__name__)


def load_pooled_data(min_nuclei: int = 1, prop_dir: str = None):
    """Load embeddings and proportions from all samples.

    Args:
        min_nuclei: Minimum nuclei per spot.
        prop_dir: If set, load proportions from CSV files in this directory
            instead of per-spot proportions.npy files.

    Returns:
        List of (embeddings, proportions, sample_name, barcode) tuples
    """
    import pandas as pd

    # Pre-load proportion CSVs if prop_dir specified
    prop_csvs = {}
    if prop_dir:
        prop_path = Path(prop_dir)
        for sample in SAMPLES:
            csv_path = prop_path / sample / f"{sample}_cell_prop_global_results.csv"
            if csv_path.exists():
                df = pd.read_csv(csv_path, index_col=0)
                available = [c for c in CELL_TYPES_FOR_PROPS if c in df.columns]
                prop_csvs[sample] = df[available]
                logger.info("  Loaded QP proportions for %s: %d spots x %d types",
                            sample, len(df), len(available))
            else:
                logger.warning("  No QP proportion CSV for %s at %s", sample, csv_path)

    pooled = []
    for sample in SAMPLES:
        emb_dir = OUTPUT_ROOT / sample / "embeddings"
        if not emb_dir.exists():
            logger.warning("No embeddings dir for %s, skipping", sample)
            continue

        n_loaded = 0
        for spot_dir in sorted(emb_dir.iterdir()):
            if not spot_dir.is_dir():
                continue
            emb_path = spot_dir / "embeddings.npy"
            if not emb_path.exists():
                continue

            emb = np.load(emb_path)
            if emb.shape[0] < min_nuclei:
                continue

            # Get proportions: from CSV override or from saved .npy
            barcode = spot_dir.name
            if prop_dir and sample in prop_csvs:
                df = prop_csvs[sample]
                if barcode not in df.index:
                    continue
                prop = df.loc[barcode].values.astype(np.float32)
            else:
                prop_path_npy = spot_dir / "proportions.npy"
                if not prop_path_npy.exists():
                    continue
                prop = np.load(prop_path_npy)

            pooled.append((emb, prop, sample, barcode))
            n_loaded += 1

        logger.info("  %s: loaded %d spots", sample, n_loaded)

    return pooled


def stratified_split(pooled, val_fraction=0.2, seed=42):
    """Split data with stratification by sample.

    Ensures each sample contributes proportionally to train and val sets.

    Returns:
        (train_data, val_data) where each is a list of (emb, prop, sample, barcode)
    """
    rng = np.random.default_rng(seed)
    by_sample = defaultdict(list)
    for item in pooled:
        by_sample[item[2]].append(item)

    train_data, val_data = [], []
    for sample, spots in sorted(by_sample.items()):
        idx = rng.permutation(len(spots))
        n_val = max(1, int(len(spots) * val_fraction))
        val_data.extend([spots[i] for i in idx[:n_val]])
        train_data.extend([spots[i] for i in idx[n_val:]])
        logger.info("  %s: %d train, %d val", sample, len(idx) - n_val, n_val)

    return train_data, val_data


def train_pooled(
    train_data,
    val_data,
    n_epochs=300,
    lr=1e-4,
    hidden_dim=256,
    entropy_weight=0.05,
    device="cpu",
):
    """Train ProportionGuidedMIL on pooled data.

    Returns:
        (best_state_dict, config, history)
    """
    import torch
    import torch.nn.functional as F
    from proportion_mil import ProportionGuidedMIL, proportion_loss, entropy_regularization
    from scipy.stats import pearsonr

    embed_dim = train_data[0][0].shape[1]
    n_cell_types = len(CELL_TYPES)

    model = ProportionGuidedMIL(
        input_dim=embed_dim,
        n_cell_types=n_cell_types,
        hidden_dim=hidden_dim,
    ).to(device)

    optimizer = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-4)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=n_epochs)

    best_val_r = -1.0
    best_state = None
    history = {
        "train_loss": [],
        "val_loss": [],
        "val_r": [],
        "val_r_per_sample": [],
        "lr": [],
    }

    # Pre-convert to tensors
    train_tensors = [
        (torch.from_numpy(e).float(), torch.from_numpy(p).float())
        for e, p, _, _ in train_data
    ]
    val_tensors = [
        (torch.from_numpy(e).float(), torch.from_numpy(p).float())
        for e, p, _, _ in val_data
    ]
    val_samples = [item[2] for item in val_data]

    # Compute inverse-frequency weights per cell type across training set
    # so rare types get proportionally more gradient signal
    all_props = np.array([p for _, p, _, _ in train_data])  # (N_spots, K)
    mean_props = all_props.mean(axis=0)  # (K,)
    # Inverse frequency: weight = 1/freq, clamped and normalized to mean=1
    inv_freq = 1.0 / (mean_props + 1e-3)
    type_weights = inv_freq / inv_freq.mean()  # normalize so mean weight = 1
    type_weights_t = torch.from_numpy(type_weights.astype(np.float32)).to(device)
    logger.info("  Type weights (inverse-freq): %s",
                {ct: f"{w:.2f}" for ct, w in zip(CELL_TYPES, type_weights)})

    for epoch in range(n_epochs):
        # --- Train ---
        model.train()
        epoch_loss = 0.0
        # Shuffle training order each epoch
        perm = np.random.permutation(len(train_tensors))
        for idx in perm:
            emb, target = train_tensors[idx]
            emb, target = emb.to(device), target.to(device)

            pred, attention = model(emb)

            # Weighted MSE: rare types contribute more to loss
            weighted_mse = (type_weights_t * (pred - target) ** 2).mean()

            # KL divergence (unweighted, for distribution matching)
            eps = 1e-8
            kl = (pred * torch.log((pred + eps) / (target + eps))).sum()

            # Attention diversity: penalize attention columns with near-zero
            # mean activation (forces model to USE all type channels)
            attn_mean = attention.mean(dim=0)  # (K,) mean attention per type
            attn_diversity = -torch.log(attn_mean + eps).mean()

            loss = (weighted_mse + 0.1 * kl
                    + entropy_weight * entropy_regularization(attention)
                    + 0.01 * attn_diversity)

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            epoch_loss += loss.item()

        scheduler.step()
        avg_loss = epoch_loss / len(train_tensors)
        history["train_loss"].append(avg_loss)
        history["lr"].append(scheduler.get_last_lr()[0])

        # --- Validate ---
        model.eval()
        val_preds, val_targets = [], []
        val_loss = 0.0
        with torch.no_grad():
            for emb, target in val_tensors:
                emb, target = emb.to(device), target.to(device)
                pred, attention = model(emb)
                val_loss += proportion_loss(pred, target).item()
                val_preds.append(pred.cpu().numpy())
                val_targets.append(target.cpu().numpy())

        val_loss /= len(val_tensors)
        val_preds_arr = np.array(val_preds)
        val_targets_arr = np.array(val_targets)

        # Overall Pearson r
        val_r = pearsonr(val_preds_arr.flatten(), val_targets_arr.flatten())[0]
        history["val_loss"].append(val_loss)
        history["val_r"].append(float(val_r))

        # Per-sample breakdown
        per_sample_r = {}
        sample_set = sorted(set(val_samples))
        for s in sample_set:
            mask = [i for i, vs in enumerate(val_samples) if vs == s]
            if len(mask) < 2:
                continue
            s_preds = val_preds_arr[mask].flatten()
            s_targets = val_targets_arr[mask].flatten()
            if np.std(s_preds) > 0 and np.std(s_targets) > 0:
                per_sample_r[s] = float(pearsonr(s_preds, s_targets)[0])
        history["val_r_per_sample"].append(per_sample_r)

        if val_r > best_val_r:
            best_val_r = val_r
            best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}

        if (epoch + 1) % 10 == 0:
            logger.info(
                "  Epoch %d/%d  loss=%.4f  val_loss=%.4f  val_r=%.4f  lr=%.2e",
                epoch + 1, n_epochs, avg_loss, val_loss, val_r,
                scheduler.get_last_lr()[0],
            )

        if (epoch + 1) % 50 == 0:
            logger.info("  Per-sample val r: %s",
                        {k: f"{v:.3f}" for k, v in per_sample_r.items()})

    config = {
        "input_dim": embed_dim,
        "n_cell_types": n_cell_types,
        "hidden_dim": hidden_dim,
    }
    logger.info("Best val Pearson r: %.4f", best_val_r)
    return best_state, config, history


def main():
    parser = argparse.ArgumentParser(description="Pooled multi-sample MIL training")
    parser.add_argument("--epochs", type=int, default=300, help="Training epochs")
    parser.add_argument("--lr", type=float, default=1e-4, help="Learning rate")
    parser.add_argument("--hidden-dim", type=int, default=256, help="MIL hidden dimension")
    parser.add_argument("--entropy-weight", type=float, default=0.05,
                        help="Entropy regularization weight (anti-collapse)")
    parser.add_argument("--device", default="cuda", help="Device (cuda or cpu)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--val-fraction", type=float, default=0.2, help="Validation fraction")
    parser.add_argument("--min-nuclei", type=int, default=1, help="Min nuclei per spot")
    parser.add_argument("--prop-dir", type=str, default=None,
                        help="Directory with cuOPT QP proportions per sample "
                             "(overrides proportions.npy in embedding dirs)")
    parser.add_argument("--output-dir", type=str, default=None,
                        help="Output directory for checkpoint (default: OUTPUT_ROOT)")
    args = parser.parse_args()

    import torch
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)

    logger.info("Loading pooled data from %d samples...", len(SAMPLES))
    pooled = load_pooled_data(min_nuclei=args.min_nuclei, prop_dir=args.prop_dir)
    logger.info("Total spots loaded: %d", len(pooled))

    if len(pooled) == 0:
        logger.error("No data found! Check that Stage 2 embeddings exist.")
        sys.exit(1)

    logger.info("Stratified train/val split (%.0f%% val)...", args.val_fraction * 100)
    train_data, val_data = stratified_split(pooled, val_fraction=args.val_fraction, seed=args.seed)
    logger.info("Train: %d spots, Val: %d spots", len(train_data), len(val_data))

    logger.info("Training pooled MIL (epochs=%d, lr=%g, entropy_weight=%g)...",
                args.epochs, args.lr, args.entropy_weight)
    best_state, config, history = train_pooled(
        train_data, val_data,
        n_epochs=args.epochs,
        lr=args.lr,
        hidden_dim=args.hidden_dim,
        entropy_weight=args.entropy_weight,
        device=args.device,
    )

    # Save checkpoint
    save_dir = Path(args.output_dir) if args.output_dir else OUTPUT_ROOT
    save_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_path = save_dir / "pooled_mil_checkpoint.pt"
    torch.save({
        "model_state_dict": best_state,
        "config": config,
    }, checkpoint_path)
    logger.info("Saved checkpoint: %s", checkpoint_path)

    # Save training history
    history_path = save_dir / "pooled_training_history.json"
    with open(history_path, "w") as f:
        json.dump(history, f, indent=2)
    logger.info("Saved training history: %s", history_path)

    # Summary
    logger.info("=== Summary ===")
    logger.info("  Total spots: %d (train=%d, val=%d)", len(pooled), len(train_data), len(val_data))
    logger.info("  Best val Pearson r: %.4f", max(history["val_r"]))
    logger.info("  Final val Pearson r: %.4f", history["val_r"][-1])

    # Check for type collapse
    import torch
    from proportion_mil import ProportionGuidedMIL
    model = ProportionGuidedMIL(**config).to(args.device)
    model.load_state_dict(best_state)
    model.eval()

    all_preds = []
    with torch.no_grad():
        for emb, prop, _, _ in val_data:
            emb_t = torch.from_numpy(emb).float().to(args.device)
            pred, _ = model(emb_t)
            all_preds.append(pred.cpu().numpy())

    all_preds = np.array(all_preds)
    mean_pred = all_preds.mean(axis=0)
    logger.info("  Mean predicted proportions: %s",
                {ct: f"{v:.3f}" for ct, v in zip(CELL_TYPES, mean_pred)})

    n_active = (mean_pred > 0.02).sum()
    logger.info("  Active cell types (mean > 2%%): %d/%d", n_active, len(CELL_TYPES))
    if n_active <= 2:
        logger.warning("  WARNING: Possible type collapse — only %d active types!", n_active)


if __name__ == "__main__":
    main()
