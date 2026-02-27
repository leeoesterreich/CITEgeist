"""Prototype learning training script (Stage 2).

This script trains the projection heads and prototypes for cell type assignment
using a pre-trained VAE encoder (from Stage 1). The training uses Sinkhorn
optimal transport with spot-level cell type proportions as supervision.

INCLUDES ANTI-COLLAPSE MEASURES:
- K-means initialization of prototypes from VAE embeddings
- Temperature warmup (start sharp, then anneal)
- Gradient clipping
- Multi-spot batching for stable gradients
- Per-epoch monitoring (projection variance, transport entropy)

Usage:
    python -m CITEgeist.model.train_prototypes \
        --vae-checkpoint /path/to/vae_final.pt \
        --patches-dir /path/to/patches \
        --proportions-file /path/to/proportions.csv \
        --output-dir /path/to/output \
        --epochs 50

Input Requirements:
    - VAE checkpoint: Trained VAE model from Stage 1
    - Patches directory: Contains spot_{spot_id}_patches.npy files
    - Proportions CSV: spot_id column + cell type columns (must sum to 1 per row)

Output:
    - prototypes_final.pt: Final trained model (heads + prototypes)
    - prototypes_checkpoint_epoch_{N}.pt: Checkpoints every 10 epochs
    - training_history.json: Training metrics per epoch
"""
import argparse
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import torch
from sklearn.cluster import KMeans
from tqdm import tqdm

# Support both package import and direct import
try:
    from .vae import VAEEncoder
    from .prototype_learning import PrototypeLearningModel
    from .shared_prototype_learning import SharedPrototypeLearningModel
    from .direct_softmax_model import DirectSoftmaxModel
except ImportError:
    from vae import VAEEncoder
    from prototype_learning import PrototypeLearningModel
    from shared_prototype_learning import SharedPrototypeLearningModel
    from direct_softmax_model import DirectSoftmaxModel


logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def load_vae_encoder(
    checkpoint_path: str,
    device: torch.device,
) -> Tuple[VAEEncoder, int]:
    """Load frozen VAE encoder from checkpoint.

    Args:
        checkpoint_path: Path to VAE checkpoint.
        device: Device to load model on.

    Returns:
        encoder: Loaded VAE encoder (frozen).
        latent_dim: Latent dimensionality from checkpoint.
    """
    checkpoint = torch.load(checkpoint_path, map_location=device)

    # Get model parameters from checkpoint
    in_channels = checkpoint.get("in_channels", 3)
    latent_dim = checkpoint.get("latent_dim", 128)

    # Initialize and load encoder
    encoder = VAEEncoder(in_channels=in_channels, latent_dim=latent_dim)

    # Extract encoder state from full VAE state
    state_dict = checkpoint["model_state_dict"]
    encoder_state = {
        k.replace("encoder.", ""): v
        for k, v in state_dict.items()
        if k.startswith("encoder.")
    }
    encoder.load_state_dict(encoder_state)

    # Freeze encoder
    encoder = encoder.to(device)
    encoder.eval()
    for p in encoder.parameters():
        p.requires_grad = False

    logger.info(f"Loaded frozen encoder: in_channels={in_channels}, latent_dim={latent_dim}")

    return encoder, latent_dim


def load_spot_patches(
    patches_dir: Path,
    spot_id: str,
    device: torch.device,
) -> Optional[torch.Tensor]:
    """Load patches for a single spot.

    Args:
        patches_dir: Directory containing patch files.
        spot_id: Spot identifier.
        device: Device to load patches on.

    Returns:
        patches: (N, C, H, W) tensor or None if file not found.
    """
    # Try both naming conventions: {spot_id}_patches.npy and spot_{spot_id}_patches.npy
    patch_file = patches_dir / f"{spot_id}_patches.npy"
    if not patch_file.exists():
        patch_file = patches_dir / f"spot_{spot_id}_patches.npy"
    if not patch_file.exists():
        return None

    patches = np.load(patch_file).astype(np.float32)
    return torch.from_numpy(patches).to(device)


def train_prototypes(
    vae_checkpoint: str,
    patches_dir: str,
    proportions_file: str,
    output_dir: str,
    n_types: int = 7,
    latent_dim: int = 128,
    projection_dim: int = 32,
    epochs: int = 50,
    lr: float = 1e-3,
    sinkhorn_temp: float = 0.1,
    sinkhorn_iters: int = 50,
    device: str = "cuda",
    checkpoint_interval: int = 10,
    min_cells_per_spot: int = 3,
    use_cosine_distance: bool = True,
    spread_weight: float = 1.0,
    diversity_weight: float = 0.1,
    ortho_weight: float = 0.5,
    entropy_weight: float = 0.1,
    min_variance: float = 0.1,
    temp_anneal: bool = True,
    temp_start: float = 0.5,
    temp_end: float = 0.1,
    warmup_epochs: int = 5,
    grad_clip: float = 1.0,
    batch_size: int = 8,
    use_kmeans_init: bool = True,
    max_init_samples: int = 10000,
    # Shared space architecture parameters
    use_shared_space: bool = True,
    repulsion_weight: float = 0.5,
    repulsion_margin: float = 0.5,
    contrastive_weight: float = 0.1,
    contrastive_margin: float = 0.3,
    # Direct softmax parameters (recommended - no Sinkhorn)
    use_direct_softmax: bool = False,
    softmax_temperature: float = 0.1,
    # Attention aggregation parameters
    use_attention: bool = False,
    use_per_class_attention: bool = False,
    attention_entropy_weight: float = 0.1,
) -> Dict[str, List[float]]:
    """Train projection heads and prototypes.

    Supports three architectures:
    - Direct softmax (recommended): No Sinkhorn, KL loss on aggregated proportions
    - Shared space: Single projection head with Sinkhorn OT
    - Per-type heads (legacy): K separate projection heads (prone to collapse)

    Args:
        vae_checkpoint: Path to trained VAE checkpoint.
        patches_dir: Directory containing spot_{spot_id}_patches.npy files.
        proportions_file: CSV with spot_id and cell type proportion columns.
        output_dir: Directory to save model and history.
        n_types: Number of cell types.
        latent_dim: VAE latent dimensionality.
        projection_dim: Projection head output dimensionality.
        epochs: Number of training epochs.
        lr: Learning rate.
        sinkhorn_temp: Sinkhorn temperature (final value if annealing).
        sinkhorn_iters: Number of Sinkhorn iterations.
        device: Device to train on.
        checkpoint_interval: Save checkpoint every N epochs.
        min_cells_per_spot: Minimum cells in spot to include in training.
        use_cosine_distance: Use cosine distance (per-type heads only).
        spread_weight: Weight for projection spread regularization (per-type only).
        diversity_weight: Weight for prototype diversity regularization (per-type only).
        ortho_weight: Weight for projection head orthogonality (per-type only).
        entropy_weight: Weight for transport plan entropy (per-type only).
        min_variance: Minimum variance floor (per-type only).
        temp_anneal: Whether to anneal temperature from high to low.
        temp_start: Starting temperature (warmup and anneal start).
        temp_end: Ending temperature (if annealing).
        warmup_epochs: Number of epochs to hold temperature at temp_start.
        grad_clip: Maximum gradient norm for clipping.
        batch_size: Number of spots to batch together.
        use_kmeans_init: Whether to initialize prototypes from K-means.
        max_init_samples: Maximum samples for K-means initialization.
        use_shared_space: Use shared embedding space (recommended, prevents collapse).
        repulsion_weight: Weight for prototype repulsion loss (shared space only).
        repulsion_margin: Minimum distance between prototypes (shared space only).
        contrastive_weight: Weight for contrastive margin loss (shared space only).
        contrastive_margin: Margin for contrastive loss (shared space only).
        use_direct_softmax: Use direct softmax model (no Sinkhorn, recommended).
        softmax_temperature: Temperature for softmax assignments (direct softmax only).

    Returns:
        History dict with loss per epoch.
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    patches_path = Path(patches_dir)

    # Setup device
    if device == "cuda" and not torch.cuda.is_available():
        logger.warning("CUDA not available, falling back to CPU")
        device = "cpu"
    device = torch.device(device)
    logger.info(f"Using device: {device}")

    # Load VAE encoder
    encoder, vae_latent_dim = load_vae_encoder(vae_checkpoint, device)
    if latent_dim != vae_latent_dim:
        logger.warning(
            f"Specified latent_dim={latent_dim} differs from VAE's {vae_latent_dim}, "
            f"using VAE's latent_dim"
        )
        latent_dim = vae_latent_dim

    # Load proportions
    proportions_df = pd.read_csv(proportions_file)
    logger.info(f"Loaded proportions: {len(proportions_df)} spots")

    # Identify cell type columns (all columns except spot_id)
    celltype_cols = [c for c in proportions_df.columns if c != "spot_id"]
    if len(celltype_cols) != n_types:
        logger.warning(
            f"Found {len(celltype_cols)} cell type columns, but n_types={n_types}. "
            f"Using {len(celltype_cols)} types."
        )
        n_types = len(celltype_cols)

    logger.info(f"Cell type columns: {celltype_cols}")

    # Validate proportions (should sum to 1)
    prop_sums = proportions_df[celltype_cols].sum(axis=1)
    if not np.allclose(prop_sums, 1.0, atol=1e-3):
        logger.warning("Some proportion rows don't sum to 1, normalizing...")
        proportions_df[celltype_cols] = proportions_df[celltype_cols].div(
            prop_sums, axis=0
        )

    # Filter spots that have patches
    valid_spots = []
    for _, row in tqdm(proportions_df.iterrows(), desc="Checking patches", total=len(proportions_df)):
        spot_id = str(row["spot_id"])
        # Try both naming conventions: {spot_id}_patches.npy and spot_{spot_id}_patches.npy
        patch_file = patches_path / f"{spot_id}_patches.npy"
        if not patch_file.exists():
            patch_file = patches_path / f"spot_{spot_id}_patches.npy"
        if patch_file.exists():
            # Check number of patches
            patches = np.load(patch_file)
            if len(patches) >= min_cells_per_spot:
                valid_spots.append(spot_id)

    logger.info(f"Found {len(valid_spots)} valid spots with >= {min_cells_per_spot} cells")

    if len(valid_spots) == 0:
        raise ValueError("No valid spots found with patches")

    # Create mapping from spot_id to proportions
    proportions_df["spot_id"] = proportions_df["spot_id"].astype(str)
    spot_to_props = {
        row["spot_id"]: torch.tensor(
            row[celltype_cols].values.astype(np.float32),
            device=device
        )
        for _, row in proportions_df.iterrows()
    }

    # Initialize model based on architecture choice
    initial_temp = temp_start if temp_anneal else sinkhorn_temp

    if use_direct_softmax:
        # RECOMMENDED: Direct softmax (no Sinkhorn)
        model = DirectSoftmaxModel(
            encoder=encoder,
            n_types=n_types,
            latent_dim=latent_dim,
            projection_dim=projection_dim,
            temperature=softmax_temperature,
            repulsion_weight=repulsion_weight,
            repulsion_margin=repulsion_margin,
            entropy_weight=entropy_weight,
            use_cosine=True,
            use_attention=use_attention,
            use_per_class_attention=use_per_class_attention,
            attention_entropy_weight=attention_entropy_weight,
        )
        model = model.to(device)
        logger.info(
            f"Initialized DirectSoftmaxModel: n_types={n_types}, "
            f"latent_dim={latent_dim}, projection_dim={projection_dim}, "
            f"temperature={softmax_temperature}"
        )
        if use_attention:
            attn_type = "per-class" if use_per_class_attention else "shared"
            logger.info(f"  Attention: {attn_type}, entropy_weight={attention_entropy_weight}")
        logger.info(
            f"  Losses: KL (proportion matching), repulsion={repulsion_weight} "
            f"(margin={repulsion_margin}), entropy={entropy_weight}"
        )
    elif use_shared_space:
        # RECOMMENDED: Shared embedding space (prevents collapse)
        model = SharedPrototypeLearningModel(
            encoder=encoder,
            n_types=n_types,
            latent_dim=latent_dim,
            projection_dim=projection_dim,
            sinkhorn_temp=initial_temp,
            sinkhorn_iters=sinkhorn_iters,
            repulsion_weight=repulsion_weight,
            repulsion_margin=repulsion_margin,
            contrastive_weight=contrastive_weight,
            contrastive_margin=contrastive_margin,
        )
        model = model.to(device)
        logger.info(
            f"Initialized SharedPrototypeLearningModel: n_types={n_types}, "
            f"latent_dim={latent_dim}, projection_dim={projection_dim}"
        )
        logger.info(
            f"  Losses: repulsion={repulsion_weight} (margin={repulsion_margin}), "
            f"contrastive={contrastive_weight} (margin={contrastive_margin})"
        )
    else:
        # LEGACY: Per-type projection heads (prone to collapse)
        model = PrototypeLearningModel(
            encoder=encoder,
            n_types=n_types,
            latent_dim=latent_dim,
            projection_dim=projection_dim,
            sinkhorn_temp=initial_temp,
            sinkhorn_iters=sinkhorn_iters,
            use_cosine_distance=use_cosine_distance,
            spread_weight=spread_weight,
            diversity_weight=diversity_weight,
            ortho_weight=ortho_weight,
            entropy_weight=entropy_weight,
            min_variance=min_variance,
        )
        model = model.to(device)
        logger.info(
            f"Initialized PrototypeLearningModel (LEGACY): n_types={n_types}, "
            f"latent_dim={latent_dim}, projection_dim={projection_dim}, "
            f"cosine_dist={use_cosine_distance}"
        )
        logger.info(
            f"  Regularization: spread={spread_weight}, diversity={diversity_weight}, "
            f"ortho={ortho_weight}, entropy={entropy_weight}, min_var={min_variance}"
        )

    if temp_anneal:
        logger.info(f"Temperature: warmup {warmup_epochs} epochs at {temp_start}, then anneal to {temp_end}")
    if grad_clip > 0:
        logger.info(f"Gradient clipping: max_norm={grad_clip}")
    logger.info(f"Multi-spot batching: batch_size={batch_size}")

    # K-means initialization of prototypes
    if use_kmeans_init:
        logger.info("Initializing prototypes from K-means clustering...")
        all_embeddings = []
        sample_count = 0
        np.random.shuffle(valid_spots)

        for spot_id in tqdm(valid_spots, desc="Collecting embeddings for K-means"):
            if sample_count >= max_init_samples:
                break

            patches = load_spot_patches(patches_path, spot_id, device)
            if patches is None:
                continue

            with torch.no_grad():
                mu, _ = encoder(patches)
            all_embeddings.append(mu.cpu().numpy())
            sample_count += len(mu)

        all_z = np.concatenate(all_embeddings, axis=0)
        logger.info(f"Collected {len(all_z)} embeddings for K-means")

        # Run K-means clustering on VAE embeddings
        kmeans = KMeans(n_clusters=n_types, n_init=10, random_state=42).fit(all_z)

        if use_shared_space:
            # For shared space: project cluster centers through single head
            with torch.no_grad():
                cluster_centers = torch.tensor(kmeans.cluster_centers_, dtype=torch.float32, device=device)
                projected_centers = model.projection(cluster_centers)  # (K, D)
                model.prototypes.prototypes.data = projected_centers
            logger.info("Initialized prototypes from K-means (shared space)")
        else:
            # For per-type heads: project each center through its head
            with torch.no_grad():
                cluster_centers = torch.tensor(kmeans.cluster_centers_, dtype=torch.float32, device=device)
                for k in range(n_types):
                    proj_center = model.heads.heads[k](cluster_centers[k].unsqueeze(0))
                    model.prototypes.prototypes.data[k] = proj_center.squeeze()
            logger.info("Initialized prototypes from K-means (per-type heads)")

    # Optimizer (only for learnable components, encoder is frozen)
    if use_direct_softmax:
        optimizer = torch.optim.Adam(
            list(model.projection.parameters()) + [model.prototypes],
            lr=lr
        )
    elif use_shared_space:
        optimizer = torch.optim.Adam(
            list(model.projection.parameters()) + list(model.prototypes.parameters()),
            lr=lr
        )
    else:
        optimizer = torch.optim.Adam(
            list(model.heads.parameters()) + list(model.prototypes.parameters()),
            lr=lr
        )

    # Training history with detailed loss components for monitoring
    if use_direct_softmax:
        history = {
            "loss": [],
            "loss_std": [],
            "kl_loss": [],
            "repulsion_loss": [],
            "entropy": [],
            "prop_mae": [],
            "confidence": [],
            "proto_min_dist": [],
            "logit_std": [],
            "temperature": [],
        }
    elif use_shared_space:
        history = {
            "loss": [],
            "loss_std": [],
            "ot_loss": [],
            "spread_loss": [],
            "diversity_loss": [],
            "ortho_loss": [],
            "entropy_loss": [],
            "proj_variance": [],
            "transport_entropy": [],
            "temperature": [],
        }

    # Training loop
    logger.info(f"Starting training for {epochs} epochs over {len(valid_spots)} spots")
    for epoch in range(epochs):
        model.train()
        epoch_losses = []
        if use_direct_softmax:
            epoch_components = {
                "kl": [], "repulsion": [], "entropy": [],
                "prop_mae": [], "confidence": [], "proto_min_dist": [], "logit_std": []
            }
        elif use_shared_space:
            epoch_components = {
                "ot": [], "spread": [], "diversity": [], "ortho": [],
                "entropy": [], "proj_variance": [], "transport_entropy": []
            }

        # Temperature schedule: warmup then anneal
        if temp_anneal:
            if epoch < warmup_epochs:
                # Warmup phase: hold at temp_start
                current_temp = temp_start
            else:
                # Annealing phase: linear decay from temp_start to temp_end
                anneal_progress = (epoch - warmup_epochs) / max(epochs - warmup_epochs - 1, 1)
                current_temp = temp_start + anneal_progress * (temp_end - temp_start)
            model.sinkhorn_temp = current_temp
        else:
            current_temp = sinkhorn_temp

        if epoch % 10 == 0 or epoch < warmup_epochs:
            logger.info(f"Epoch {epoch+1}: Temperature = {current_temp:.4f}")

        # Shuffle spots each epoch
        np.random.shuffle(valid_spots)

        # Multi-spot batching
        pbar = tqdm(range(0, len(valid_spots), batch_size), desc=f"Epoch {epoch+1}/{epochs}")
        for batch_start in pbar:
            batch_spots = valid_spots[batch_start:batch_start + batch_size]

            # Accumulate gradients over batch
            optimizer.zero_grad()
            batch_loss = 0.0
            batch_nuclei = 0
            batch_valid_spots = 0

            for spot_id in batch_spots:
                # Load patches for this spot
                patches = load_spot_patches(patches_path, spot_id, device)
                if patches is None:
                    continue

                # Get proportions
                proportions = spot_to_props.get(spot_id)
                if proportions is None:
                    continue

                # Skip if too few cells
                n_cells = len(patches)
                if n_cells < min_cells_per_spot:
                    continue

                # Forward pass with component tracking
                components, _ = model(patches, proportions, return_components=True)

                # Weight loss by nuclei count
                weighted_loss = components["total"] * n_cells
                weighted_loss.backward()

                batch_loss += components["total"].item() * n_cells
                batch_nuclei += n_cells
                batch_valid_spots += 1

                # Track components (different keys for each architecture)
                for key in epoch_components.keys():
                    if key in components:
                        val = components[key]
                        epoch_components[key].append(val.item() if torch.is_tensor(val) else val)

            # Gradient clipping
            if grad_clip > 0 and batch_nuclei > 0:
                if use_direct_softmax:
                    params = list(model.projection.parameters()) + [model.prototypes]
                elif use_shared_space:
                    params = list(model.projection.parameters()) + list(model.prototypes.parameters())
                else:
                    params = list(model.heads.parameters()) + list(model.prototypes.parameters())
                torch.nn.utils.clip_grad_norm_(params, max_norm=grad_clip)

            # Optimizer step
            if batch_nuclei > 0:
                optimizer.step()
                epoch_losses.append(batch_loss / batch_nuclei)

            # Update progress bar
            if epoch_losses:
                pbar.set_postfix({"loss": f"{epoch_losses[-1]:.4f}"})

        # Epoch statistics
        avg_loss = np.mean(epoch_losses) if epoch_losses else 0.0
        std_loss = np.std(epoch_losses) if epoch_losses else 0.0

        history["loss"].append(avg_loss)
        history["loss_std"].append(std_loss)
        history["temperature"].append(current_temp)

        # Track loss components
        for key in epoch_components.keys():
            if epoch_components[key]:
                mean_val = np.mean(epoch_components[key])
            else:
                mean_val = 0.0

            # Map to history key
            if use_direct_softmax:
                if key in ["prop_mae", "confidence", "proto_min_dist", "logit_std", "entropy"]:
                    history_key = key
                else:
                    history_key = f"{key}_loss"
            elif use_shared_space:
                if key in ["dist_mean", "dist_std", "proto_min_dist", "transport_entropy"]:
                    history_key = key
                else:
                    history_key = f"{key}_loss"
            else:
                if key in ["proj_variance", "transport_entropy"]:
                    history_key = key
                else:
                    history_key = f"{key}_loss"

            history[history_key].append(mean_val)

        # Log with collapse detection
        logger.info(
            f"Epoch {epoch+1}/{epochs}: "
            f"loss={avg_loss:.4f} +/- {std_loss:.4f} "
            f"({len(epoch_losses)} batches)"
        )

        if use_direct_softmax:
            prop_mae = history["prop_mae"][-1]
            confidence = history["confidence"][-1]
            proto_min = history["proto_min_dist"][-1]
            logit_std = history["logit_std"][-1]

            logger.info(
                f"  Components: KL={history['kl_loss'][-1]:.4f}, "
                f"repulsion={history['repulsion_loss'][-1]:.4f}, "
                f"entropy={history['entropy'][-1]:.4f}"
            )
            logger.info(
                f"  Monitoring: prop_MAE={prop_mae:.4f}, confidence={confidence:.3f}, "
                f"proto_min={proto_min:.3f}, logit_std={logit_std:.3f}"
            )

            # Collapse warnings for direct softmax
            if logit_std < 0.1:
                logger.warning(f"  WARNING: Logit std {logit_std:.4f} < 0.1 - logits collapsing!")
            if proto_min < 0.1:
                logger.warning(f"  WARNING: Proto min dist {proto_min:.4f} < 0.1 - prototypes collapsing!")
            if confidence < 0.2:
                logger.warning(f"  WARNING: Confidence {confidence:.3f} < 0.2 - uncertain assignments!")

        elif use_shared_space:
            trans_ent = history["transport_entropy"][-1]
            dist_mean = history["dist_mean"][-1]
            dist_std = history["dist_std"][-1]
            proto_min = history["proto_min_dist"][-1]

            logger.info(
                f"  Components: OT={history['ot_loss'][-1]:.4f}, "
                f"repulsion={history['repulsion_loss'][-1]:.4f}, "
                f"contrastive={history['contrastive_loss'][-1]:.4f}"
            )
            logger.info(
                f"  Monitoring: dist={dist_mean:.3f}+/-{dist_std:.3f}, "
                f"proto_min_dist={proto_min:.3f}, entropy={trans_ent:.2f}"
            )

            # Collapse warnings for shared space
            if dist_std < 0.05:
                logger.warning(f"  WARNING: Distance std {dist_std:.4f} < 0.05 - distances collapsing!")
            if proto_min < 0.1:
                logger.warning(f"  WARNING: Proto min dist {proto_min:.4f} < 0.1 - prototypes collapsing!")
            if history['ot_loss'][-1] < 0.01 and epoch >= warmup_epochs:
                logger.warning(f"  WARNING: OT loss < 0.01 - possible trivial solution!")
        else:
            proj_var = history["proj_variance"][-1]
            trans_ent = history["transport_entropy"][-1]

            logger.info(
                f"  Components: OT={history['ot_loss'][-1]:.4f}, "
                f"spread={history['spread_loss'][-1]:.4f}, "
                f"ortho={history['ortho_loss'][-1]:.4f}, "
                f"entropy={history['entropy_loss'][-1]:.4f}"
            )
            logger.info(
                f"  Monitoring: proj_variance={proj_var:.4f}, transport_entropy={trans_ent:.4f}"
            )

            # Collapse warnings for per-type heads
            if proj_var < 0.01:
                logger.warning(f"  WARNING: Projection variance {proj_var:.4f} < 0.01 - heads may be collapsing!")
            if trans_ent < 0.5 and epoch >= warmup_epochs:
                logger.warning(f"  WARNING: Transport entropy {trans_ent:.4f} < 0.5 - degenerate assignments!")

        # Save checkpoint
        if (epoch + 1) % checkpoint_interval == 0:
            checkpoint_path = output_path / f"prototypes_checkpoint_epoch_{epoch+1}.pt"
            checkpoint_data = {
                "epoch": epoch + 1,
                "optimizer_state_dict": optimizer.state_dict(),
                "history": history,
                "n_types": n_types,
                "latent_dim": latent_dim,
                "projection_dim": projection_dim,
                "celltype_cols": celltype_cols,
                "use_direct_softmax": use_direct_softmax,
                "use_shared_space": use_shared_space,
            }
            if use_direct_softmax:
                checkpoint_data["projection_state_dict"] = model.projection.state_dict()
                checkpoint_data["prototypes"] = model.prototypes.data
                checkpoint_data["temperature"] = softmax_temperature
                checkpoint_data["repulsion_weight"] = repulsion_weight
                checkpoint_data["repulsion_margin"] = repulsion_margin
                checkpoint_data["entropy_weight"] = entropy_weight
            elif use_shared_space:
                checkpoint_data["prototypes_state_dict"] = model.prototypes.state_dict()
                checkpoint_data["sinkhorn_temp"] = current_temp
                checkpoint_data["sinkhorn_iters"] = sinkhorn_iters
                checkpoint_data["projection_state_dict"] = model.projection.state_dict()
                checkpoint_data["repulsion_weight"] = repulsion_weight
                checkpoint_data["repulsion_margin"] = repulsion_margin
                checkpoint_data["contrastive_weight"] = contrastive_weight
                checkpoint_data["contrastive_margin"] = contrastive_margin
            else:
                checkpoint_data["heads_state_dict"] = model.heads.state_dict()
                checkpoint_data["use_cosine_distance"] = use_cosine_distance
                checkpoint_data["spread_weight"] = spread_weight
                checkpoint_data["diversity_weight"] = diversity_weight
                checkpoint_data["ortho_weight"] = ortho_weight
                checkpoint_data["entropy_weight"] = entropy_weight
                checkpoint_data["min_variance"] = min_variance
            torch.save(checkpoint_data, checkpoint_path)
            logger.info(f"Saved checkpoint: {checkpoint_path}")

    # Save final model
    final_path = output_path / "prototypes_final.pt"
    final_data = {
        "epoch": epochs,
        "n_types": n_types,
        "latent_dim": latent_dim,
        "projection_dim": projection_dim,
        "celltype_cols": celltype_cols,
        "use_direct_softmax": use_direct_softmax,
        "use_shared_space": use_shared_space,
        "history": history,
    }
    if use_direct_softmax:
        final_data["projection_state_dict"] = model.projection.state_dict()
        final_data["prototypes"] = model.prototypes.data
        final_data["temperature"] = softmax_temperature
        final_data["repulsion_weight"] = repulsion_weight
        final_data["repulsion_margin"] = repulsion_margin
        final_data["entropy_weight"] = entropy_weight
    elif use_shared_space:
        final_data["prototypes_state_dict"] = model.prototypes.state_dict()
        final_data["sinkhorn_temp"] = temp_end if temp_anneal else sinkhorn_temp
        final_data["sinkhorn_iters"] = sinkhorn_iters
        final_data["projection_state_dict"] = model.projection.state_dict()
        final_data["repulsion_weight"] = repulsion_weight
        final_data["repulsion_margin"] = repulsion_margin
        final_data["contrastive_weight"] = contrastive_weight
        final_data["contrastive_margin"] = contrastive_margin
    else:
        final_data["heads_state_dict"] = model.heads.state_dict()
        final_data["use_cosine_distance"] = use_cosine_distance
        final_data["spread_weight"] = spread_weight
        final_data["diversity_weight"] = diversity_weight
        final_data["ortho_weight"] = ortho_weight
        final_data["entropy_weight"] = entropy_weight
        final_data["min_variance"] = min_variance
    torch.save(final_data, final_path)
    logger.info(f"Saved final model: {final_path}")

    # Save history
    history_path = output_path / "training_history.json"
    with open(history_path, "w") as f:
        json.dump(history, f, indent=2)
    logger.info(f"Saved training history: {history_path}")

    return history


def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Train projection heads and prototypes for cell type assignment",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--vae-checkpoint",
        type=str,
        required=True,
        help="Path to trained VAE checkpoint (vae_final.pt)"
    )
    parser.add_argument(
        "--patches-dir",
        type=str,
        required=True,
        help="Directory containing spot_{spot_id}_patches.npy files"
    )
    parser.add_argument(
        "--proportions-file",
        type=str,
        required=True,
        help="CSV file with spot_id and cell type proportion columns"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        required=True,
        help="Directory to save model and history"
    )
    parser.add_argument(
        "--n-types",
        type=int,
        default=7,
        help="Number of cell types"
    )
    parser.add_argument(
        "--latent-dim",
        type=int,
        default=128,
        help="VAE latent dimensionality"
    )
    parser.add_argument(
        "--projection-dim",
        type=int,
        default=32,
        help="Projection head output dimensionality"
    )
    parser.add_argument(
        "--epochs",
        type=int,
        default=50,
        help="Number of training epochs"
    )
    parser.add_argument(
        "--lr",
        type=float,
        default=5e-4,
        help="Learning rate"
    )
    parser.add_argument(
        "--sinkhorn-temp",
        type=float,
        default=0.1,
        help="Sinkhorn temperature (lower = sharper)"
    )
    parser.add_argument(
        "--sinkhorn-iters",
        type=int,
        default=50,
        help="Number of Sinkhorn iterations"
    )
    parser.add_argument(
        "--device",
        type=str,
        default="cuda",
        choices=["cuda", "cpu"],
        help="Device to train on"
    )
    parser.add_argument(
        "--checkpoint-interval",
        type=int,
        default=10,
        help="Save checkpoint every N epochs"
    )
    parser.add_argument(
        "--min-cells-per-spot",
        type=int,
        default=3,
        help="Minimum cells in spot to include in training"
    )
    # Anti-collapse regularization arguments
    parser.add_argument(
        "--use-cosine-distance",
        action="store_true",
        default=True,
        help="Use cosine distance (prevents collapse, recommended)"
    )
    parser.add_argument(
        "--no-cosine-distance",
        action="store_false",
        dest="use_cosine_distance",
        help="Use Euclidean distance instead of cosine"
    )
    parser.add_argument(
        "--spread-weight",
        type=float,
        default=1.0,
        help="Weight for projection spread regularization (variance floor hinge)"
    )
    parser.add_argument(
        "--diversity-weight",
        type=float,
        default=0.1,
        help="Weight for prototype diversity regularization"
    )
    parser.add_argument(
        "--ortho-weight",
        type=float,
        default=0.5,
        help="Weight for projection head orthogonality regularization"
    )
    parser.add_argument(
        "--entropy-weight",
        type=float,
        default=0.1,
        help="Weight for transport plan entropy regularization"
    )
    parser.add_argument(
        "--min-variance",
        type=float,
        default=0.1,
        help="Minimum variance floor for spread hinge loss"
    )
    parser.add_argument(
        "--temp-anneal",
        action="store_true",
        default=True,
        help="Anneal temperature from high to low (recommended)"
    )
    parser.add_argument(
        "--no-temp-anneal",
        action="store_false",
        dest="temp_anneal",
        help="Use fixed temperature throughout training"
    )
    parser.add_argument(
        "--temp-start",
        type=float,
        default=0.5,
        help="Starting temperature for warmup and annealing"
    )
    parser.add_argument(
        "--temp-end",
        type=float,
        default=0.1,
        help="Ending temperature for annealing"
    )
    parser.add_argument(
        "--warmup-epochs",
        type=int,
        default=5,
        help="Number of epochs to hold temperature at temp-start before annealing"
    )
    parser.add_argument(
        "--grad-clip",
        type=float,
        default=1.0,
        help="Maximum gradient norm for clipping (0 to disable)"
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=8,
        help="Number of spots to batch together"
    )
    parser.add_argument(
        "--use-kmeans-init",
        action="store_true",
        default=True,
        help="Initialize prototypes from K-means clustering (recommended)"
    )
    parser.add_argument(
        "--no-kmeans-init",
        action="store_false",
        dest="use_kmeans_init",
        help="Use random prototype initialization instead of K-means"
    )
    parser.add_argument(
        "--max-init-samples",
        type=int,
        default=10000,
        help="Maximum samples for K-means initialization"
    )
    # Shared space architecture arguments
    parser.add_argument(
        "--use-shared-space",
        action="store_true",
        default=True,
        help="Use shared embedding space architecture (recommended, prevents collapse)"
    )
    parser.add_argument(
        "--no-shared-space",
        action="store_false",
        dest="use_shared_space",
        help="Use legacy per-type heads architecture (prone to collapse)"
    )
    parser.add_argument(
        "--repulsion-weight",
        type=float,
        default=0.5,
        help="Weight for prototype repulsion loss (shared space only)"
    )
    parser.add_argument(
        "--repulsion-margin",
        type=float,
        default=0.5,
        help="Minimum distance between prototypes (shared space only)"
    )
    parser.add_argument(
        "--contrastive-weight",
        type=float,
        default=0.1,
        help="Weight for contrastive margin loss (shared space only)"
    )
    parser.add_argument(
        "--contrastive-margin",
        type=float,
        default=0.3,
        help="Margin for contrastive loss (shared space only)"
    )
    # Direct softmax arguments (recommended - no Sinkhorn)
    parser.add_argument(
        "--use-direct-softmax",
        action="store_true",
        default=False,
        help="Use direct softmax model (no Sinkhorn, recommended)"
    )
    parser.add_argument(
        "--softmax-temperature",
        type=float,
        default=0.1,
        help="Temperature for softmax assignments (direct softmax only)"
    )
    # Attention aggregation parameters
    parser.add_argument(
        "--use-attention",
        action="store_true",
        help="Use attention-weighted MIL aggregation instead of mean pooling"
    )
    parser.add_argument(
        "--use-per-class-attention",
        action="store_true",
        help="Use per-class (MoE) attention heads"
    )
    parser.add_argument(
        "--attention-entropy-weight",
        type=float,
        default=0.1,
        help="Weight for attention entropy regularization"
    )

    args = parser.parse_args()

    train_prototypes(
        vae_checkpoint=args.vae_checkpoint,
        patches_dir=args.patches_dir,
        proportions_file=args.proportions_file,
        output_dir=args.output_dir,
        n_types=args.n_types,
        latent_dim=args.latent_dim,
        projection_dim=args.projection_dim,
        epochs=args.epochs,
        lr=args.lr,
        sinkhorn_temp=args.sinkhorn_temp,
        sinkhorn_iters=args.sinkhorn_iters,
        device=args.device,
        checkpoint_interval=args.checkpoint_interval,
        min_cells_per_spot=args.min_cells_per_spot,
        use_cosine_distance=args.use_cosine_distance,
        spread_weight=args.spread_weight,
        diversity_weight=args.diversity_weight,
        ortho_weight=args.ortho_weight,
        entropy_weight=args.entropy_weight,
        min_variance=args.min_variance,
        temp_anneal=args.temp_anneal,
        temp_start=args.temp_start,
        temp_end=args.temp_end,
        warmup_epochs=args.warmup_epochs,
        grad_clip=args.grad_clip,
        batch_size=args.batch_size,
        use_kmeans_init=args.use_kmeans_init,
        max_init_samples=args.max_init_samples,
        # Shared space arguments
        use_shared_space=args.use_shared_space,
        repulsion_weight=args.repulsion_weight,
        repulsion_margin=args.repulsion_margin,
        contrastive_weight=args.contrastive_weight,
        contrastive_margin=args.contrastive_margin,
        # Direct softmax arguments
        use_direct_softmax=args.use_direct_softmax,
        softmax_temperature=args.softmax_temperature,
        # Attention aggregation arguments
        use_attention=args.use_attention,
        use_per_class_attention=args.use_per_class_attention,
        attention_entropy_weight=args.attention_entropy_weight,
    )


if __name__ == "__main__":
    main()
