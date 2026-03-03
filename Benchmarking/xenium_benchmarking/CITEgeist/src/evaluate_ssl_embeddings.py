"""Evaluate SSL embeddings quality and classification performance.

Computes:
- Silhouette score
- t-SNE visualization
- XGBoost classification accuracy

Usage:
    python evaluate_ssl_embeddings.py \
        --embeddings-dir output/ssl_embeddings \
        --output-dir output/ssl_evaluation
"""
import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_score, classification_report, accuracy_score
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.model_selection import train_test_split

REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

try:
    import xgboost as xgb
    HAS_XGB = True
except ImportError:
    from sklearn.ensemble import GradientBoostingClassifier
    HAS_XGB = False

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

PSEUDOVISIUM_DIR = REPO_ROOT / "Benchmarking/xenium_pseudovisium/data_protein_gt"
XENIUM_DIR = Path("/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma")

CELL_TYPES = ["B cells", "CD4+ T cells", "CD8+ T cells", "Macrophages", "Endothelial", "Epithelial", "Fibroblasts"]


def load_ground_truth():
    """Load cell type ground truth."""
    gt_df = pd.read_csv(PSEUDOVISIUM_DIR / "cell_type_assignments.csv", index_col=0)
    gt_df = gt_df[gt_df["cell_type"].isin(CELL_TYPES)]
    return gt_df


def evaluate_embeddings(
    embeddings_dir: str,
    output_dir: str,
    model_names: list = None,
):
    """Evaluate all embedding files in directory."""
    embeddings_path = Path(embeddings_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Load ground truth
    gt_df = load_ground_truth()
    labels = gt_df["cell_type"].values

    if model_names is None:
        model_names = ["mae", "dino", "vae"]

    results = {}

    for model_name in model_names:
        emb_file = embeddings_path / f"{model_name}_embeddings.npy"
        if not emb_file.exists():
            logger.warning(f"Embeddings not found: {emb_file}")
            continue

        logger.info(f"Evaluating {model_name}...")
        embeddings = np.load(emb_file)

        # Match to GT (assume same order)
        n_samples = min(len(embeddings), len(labels))
        embeddings = embeddings[:n_samples]
        sample_labels = labels[:n_samples]

        # Standardize
        scaler = StandardScaler()
        embeddings_scaled = scaler.fit_transform(embeddings)
        embeddings_scaled = np.nan_to_num(embeddings_scaled)

        # Silhouette score
        try:
            sil_score = silhouette_score(embeddings_scaled, sample_labels)
        except Exception as e:
            logger.warning(f"Silhouette failed: {e}")
            sil_score = -1.0

        logger.info(f"  Silhouette score: {sil_score:.4f}")

        # XGBoost classification
        le = LabelEncoder()
        y = le.fit_transform(sample_labels)

        X_train, X_test, y_train, y_test = train_test_split(
            embeddings_scaled, y, test_size=0.2, stratify=y, random_state=42
        )

        if HAS_XGB:
            clf = xgb.XGBClassifier(n_estimators=100, max_depth=6, random_state=42, n_jobs=-1)
        else:
            clf = GradientBoostingClassifier(n_estimators=100, max_depth=6, random_state=42)

        clf.fit(X_train, y_train)
        y_pred = clf.predict(X_test)
        accuracy = accuracy_score(y_test, y_pred)

        logger.info(f"  Classification accuracy: {accuracy:.4f}")

        # Per-class accuracy
        report = classification_report(y_test, y_pred, target_names=le.classes_, output_dict=True)

        results[model_name] = {
            "silhouette": float(sil_score),
            "accuracy": float(accuracy),
            "n_samples": n_samples,
            "embed_dim": embeddings.shape[1],
            "per_class": {ct: report[ct]["recall"] for ct in le.classes_ if ct in report},
        }

        # t-SNE visualization
        logger.info(f"  Computing t-SNE...")
        tsne = TSNE(n_components=2, random_state=42, perplexity=30)
        emb_2d = tsne.fit_transform(embeddings_scaled[:5000])  # Subsample for speed

        plt.figure(figsize=(10, 8))
        for ct in CELL_TYPES:
            mask = sample_labels[:5000] == ct
            if mask.sum() > 0:
                plt.scatter(emb_2d[mask, 0], emb_2d[mask, 1], label=ct, alpha=0.5, s=10)
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.title(f"{model_name.upper()} t-SNE (silhouette={sil_score:.3f}, acc={accuracy:.3f})")
        plt.tight_layout()
        plt.savefig(output_path / f"{model_name}_tsne.png", dpi=150)
        plt.close()

    # Save results
    results_file = output_path / "evaluation_results.json"
    with open(results_file, "w") as f:
        json.dump(results, f, indent=2)
    logger.info(f"Saved results to {results_file}")

    # Print summary
    logger.info("\n" + "="*60)
    logger.info("SUMMARY")
    logger.info("="*60)
    for model, res in results.items():
        logger.info(f"{model}: silhouette={res['silhouette']:.4f}, accuracy={res['accuracy']:.4f}")

    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--embeddings-dir", type=str, required=True)
    parser.add_argument("--output-dir", type=str, required=True)
    parser.add_argument("--models", nargs="+", default=["mae", "dino", "vae"])

    args = parser.parse_args()

    evaluate_embeddings(
        embeddings_dir=args.embeddings_dir,
        output_dir=args.output_dir,
        model_names=args.models,
    )


if __name__ == "__main__":
    main()
