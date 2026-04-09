#!/usr/bin/env python3
"""
Constrained single-cell assignment module.

Assigns nuclei to cell types given spot-level counts from protein deconvolution.
Supports three assignment modes:
- random: Constrained random shuffle within counts (baseline, ~22% accuracy)
- morphology: 7-class morphology-based Hungarian with Gaussian likelihood (~46% accuracy)
- xgboost: 7-class assignment using XGBoost + VAE + comprehensive morphology (~50%+ accuracy)

The key insight is that global morphology classification fails (~25% accuracy),
but constrained assignment given spot counts achieves significant gains by
converting "what type is this cell?" to "given counts, which assignment?".
"""

import logging
import pickle
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.optimize import linear_sum_assignment
from scipy.stats import multivariate_normal
from sklearn.preprocessing import StandardScaler

from .morphology_features import largest_remainder_discretize

logger = logging.getLogger(__name__)

# Coarse label mapping for 'groups' mode
COARSE_MAPPING = {
    "B cells": "Immune",
    "CD4+ T cells": "Immune",
    "CD8+ T cells": "Immune",
    "Macrophages": "Immune",
    "Fibroblasts": "Stromal",
    "Endothelial": "Stromal",
    "Epithelial": "Epithelial",
}

COARSE_LABELS = ["Immune", "Stromal", "Epithelial"]

# Alias for API consistency
proportions_to_counts = largest_remainder_discretize


def extract_patch_features(patch: np.ndarray) -> np.ndarray:
    """
    Extract morphological features from a 2-channel patch (DAPI + boundary).

    Args:
        patch: (2, H, W) array with DAPI channel 0, boundary channel 1

    Returns:
        Feature vector (13 features)
    """
    dapi = patch[0]
    boundary = patch[1]

    # Basic intensity features
    features = [
        dapi.mean(),
        dapi.std(),
        dapi.max(),
        boundary.mean(),
        boundary.std(),
        boundary.max(),
    ]

    # Nuclear area proxy (pixels above mean + std)
    nuc_mask = dapi > (dapi.mean() + dapi.std())
    nuc_area = nuc_mask.sum()
    features.append(nuc_area)

    # Shape features from nuclear mask
    if nuc_area > 10:
        rows = np.any(nuc_mask, axis=1)
        cols = np.any(nuc_mask, axis=0)
        if rows.any() and cols.any():
            rmin, rmax = np.where(rows)[0][[0, -1]]
            cmin, cmax = np.where(cols)[0][[0, -1]]
            height = rmax - rmin + 1
            width = cmax - cmin + 1
            aspect_ratio = width / max(height, 1)
            extent = nuc_area / max(height * width, 1)
            features.extend([aspect_ratio, extent, height, width])
        else:
            features.extend([1.0, 0.5, 48, 48])
    else:
        features.extend([1.0, 0.5, 48, 48])

    # Boundary/membrane edge signal
    boundary_edge = boundary.copy()
    boundary_edge[10:-10, 10:-10] = 0
    edge_signal = boundary_edge.mean()
    features.append(edge_signal)

    # Channel correlation
    corr = np.corrcoef(dapi.flatten(), boundary.flatten())[0, 1]
    if np.isnan(corr):
        corr = 0
    features.append(corr)

    return np.array(features, dtype=np.float32)


def collapse_to_coarse(
    fine_counts: np.ndarray,
    celltype_names: List[str],
) -> np.ndarray:
    """Collapse fine-grained counts to coarse categories."""
    coarse_counts = np.zeros(3)  # Immune, Stromal, Epithelial

    for i, ct in enumerate(celltype_names):
        coarse = COARSE_MAPPING.get(ct)
        if coarse == "Immune":
            coarse_counts[0] += fine_counts[i]
        elif coarse == "Stromal":
            coarse_counts[1] += fine_counts[i]
        elif coarse == "Epithelial":
            coarse_counts[2] += fine_counts[i]

    return coarse_counts


def expand_coarse_to_fine(
    coarse_assignments: np.ndarray,
    fine_counts: np.ndarray,
    celltype_names: List[str],
) -> np.ndarray:
    """
    Expand coarse assignments to fine-grained types via random sub-assignment.

    Args:
        coarse_assignments: (N_nuclei,) with values 0=Immune, 1=Stromal, 2=Epithelial
        fine_counts: (N_types,) integer counts per fine type
        celltype_names: List of fine-grained cell type names

    Returns:
        (N_nuclei,) fine-grained assignments
    """
    n_nuclei = len(coarse_assignments)
    fine_assignments = np.zeros(n_nuclei, dtype=int)

    # Group fine types by coarse category
    coarse_to_fine = {0: [], 1: [], 2: []}  # Immune, Stromal, Epithelial
    for i, ct in enumerate(celltype_names):
        coarse = COARSE_MAPPING.get(ct)
        if coarse == "Immune":
            coarse_to_fine[0].append((i, fine_counts[i]))
        elif coarse == "Stromal":
            coarse_to_fine[1].append((i, fine_counts[i]))
        elif coarse == "Epithelial":
            coarse_to_fine[2].append((i, fine_counts[i]))

    # For each coarse category, randomly assign to fine types
    for coarse_idx in range(3):
        nuclei_mask = coarse_assignments == coarse_idx
        nuclei_indices = np.where(nuclei_mask)[0]
        n_in_coarse = len(nuclei_indices)

        if n_in_coarse == 0:
            continue

        # Build pool of fine type labels
        fine_pool = []
        for fine_idx, count in coarse_to_fine[coarse_idx]:
            fine_pool.extend([fine_idx] * int(count))

        # Handle mismatch
        if len(fine_pool) < n_in_coarse:
            if coarse_to_fine[coarse_idx]:
                most_common = max(coarse_to_fine[coarse_idx], key=lambda x: x[1])[0]
                fine_pool.extend([most_common] * (n_in_coarse - len(fine_pool)))
        fine_pool = fine_pool[:n_in_coarse]

        # Shuffle and assign
        np.random.shuffle(fine_pool)
        fine_assignments[nuclei_indices] = fine_pool

    return fine_assignments


def hungarian_assign(log_likes: np.ndarray, counts: np.ndarray) -> np.ndarray:
    """
    Assign samples to types using Hungarian algorithm with count constraints.

    Args:
        log_likes: (N, K) log-likelihoods for each sample and type
        counts: (K,) integer counts for each type

    Returns:
        (N,) assignments - type index for each sample
    """
    N = log_likes.shape[0]
    K = len(counts)

    # Adjust counts to sum to N
    counts = counts.copy().astype(int)
    while counts.sum() < N:
        counts[np.argmax(counts)] += 1
    while counts.sum() > N:
        idx = np.where(counts > 0)[0]
        counts[idx[np.argmin(counts[idx])]] -= 1

    # Build expanded cost matrix
    expanded_cols = []
    col_to_type = []

    for k in range(K):
        for _ in range(int(counts[k])):
            # Negative because we maximize likelihood but minimize cost
            expanded_cols.append(-log_likes[:, k])
            col_to_type.append(k)

    if len(expanded_cols) == 0:
        return np.zeros(N, dtype=int)

    cost_matrix = np.column_stack(expanded_cols)

    # Solve assignment problem
    row_ind, col_ind = linear_sum_assignment(cost_matrix)

    # Map column indices back to type indices
    assignments = np.array([col_to_type[c] for c in col_ind])

    return assignments


def random_assign(counts: np.ndarray, n_samples: int) -> np.ndarray:
    """
    Random assignment respecting count constraints.

    Args:
        counts: (K,) integer counts for each type
        n_samples: total number of samples

    Returns:
        (N,) random assignments
    """
    counts = counts.copy().astype(int)

    # Adjust counts to sum to n_samples
    while counts.sum() < n_samples:
        counts[np.argmax(counts)] += 1
    while counts.sum() > n_samples:
        idx = np.where(counts > 0)[0]
        counts[idx[np.argmin(counts[idx])]] -= 1

    assignments = []
    for k, c in enumerate(counts):
        assignments.extend([k] * int(c))

    np.random.shuffle(assignments)
    return np.array(assignments)


class ConstrainedAssignment:
    """
    Assigns nuclei to cell types given spot-level counts.

    Supports three modes:
    - 'random': Constrained random shuffle (baseline, ~22% accuracy)
    - 'morphology': Morphology-based Hungarian with Gaussian likelihood (~46% accuracy)
    - 'xgboost': XGBoost classifier with VAE + comprehensive morphology (~50%+ accuracy)

    Example (morphology mode):
        assigner = ConstrainedAssignment(mode='morphology')
        assigner.fit(patches_dir, proportions_df)
        assignments = assigner.assign_spot(nuclei_features, cell_counts)

    Example (xgboost mode):
        assigner = ConstrainedAssignment(mode='xgboost')
        assigner.load_xgboost_model(model_path, vae_checkpoint_path)
        assignments = assigner.assign_spot(nuclei_patches, cell_counts)
    """

    def __init__(self, mode: str = "random", seed: int = 42):
        """
        Initialize assigner.

        Args:
            mode: 'random' | 'morphology' | 'xgboost'
            seed: Random seed for reproducibility
        """
        if mode not in ("random", "morphology", "xgboost"):
            raise ValueError(f"Unknown mode: {mode}. Use 'random', 'morphology', or 'xgboost'")

        self.mode = mode
        self.seed = seed
        np.random.seed(seed)

        # Learned distributions (for 'morphology' mode)
        self.fine_params: Optional[Dict[int, Tuple[np.ndarray, np.ndarray]]] = None
        self.fine_scaler: Optional[StandardScaler] = None
        self.celltype_names: Optional[List[str]] = None

        # XGBoost model components (for 'xgboost' mode)
        self.xgb_model = None
        self.xgb_scaler: Optional[StandardScaler] = None
        self.xgb_label_encoder = None
        self.vae_encoder = None
        self.vae_device = None

        self._fitted = False

    def fit(
        self,
        patches_dir: Path,
        proportions_df: pd.DataFrame,
        max_spots: int = 2000,
        purity_threshold: float = 0.5,
        layout: str = "benchmark",
    ) -> "ConstrainedAssignment":
        """
        Learn class-conditional distributions from high-purity spots.

        Args:
            patches_dir: Directory containing patch files.
            proportions_df: DataFrame with celltype proportion columns.
                For "benchmark" layout: must have spot_id, region columns.
                For "citegeist" layout: must have spot_id column (index or column),
                    celltype columns only (no region column needed).
            max_spots: Maximum spots to use for training
            purity_threshold: Minimum proportion for "high-purity" spots
            layout: "benchmark" (region_*/spot_*_patches.npy) or
                    "citegeist" ({spot_id}_patches.npy flat directory)

        Returns:
            self (for chaining)
        """
        if layout not in ("benchmark", "citegeist"):
            raise ValueError(f"Unknown layout: {layout}. Use 'benchmark' or 'citegeist'.")

        if self.mode == "random":
            logger.info("Random mode - no fitting required")
            self._fitted = True
            return self

        # Identify celltype columns
        non_celltype_cols = {"spot_id", "region"}
        celltype_cols = [
            c for c in proportions_df.columns if c not in non_celltype_cols
        ]
        self.celltype_names = celltype_cols
        n_types = len(celltype_cols)

        logger.info(f"Fitting {self.mode} mode with {n_types} cell types")

        # Collect features grouped by type
        type_features: Dict[int, List] = {i: [] for i in range(n_types)}

        # Find patch files based on layout
        patches_dir = Path(patches_dir)

        if layout == "citegeist":
            spot_files = sorted(patches_dir.glob("*_patches.npy"))[:max_spots]
            logger.info(f"Found {len(spot_files)} patch files (citegeist layout)")

            # Build spot_id lookup from proportions
            if "spot_id" in proportions_df.columns:
                prop_lookup = proportions_df.set_index("spot_id")
            else:
                prop_lookup = proportions_df  # assume index is spot_id

            for patch_file in spot_files:
                spot_id = patch_file.stem.replace("_patches", "")

                if spot_id not in prop_lookup.index:
                    continue

                props = prop_lookup.loc[spot_id, celltype_cols].values.astype(float)
                dominant = int(np.argmax(props))
                confidence = props[dominant]

                if confidence < purity_threshold:
                    continue

                try:
                    patches = np.load(patch_file).astype(np.float32)
                except (IOError, ValueError) as e:
                    logger.warning(f"Failed to load {patch_file}: {e}")
                    continue

                for patch in patches[:10]:
                    try:
                        feats = extract_patch_features(patch)
                        if not np.any(np.isnan(feats)):
                            type_features[dominant].append(feats)
                    except (ValueError, IndexError):
                        continue
        else:
            # benchmark layout: region_*/spot_*_patches.npy
            spot_files = []
            for region_dir in sorted(patches_dir.glob("region_*")):
                spot_files.extend(list(region_dir.glob("spot_*_patches.npy")))
            spot_files = spot_files[:max_spots]
            logger.info(f"Found {len(spot_files)} patch files (benchmark layout)")

            for patch_file in spot_files:
                region = int(patch_file.parent.name.split("_")[1])
                stem = patch_file.stem.replace("_patches", "")
                if stem.startswith("spot_"):
                    spot_id = stem[5:]
                else:
                    spot_id = stem

                spot_row = proportions_df[
                    (proportions_df["spot_id"] == spot_id)
                    & (proportions_df["region"] == region)
                ]
                if len(spot_row) == 0:
                    continue

                props = spot_row[celltype_cols].values[0]
                dominant = int(np.argmax(props))
                confidence = props[dominant]

                if confidence < purity_threshold:
                    continue

                try:
                    patches = np.load(patch_file).astype(np.float32)
                except (IOError, ValueError) as e:
                    logger.warning(f"Failed to load {patch_file}: {e}")
                    continue

                for patch in patches[:10]:
                    try:
                        feats = extract_patch_features(patch)
                        if not np.any(np.isnan(feats)):
                            type_features[dominant].append(feats)
                    except (ValueError, IndexError):
                        continue

        # Fit scaler on all features
        all_features = []
        for i in type_features:
            all_features.extend(type_features[i])

        if len(all_features) == 0:
            raise ValueError("No features collected! Check patches_dir and proportions_df.")

        all_features = np.array(all_features)
        scaler = StandardScaler()
        scaler.fit(all_features)

        # Fit Gaussian for each type
        class_params = {}
        n_feats = all_features.shape[1]

        for i, label in enumerate(celltype_cols):
            n_samples = len(type_features[i])
            logger.info(f"  {label}: {n_samples} samples")

            if n_samples > 20:
                feats = scaler.transform(np.array(type_features[i]))
                mean = feats.mean(axis=0)
                cov = np.cov(feats.T) + 1e-6 * np.eye(n_feats)
                class_params[i] = (mean, cov)
            else:
                logger.warning(f"  {label}: insufficient samples, using fallback")
                class_params[i] = (np.zeros(n_feats), np.eye(n_feats))

        # Store results
        self.fine_params = class_params
        self.fine_scaler = scaler

        self._fitted = True
        return self

    def load_xgboost_model(
        self,
        model_path: Union[str, Path],
        vae_checkpoint_path: Optional[Union[str, Path]] = None,
        device: str = "cuda",
    ) -> "ConstrainedAssignment":
        """
        Load a pre-trained XGBoost model for assignment.

        Args:
            model_path: Path to saved model (.pkl file containing model, scaler, label_encoder)
            vae_checkpoint_path: Path to VAE checkpoint (.pt file) for feature extraction.
                                 If None, assumes features are pre-extracted.
            device: Device for VAE inference ('cuda' or 'cpu')

        Returns:
            self (for chaining)
        """
        import torch

        model_path = Path(model_path)
        if not model_path.exists():
            raise FileNotFoundError(f"Model not found: {model_path}")

        # Load the model bundle
        with open(model_path, "rb") as f:
            bundle = pickle.load(f)

        self.xgb_model = bundle["model"]
        self.xgb_scaler = bundle["scaler"]
        self.xgb_label_encoder = bundle["label_encoder"]
        self.celltype_names = list(self.xgb_label_encoder.classes_)

        logger.info(f"Loaded XGBoost model with {len(self.celltype_names)} classes")

        # Load VAE encoder if provided
        if vae_checkpoint_path is not None:
            vae_checkpoint_path = Path(vae_checkpoint_path)
            if not vae_checkpoint_path.exists():
                raise FileNotFoundError(f"VAE checkpoint not found: {vae_checkpoint_path}")

            from .vae import VAEEncoder

            self.vae_device = torch.device(device if torch.cuda.is_available() else "cpu")
            checkpoint = torch.load(vae_checkpoint_path, map_location=self.vae_device)

            in_channels = checkpoint.get("in_channels", 2)
            latent_dim = checkpoint.get("latent_dim", 128)

            encoder = VAEEncoder(in_channels=in_channels, latent_dim=latent_dim)
            state_dict = checkpoint["model_state_dict"]
            encoder_state = {
                k.replace("encoder.", ""): v
                for k, v in state_dict.items()
                if k.startswith("encoder.")
            }
            encoder.load_state_dict(encoder_state)
            encoder = encoder.to(self.vae_device)
            encoder.eval()

            self.vae_encoder = encoder
            logger.info(f"Loaded VAE encoder: latent_dim={latent_dim}")

        self._fitted = True
        return self

    def _extract_xgboost_features(
        self,
        patches: np.ndarray,
        nucleus_masks: Optional[np.ndarray] = None,
        cell_masks: Optional[np.ndarray] = None,
    ) -> np.ndarray:
        """
        Extract combined VAE + morphology features for XGBoost prediction.

        Args:
            patches: (N, 2, H, W) nucleus patches
            nucleus_masks: Optional (N, H, W) nucleus masks for morphology
            cell_masks: Optional (N, H, W) cell masks for morphology

        Returns:
            (N, D) combined feature array
        """
        import torch
        from .comprehensive_morphology import extract_from_patch

        n_samples = len(patches)
        features_list = []

        # Extract VAE embeddings if encoder is loaded
        if self.vae_encoder is not None:
            with torch.no_grad():
                batch_size = 256
                vae_embeddings = []
                for i in range(0, n_samples, batch_size):
                    batch = torch.from_numpy(patches[i:i+batch_size]).float().to(self.vae_device)
                    mu, _ = self.vae_encoder(batch)
                    vae_embeddings.append(mu.cpu().numpy())
                vae_embeddings = np.concatenate(vae_embeddings, axis=0)
            features_list.append(vae_embeddings)

        # Extract morphology features
        morph_features = []
        for i in range(n_samples):
            patch = patches[i]
            nuc_mask = nucleus_masks[i] if nucleus_masks is not None else None
            cell_mask = cell_masks[i] if cell_masks is not None else None
            morph_feat = extract_from_patch(patch, nuc_mask, cell_mask)
            morph_features.append(morph_feat)
        morph_features = np.array(morph_features)
        features_list.append(morph_features)

        # Combine features
        combined = np.hstack(features_list)

        # Scale features
        if self.xgb_scaler is not None:
            combined = self.xgb_scaler.transform(combined)

        # Handle NaN/Inf
        combined = np.nan_to_num(combined, nan=0.0, posinf=0.0, neginf=0.0)

        return combined

    def _compute_xgboost_log_probs(self, features: np.ndarray) -> np.ndarray:
        """
        Get log-probabilities from XGBoost model.

        Args:
            features: (N, D) scaled feature array

        Returns:
            (N, K) log-probability matrix
        """
        # Get probability predictions
        probs = self.xgb_model.predict_proba(features)

        # Convert to log-probs (add small epsilon to avoid log(0))
        log_probs = np.log(probs + 1e-10)

        return log_probs

    def _compute_log_likelihoods(
        self,
        features: np.ndarray,
        params: Dict[int, Tuple[np.ndarray, np.ndarray]],
        scaler: StandardScaler,
    ) -> np.ndarray:
        """Compute log-likelihood of features for each type."""
        features_scaled = scaler.transform(features)
        n_nuclei = len(features)
        n_types = len(params)

        log_likes = np.zeros((n_nuclei, n_types))

        for k, (mean, cov) in params.items():
            try:
                rv = multivariate_normal(mean=mean, cov=cov, allow_singular=True)
                log_likes[:, k] = rv.logpdf(features_scaled)
            except Exception:
                log_likes[:, k] = -100

        return log_likes

    def assign_spot(
        self,
        nuclei_features: np.ndarray,
        cell_counts: np.ndarray,
        celltype_names: Optional[List[str]] = None,
        patches: Optional[np.ndarray] = None,
        nucleus_masks: Optional[np.ndarray] = None,
        cell_masks: Optional[np.ndarray] = None,
    ) -> np.ndarray:
        """
        Assign nuclei to cell types for a single spot.

        Args:
            nuclei_features: (N_nuclei, N_features) feature matrix.
                For 'random' and 'morphology' modes this contains the features.
                For 'xgboost' mode this is ignored (use patches instead).
            cell_counts: (N_types,) integer counts per type
            celltype_names: List of cell type names (unused, kept for API compatibility)
            patches: Optional (N_nuclei, 2, H, W) nucleus patches for xgboost mode.
                Required when mode='xgboost'.
            nucleus_masks: Optional (N_nuclei, H, W) nucleus masks for xgboost morphology.
            cell_masks: Optional (N_nuclei, H, W) cell masks for xgboost morphology.

        Returns:
            (N_nuclei,) array of cell type indices
        """
        if not self._fitted and self.mode != "random":
            raise RuntimeError("Must call fit() or load_xgboost_model() before assign_spot()")

        # Determine number of nuclei
        if self.mode == "xgboost":
            if patches is None:
                raise ValueError("patches required for xgboost mode")
            n_nuclei = len(patches)
        else:
            n_nuclei = len(nuclei_features)

        if n_nuclei == 0:
            return np.array([], dtype=int)

        if self.mode == "random":
            return random_assign(cell_counts, n_nuclei)

        elif self.mode == "xgboost":
            # Extract combined VAE + morphology features
            features = self._extract_xgboost_features(patches, nucleus_masks, cell_masks)
            # Get log-probabilities from XGBoost
            log_probs = self._compute_xgboost_log_probs(features)
            # Use Hungarian assignment with XGBoost probabilities
            return hungarian_assign(log_probs, cell_counts)

        else:  # morphology
            log_likes = self._compute_log_likelihoods(
                nuclei_features, self.fine_params, self.fine_scaler
            )
            return hungarian_assign(log_likes, cell_counts)

    def create_single_cell_adata(
        self,
        spot_adata: sc.AnnData,
        assignments: Dict[str, np.ndarray],
        nuclei_centroids: Dict[str, np.ndarray],
        deconvolved_gex: pd.DataFrame,
        celltype_names: List[str],
    ) -> sc.AnnData:
        """
        Create single-cell AnnData from spot-level data and assignments.

        Args:
            spot_adata: Original spot-level AnnData
            assignments: Dict[spot_id] -> (N_nuclei,) assignment indices
            nuclei_centroids: Dict[spot_id] -> (N_nuclei, 2) coordinates
            deconvolved_gex: DataFrame with index 'barcode:::celltype', columns = genes
            celltype_names: List of cell type names

        Returns:
            Single-cell AnnData with X = deconvolved GEX per cell
        """
        cell_data = []
        cell_gex = []
        gene_names = deconvolved_gex.columns.tolist()

        for spot_id, spot_assignments in assignments.items():
            centroids = nuclei_centroids.get(spot_id, np.array([]))
            n_nuclei = len(spot_assignments)

            if len(centroids) != n_nuclei:
                logger.warning(
                    f"Spot {spot_id}: centroid/assignment mismatch "
                    f"({len(centroids)} vs {n_nuclei})"
                )
                centroids = np.zeros((n_nuclei, 2))

            # Count cells per type in this spot
            type_counts = np.bincount(spot_assignments, minlength=len(celltype_names))

            for i, (assignment, centroid) in enumerate(zip(spot_assignments, centroids)):
                celltype = celltype_names[assignment]
                coarse = COARSE_MAPPING.get(celltype, "Unknown")

                # Get deconvolved GEX for this cell type in this spot
                gex_key = f"{spot_id}:::{celltype}"
                if gex_key in deconvolved_gex.index:
                    gex = deconvolved_gex.loc[gex_key].values.copy()
                    # Divide by number of cells of this type
                    n_of_type = type_counts[assignment]
                    if n_of_type > 0:
                        gex = gex / n_of_type
                else:
                    gex = np.zeros(len(gene_names))

                cell_data.append({
                    "cell_id": f"{spot_id}_cell_{i}",
                    "spot_id": spot_id,
                    "cell_type": celltype,
                    "cell_type_coarse": coarse,
                    "cell_type_idx": assignment,
                    "assignment_mode": self.mode,
                    "x_centroid": centroid[0],
                    "y_centroid": centroid[1],
                })
                cell_gex.append(gex)

        # Build AnnData
        obs_df = pd.DataFrame(cell_data)
        obs_df.index = obs_df["cell_id"]

        X = np.array(cell_gex)
        var_df = pd.DataFrame(index=gene_names)

        adata = sc.AnnData(X=X, obs=obs_df, var=var_df)

        # Add spatial coordinates
        adata.obsm["spatial"] = obs_df[["x_centroid", "y_centroid"]].values

        # Metadata
        adata.uns["assignment_mode"] = self.mode
        adata.uns["celltype_names"] = celltype_names
        adata.uns["coarse_mapping"] = COARSE_MAPPING
        adata.uns["n_spots"] = len(assignments)
        adata.uns["n_cells"] = len(cell_data)

        logger.info(
            f"Created single-cell AnnData: {adata.n_obs} cells, {adata.n_vars} genes"
        )

        return adata
