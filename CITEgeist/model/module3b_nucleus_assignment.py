"""Module 3b: Per-Nucleus Cell Type Assignment.

Assigns individual nuclei to cell types using spot-level proportions.

Three assignment methods are available:
1. Random (default): Random assignment within spots, respecting cell type counts.
2. Morphology-guided: Uses nucleus/cell morphology features with soft-label classifier
   and Hungarian matching. Provides minimal improvement (~1% accuracy gain).
3. VAE-Sinkhorn: Uses deep learning (VAE encoder + learned prototypes) with
   optimal transport for assignment. Requires pre-trained models.
"""
import numpy as np
import pandas as pd
from dataclasses import dataclass
from typing import Dict, List, Optional

from .morphology_features import extract_nucleus_features, extract_cell_features, largest_remainder_discretize
from .soft_label_classifier import SoftLabelClassifier
from .hungarian_assignment import assign_nuclei_to_types


@dataclass
class NucleusAssignmentResult:
    """Result of nucleus assignment."""
    assignments: Dict[int, str]  # nucleus_id -> cell_type name
    morphology_features: Optional[pd.DataFrame]  # features for all nuclei (None if random)
    classifier: Optional[SoftLabelClassifier]  # trained classifier (None if random)
    assignment_probs: Optional[pd.DataFrame]  # probability matrix (None if random)
    method: str = "random"  # "random" or "morphology"


def random_assign_nuclei_to_types(
    nucleus_ids: np.ndarray,
    counts: np.ndarray,
    rng: Optional[np.random.Generator] = None,
) -> Dict[int, int]:
    """
    Randomly assign nuclei to cell types respecting count constraints.

    Args:
        nucleus_ids: (n_nuclei,) unique identifier for each nucleus
        counts: (n_types,) integer counts per cell type
        rng: Optional random generator for reproducibility

    Returns:
        Dict mapping nucleus_id -> cell_type_index
    """
    if rng is None:
        rng = np.random.default_rng()

    n_nuclei = len(nucleus_ids)
    n_slots = int(counts.sum())

    # Build slot list from counts
    slots = []
    for type_idx, count in enumerate(counts):
        slots.extend([type_idx] * int(count))

    # Shuffle slots randomly
    rng.shuffle(slots)

    # Assign nuclei to shuffled slots
    assignments = {}
    for i, nid in enumerate(nucleus_ids):
        if i < len(slots):
            assignments[int(nid)] = slots[i]
        else:
            # More nuclei than slots - assign uniformly at random
            assignments[int(nid)] = rng.integers(0, len(counts))

    return assignments


def run_nucleus_assignment(
    mask: np.ndarray,
    nuclei_spot_map: pd.DataFrame,
    proportions: pd.DataFrame,
    nuclei_counts: pd.Series,
    cell_types: List[str],
    cell_mask: Optional[np.ndarray] = None,
    use_morphology: bool = False,
    random_seed: Optional[int] = None,
) -> NucleusAssignmentResult:
    """
    Run full nucleus assignment pipeline.

    Args:
        mask: Cellpose label mask (H, W) with nucleus labels
        nuclei_spot_map: DataFrame with 'nucleus_id' and 'spot_id' columns
        proportions: DataFrame with 'spot_id' and one column per cell type
        nuclei_counts: Series mapping spot_id -> nuclei count
        cell_types: List of cell type names (column names in proportions)
        cell_mask: Optional cell mask from watershed (same labels as nucleus mask).
                   If provided and use_morphology=True, extracts cell-level features.
        use_morphology: If True, use morphology-guided Hungarian assignment.
                        If False (default), use random assignment within spots.
                        Note: Morphology provides only ~1% accuracy improvement
                        over random in benchmarks.
        random_seed: Random seed for reproducibility (only used when use_morphology=False)

    Returns:
        NucleusAssignmentResult with assignments and metadata
    """
    # Set up proportions lookup
    prop_cols = cell_types
    spot_props = proportions.set_index('spot_id')[prop_cols]

    # Random assignment path (default)
    if not use_morphology:
        rng = np.random.default_rng(random_seed)
        all_assignments = {}

        for spot_id in nuclei_spot_map['spot_id'].unique():
            spot_nuclei = nuclei_spot_map[nuclei_spot_map['spot_id'] == spot_id]
            nucleus_ids = spot_nuclei['nucleus_id'].values

            # Get proportions and nuclei count for this spot
            spot_proportions = spot_props.loc[spot_id].values
            n_nuclei = int(nuclei_counts.get(spot_id, len(spot_nuclei)))

            # Convert proportions to counts
            counts = largest_remainder_discretize(spot_proportions, n_nuclei)

            # Random assignment
            type_assignments = random_assign_nuclei_to_types(nucleus_ids, counts, rng)

            # Convert type indices to names
            for nid, type_idx in type_assignments.items():
                all_assignments[nid] = cell_types[type_idx]

        return NucleusAssignmentResult(
            assignments=all_assignments,
            morphology_features=None,
            classifier=None,
            assignment_probs=None,
            method="random",
        )

    # Morphology-guided assignment path
    # Step 1: Extract morphology features
    if cell_mask is not None:
        # Use cell-level features (12 features)
        morph_df = extract_cell_features(mask, cell_mask)
        morph_df = morph_df.rename(columns={'cell_id': 'nucleus_id'})
        feature_cols = [
            'nucleus_area', 'nucleus_circularity', 'nucleus_eccentricity',
            'nucleus_solidity', 'nucleus_aspect_ratio',
            'cell_area', 'cell_circularity', 'cell_eccentricity',
            'cell_solidity', 'cell_aspect_ratio',
            'nc_ratio', 'cytoplasm_area',
        ]
    else:
        # Use nuclear features only (7 features)
        morph_df = extract_nucleus_features(mask)
        feature_cols = ['area', 'circularity', 'eccentricity', 'solidity',
                        'major_axis_length', 'minor_axis_length', 'aspect_ratio']

    # Merge with spot assignments
    morph_df = morph_df.merge(nuclei_spot_map, on='nucleus_id', how='inner')

    # Step 2: Build training data (soft labels from spot proportions)
    # Get soft labels for each nucleus
    y_soft = morph_df['spot_id'].map(lambda s: spot_props.loc[s].values).values
    y_soft = np.vstack(y_soft)  # (n_nuclei, n_types)

    # Feature matrix
    feature_cols = [c for c in feature_cols if c in morph_df.columns]
    X = morph_df[feature_cols].values

    # Step 3: Train classifier
    clf = SoftLabelClassifier(n_cell_types=len(cell_types))
    clf.fit(X, y_soft)

    # Step 4: Predict probabilities
    probs = clf.predict_proba(X)
    probs_df = pd.DataFrame(probs, columns=cell_types)
    probs_df['nucleus_id'] = morph_df['nucleus_id'].values
    probs_df['spot_id'] = morph_df['spot_id'].values

    # Step 5: Per-spot Hungarian assignment
    all_assignments = {}

    for spot_id in morph_df['spot_id'].unique():
        spot_mask = morph_df['spot_id'] == spot_id
        spot_nuclei = morph_df[spot_mask]
        spot_probs = probs[spot_mask]

        # Get proportions and nuclei count for this spot
        spot_proportions = spot_props.loc[spot_id].values
        n_nuclei = int(nuclei_counts.get(spot_id, len(spot_nuclei)))

        # Convert proportions to counts
        counts = largest_remainder_discretize(spot_proportions, n_nuclei)

        # Run Hungarian assignment
        nucleus_ids = spot_nuclei['nucleus_id'].values
        type_assignments = assign_nuclei_to_types(spot_probs, counts, nucleus_ids)

        # Convert type indices to names
        for nid, type_idx in type_assignments.items():
            all_assignments[nid] = cell_types[type_idx]

    return NucleusAssignmentResult(
        assignments=all_assignments,
        morphology_features=morph_df,
        classifier=clf,
        assignment_probs=probs_df,
        method="morphology",
    )


def run_nucleus_assignment_vae(
    image: np.ndarray,
    nuclei_df: pd.DataFrame,
    proportions: pd.DataFrame,
    cell_types: List[str],
    vae_checkpoint: str,
    prototype_checkpoint: str,
    device: str = "cuda",
    patch_expansion: float = 0.75,
    patch_size: int = 96,
    latent_dim: int = 128,
    projection_dim: int = 32,
    batch_size: int = 64,
) -> NucleusAssignmentResult:
    """
    Run VAE-Sinkhorn nucleus assignment.

    Uses a trained VAE encoder and prototype learning model to assign
    nuclei to cell types via optimal transport with spot-level proportion
    constraints.

    Args:
        image: (C, H, W) multi-channel image (e.g., DAPI + morphology channels)
        nuclei_df: DataFrame with columns:
            - 'nucleus_id': unique identifier for each nucleus
            - 'spot_id': spot assignment for each nucleus
            - 'bbox_x_min', 'bbox_y_min', 'bbox_x_max', 'bbox_y_max': bounding boxes
        proportions: DataFrame with 'spot_id' column and one column per cell type
        cell_types: List of cell type names (column names in proportions)
        vae_checkpoint: Path to trained VAE checkpoint (.pt file)
        prototype_checkpoint: Path to trained prototype model checkpoint (.pt file)
        device: Device to run inference on ('cuda' or 'cpu')
        patch_expansion: Fraction to expand bbox in each direction (default 0.75)
        patch_size: Patch size for VAE input (default 96)
        latent_dim: VAE latent dimension (must match checkpoint, default 128)
        projection_dim: Projection head output dimension (must match checkpoint, default 32)
        batch_size: Batch size for patch encoding (default 64)

    Returns:
        NucleusAssignmentResult with:
            - assignments: Dict mapping nucleus_id -> cell type name
            - morphology_features: None (not used in VAE method)
            - classifier: None (not used in VAE method)
            - assignment_probs: DataFrame with nucleus_id and confidence scores
            - method: "vae_sinkhorn"

    Raises:
        ImportError: If torch is not available
        FileNotFoundError: If checkpoint files don't exist
        RuntimeError: If checkpoint architecture doesn't match parameters

    Example:
        >>> result = run_nucleus_assignment_vae(
        ...     image=morphology_image,
        ...     nuclei_df=nuclei_with_bboxes,
        ...     proportions=spot_proportions,
        ...     cell_types=['Epithelial', 'Stromal', 'Immune'],
        ...     vae_checkpoint='models/vae_stage1.pt',
        ...     prototype_checkpoint='models/prototypes_stage2.pt',
        ... )
        >>> print(result.assignments[123])  # Cell type for nucleus 123
        'Epithelial'
    """
    # Lazy imports to avoid breaking existing code if torch is not available
    try:
        import torch
    except ImportError as e:
        raise ImportError(
            "PyTorch is required for VAE-Sinkhorn assignment. "
            "Install with: pip install torch"
        ) from e

    try:
        from .vae import VAE
        from .prototype_learning import PrototypeLearningModel
        from .patch_extraction import extract_patch
    except ImportError:
        # Support direct import for testing
        from vae import VAE
        from prototype_learning import PrototypeLearningModel
        from patch_extraction import extract_patch

    import os

    # Validate checkpoint files exist
    if not os.path.exists(vae_checkpoint):
        raise FileNotFoundError(f"VAE checkpoint not found: {vae_checkpoint}")
    if not os.path.exists(prototype_checkpoint):
        raise FileNotFoundError(f"Prototype checkpoint not found: {prototype_checkpoint}")

    # Validate required columns in nuclei_df
    required_cols = ['nucleus_id', 'spot_id', 'bbox_x_min', 'bbox_y_min', 'bbox_x_max', 'bbox_y_max']
    missing_cols = [c for c in required_cols if c not in nuclei_df.columns]
    if missing_cols:
        raise ValueError(f"nuclei_df missing required columns: {missing_cols}")

    n_types = len(cell_types)
    in_channels = image.shape[0]

    # Load VAE and extract encoder
    vae = VAE(in_channels=in_channels, latent_dim=latent_dim)
    vae_state = torch.load(vae_checkpoint, map_location=device)
    vae.load_state_dict(vae_state)
    vae.eval()
    encoder = vae.encoder

    # Load prototype model
    model = PrototypeLearningModel(
        encoder=encoder,
        n_types=n_types,
        latent_dim=latent_dim,
        projection_dim=projection_dim,
    )
    proto_state = torch.load(prototype_checkpoint, map_location=device)
    model.load_state_dict(proto_state, strict=False)  # strict=False since encoder is already loaded
    model.to(device)
    model.eval()

    # Set up proportions lookup by spot_id
    prop_cols = cell_types
    spot_props = proportions.set_index('spot_id')[prop_cols]

    # Process all spots
    all_assignments = {}
    all_confidences = []

    unique_spots = nuclei_df['spot_id'].unique()

    for spot_id in unique_spots:
        # Get nuclei for this spot
        spot_nuclei = nuclei_df[nuclei_df['spot_id'] == spot_id].copy()
        nucleus_ids = spot_nuclei['nucleus_id'].values

        if len(nucleus_ids) == 0:
            continue

        # Get proportions for this spot
        if spot_id not in spot_props.index:
            # Skip spots without proportions
            continue
        spot_proportions = spot_props.loc[spot_id].values.astype(np.float32)

        # Normalize proportions to sum to 1 (handle edge cases)
        prop_sum = spot_proportions.sum()
        if prop_sum <= 0:
            # No cell types expected, skip
            continue
        spot_proportions = spot_proportions / prop_sum

        # Extract patches for all nuclei in this spot
        patches_list = []
        for _, row in spot_nuclei.iterrows():
            bbox = (
                int(row['bbox_x_min']),
                int(row['bbox_y_min']),
                int(row['bbox_x_max']),
                int(row['bbox_y_max']),
            )
            patch = extract_patch(image, bbox, expansion=patch_expansion, output_size=patch_size)
            patches_list.append(patch)

        patches = np.stack(patches_list, axis=0)  # (N, C, H, W)
        patches_tensor = torch.tensor(patches, dtype=torch.float32, device=device)
        proportions_tensor = torch.tensor(spot_proportions, dtype=torch.float32, device=device)

        # Run assignment with batching for memory efficiency
        N = patches_tensor.shape[0]
        if N > batch_size:
            # For large spots, we still need to process all patches together for Sinkhorn
            # But encode in batches
            all_mu = []
            with torch.no_grad():
                for i in range(0, N, batch_size):
                    batch_patches = patches_tensor[i:i+batch_size]
                    mu, _ = model.encoder(batch_patches)
                    all_mu.append(mu)
                mu = torch.cat(all_mu, dim=0)

                # Project and compute distances
                projected = model.heads(mu)
                proto = model.prototypes()
                from .projection_heads import compute_distances
                distances = compute_distances(projected, proto)

                # Sinkhorn OT
                from .sinkhorn import sinkhorn
                row_marginal = torch.ones(N, device=device) / N
                transport_plan = sinkhorn(
                    distances,
                    row_marginal,
                    proportions_tensor,
                    temperature=0.05,  # Sharp for inference
                    n_iters=100,
                )

                # Hard assignment
                assignments = transport_plan.argmax(dim=1)
                confidence = transport_plan.max(dim=1).values
        else:
            # Use model.assign() directly for smaller spots
            with torch.no_grad():
                assignments, confidence = model.assign(patches_tensor, proportions_tensor)

        # Store results
        assignments_np = assignments.cpu().numpy()
        confidence_np = confidence.cpu().numpy()

        for i, nid in enumerate(nucleus_ids):
            type_idx = int(assignments_np[i])
            all_assignments[int(nid)] = cell_types[type_idx]
            all_confidences.append({
                'nucleus_id': int(nid),
                'spot_id': spot_id,
                'assigned_type': cell_types[type_idx],
                'confidence': float(confidence_np[i]),
            })

    # Create assignment probabilities DataFrame
    if all_confidences:
        assignment_probs = pd.DataFrame(all_confidences)
    else:
        assignment_probs = pd.DataFrame(columns=['nucleus_id', 'spot_id', 'assigned_type', 'confidence'])

    return NucleusAssignmentResult(
        assignments=all_assignments,
        morphology_features=None,
        classifier=None,
        assignment_probs=assignment_probs,
        method="vae_sinkhorn",
    )
