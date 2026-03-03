"""Training pipeline for proportion-guided MIL.

Trains MIL head on pre-extracted ViT embeddings with spot-level
proportion supervision.
"""
import sys
from pathlib import Path
import importlib.util
import numpy as np
import torch
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
from typing import List, Tuple, Dict, Optional
import logging


def _pearsonr(x: np.ndarray, y: np.ndarray) -> tuple:
    """Compute Pearson correlation coefficient using numpy.

    Args:
        x: First array
        y: Second array

    Returns:
        (r, p) tuple where r is correlation and p is two-sided p-value
    """
    n = len(x)
    if n < 2:
        return np.nan, np.nan

    # Center the data
    x_centered = x - np.mean(x)
    y_centered = y - np.mean(y)

    # Correlation coefficient
    numerator = np.sum(x_centered * y_centered)
    denominator = np.sqrt(np.sum(x_centered**2) * np.sum(y_centered**2))

    if denominator == 0:
        return np.nan, np.nan

    r = numerator / denominator

    # P-value calculation (two-tailed t-test)
    if abs(r) == 1.0:
        p = 0.0
    else:
        t_stat = r * np.sqrt((n - 2) / (1 - r**2))
        # Use numerical approximation for t-distribution CDF
        # For simplicity, return approximate p-value
        p = 2 * (1 - _t_cdf(abs(t_stat), n - 2))

    return float(r), float(p)


def _t_cdf(t: float, df: int) -> float:
    """Approximate t-distribution CDF using normal approximation for large df.

    Args:
        t: t-statistic
        df: Degrees of freedom

    Returns:
        Approximate CDF value
    """
    # For df > 30, t-distribution is approximately normal
    if df > 30:
        # Normal CDF approximation
        return 0.5 * (1 + np.tanh(t * 0.797884560802865))  # approx erf(t/sqrt(2))
    else:
        # Use a better approximation for smaller df
        # Hill approximation
        x = df / (df + t * t)
        return 0.5 + 0.5 * np.sign(t) * (1 - np.power(x, df / 2))

# Direct import to avoid __init__.py issues
REPO_ROOT = Path(__file__).parent.parent.parent.parent
_spec = importlib.util.spec_from_file_location(
    'proportion_mil',
    REPO_ROOT / 'CITEgeist/model/proportion_mil.py'
)
_mil_module = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mil_module)

ProportionGuidedMIL = _mil_module.ProportionGuidedMIL
proportion_loss = _mil_module.proportion_loss
entropy_regularization = _mil_module.entropy_regularization

logger = logging.getLogger(__name__)


class SpotDataset(Dataset):
    """Dataset of spot embeddings and proportions.

    Expects directory structure:
        data_dir/
            spot_0/
                embeddings.npy  # (N, embed_dim)
                proportions.npy  # (K,)
            spot_1/
            ...

    Attributes:
        data_dir: Directory containing spot subdirectories
        n_cell_types: Number of cell types
        spots: List of valid spot directories
    """

    def __init__(
        self,
        data_dir: Path,
        n_cell_types: int,
        min_nuclei: int = 5,
    ):
        """Initialize dataset.

        Args:
            data_dir: Directory containing spot subdirectories
            n_cell_types: Number of cell types
            min_nuclei: Minimum nuclei per spot to include
        """
        self.data_dir = Path(data_dir)
        self.n_cell_types = n_cell_types

        # Find valid spots
        self.spots = []
        for spot_dir in sorted(self.data_dir.glob("spot_*")):
            emb_path = spot_dir / "embeddings.npy"
            prop_path = spot_dir / "proportions.npy"

            if emb_path.exists() and prop_path.exists():
                emb = np.load(emb_path)
                if len(emb) >= min_nuclei:
                    self.spots.append(spot_dir)

        logger.info(f"Found {len(self.spots)} valid spots in {data_dir}")

    def __len__(self) -> int:
        """Return number of spots in dataset."""
        return len(self.spots)

    def __getitem__(self, idx: int) -> Tuple[torch.Tensor, torch.Tensor]:
        """Get embeddings and proportions for a spot.

        Args:
            idx: Spot index

        Returns:
            embeddings: (N, embed_dim) tensor of nucleus embeddings
            proportions: (K,) tensor of cell type proportions
        """
        spot_dir = self.spots[idx]

        embeddings = np.load(spot_dir / "embeddings.npy")
        proportions = np.load(spot_dir / "proportions.npy")

        return (
            torch.from_numpy(embeddings).float(),
            torch.from_numpy(proportions).float(),
        )


def collate_spots(
    batch: List[Tuple[torch.Tensor, torch.Tensor]]
) -> List[Tuple[torch.Tensor, torch.Tensor]]:
    """Custom collate function - each spot is processed independently.

    Since spots have varying numbers of nuclei, we cannot stack them.
    Return as list for iteration.

    Args:
        batch: List of (embeddings, proportions) tuples

    Returns:
        Same list (no stacking)
    """
    return batch


def train_epoch(
    model: 'ProportionGuidedMIL',
    data: List[Tuple[torch.Tensor, torch.Tensor]],
    optimizer: torch.optim.Optimizer,
    device: str = 'cpu',
    entropy_weight: float = 0.01,
) -> float:
    """Train for one epoch.

    Args:
        model: MIL model
        data: List of (embeddings, proportions) tuples
        optimizer: Optimizer
        device: Device to use ('cpu' or 'cuda')
        entropy_weight: Weight for entropy regularization

    Returns:
        Average loss over epoch
    """
    model.train()
    total_loss = 0.0

    for embeddings, gt_proportions in data:
        embeddings = embeddings.to(device)
        gt_proportions = gt_proportions.to(device)

        # Forward pass
        pred_proportions, attention = model(embeddings)

        # Compute loss
        loss = proportion_loss(pred_proportions, gt_proportions)
        loss += entropy_weight * entropy_regularization(attention)

        # Backward pass
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        total_loss += loss.item()

    return total_loss / len(data)


@torch.no_grad()
def evaluate(
    model: 'ProportionGuidedMIL',
    data: List[Tuple[torch.Tensor, torch.Tensor]],
    device: str = 'cpu',
) -> Dict[str, float]:
    """Evaluate model on data.

    Args:
        model: MIL model
        data: List of (embeddings, proportions) tuples
        device: Device to use

    Returns:
        Dictionary containing:
            - loss: Average loss
            - pearson_r: Pearson correlation coefficient
            - pearson_p: P-value for correlation
    """
    model.eval()

    total_loss = 0.0
    all_pred = []
    all_gt = []

    for embeddings, gt_proportions in data:
        embeddings = embeddings.to(device)
        gt_proportions = gt_proportions.to(device)

        pred_proportions, _ = model(embeddings)

        loss = proportion_loss(pred_proportions, gt_proportions)
        total_loss += loss.item()

        all_pred.append(pred_proportions.cpu().numpy())
        all_gt.append(gt_proportions.cpu().numpy())

    # Flatten all predictions and ground truths
    all_pred = np.concatenate([p.flatten() for p in all_pred])
    all_gt = np.concatenate([g.flatten() for g in all_gt])

    # Compute Pearson correlation
    r, p = _pearsonr(all_pred, all_gt)

    return {
        'loss': total_loss / len(data),
        'pearson_r': r,
        'pearson_p': p,
    }


def train(
    model: 'ProportionGuidedMIL',
    train_data: List[Tuple[torch.Tensor, torch.Tensor]],
    val_data: List[Tuple[torch.Tensor, torch.Tensor]],
    n_epochs: int = 100,
    lr: float = 1e-4,
    device: str = 'cpu',
    save_path: Optional[Path] = None,
    entropy_weight: float = 0.01,
    log_interval: int = 10,
) -> Dict[str, List[float]]:
    """Full training loop with validation and checkpointing.

    Args:
        model: MIL model to train
        train_data: Training data (list of spot tuples)
        val_data: Validation data (list of spot tuples)
        n_epochs: Number of training epochs
        lr: Learning rate
        device: Device to train on
        save_path: Path to save best model checkpoint
        entropy_weight: Weight for entropy regularization
        log_interval: Log metrics every N epochs

    Returns:
        Training history dictionary containing:
            - train_loss: List of training losses per epoch
            - val_loss: List of validation losses per epoch
            - val_r: List of validation Pearson r per epoch
    """
    model = model.to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-4)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, n_epochs)

    history = {
        'train_loss': [],
        'val_loss': [],
        'val_r': [],
    }
    best_r = -1

    for epoch in range(n_epochs):
        # Train
        train_loss = train_epoch(
            model, train_data, optimizer, device, entropy_weight
        )

        # Evaluate
        val_metrics = evaluate(model, val_data, device)

        # Update learning rate
        scheduler.step()

        # Record history
        history['train_loss'].append(train_loss)
        history['val_loss'].append(val_metrics['loss'])
        history['val_r'].append(val_metrics['pearson_r'])

        # Save best model
        if val_metrics['pearson_r'] > best_r:
            best_r = val_metrics['pearson_r']
            if save_path:
                torch.save(model.state_dict(), save_path)
                logger.debug(f"Saved best model with r={best_r:.4f}")

        # Log progress
        if epoch % log_interval == 0:
            logger.info(
                f"Epoch {epoch}: train_loss={train_loss:.4f}, "
                f"val_loss={val_metrics['loss']:.4f}, "
                f"val_r={val_metrics['pearson_r']:.4f}"
            )

    logger.info(f"Training complete. Best validation r={best_r:.4f}")

    return history


def create_data_loaders(
    dataset: SpotDataset,
    train_frac: float = 0.8,
    batch_size: int = 32,
    seed: int = 42,
) -> Tuple[DataLoader, DataLoader]:
    """Create train/val data loaders from dataset.

    Args:
        dataset: SpotDataset instance
        train_frac: Fraction of data for training
        batch_size: Batch size (number of spots per batch)
        seed: Random seed for split

    Returns:
        train_loader, val_loader
    """
    # Split indices
    n_samples = len(dataset)
    indices = np.arange(n_samples)

    np.random.seed(seed)
    np.random.shuffle(indices)

    n_train = int(n_samples * train_frac)
    train_indices = indices[:n_train]
    val_indices = indices[n_train:]

    # Create subset datasets
    train_subset = torch.utils.data.Subset(dataset, train_indices)
    val_subset = torch.utils.data.Subset(dataset, val_indices)

    # Create data loaders
    train_loader = DataLoader(
        train_subset,
        batch_size=batch_size,
        shuffle=True,
        collate_fn=collate_spots,
    )
    val_loader = DataLoader(
        val_subset,
        batch_size=batch_size,
        shuffle=False,
        collate_fn=collate_spots,
    )

    return train_loader, val_loader


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Train MIL model')
    parser.add_argument(
        '--data_dir',
        type=Path,
        required=True,
        help='Directory containing spot data',
    )
    parser.add_argument(
        '--n_cell_types',
        type=int,
        required=True,
        help='Number of cell types',
    )
    parser.add_argument(
        '--output_dir',
        type=Path,
        required=True,
        help='Output directory for model and logs',
    )
    parser.add_argument(
        '--n_epochs',
        type=int,
        default=100,
        help='Number of training epochs',
    )
    parser.add_argument(
        '--lr',
        type=float,
        default=1e-4,
        help='Learning rate',
    )
    parser.add_argument(
        '--hidden_dim',
        type=int,
        default=256,
        help='Hidden dimension in MIL model',
    )
    parser.add_argument(
        '--device',
        type=str,
        default='cuda' if torch.cuda.is_available() else 'cpu',
        help='Device to train on',
    )

    args = parser.parse_args()

    # Setup logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
    )

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Load dataset
    logger.info(f"Loading data from {args.data_dir}")
    dataset = SpotDataset(
        args.data_dir,
        n_cell_types=args.n_cell_types,
    )

    # Create data loaders
    train_loader, val_loader = create_data_loaders(dataset)

    # Flatten loaders for train function
    train_data = [item for batch in train_loader for item in batch]
    val_data = [item for batch in val_loader for item in batch]

    logger.info(f"Train: {len(train_data)} spots, Val: {len(val_data)} spots")

    # Create model
    model = ProportionGuidedMIL(
        input_dim=768,
        n_cell_types=args.n_cell_types,
        hidden_dim=args.hidden_dim,
    )

    # Train
    save_path = args.output_dir / 'best_model.pt'
    history = train(
        model,
        train_data,
        val_data,
        n_epochs=args.n_epochs,
        lr=args.lr,
        device=args.device,
        save_path=save_path,
    )

    # Save history
    import json
    history_path = args.output_dir / 'training_history.json'
    with open(history_path, 'w') as f:
        json.dump(history, f, indent=2)

    logger.info(f"Training complete. Model saved to {save_path}")
