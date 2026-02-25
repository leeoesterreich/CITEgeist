"""Differentiable Sinkhorn optimal transport."""
import torch


def sinkhorn(
    cost: torch.Tensor,
    row_marginal: torch.Tensor,
    col_marginal: torch.Tensor,
    temperature: float = 0.1,
    n_iters: int = 50,
    eps: float = 1e-8,
) -> torch.Tensor:
    """
    Compute optimal transport plan via Sinkhorn iterations.

    Args:
        cost: (N, K) cost matrix (e.g., distances to prototypes)
        row_marginal: (N,) target row sums (typically uniform 1/N)
        col_marginal: (K,) target column sums (spot proportions)
        temperature: Lower = sharper assignments (default 0.1)
        n_iters: Number of Sinkhorn iterations (default 50)
        eps: Small value for numerical stability

    Returns:
        P: (N, K) transport plan (soft assignment matrix)
    """
    # Initialize kernel from cost
    K = torch.exp(-cost / temperature)

    # Normalize marginals
    row_marginal = row_marginal / row_marginal.sum()
    col_marginal = col_marginal / col_marginal.sum()

    # Sinkhorn iterations
    for _ in range(n_iters):
        # Row normalization
        row_sum = K.sum(dim=1, keepdim=True) + eps
        K = K / row_sum * row_marginal.unsqueeze(1)

        # Column normalization
        col_sum = K.sum(dim=0, keepdim=True) + eps
        K = K / col_sum * col_marginal.unsqueeze(0)

    return K
