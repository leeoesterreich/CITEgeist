"""
Checkpoint management for saving and resuming optimization state.
"""

import logging
from pathlib import Path

import numpy as np


class CheckpointManager:
    """Manages loading and saving of optimization checkpoints."""

    def __init__(self, output_dir, sample_name):
        """
        Initialize checkpoint manager.

        Args:
            output_dir (str): Directory for checkpoint files
            sample_name (str): Unique identifier for this sample
        """
        self.output_dir = Path(output_dir)
        self.sample_name = sample_name
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def check_complete_run(self, N, T, M):
        """
        Check if a complete run exists.

        Args:
            N (int): Number of spots
            T (int): Number of cell types
            M (int): Number of genes/features

        Returns:
            dict or None: Dictionary of spot profiles if complete run exists, else None
        """
        complete_file = self.output_dir / f"{self.sample_name}_gene_expression_complete.npz"

        if complete_file.exists():
            try:
                complete_data = np.load(complete_file)
                if "profiles" in complete_data and "completed_spots" in complete_data:
                    profiles = complete_data["profiles"]
                    if profiles.shape == (N, T, M):
                        return {i: profiles[i] for i in range(N)}
            except (OSError, ValueError) as e:

                logging.error("Error loading complete file: %s", e)
                self._cleanup_corrupted_files()
        return None

    def load_latest_checkpoint(self, N, T, M):
        """
        Load the latest valid checkpoint.

        Args:
            N (int): Number of spots
            T (int): Number of cell types
            M (int): Number of genes/features

        Returns:
            tuple: (completed_spots set, spotwise_profiles dict)
        """
        checkpoints = list(self.output_dir.glob(f"{self.sample_name}_gene_expression_checkpoint_*.npz"))

        if not checkpoints:
            return set(), {}

        checkpoint_numbers = [int(f.stem.split("_")[-1]) for f in checkpoints]
        latest_number = max(checkpoint_numbers)
        latest_checkpoint = self.output_dir / f"{self.sample_name}_gene_expression_checkpoint_{latest_number}.npz"

        try:
            checkpoint_data = np.load(latest_checkpoint)
            if "profiles" in checkpoint_data and "completed_spots" in checkpoint_data:
                profiles = checkpoint_data["profiles"]
                completed_spots = set(checkpoint_data["completed_spots"].tolist())

                if profiles.shape == (N, T, M):
                    spotwise_profiles = {i: profiles[i] for i in completed_spots if not np.any(np.isnan(profiles[i]))}
                    completed_spots = set(spotwise_profiles.keys())
                    logging.info("Loaded %s valid profiles from checkpoint", len(completed_spots))
                    return completed_spots, spotwise_profiles

        except (OSError, ValueError) as e:

            logging.error("Error loading checkpoint: %s", e)
            self._cleanup_corrupted_files()

        return set(), {}

    def save_checkpoint(
        self, completed_spots, spotwise_profiles, N, T, M
    ):  # pylint: disable=too-many-positional-arguments
        """
        Save current progress as checkpoint.

        Args:
            completed_spots (set): Set of completed spot indices
            spotwise_profiles (dict): Dictionary of spot profiles
            N (int): Number of spots
            T (int): Number of cell types
            M (int): Number of genes/features
        """
        try:
            n_completed = len(completed_spots)
            checkpoint_path = self.output_dir / f"{self.sample_name}_gene_expression_checkpoint_{n_completed}.npz"

            # Always allocate the full (N, T, M) array so loaders can verify shape.
            profiles_array = np.full((N, T, M), np.nan)

            for spot_idx, profile in spotwise_profiles.items():
                profiles_array[spot_idx] = profile

            # Save checkpoint
            np.savez_compressed(
                checkpoint_path,
                profiles=profiles_array,
                completed_spots=np.array(list(completed_spots)),
                n_completed=n_completed,
            )

            # Cleanup old checkpoints
            self._cleanup_old_checkpoints(checkpoint_path)

            logging.info("Saved checkpoint after %s completed spots", n_completed)

        except (OSError, ValueError) as e:

            logging.error("Failed to save checkpoint: %s", e)

    def save_final_results(
        self, spotwise_profiles, completed_spots, N, T, M
    ):  # pylint: disable=too-many-positional-arguments
        """
        Save final results.

        Args:
            spotwise_profiles (dict): Dictionary of spot profiles
            completed_spots (set): Set of completed spot indices
            N (int): Number of spots
            T (int): Number of cell types
            M (int): Number of genes/features
        """
        final_path = self.output_dir / f"{self.sample_name}_gene_expression_complete.npz"
        final_profiles = np.full((N, T, M), np.nan)

        for spot_idx, profile in spotwise_profiles.items():
            final_profiles[spot_idx] = profile

        np.savez_compressed(final_path, profiles=final_profiles, completed_spots=np.array(list(completed_spots)))
        logging.info("Saved final results with %s completed spots", len(completed_spots))

    def _cleanup_corrupted_files(self):
        """Remove all checkpoint files if corruption is detected."""
        for file in self.output_dir.glob(f"{self.sample_name}_gene_expression*.npz"):
            try:
                file.unlink()
                logging.info("Deleted corrupted checkpoint: %s", file)
            except (OSError, ValueError) as e:

                logging.warning("Failed to delete %s: %s", file, e)

    def _cleanup_old_checkpoints(self, current_checkpoint):
        """Remove old checkpoints, keeping only the latest."""
        for checkpoint in self.output_dir.glob(f"{self.sample_name}_gene_expression_checkpoint_*.npz"):
            if checkpoint != current_checkpoint:
                try:
                    checkpoint.unlink()
                except (OSError, ValueError) as e:

                    logging.warning("Failed to delete old checkpoint %s: %s", checkpoint, e)
