import os
import logging
import numpy as np
from pathlib import Path

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
                if 'profiles' in complete_data and 'completed_spots' in complete_data:
                    profiles = complete_data['profiles']
                    if profiles.shape == (N, T, M):
                        return {i: profiles[i] for i in range(N)}
            except Exception as e:
                logging.error(f"Error loading complete file: {e}")
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
        checkpoints = list(self.output_dir.glob(
            f"{self.sample_name}_gene_expression_checkpoint_*.npz"
        ))
        
        if not checkpoints:
            return set(), {}
            
        checkpoint_numbers = [
            int(f.stem.split('_')[-1]) 
            for f in checkpoints
        ]
        latest_number = max(checkpoint_numbers)
        latest_checkpoint = self.output_dir / f"{self.sample_name}_gene_expression_checkpoint_{latest_number}.npz"
        
        try:
            checkpoint_data = np.load(latest_checkpoint)
            if 'profiles' in checkpoint_data and 'completed_spots' in checkpoint_data:
                profiles = checkpoint_data['profiles']
                completed_spots = set(checkpoint_data['completed_spots'].tolist())
                
                if profiles.shape == (N, T, M):
                    spotwise_profiles = {
                        i: profiles[i] 
                        for i in completed_spots 
                        if not np.any(np.isnan(profiles[i]))
                    }
                    completed_spots = set(spotwise_profiles.keys())
                    logging.info(f"Loaded {len(completed_spots)} valid profiles from checkpoint")
                    return completed_spots, spotwise_profiles
                    
        except Exception as e:
            logging.error(f"Error loading checkpoint: {e}")
            self._cleanup_corrupted_files()
            
        return set(), {}

    def save_checkpoint(self, completed_spots, spotwise_profiles, N, T, M):
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
            
            # Convert dictionary to numpy array
            max_spot = max(spotwise_profiles.keys())
            profiles_array = np.full((max_spot + 1, T, M), np.nan)
            
            for spot_idx, profile in spotwise_profiles.items():
                profiles_array[spot_idx] = profile
            
            # Save checkpoint
            np.savez_compressed(
                checkpoint_path,
                profiles=profiles_array,
                completed_spots=np.array(list(completed_spots)),
                n_completed=n_completed
            )
            
            # Cleanup old checkpoints
            self._cleanup_old_checkpoints(checkpoint_path)
            
            logging.info(f"Saved checkpoint after {n_completed} completed spots")
            
        except Exception as e:
            logging.error(f"Failed to save checkpoint: {e}")

    def save_final_results(self, spotwise_profiles, completed_spots, N, T, M):
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
        max_spot = max(spotwise_profiles.keys())
        final_profiles = np.full((max_spot + 1, T, M), np.nan)
        
        for spot_idx, profile in spotwise_profiles.items():
            final_profiles[spot_idx] = profile
            
        np.savez_compressed(
            final_path, 
            profiles=final_profiles, 
            completed_spots=np.array(list(completed_spots))
        )
        logging.info(f"Saved final results with {len(completed_spots)} completed spots")

    def _cleanup_corrupted_files(self):
        """Remove all checkpoint files if corruption is detected."""
        for file in self.output_dir.glob(f"{self.sample_name}_gene_expression*.npz"):
            try:
                file.unlink()
                logging.info(f"Deleted corrupted checkpoint: {file}")
            except Exception as e:
                logging.warning(f"Failed to delete {file}: {e}")

    def _cleanup_old_checkpoints(self, current_checkpoint):
        """Remove old checkpoints, keeping only the latest."""
        for checkpoint in self.output_dir.glob(f"{self.sample_name}_gene_expression_checkpoint_*.npz"):
            if checkpoint != current_checkpoint:
                try:
                    checkpoint.unlink()
                except Exception as e:
                    logging.warning(f"Failed to delete old checkpoint {checkpoint}: {e}") 