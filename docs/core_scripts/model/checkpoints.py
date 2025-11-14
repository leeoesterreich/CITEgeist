"""
Checkpoint management for CITEgeist optimization processes.
"""

import os
import pickle
import logging
from typing import Dict, Any, Optional, Set, Tuple


class CheckpointManager:
    """
    Manages checkpointing for long-running optimization processes.
    """
    
    def __init__(self, output_dir: str, sample_name: str):
        """
        Initialize checkpoint manager.
        
        Args:
            output_dir: Directory to save checkpoints
            sample_name: Name of the sample being processed
        """
        self.output_dir = output_dir
        self.sample_name = sample_name
        self.checkpoint_dir = os.path.join(output_dir, "checkpoints")
        os.makedirs(self.checkpoint_dir, exist_ok=True)
        
    def check_complete_run(self, N: int, T: int, M: int) -> Optional[Dict[int, Any]]:
        """
        Check if a complete run already exists.
        
        Args:
            N: Number of spots
            T: Number of cell types  
            M: Number of genes
            
        Returns:
            Complete results if found, None otherwise
        """
        final_file = os.path.join(self.checkpoint_dir, f"{self.sample_name}_final_results.pkl")
        if os.path.exists(final_file):
            try:
                with open(final_file, 'rb') as f:
                    results = pickle.load(f)
                if len(results) == N:
                    logging.info(f"Found complete run with {len(results)} spots")
                    return results
            except Exception as e:
                logging.warning(f"Failed to load final results: {e}")
        return None
        
    def load_latest_checkpoint(self, N: int, T: int, M: int) -> Tuple[Set[int], Dict[int, Any]]:
        """
        Load the latest checkpoint if available.
        
        Args:
            N: Number of spots
            T: Number of cell types
            M: Number of genes
            
        Returns:
            Tuple of (completed_spots, spotwise_profiles)
        """
        checkpoint_files = [f for f in os.listdir(self.checkpoint_dir) 
                          if f.startswith(f"{self.sample_name}_checkpoint_") and f.endswith(".pkl")]
        
        if not checkpoint_files:
            return set(), {}
            
        # Get the latest checkpoint
        latest_file = sorted(checkpoint_files)[-1]
        checkpoint_path = os.path.join(self.checkpoint_dir, latest_file)
        
        try:
            with open(checkpoint_path, 'rb') as f:
                data = pickle.load(f)
            completed_spots = data.get('completed_spots', set())
            spotwise_profiles = data.get('spotwise_profiles', {})
            logging.info(f"Loaded checkpoint with {len(completed_spots)} completed spots")
            return completed_spots, spotwise_profiles
        except Exception as e:
            logging.warning(f"Failed to load checkpoint: {e}")
            return set(), {}
            
    def save_checkpoint(self, completed_spots: Set[int], spotwise_profiles: Dict[int, Any], 
                       N: int, T: int, M: int):
        """
        Save a checkpoint.
        
        Args:
            completed_spots: Set of completed spot indices
            spotwise_profiles: Dictionary of spot profiles
            N: Number of spots
            T: Number of cell types
            M: Number of genes
        """
        checkpoint_data = {
            'completed_spots': completed_spots,
            'spotwise_profiles': spotwise_profiles,
            'N': N,
            'T': T,
            'M': M
        }
        
        checkpoint_file = os.path.join(self.checkpoint_dir, f"{self.sample_name}_checkpoint_{len(completed_spots)}.pkl")
        try:
            with open(checkpoint_file, 'wb') as f:
                pickle.dump(checkpoint_data, f)
            logging.info(f"Saved checkpoint with {len(completed_spots)} completed spots")
        except Exception as e:
            logging.error(f"Failed to save checkpoint: {e}")
            
    def save_final_results(self, spotwise_profiles: Dict[int, Any], completed_spots: Set[int],
                          N: int, T: int, M: int):
        """
        Save final results.
        
        Args:
            spotwise_profiles: Dictionary of spot profiles
            completed_spots: Set of completed spot indices
            N: Number of spots
            T: Number of cell types
            M: Number of genes
        """
        final_file = os.path.join(self.checkpoint_dir, f"{self.sample_name}_final_results.pkl")
        try:
            with open(final_file, 'wb') as f:
                pickle.dump(spotwise_profiles, f)
            logging.info(f"Saved final results with {len(spotwise_profiles)} spot profiles")
        except Exception as e:
            logging.error(f"Failed to save final results: {e}")
