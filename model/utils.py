import os
import gc
import logging
import pandas as pd
import numpy as np
import scanpy as sc

def validate_cell_profile_dict(cell_profile_dict):
    """
    Validate the structure of a cell profile dictionary.
    """
    if not isinstance(cell_profile_dict, dict):
        return False
    return all(isinstance(k, str) and isinstance(v, dict) for k, v in cell_profile_dict.items())

def save_results_to_output(results, filepath):
    """
    Save results as a CSV file.
    """
    df = pd.DataFrame(results)
    df.to_csv(filepath)

def cleanup_memory():
    """
    Force garbage collection to free memory.
    """
    gc.collect()

def setup_logging(output_folder, sample_name):
    """
    Set up dynamic logging.
    """
    log_file = os.path.join(output_folder, f"{sample_name}_CITEgeist.log")
    logging.basicConfig(
        filename=log_file,
        filemode='w',
        format='%(asctime)s - %(levelname)s - %(message)s',
        level=logging.INFO
    )
    logging.info('Logging initialized.')

### 📏 **Spatial Neighbor Functions**



def find_fixed_radius_neighbors(spot_index, adata, radius=50):
    """
    Find neighbors within a fixed radius of a given spot.

    Args:
        spot_index (int): Index of the central spot in the AnnData object.
        adata (AnnData): Spatial transcriptomics dataset with obsm['spatial'] containing spot coordinates.
        radius (float): Fixed radius to identify neighboring spots.

    Returns:
        tuple: (central_spot_name, list of neighbor spot names)
    """
    coordinates = adata.obsm['spatial']
    central_coord = coordinates[spot_index]
    
    # Identify all spots within the given radius, excluding the central spot
    neighbors = [idx for idx, coord in enumerate(coordinates)
                 if idx != spot_index and np.linalg.norm(coord - central_coord) <= radius]
    
    # Convert indices to spot names
    neighbor_names = adata.obs_names[neighbors].tolist()
    central_spot_name = adata.obs_names[spot_index]
    
    return central_spot_name, neighbor_names, spot_index, neighbors


def get_neighbors_with_fixed_radius(spot_index, adata, radius=50, include_center=True): # In Visium, 2 rings gives you the adjacent 6 units
    """
    Get indices of neighboring spots based on a fixed radius around the central spot.

    Parameters:
    - spot_index (int): The index of the central spot.
    - adata (AnnData): The AnnData object with spatial coordinates in obsm['spatial'].
    - radius (float): Fixed radius for finding neighbors.
    - include_center (bool): Whether to include the central spot in the neighbor list.

    Returns:
    - List of indices representing the central spot and its neighbors.
    """
    # Find neighbors within the given radius
    central_spot_names, neighbor_spots_names, spot_index, neighbors = find_fixed_radius_neighbors(spot_index, adata, radius)
    
    # Optionally include the central spot itself
    if include_center:
        neighbors = [spot_index] + neighbors
    
    logging.debug(f"Total neighbors for spot {spot_index} within radius {radius}: {neighbors}")
    return neighbors

def plot_neighbors_with_fixed_radius(adata, radius=50, num_spots=5):
    """
    Plot neighbors for multiple random central spots using `sc.pl.spatial`.

    Args:
        adata (AnnData): Spatial transcriptomics dataset with obsm['spatial'].
        radius (float): Fixed radius to identify neighboring spots.
        num_spots (int): Number of random spots to visualize.

    Returns:
        None: Displays a series of spatial plots showing neighbors.
    """
    import random
    
    # Select random spots
    random_spots = random.sample(range(adata.shape[0]), min(num_spots, adata.shape[0]))
    
    for spot_index in random_spots:
        # Find neighbors within the given radius
        central_spot_names, neighbor_spots_names, spot_index, neighbors = find_fixed_radius_neighbors(spot_index, adata, radius)

        # Create a temporary column to highlight spots
        adata.obs['highlight'] = 'Other spots'
        adata.obs.loc[neighbor_spots_names, 'highlight'] = 'Neighbor'
        adata.obs.loc[central_spot_names, 'highlight'] = 'Central spot'
        
        # Plot using `sc.pl.spatial`
        sc.pl.spatial(
            adata,
            color='highlight',
            title=f"Neighbors within {radius} units for Spot {central_spot_names}",
            spot_size=50,
            frameon=False
        )
        
        # Clean up temporary column after each plot
        adata.obs.drop(columns=['highlight'], inplace=True)
        
def assert_neighborhood_size(adata, cell_profile_dict, radius=50, num_spots=5):
    """

    """
    import random
    
    # Select random spots
    random_spots = random.sample(range(adata.shape[0]), min(num_spots, adata.shape[0]))
    
    neighborhood_sizes = []
    
    for spot_index in random_spots:
        # Find neighbors within the given radius
        central_spot_names, neighbor_spots_names, spot_index, neighbors = find_fixed_radius_neighbors(spot_index, adata, radius)

    central_spot_names = list(central_spot_names) if not isinstance(central_spot_names, list) else central_spot_names
    neighbor_spots_names = list(neighbor_spots_names) if not isinstance(neighbor_spots_names, list) else neighbor_spots_names

    neighborhood_size = len(central_spot_names + neighbor_spots_names)
    
    
    assert all(x <= len(cell_profile_dict) for x in neighborhood_sizes), f"Some neighborhood values in the list are less than {len(cell_profile_dict)} celltypes being deconvoluted"
