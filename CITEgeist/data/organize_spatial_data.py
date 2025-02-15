import os
import argparse
import re
from pathlib import Path
import logging
import shutil
import gzip
from typing import Dict, List, Set, Tuple

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# Default main folder path
DEFAULT_FOLDER: Path = Path("/bgfs/alee/LO_LAB/General/Lab_Data/20250210_CITEGeistPublicData_GEO_Alex/processed_files/")

def validate_directory_structure(sample_folder: Path) -> bool:
    """Validate that a sample folder has the correct directory structure.
    
    Args:
        sample_folder: Path to the sample folder
        
    Returns:
        bool: True if structure is valid, False otherwise
    """
    required_paths = [
        sample_folder / "outs",
        sample_folder / "outs" / "filtered_feature_bc_matrix",
        sample_folder / "outs" / "spatial"
    ]
    
    return all(path.is_dir() for path in required_paths)

def fix_directory_structure(sample_folder: Path) -> None:
    """Fix the directory structure if it's not correct.
    
    Args:
        sample_folder: Path to the sample folder
    """
    if not validate_directory_structure(sample_folder):
        # Create the correct structure
        matrix_dir = sample_folder / "outs" / "filtered_feature_bc_matrix"
        spatial_dir = sample_folder / "outs" / "spatial"
        
        for dir_path in [matrix_dir, spatial_dir]:
            dir_path.mkdir(parents=True, exist_ok=True)
        
        # Move any misplaced files to their correct locations
        for file in sample_folder.rglob("*"):
            if file.is_file():
                dest_path = get_destination_path(file, matrix_dir, spatial_dir, sample_folder / "outs")
                if dest_path and dest_path != file:
                    dest_path.parent.mkdir(parents=True, exist_ok=True)
                    shutil.move(str(file), str(dest_path))

def create_directory_structure(sample_folder: Path) -> tuple[Path, Path]:
    """Create the required directory structure for a sample.
    
    Args:
        sample_folder: Path to the sample folder
        
    Returns:
        tuple[Path, Path]: Paths to the filtered_feature_bc_matrix and spatial directories
    """
    outs_dir = sample_folder / "outs"
    matrix_dir = outs_dir / "filtered_feature_bc_matrix"
    spatial_dir = outs_dir / "spatial"
    
    for dir_path in [outs_dir, matrix_dir, spatial_dir]:
        dir_path.mkdir(parents=True, exist_ok=True)
        
    return matrix_dir, spatial_dir

def get_destination_path(file_path: Path, matrix_dir: Path, spatial_dir: Path, outs_dir: Path) -> tuple[Path, str] | None:
    """Determine the destination path and new name for a file.
    
    Args:
        file_path: Original file path
        matrix_dir: Path to filtered_feature_bc_matrix directory
        spatial_dir: Path to spatial directory
        outs_dir: Path to outs directory
        
    Returns:
        tuple[Path, str] | None: (Destination directory, new filename) or None if no match
    """
    # Define file mappings with their destinations
    file_mappings = {
        # Matrix files
        'barcodes.tsv.gz': (matrix_dir, 'barcodes.tsv.gz'),
        'features.tsv.gz': (matrix_dir, 'features.tsv.gz'),
        'matrix.mtx.gz': (matrix_dir, 'matrix.mtx.gz'),
        
        # Spatial files
        'aligned_fiducials.jpg.gz': (spatial_dir, 'aligned_fiducials.jpg'),
        'aligned_tissue_image.jpg.gz': (spatial_dir, 'aligned_tissue_image.jpg'),
        'cytassist_image.tiff.gz': (spatial_dir, 'cytassist_image.tiff'),
        'detected_tissue_image.jpg.gz': (spatial_dir, 'detected_tissue_image.jpg'),
        'scalefactors_json.json.gz': (spatial_dir, 'scalefactors_json.json'),
        'spatial_enrichment.csv.gz': (spatial_dir, 'spatial_enrichment.csv'),
        'tissue_hires_image.png.gz': (spatial_dir, 'tissue_hires_image.png'),
        'tissue_lowres_image.png.gz': (spatial_dir, 'tissue_lowres_image.png'),
        'tissue_positions.csv.gz': (spatial_dir, 'tissue_positions.csv'),
        
        # Outs files
        'filtered_feature_bc_matrix.h5': (outs_dir, 'filtered_feature_bc_matrix.h5'),
        'isotype_normalization_factors.csv.gz': (outs_dir, 'isotype_normalization_factors.csv'),
        'molecule_info.h5': (outs_dir, 'molecule_info.h5')
    }
    
    # Extract the part after sample ID by finding the last underscore before known suffixes
    name_parts = file_path.name.split('_')
    for i in range(len(name_parts)-1, -1, -1):
        suffix = '_'.join(name_parts[i:])
        if suffix in file_mappings:
            return file_mappings[suffix]
    
    return None

def organize_files(folder_path: Path) -> Dict[str, int]:
    """Organize files into the proper directory structure.
    
    Args:
        folder_path: Path to the folder containing the files to organize
        
    Returns:
        Dict[str, int]: Statistics about processed files
    """
    stats = {
        'total_files': 0,
        'successful_moves': 0,
        'failed_moves': 0,
        'skipped_files': 0
    }
    
    # Define the exact list of valid sample IDs
    VALID_SAMPLE_IDS = {
        'HCC22-088-P1-S1',
        'HCC22-088-P1-S2',
        'HCC22-088-P2-S1',
        'HCC22-088-P2-S2',
        'HCC22-088-P3-S1_A',
        'HCC22-088-P3-S2',
        'HCC22-088-P4-S1',
        'HCC22-088-P4-S2',
        'HCC22-088-P4-S2_1i_rep',
        'HCC22-088-P5-S1',
        'HCC22-088-P5-S2',
        'HCC22-088-P5-S2_F_rep',
        'HCC22-088-P6-S1',
        'HCC22-088-P6-S2_D'
    }
    
    try:
        # Get all files
        files = [f for f in folder_path.iterdir() if f.is_file()]
        stats['total_files'] = len(files)
        
        # Process each sample
        for sample_id in VALID_SAMPLE_IDS:
            sample_folder = folder_path / sample_id
            outs_dir = sample_folder / "outs"
            matrix_dir, spatial_dir = create_directory_structure(sample_folder)
            logging.info(f"Processing sample: {sample_id}")
            
            # Find files that exactly match this sample
            for file in files:
                # Remove GSM prefix for matching
                without_gsm = re.sub(r'^GSM\d+_', '', file.name)
                # Check for exact match at start of filename (after GSM removal)
                if without_gsm.startswith(f"{sample_id}_"):
                    dest_info = get_destination_path(file, matrix_dir, spatial_dir, outs_dir)
                    if dest_info:
                        dest_dir, new_name = dest_info
                        try:
                            # Create parent directories if they don't exist
                            dest_dir.mkdir(parents=True, exist_ok=True)
                            dest_path = dest_dir / new_name
                            
                            if file.name.endswith('.gz') and not new_name.endswith('.gz'):
                                # Decompress while moving
                                with gzip.open(file, 'rb') as f_in:
                                    with dest_path.open('wb') as f_out:
                                        shutil.copyfileobj(f_in, f_out)
                                file.unlink()
                            else:
                                # Simple move
                                shutil.move(str(file), str(dest_path))
                            logging.info(f"Moved {file.name} to {dest_path}")
                            stats['successful_moves'] += 1
                        except Exception as e:
                            logging.error(f"Failed to move {file.name} to {dest_path}: {str(e)}")
                            stats['failed_moves'] += 1
                    else:
                        logging.warning(f"Skipping file {file.name} - unknown destination")
                        stats['skipped_files'] += 1
                
    except Exception as e:
        logging.error(f"Error while organizing files: {str(e)}")
        raise
        
    return stats

def test_pattern_matching():
    """Test function to verify sample ID extraction and file suffix matching."""
    test_files = [
        'GSM8789203_HCC22-088-P1-S1_aligned_fiducials.jpg.gz',
        'GSM8789207_HCC22-088-P3-S1_A_aligned_fiducials.jpg.gz',
        'GSM8789211_HCC22-088-P4-S2_1i_rep_aligned_fiducials.jpg.gz',
        'GSM8789214_HCC22-088-P5-S2_F_rep_tissue_positions.csv.gz',
        'GSM8789216_HCC22-088-P6-S2_D_matrix.mtx.gz'
    ]
    
    for test_file in test_files:
        # Remove GSM prefix
        without_gsm = re.sub(r'^GSM\d+_', '', test_file)
        
        # Find the sample ID
        match = re.match(r'(HCC22-\d{3}-P\d+-S\d+(?:_[A-Za-z0-9]+(?:_rep)?)?(?=_(?:aligned|barcodes|cytassist|detected|features|filtered|isotype|matrix|molecule|scalefactors|spatial|tissue)))', without_gsm)
        if match:
            sample_id = match.group(1)
            logging.info(f"\nTest file: {test_file}")
            logging.info(f"Extracted sample ID: {sample_id}")
            
            # Test destination path
            test_path = Path(test_file)
            matrix_dir = Path('test/matrix')
            spatial_dir = Path('test/spatial')
            outs_dir = Path('test/outs')
            dest_info = get_destination_path(test_path, matrix_dir, spatial_dir, outs_dir)
            if dest_info:
                dest_dir, new_name = dest_info
                logging.info(f"Destination directory: {dest_dir}")
                logging.info(f"New filename: {new_name}")
            else:
                logging.error(f"Failed to get destination path for {test_file}")
        else:
            logging.error(f"Failed to extract sample ID from {test_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Organize files into proper directory structure.')
    parser.add_argument('--folder', type=Path, default=DEFAULT_FOLDER,
                      help='Path to the main folder to process (default: %(default)s)')
    parser.add_argument('--test', action='store_true',
                      help='Run pattern matching tests')
    args = parser.parse_args()
    
    if args.test:
        test_pattern_matching()
    else:
        confirmation = input(f"Proceed with organizing files in {args.folder}? (yes/no): ").lower()
        if confirmation == 'yes':
            stats = organize_files(args.folder)
            logging.info("\nProcessing Summary:")
            logging.info(f"Total files processed: {stats['total_files']}")
            logging.info(f"Successfully moved: {stats['successful_moves']}")
            logging.info(f"Failed to move: {stats['failed_moves']}")
            logging.info(f"Skipped files: {stats['skipped_files']}")
            
            if stats['failed_moves'] == 0 and stats['skipped_files'] == 0:
                logging.info("All files were successfully processed and placed in the correct structure!")
            else:
                logging.warning("Some files were not processed successfully. Please check the logs above for details.")
        else:
            logging.info("Operation cancelled by user")
