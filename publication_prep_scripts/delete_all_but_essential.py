import os
from pathlib import Path
from typing import Set, List
import logging
import shutil

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# Main folder path
main_folder: Path = Path("/bgfs/alee/LO_LAB/General/Lab_Data/20250210_CITEGeistPublicData_GEO_Alex/processed_files/")

# List of acceptable files to keep
acceptable_files: Set[str] = {
    'barcodes.tsv.gz',
    'features.tsv.gz',
    'matrix.mtx.gz',
    'aligned_fiducials.jpg',
    'aligned_tissue_image.jpg',
    'cytassist_image.tiff',
    'detected_tissue_image.jpg',
    'scalefactors_json.json',
    'spatial_enrichment.csv',
    'tissue_hires_image.png',
    'tissue_lowres_image.png',
    'tissue_positions.csv',
    'filtered_feature_bc_matrix.h5',
    'isotype_normalization_factors.csv',
    'molecule_info.h5'
}

def rename_essential_file(file_path: Path, parent_folder_name: str) -> None:
    """
    Rename an essential file by adding its parent folder name as prefix.
    
    Args:
        file_path: Path to the file to rename
        parent_folder_name: Name of the parent folder to use as prefix
        
    Raises:
        OSError: If there are issues renaming the file
    """
    try:
        new_name = f"{parent_folder_name}_{file_path.name}"
        new_path = file_path.parent / new_name
        
        # Skip if the prefix is already there
        if file_path.name.startswith(f"{parent_folder_name}_"):
            logging.info(f"Skipping {file_path.name} as prefix already exists")
            return
            
        file_path.rename(new_path)
        logging.info(f"Renamed {file_path.name} to {new_name}")
    except Exception as e:
        logging.error(f"Failed to rename file {file_path}: {str(e)}")

def delete_non_essential_files(folder_path: Path) -> None:
    """
    Delete all files in the given folder except those listed in acceptable_files.
    First deletes non-essential files, then renames remaining essential files to include their parent folder name as prefix.
    Also deletes any folder named 'SPATIAL_RNA_COUNTER_CS'.
    
    Args:
        folder_path: Path to the folder to process
        
    Raises:
        OSError: If there are issues accessing files or folders
    """
    try:
        # Ensure the folder exists
        if not folder_path.exists():
            logging.error(f"Folder does not exist: {folder_path}")
            return

        # First pass: delete non-essential files and directories
        for root, dirs, files in os.walk(folder_path, topdown=False):
            root_path = Path(root)
            
            # Process files
            for file in files:
                file_path = root_path / file
                # Check if any acceptable file appears at the end of the filename
                is_acceptable = any(file.endswith(acceptable_file) for acceptable_file in acceptable_files)
                if not is_acceptable and file not in acceptable_files:
                    try:
                        file_path.unlink()
                        logging.info(f"Deleted file: {file_path}")
                    except Exception as e:
                        logging.error(f"Failed to delete file {file_path}: {str(e)}")
            
            # Process directories
            for dir_name in dirs:
                dir_path = root_path / dir_name
                try:
                    # Delete SPATIAL_RNA_COUNTER_CS folder and its contents
                    if dir_name == 'SPATIAL_RNA_COUNTER_CS':
                        shutil.rmtree(dir_path)
                        logging.info(f"Deleted SPATIAL_RNA_COUNTER_CS directory: {dir_path}")
                        continue
                        
                    # Delete raw_feature_bc_matrix folder and its contents
                    if dir_name == 'raw_feature_bc_matrix':
                        shutil.rmtree(dir_path)
                        logging.info(f"Deleted raw_feature_bc_matrix directory: {dir_path}")
                        continue
                        
                    if not any(dir_path.iterdir()):  # Check if directory is empty
                        dir_path.rmdir()
                        logging.info(f"Deleted empty directory: {dir_path}")
                except Exception as e:
                    logging.error(f"Failed to delete directory {dir_path}: {str(e)}")

        # Second pass: rename remaining essential files
        for root, _, files in os.walk(folder_path):
            root_path = Path(root)
            # Get the immediate parent folder name
            parent_folder_name = root_path.relative_to(folder_path).parts[0] if len(root_path.relative_to(folder_path).parts) > 0 else ""
            
            for file in files:
                if file in acceptable_files:
                    file_path = root_path / file
                    if parent_folder_name:  # Only rename if there's a parent folder
                        rename_essential_file(file_path, parent_folder_name)

    except Exception as e:
        logging.error(f"An error occurred while processing the folder: {str(e)}")
        raise

def remove_prefix_from_essential_files(folder_path: Path) -> None:
    """
    Remove the folder prefix from all essential files that were previously renamed.
    
    Args:
        folder_path: Path to the folder to process
        
    Raises:
        OSError: If there are issues accessing or renaming files
    """
    try:
        # Ensure the folder exists
        if not folder_path.exists():
            logging.error(f"Folder does not exist: {folder_path}")
            return

        for root, _, files in os.walk(folder_path):
            root_path = Path(root)
            
            for file in files:
                file_path = root_path / file
                # Check if this is a renamed essential file
                original_name = file.split('_', 1)[-1] if '_' in file else file
                if original_name in acceptable_files and '_' in file:
                    try:
                        new_path = file_path.parent / original_name
                        file_path.rename(new_path)
                        logging.info(f"Removed prefix from {file_path.name} -> {original_name}")
                    except Exception as e:
                        logging.error(f"Failed to remove prefix from {file_path}: {str(e)}")

    except Exception as e:
        logging.error(f"An error occurred while removing prefixes: {str(e)}")
        raise

if __name__ == "__main__":
    print(f"Choose operation:")
    print("1. Delete non-essential files and add prefixes")
    print("2. Remove prefixes from essential files")
    
    choice = input("Enter choice (1 or 2): ").strip()
    
    if choice == "1":
        print(f"This will delete all files except the following from {main_folder}:")
        for file in sorted(acceptable_files):
            print(f"  - {file}")
        
        confirmation = input("Do you want to proceed? (yes/no): ").lower()
        if confirmation == 'yes':
            delete_non_essential_files(main_folder)
            logging.info("Cleanup completed successfully")
        else:
            logging.info("Operation cancelled by user")
            
    elif choice == "2":
        confirmation = input("Remove prefixes from essential files? (yes/no): ").lower()
        if confirmation == 'yes':
            remove_prefix_from_essential_files(main_folder)
            logging.info("Prefix removal completed successfully")
        else:
            logging.info("Operation cancelled by user")
    
    else:
        logging.error("Invalid choice. Please run again and select 1 or 2.")
