



def main():

    cell_profiles = {
        "Cancer Cells": {
            "Major": ["EPCAM-1"],
            "Minor": ["SDC1-1",  "KRT5-1"]  # CD138 - possible cancer stem cell marker
        },
        "Macrophages": {
            "Major": ["CD68-1" ],  # General macrophage and M2-polarized macrophages
            "Minor": ["CD14-1"]  # Monocyte/macrophage lineage marker
        },
        "CD4 T Cells": {
            "Major": ["CD3E-1", "CD4-1"],  # General, Helper, and Cytotoxic T cells
        },
        "CD8 T Cells": {
            "Major": ["CD3E-1", "CD8A-1"],  # General, Helper, and Cytotoxic T cells
        },
        "B Cells": {
            "Major": ["MS4A1-1", "CD19-1"],  # General B cell markers and developmental marker
        },
        "Endothelial Cells": {
            "Major": ["PECAM1-1"],  # CD31 - endothelial cell marker
        },
        "Fibroblasts": {
            "Major": ["ACTA2-1"],  # α-SMA - myofibroblast marker, indicates activated stroma
        }
    }

    import scanpy as sc
    import importlib
    import model
    importlib.reload(model)
    from model import CitegeistModel  # Now directly import from the module
    import argparse
    # Reload the model module if you're actively developing

    parser = argparse.ArgumentParser(description='Run CITEgeist on a single sample.')
    parser.add_argument('--path', type=str, required=True, help='Path to the sample')
    parser.add_argument('--radius', type=float, required=True, help='Radius for neighbor detection')
    parser.add_argument('--min_counts', type=float, default=100, required=True, 
                       help='Minimum counts for filtering')
    parser.add_argument('--output_folder', type=str, default='output', help='Output folder')

    args = parser.parse_args()


    adata = sc.read_visium(args.path, count_file='filtered_feature_bc_matrix.h5', load_images=True, gex_only=False)
    sample_name = args.path.split('/')[-2]
    print(sample_name)

        # Initialize the model
    model = CitegeistModel(sample_name=sample_name, adata=adata, output_folder='output')

    # Load cell profile dictionary
    model.load_cell_profile_dict(cell_profiles)

    # Preprocess and run models
    # Split into gene expression and antibody capture datasets
    model.split_adata()

    model.filter_gex(nonzero_percentage=0.01, mean_expression_threshold=1.1, min_counts=100)

    model.copy_gex_to_protein_adata()

    # Preprocess datasets
    model.preprocess_gex()
    model.preprocess_antibody()
    print(model)

    # Register Gurobi license
    model.register_gurobi("/ihome/crc/install/gurobi/gurobi1102/linux64/lic/gurobi.lic")

    global_cell_type_proportions_df, finetuned_cell_type_proportions_df = model.run_cell_proportion_model(radius=args.radius)

    # Plot cell proportions (Append cell proportions) 
    model.append_proportions_to_adata(key='global')

    # Run gene expression model (if needed)
    # model.run_gex_model()
    prop_adata = model.get_adata()

    # Plot cell proportions (Append cell proportions) 
    model.append_proportions_to_adata(key='finetuned')

    prop_adata = model.get_adata()

    sc.pl.spatial(prop_adata, color = "Cancer Cells")


    pass1_results = model.run_cell_expression_pass1(
                radius=args.radius, 
                max_workers=None, 
                checkpoint_interval=100, 
                output_dir="checkpoints", 
                rerun=True
            )

    # Plot cell proportions (Append cell proportions) 
    model.append_gex_to_adata(pass_number=1)


    prop_gex_adata = model.get_adata()

    print(prop_gex_adata)




if __name__ == "__main__":
    main()
