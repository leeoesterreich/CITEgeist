
import scanpy as sc
import anndata as ad


def filter_spots_by_celltype_proportion(adata, celltype, proportion_threshold):
    """
    Filter spots by celltype proportion based on if their .obs[celltype] is greater than the proportion_threshold.
    """
    return adata[adata.obs[celltype] > proportion_threshold].copy()


def expand_prop_gex_adata(prop_gex_adata, celltype_cutoff=0.2, celltype_profile_dict=None):
    """
    Expands celltype level spots to have multiple spots for each celltype, but retain spatial information
    for COMMMOT analysis.
    """
    
    if celltype_profile_dict is None:
        raise ValueError("celltype_profile_dict is required for this function.")
    celltypes = celltype_profile_dict.keys()

    all_adatas = []
    
    for celltype in celltypes:
        filtered_adata = filter_spots_by_celltype_proportion(prop_gex_adata, celltype, celltype_cutoff)
        celltype_label = celltype.replace(" ", "_")
        layer_name = f"{celltype_label}_genes_pass1"
        filtered_adata.X = filtered_adata.layers[layer_name]
        filtered_adata.obs['celltype'] = celltype
        all_adatas.append(filtered_adata)
    # Extract the corresponding layer for the filtered_adata, including the spatial information
    combined_adata = ad.concat(all_adatas, join="outer", keys=celltypes, index_unique = "_", merge="same", uns_merge = "first")
    return combined_adata

def get_celltype_expression_data(adata, celltype):
    """
    Get expression data for a specific celltype.
    """
    return adata[adata.obs['celltype'] == celltype].copy()


