"""
Train Cell2Location reference model on Wu 2021 breast cancer atlas.

One-time GPU job. Saves trained model + signatures for use by run_benchmark.py.
"""

import argparse
import logging
from pathlib import Path

import scanpy as sc

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref-path", required=True, help="Path to Wu 2021 reference.h5ad")
    parser.add_argument("--output-dir", required=True, help="Output directory for trained model")
    parser.add_argument("--max-epochs", type=int, default=250, help="Reference training epochs")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    import cell2location
    from cell2location.models import RegressionModel

    logger.info(f"Loading reference from {args.ref_path}")
    adata_ref = sc.read_h5ad(args.ref_path)

    sc.pp.filter_genes(adata_ref, min_cells=10)
    logger.info(f"After filtering: {adata_ref.n_obs} cells x {adata_ref.n_vars} genes")

    RegressionModel.setup_anndata(adata_ref, labels_key="cell_type")
    model = RegressionModel(adata_ref)
    model.train(max_epochs=args.max_epochs, accelerator="gpu")

    adata_ref = model.export_posterior(adata_ref, use_quantiles=True)

    model.save(str(output_dir / "model"), overwrite=True)

    # Save in format expected by existing run_benchmark.py
    # run_benchmark.py looks for: reference_trained.h5ad + cell_type_signatures.csv
    adata_ref.write(output_dir / "reference_trained.h5ad")

    vals = adata_ref.varm["q05_per_cluster_mu_fg"].values
    inf_aver = pd.DataFrame(vals, index=adata_ref.var_names, columns=adata_ref.uns["mod"]["factor_names"])
    inf_aver.to_csv(output_dir / "cell_type_signatures.csv")

    logger.info(f"Reference model saved to {output_dir}")


if __name__ == "__main__":
    main()
