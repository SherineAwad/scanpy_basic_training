import scanpy as sc
import importlib_metadata
import sys
import numpy as np
import scipy.sparse as sp

sys.modules['importlib.metadata'] = importlib_metadata


def main():

    # -----------------------------
    # INPUT FILES
    # -----------------------------
    samples = {
        "SITTB2": "SITTB2_filtered_feature_bc_matrix.h5",
        "SITTB10": "SITTB10_filtered_feature_bc_matrix.h5"
    }

    output_file = "digittip.h5ad"
    figure_prefix = "digittip"

    # -----------------------------
    # QC PARAMETERS
    # -----------------------------
    min_genes = 800
    max_genes = 8000

    min_counts = 1200
    max_counts = 50000

    max_mt = 10

    min_cells_gene = 3
    min_genes_cell = 100

    # -----------------------------
    # READ AND COMBINE SAMPLES
    # -----------------------------
    adatas = {}

    for sample_id, filename in samples.items():

        print(f"Reading {filename}...")

        adata = sc.read_10x_h5(filename)
        adata.var_names_make_unique()
        adata.obs['sample'] = sample_id

        adatas[sample_id] = adata

    combined_adata = sc.concat(
        adatas.values(),
        label='sample',
        keys=adatas.keys()
    )

    # -----------------------------
    # QC METRICS
    # -----------------------------
    combined_adata.var["mt"] = combined_adata.var_names.str.startswith("mt-")

    sc.pp.calculate_qc_metrics(
        combined_adata,
        qc_vars=["mt"],
        inplace=True,
        log1p=True
    )

    # -----------------------------
    # BEFORE QC PLOT
    # -----------------------------
    sc.pl.violin(
        combined_adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        save=f"_{figure_prefix}_preQC.png"
    )

    # -----------------------------
    # CELL FILTERING
    # -----------------------------
    combined_adata = combined_adata[
        (combined_adata.obs['n_genes_by_counts'] > min_genes) &
        (combined_adata.obs['n_genes_by_counts'] < max_genes) &
        (combined_adata.obs['total_counts'] > min_counts) &
        (combined_adata.obs['total_counts'] < max_counts) &
        (combined_adata.obs['pct_counts_mt'] < max_mt),
        :
    ]

    sc.pp.filter_cells(combined_adata, min_genes=min_genes_cell)
    sc.pp.filter_genes(combined_adata, min_cells=min_cells_gene)

    # -----------------------------
    # REMOVE GENES WITH ZERO EXPRESSION
    # -----------------------------
    if sp.issparse(combined_adata.X):
        gene_sums = np.array(
            combined_adata.X.sum(axis=0)
        ).flatten()
    else:
        gene_sums = combined_adata.X.sum(axis=0)

    combined_adata = combined_adata[:, gene_sums > 0].copy()

    # -----------------------------
    # AFTER QC PLOT
    # -----------------------------
    sc.pl.violin(
        combined_adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        save=f"_{figure_prefix}_AfterQC.png"
    )

    # -----------------------------
    # SAVE OUTPUT
    # -----------------------------
    combined_adata.write(output_file)

    print(f"Saved combined AnnData object: {output_file}")
    print(f"Final dimensions: {combined_adata.shape}")


if __name__ == "__main__":
    main()

