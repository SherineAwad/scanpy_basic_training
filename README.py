# scRNA-seq: Quality Control and Preprocessing

## Quality Control (QC)

```python
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"], jitter=0.4, multi_panel=True)
```

- `n_genes_by_counts`: Number of different genes detected per cell.
  - Below 200–500: Potentially poor-quality cells.
  - Above 5,000–6,000: Investigate possible doublets.
- `total_counts`: Total RNA counts (UMIs) detected per cell.
  - Below 500–1,000: Potentially poor-quality cells.
  - Above 20,000–30,000: Investigate possible doublets.
- `pct_counts_mt`: Percentage of counts originating from mitochondrial genes.
  - Below 5%: Generally good.
  - 5–10%: Generally acceptable.
  - 10–20%: Investigate possible cell stress or damage.
  - Above 20–30%: Potentially poor-quality cells.

*Note: These are general guidelines, not universal filtering thresholds.*

```python
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs["n_genes_by_counts"] < 6000].copy()
adata = adata[adata.obs["total_counts"] < 30000].copy()
adata = adata[adata.obs["pct_counts_mt"] < 20].copy()
```

## Normalisation (`normalize_total`)

```python
sc.pp.normalize_total(adata, target_sum=1e4)
```

Adjusts the total RNA counts of each cell to a common target (10,000 here), making cells comparable despite differences in sequencing depth.

## Log Transformation (`log1p`)

```python
sc.pp.log1p(adata)
```

Applies a logarithmic transformation to compress the expression value range, reducing the influence of highly expressed genes.

## Scaling (`scale`)

```python
sc.pp.scale(adata, max_value=10)
```

Centres and standardises gene expression values across cells, giving genes more comparable numerical ranges. This prevents highly variable genes from dominating PCA simply because of their larger values or variability. It does not reduce their actual biological expression.

## Dimensionality Reduction (PCA)

```python
sc.tl.pca(adata, n_comps=50, svd_solver="arpack")
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)
```

Simplifies data containing thousands of genes into a smaller number of features (principal components), retaining major patterns of variation between cells.

## Neighbour Graph (`neighbors`)

```python
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
```

Identifies cells with similar expression profiles based on their PCA coordinates, establishing relationships between cells.

## UMAP (`umap`)

```python
sc.tl.umap(adata)
sc.pl.umap(adata, color=["sample"], frameon=False)
```

Projects high-dimensional cell relationships into a two-dimensional representation, allowing visualisation of groups of transcriptionally similar cells.

## Global UMAP

```python
sc.pl.umap(adata, color=["sample", "batch"], frameon=False, ncols=2)
```

Colours cells according to their sample of origin, helping assess sample distribution, overlap and potential sample-specific or batch effects.

QC and Filtering → Normalisation → Log Transformation → Scaling → PCA → Neighbour Graph → UMAP → Visualisation
