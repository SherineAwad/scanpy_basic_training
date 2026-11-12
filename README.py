# scRNA-seq: Quality Control and Preprocessing

### Quality Control (QC)

- **`n_genes_by_counts`**: Number of different genes detected per cell. Very low values (below 200–500) may indicate poor-quality cells; unusually high values (above 5,000–6,000) may indicate doublets.
- **`total_counts`**: Total RNA counts (UMIs) detected per cell. Very low values (below 500–1,000) may indicate poor-quality cells; unusually high values (above 20,000–30,000) may indicate doublets.
- **`pct_counts_mt`**: Percentage of counts originating from mitochondrial genes. High values (above 10–20%) may indicate stressed or damaged cells.

*These are general guidelines, not universal filtering thresholds.*

### Normalisation (`normalize_total`)
Adjusts the total RNA counts of each cell to a common target (10,000 here), making cells comparable despite differences in sequencing depth.

### Log Transformation (`log1p`)
Applies a logarithmic transformation to compress the expression value range, reducing the influence of highly expressed genes.

### Scaling (`scale`)
Centres and standardises gene expression values across cells, giving genes more comparable numerical ranges. This prevents highly variable genes from dominating PCA simply because of their larger values or variability. It does not reduce their actual biological expression.

### Dimensionality Reduction (PCA)
Simplifies data containing thousands of genes into a smaller number of features (principal components), retaining major patterns of variation between cells.

### Neighbour Graph (`neighbors`)
Identifies cells with similar expression profiles based on their PCA coordinates, establishing relationships between cells.

### UMAP (`umap`)
Projects the high-dimensional cell relationships into a two-dimensional representation, allowing visualisation of groups of transcriptionally similar cells.

### 8. Global UMAP
Colours cells according to their sample of origin, helping assess sample distribution, overlap and potential sample-specific or batch effects.

