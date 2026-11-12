# scRNA-seq: Quality Control and Preprocessing

## 1. Quality Control (QC)

- **n_genes_by_counts:** Number of different genes detected per cell.
  - Below 200–500: Potentially poor-quality cells.
  - Above 5,000–6,000: Investigate possible doublets.

- **total_counts:** Total RNA counts (UMIs) detected per cell.
  - Below 500–1,000: Potentially poor-quality cells.
  - Above 20,000–30,000: Investigate possible doublets.

- **pct_counts_mt:** Percentage of counts originating from mitochondrial genes.
  - Below 5%: Generally good.
  - 5–10%: Generally acceptable.
  - 10–20%: Investigate possible cell stress or damage.
  - Above 20–30%: Potentially poor-quality cells.

*Note: These are general guidelines, not universal filtering thresholds.*

## 2. Normalisation (normalize_total)

Adjusts the total RNA counts of each cell to a common target (10,000 here), making cells comparable despite differences in sequencing depth.

## 3. Log Transformation (log1p)

Applies a logarithmic transformation to compress the expression value range, reducing the influence of highly expressed genes.

## 4. Scaling (scale)

Centres and standardises gene expression values across cells, giving genes more comparable numerical ranges. This prevents highly variable genes from dominating PCA simply because of their larger values or variability. It does not reduce their actual biological expression.

## 5. Dimensionality Reduction (PCA)

Simplifies data containing thousands of genes into a smaller number of features (principal components), retaining major patterns of variation between cells.

## 6. Neighbour Graph (neighbors)

Identifies cells with similar expression profiles based on their PCA coordinates, establishing relationships between cells.

## 7. UMAP (umap)

Projects high-dimensional cell relationships into a two-dimensional representation, allowing visualisation of groups of transcriptionally similar cells.

## 8. Global UMAP

Colours cells according to their sample of origin, helping assess sample distribution, overlap and potential sample-specific or batch effects.



QC and Filtering → Normalisation → Log Transformation → Scaling → PCA → Neighbour Graph → UMAP → Visualisation
