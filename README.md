# BENG 183 Final Project

## scRNA-seq and CRISPR Analysis Workflow for Malignant Melanoma Cells

Author: Nina Marie Zinke

---

## Purpose of Workflow and scRNA-seq Integration

In this project, we analyze single-cell RNA-seq and CRISPR data that has been generated from 1,000 A375 malignant melanoma cells. Each are stably expressing sCas9 and are transduced with either a non-targeting sgRNA or an sgRNA that is targeting Rab1a. Rab1a is a gene involved in vesicle trafficking and melanoma cell survival. The dataset includes paired Gene Expression as well as CRISPR Guide Capture libraries that have been generated using Chromium GEM-X Single Cell 5’ v3 chemistry, sequenced on an Illumina NovaSeq6000 (~80,000 read pairs per cell).

> **Note on this repository:** The code in this repo (`analysis.Rmd`) runs the **un-integrated** Seurat workflow: preprocessing, clustering, UMAP, pseudobulk comparison, differential expression, and marker visualization. Seurat v5 CCA integration, described in detail below, is **not performed** in the current analysis. It is included as an optional block (`RUN_INTEGRATION <- FALSE`) in the code and is **planned to be added in a later update**. The integration sections of the original project handout are kept below so the intended full workflow is documented.

### Rationale for Using Seurat

We chose to use Seurat as the main primary analysis framework for this project because it offers a very comprehensive and well-validated workflow for single-cell analysis. Seurat is widely used for data analysis and offers intuitive visualization tools, including UMAP, t-SNE, and feature plots, that allow us to interpret the clustering results and gene expression patterns effectively. Seurat’s broad use in the single-cell research community means that its methods are widely benchmarked, peer-reviewed and supported by extensive documentation.

### Integration of scRNA-seq data

Because perturbation experiments are often introducing both biological and technical differences across sgRNA conditions, the process of dataset integration becomes a very important step. Dataset integration is the process of combining and harmonizing single-cell RNA sequencing datasets. It can be done across experimental batches, donors, or conditions.

* Overcome batch effects
  * Datasets from different experiments often have technical or biological variations (batch effects) that can obscure the true biological differences
* Increase statistical power
  * Combining data increases the sample size, which improves the ability to detect rare cell populations that might not be present or detectable in any single dataset
* Enable cross-study comparisons
  * Allows for the systematic comparison of cell types across multiple samples, conditions, and even different studies
* Integrate multi-omic data
  * Can combine scRNA-seq with other single-cell data types

**Status in this repo:** integration could be added to this analysis and is planned as a future addition. The current results are from the un-integrated analysis.

---

## Integration Goals

The following tutorial is designed to give you an overview of the kinds of comparative analyses on complex cell types that are possible using the Seurat integration procedure. Here, we address a few key goals:

* Identify cell subpopulations that are present in both datasets
* Obtain cell type markers that are conserved in both nontargeted (control) and targeted cells
* Compare the datasets to observe effects of the removal of the Rab1a gene

---

## Repository Structure and How to Run

```
BENG183_project/
|-- BENG183_project.Rproj
|-- data/
|   |-- <...>_filtered_feature_bc_matrix.h5
|   `-- crispr_analysis/protospacer_calls_per_cell.csv
|-- analysis.Rmd            <- analysis code
`-- results/                <- created automatically (figures, DE table, Seurat object)
```

1. Download the dataset (the filtered feature-barcode matrix `.h5` and `protospacer_calls_per_cell.csv`) from the 10x Genomics website and place the files in `data/` as shown above.
2. Open `BENG183_project.Rproj` in RStudio.
3. In the "Locate input files" chunk of `analysis.Rmd`, change the `base` path to the location of the project folder on your machine.
4. Run the notebook top to bottom. Missing packages are installed automatically.

**Requirements:** R, Seurat v5 (the original run used 5.3.1), SeuratObject, hdf5r, ggplot2, ggrepel, patchwork, dplyr. (`metap` is only needed for the optional integration block.)

---

## Setup of Seurat Objects

* Import library(hdf5r)
* Install dataset
* Load dataset (Read10X\_h5(“”), read.csv(“”))

The Gene Expression matrix is used to build the Seurat object (`CreateSeuratObject`, `min.cells = 3`, `min.features = 200`, giving roughly 22,238 features x 1,163 cells). The CRISPR Guide Capture assignments are read from the protospacer calls file and added as metadata.

---

## Performing Analysis Without Integration

If we cluster the dataset without integration, we would typically observe:

* Separation that is driven largely by sgRNA identity
* A375 cells that are clustering per perturbation rather than by intrinsic biological states
* Some difficulty when trying to identify shared melanoma programs

**This is the analysis performed by the code in this repo.** Clustering and UMAP are run directly on the PCA of the un-integrated object.

---

## Performing Integration

*(Not performed in the current code. Planned as a later addition. The details below are from the original project handout.)*

Using Seurat v5 CCA-based integration, we align the Non\_Target-1\_MS and Rab1a-sgRNA cells into a shared low-dimensional space. The goal is not to completely get rid of the true perturbation effects, but to:

* Identify conserved malignant cell states
* Provide a fair basis for differential expression testing

Preprocessing steps before integration

* Add CRISPR metadata:
  * a375$num\_features <- crispr\_data$num\_features
  * a375$feature\_call <- crispr\_data$feature\_call
  * a375$num\_umis <- crispr\_data$num\_umis
* Normalization - NormalizeData(a375)
* Variable features - FindVariableFeatures(a375)
* Scaling - ScaleData(a375)
* PCA - RunPCA(a375, npcs = 30)

Integrate dataset (or select features for integration if there are multiple samples)

* IntegrateLayers()
  * Aligns layers or modalities inside ONE object
  * Method returns method returns a dimensional reduction (i.e. integrated.cca) which can be used for visualization and unsupervised clustering analysis
* For multiple Seurat objects
  * FindIntegrationAnchors()
    * Identify matching cells or “anchor points” between datasets
  * IntegrateData()
    * Merges separate Seurat objects into a batch-corrected assay: crispr\_data <- crispr\_data[colnames(a375), ]

**Where this stands in the code:** `analysis.Rmd` contains an optional section ("Seurat v5 CCA integration") that splits the object by `feature_call` (RAB1A-2 vs Non\_Target-1\_MS), runs `IntegrateLayers()` with `CCAIntegration`, re-clusters on `integrated.cca`, and runs `FindConservedMarkers()`. It is switched off by default (`RUN_INTEGRATION <- FALSE`) and has not yet been tested or validated on this dataset, so none of the results below come from it.

---

## Downstream Analysis

* FindNeighbors(a375, dims = 1:20)
  * Builds a k-nearest neighbor graph based on the chosen dimensional reduction
* FindClusters(a375, resolution = 0.5)
  * Number of nodes/cells: 1163
    * Number of edges/connections between cells: 38856
    * Maximum modularity in 10 random starts: 0.8187 // strong clusters
    * Number of communities/clusters: 7
  * Performs unsupervised clustering on the neighbor graph
  * Uses Louvain / Leiden algorithm internally
  * Assigns each cell to a cluster
* RunUMAP(a375, dims = 1:20)
  * 2Dmap with embedding parameters a = 0.9922 b = 1.112
  * UMAP provides a low-dimensional representation
  * Useful for qualitatively seeing whether Rab1a perturbation shifts cell-state distributions or creates new subpopulations
  * FeaturePlot(a375, group.by = nCount\_RNA) UMAP: Total RNA content per cell - high or low
  * DimPlot(a375, group.by = num\_umis) UMAP: Clusters found by Seurat - check clusters
  * DimPlot(a375, group.by = feature\_call(RAB1A-2, Non\_Target-1\_MS, Non\_Target-1\_MS|RAB1A-2, Unknown) UMAP: locates the different conditioned cells
  * FeaturePlot(a375, group.by = num\_umi) UMAP: UMI per cell for sequencing depth, data quality and processing

In the current code, the neighbor graph and UMAP are built from the un-integrated PCA. In the planned integrated version, a single reduction (integrated.cca) would drive the clustering and UMAP visualization.

### UMAP results

> **[UMAP, total RNA content per cell]**
> <img width="1800" height="1500" alt="umap_nCount_RNA" src="https://github.com/user-attachments/assets/1f9b312f-75c1-4ca3-b697-6565bfa28105" />
 This is the UMAP colored by `nCount_RNA` (total RNA content per cell, high or low).

> **[UMAP, Seurat clusters]**
> <img width="1800" height="1500" alt="umap_clusters" src="https://github.com/user-attachments/assets/2f559e2d-7629-4cb4-8a8d-c1b1903981e3" />
 This is the UMAP colored by the clusters found by Seurat (7 clusters at resolution 0.5).

> **UMAP, guide assignment]**
> <img width="2100" height="1500" alt="umap_feature_call" src="https://github.com/user-attachments/assets/98bc411e-56da-4835-9da0-ecf1e2163e6a" />
 This is the UMAP colored by `feature_call` (RAB1A-2, Non\_Target-1\_MS, Non\_Target-1\_MS|RAB1A-2, Unknown), showing where the different conditioned cells sit.

> **[UMAP, UMIs per cell]**
> <img width="1800" height="1500" alt="umap_num_umis" src="https://github.com/user-attachments/assets/7cdd2b5c-79a9-425e-9465-4f3caa418d86" />
 This is the UMAP colored by `num_umis` (guide UMIs per cell), for sequencing depth, data quality and processing.

> **[PLACEHOLDER (future): integrated UMAP]**
> Once CCA integration is added, place `results/umap_integrated_cca.png` here (UMAP on `integrated.cca`, colored by `feature_call` and `seurat_clusters`).

---

## Identifying Conserved Cell Type Markers

Even within a homogeneous A375 melanoma cell line, cells can exist in multiple transcriptional states (e.g. stress-like, proliferative, EMT-like states). We use FindConservedMarkers() in order to identify genes that consistently mark each state regardless of sgRNA condition.

*(FindConservedMarkers() depends on the integration step above, so it is part of the planned addition. The current code uses FindMarkers() between the two guide groups.)*

Why this matters in CRISPR experiments

* Ensures cluster annotations are biologically meaningful and not driven by perturbation
* Helps distinguish “true” melanoma cell states from perturbation-induced artifacts
* Helps confident identification of Rab1a-sensitive vs. Rab1a-insensitive states

### Concept

For conserved markers you want genes that robustly mark a cell cluster regardless of experimental condition/donor/batch. Seurat’s FindConservedMarkers() performs per-group tests and then meta-analyzes p-values.

### FindMarkers()

* genes\_1 <- FindMarkers(a375, ident.1 = "RAB1A-2", ident.2 = "Non\_Target-1\_MS")
* Runs one differential expression test between two sets of cells you specify (ident.1 vs ident.2)
* Treats all cells in ident.1 as one group and all in ident.2 as one group.
* TOP 6 markers were adjusted based on p-value (Benjamini–Hochberg correction), diff expressed genes show both the effectiveness of the CRISPR perturbation and downstream transcriptional effects (like possibly changes in signal cascades and sequences itself)
  * **"GAGE2A", "RAB1A", "TMPRSS15", "RPLP0", "RPSA", "SERPINE2"**

### FindConservedMarkers()

* Runs the specified test (default = Wilcoxon) within each level of a grouping variable (e.g. condition)
* Combines the p-values across groups with meta-analysis (Seurat uses MetaDE methods)
* Reports per-group statistics and a combined p-value.

FindMarkers() is used when you want DE between two groups in one dataset, FindConservedMarkers() is used across multiple groups

### Scatterplot

* Showed the genes expression for the RAB1A-2 and Non\_Target-1\_MS, and labeled the genes from the markers found above
* Code:

```r
p1 <- CellScatter( aggregate_crispr, "RAB1A-2", "Non_Target-1_MS", highlight = genes.to.label )
p2 <- LabelPoints( plot = p1, points = genes.to.label, repel = TRUE )
p1+p2
```

*(In `analysis.Rmd` the pseudobulk object is named `agg`, and the groups are written `"RAB1A-2"` and `"Non-Target-1-MS"`, because `AggregateExpression` converts underscores to dashes.)*

> **pseudobulk scatter plot]**
> <img width="3000" height="1500" alt="pseudobulk_scatter" src="https://github.com/user-attachments/assets/0e06ce7d-960b-4cb8-95af-ed2d5851d195" />
 This is the scatter plot of aggregated gene expression, RAB1A-2 vs Non\_Target-1\_MS, with the six marker genes labeled.

### Functions and Interpretation

To get more information on the markers, made the featureplots graphing each expression levels in different clusters, for each marker:

**GAGE2A**: red in RAB1A-2 cluster → upregulation in expression
**RAB1A**: mostly grey in RAB1A-2 cluster, so CRISPR worked
**TMPRSS15**: mildly red in RAB1A-2 cluster
**RPLP0**: red in most clusters, so it is acting as a control marker
**RPSA**: red in most clusters, so it is acting as a control marker
**SERPINE2**: mildly red in RAB1A-2 cluster, and brighter red in Non\_Target-1\_MS

> **[marker feature plots]**
> <img width="2700" height="3000" alt="featureplot_markers" src="https://github.com/user-attachments/assets/110da7be-7524-46e4-93e7-cd2ace71839a" />
 This is the panel of UMAP feature plots for GAGE2A, RAB1A, TMPRSS15, RPLP0, RPSA and SERPINE2.

### Volcano Plot

Code:

```r
volcano <- ggplot(de, aes(x = avg_log2FC, y = log10_pval, color = significant, label = label)) +
  geom_point(alpha = .5, size = 3) + scale_color_manual(values = c("red", "blue")) +
  geom_text_repel(show.legend = FALSE) + theme_minimal(base_size = 12) +
  xlab("Log2 Fold Change") + ylab("-Log10 P-value") +
  ggtitle("Volcano Plot: RAB1A-2 vs Non_Target-1_MS")
print(volcano)
```

*(In `analysis.Rmd`, `log10_pval` is calculated as -log10 of the adjusted p-value, and a gene is marked significant when adjusted p < 0.05 and |avg\_log2FC| > 1.)*

> **[volcano plot]**
> <img width="2100" height="1800" alt="volcano_RAB1A_vs_NT" src="https://github.com/user-attachments/assets/f0346b28-ccfe-4339-abdc-0fedc86892f4" />
 This is the volcano plot of differential expression, RAB1A-2 vs Non\_Target-1\_MS, with the six marker genes labeled.

Results:

* Upregulated genes - GAGE2A, TMPRSS15
* Downregulated genes - RAB1A, SERPINE2
* No change - RPLP0, RPSA

The full differential expression table is saved to `results/DE_RAB1A_vs_NonTarget.csv`.

---

## Identifying Differential Expressed Genes Across Conditions (Perform DE on the ORIGINAL ASSAY)

The purpose of differential expression (DE) analysis is a key step in identifying genes with significantly different expression levels between two or more conditions, such as a treated vs. untreated group. In this case it's the non-targeting control vs. Rab1a-sgRNA-perturbed cells. It includes using read count data input and performing statistical analysis to discover quantitative changes in expression. For the particular dataset used in this example, we are interested in determining how the transcriptional system of A375 melanoma cells changes when Rab1a is perturbed.

Before performing the DE testing, we first worked within the integrated Seurat object to make sure that any comparisons between conditions are not messed up by batch effects etc. After integration, we cluster the combined dataset and compute a shared UMAP in order to identify transcriptionally distinct melanoma cell states. Working with these integrated clusters is very essential because it enables within-cell-state differential expression. This is more interpretable than comparing all cells globally. Global comparisons can be dominated by potential changes in proportions of the cluster, rather than changes that have to do with the gene expression itself.

*(In the current code, DE is run on the un-integrated object, comparing guide groups globally. Within-cluster DE on integrated clusters is part of the planned integration addition.)*

Preprocessing before DE

* Do standard QC and filtering (mito%, nFeature\_RNA, nCount\_RNA).
* Normalize per-sample ex: SCTransform
* Run PCA per sample

Aggregate counts per (sample × celltype)

* AggregateExpression(a375, group.by = "feature\_call")
  * Aggregates all cells within the same cell type
  * Maximize statistical analysis by using all cells

Define groups

* After CCA-based integration, we use FindNeighbors() and FindClusters() in order to identify biologically meaningful clusters that are shared between control and Rab1a-targeted cells

FindMarkers()

* Find the genes that are different between control and treated groups (DE features)
* The bulk of Seurat’s differential expression features can be accessed through the FindMarkers() function
* By default, Seurat performs differential expression (DE) testing based on the non-parametric Wilcoxon rank sum test
* To test for DE genes between two specific groups of cells, specify the ident.1 and ident.2 parameters

Visualization

* Feature plots
* Volcano plots

---

## Discussion

RAB1A and the secreted factor SERPINE2 are downregulated, indicating that the CRISPR silencing was effective. Upregulation of stress/dedifferentiation markers GAGE2A and TMPRSS15 is most likely compensation of the cell in response to downregulation of the other genes. Since neither RPLP0 and RPSA changed significantly in the DE results, they act as controls in this context and are generally considered to be housekeeping genes.

---

## Planned Future Work

* Run and validate Seurat v5 CCA integration (`IntegrateLayers()` with `CCAIntegration`) on the RAB1A-2 and Non\_Target-1\_MS cells, as described in the integration sections above
* Cluster and visualize on the `integrated.cca` reduction
* Identify conserved markers across sgRNA conditions with `FindConservedMarkers()`
* Repeat differential expression within integrated clusters

---

## Output Files

All outputs are written to `results/`:

| File | Contents |
|---|---|
| `umap_nCount_RNA.png` | UMAP, total RNA per cell |
| `umap_clusters.png` | UMAP, Seurat clusters |
| `umap_feature_call.png` | UMAP, guide assignment |
| `umap_num_umis.png` | UMAP, guide UMIs per cell |
| `pseudobulk_scatter.png` | Pseudobulk scatter, RAB1A-2 vs Non\_Target-1\_MS |
| `DE_RAB1A_vs_NonTarget.csv` | Differential expression table |
| `volcano_RAB1A_vs_NT.png` | Volcano plot |
| `featureplot_markers.png` | Marker gene feature plots |
| `a375_seurat_object.rds` | Saved Seurat object |
