# AI4BIO-Dataset
datasets/: final results + processed matrices scripts/: the reproducible pipeline README.md: methodological narrative
Overview of the AI4BIO transcriptomic analysis pipeline:
Publicly available transcriptomic datasets were collected from the Gene Expression Omnibus (GEO), including both RNA-seq and microarray experiments. The overall objective was to generate comparable gene-level differential expression signatures across heterogeneous platforms for downstream integrative analyses.

Data acquisition and expression matrix construction:
For each dataset, raw or processed expression data were downloaded from GEO and organized into a gene-by-sample expression matrix. For microarray datasets, probe identifiers were mapped to official gene symbols using platform-specific annotation files. When multiple probes mapped to the same gene, expression values were collapsed at the gene level using the median value across probes.
In two datasets (GSE51808 and GSE40628), residual duplicated gene symbols remained after initial collapsing due to annotation inconsistencies; these duplicates were subsequently identified and resolved to ensure one row per gene.

Outlier detection and sample quality control:
Multidimensional projection (MDP) analysis was performed on each dataset to assess sample clustering and detect potential outliers. Samples were removed only when clear separation from the main cluster suggested technical artifacts rather than biological variability. In the absence of strong evidence for technical outliers, all samples were retained to preserve statistical power.

Filtering for protein-coding genes:
To improve biological interpretability and cross-dataset comparability, expression matrices were filtered to retain only protein-coding genes based on current gene annotations. This step reduced noise from poorly characterized or low-confidence transcripts while maintaining sufficient coverage of relevant biological pathways.

Handling missing values:
Genes with excessive missing values were excluded prior to differential expression analysis. Specifically, genes with more than 20% missing values across samples were removed to avoid unstable model estimates. Remaining missing values were imputed using gene-wise mean imputation, which is appropriate given the low proportion of missingness after filtering.

Differential expression analysis:
Differential expression was assessed according to the experimental design of each dataset:
RNA-seq datasets
GSE157240: limma
GSE279208: limma-voom with edgeR normalization
Microarray datasets
GSE96656, GSE28405, GSE51808, GSE18090, GSE25226, GSE40628: limma
Linear models were fitted using limma, followed by empirical Bayes moderation of variance estimates. All genes passing preprocessing steps were retained in the differential expression output, regardless of statistical significance, to support downstream ranking-based analyses.

Gene ranking and score computation:
Rather than applying hard thresholds on adjusted p-values and log-fold change—which substantially reduced gene set sizes in some datasets—a ranking-based approach was adopted. A composite score was computed for each gene as:
score = |logFC| × (−log10(P.Value))
Genes were ranked in descending order of this score. For each dataset, two outputs were generated:
- A complete ranked gene list
- A subset containing the top 1,000 genes with the highest scores
This strategy preserved consistent gene set sizes across datasets while prioritizing genes with strong and statistically supported expression changes.
