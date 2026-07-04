# Host Transcriptomics of Cutaneous Leishmaniasis

This repository contains an end-to-end bioinformatics pipeline for analyzing human host transcriptomic profiles from individuals infected with Cutaneous Leishmaniasis (**CL**) compared to healthy controls (**HS**). The workflow integrates raw quality control, pseudoalignment, robust normalization, differential gene expression (DGE), and pathway-level functional enrichment analysis.

---

# 1. Project Directory Structure

The workspace structure is optimized for computational reproducibility:

```text
rna_seq_CL/
├── data/                       # Reference cDNA indices & sample FASTQ sequences
├── results/                    # Pipeline output directories
│   ├── fastqc/                 # Individual sample raw read quality checks
│   ├── mappedReads/            # Kallisto quantification abundance files
│   ├── multiqc/                # Aggregated quality report
│   ├── plots/                  # Visualizations
│   └── tables/                 # Differential expression & enrichment tables
├── scripts/
│   ├── script.sh               # Quality control & Kallisto quantification
│   ├── R_Workflow.R            # Normalization, PCA & Differential Expression
│   └── enrichment_analysis.R   # ORA and GSEA analysis
├── README.md
└── .gitignore
```

---

# 2. Methodology & Analytical Results

## Data Preprocessing & Pipeline Normalization

**Scripts**

- `scripts/script.sh`
- `scripts/R_Workflow.R`

### Approach

Subsampled single-end FASTQ reads were assessed using **FastQC** and quantified with **Kallisto** against the **GRCh38 human cDNA transcriptome**. Lowly expressed genes were removed using a filtering threshold of:

- CPM > 1
- Present in at least 5 samples

Library-size differences were corrected using **Trimmed Mean of M-values (TMM)** normalization implemented in **edgeR**, producing normalized log₂ CPM values suitable for downstream statistical analyses.

### Results

Sequential filtering and TMM normalization successfully aligned sample expression distributions, improving cross-sample comparability while preserving biological variation.

**Figure 1**

`results/plots/CPM_violin_TMM_filtering_comparison.jpg`

Comparison of log₂ CPM distributions across:

- **A:** Raw counts
- **B:** Filtered counts
- **C:** TMM-normalized expression

---

## Dimensionality Reduction & Variance Breakdown

**Script**

- `scripts/R_Workflow.R`

### Approach

Principal Component Analysis (PCA) was performed on the TMM-normalized expression matrix to evaluate overall transcriptomic similarity among samples.

### Results

Samples separated almost entirely according to disease status.

- **PC1 explained 55.6% of total variance**
- Clear separation between Cutaneous Leishmaniasis and Healthy controls demonstrates a strong disease-associated transcriptional signature.

**Figure 2**

`results/plots/PCA_TMM_disease_VS_healthy.png`

PCA of normalized host transcriptomes.

---

## Differential Gene Expression (DGE)

**Script**

- `scripts/R_Workflow.R`

### Approach

Differential expression was performed using the **limma-voom** framework.

Analysis included:

- Voom transformation
- Linear modeling
- Empirical Bayes moderation (eBayes)

Genes were considered significant when satisfying:

- Adjusted *p* < 0.01
- |log₂ Fold Change| > 1

### Results

Thousands of genes showed significant differential expression between infected individuals and healthy controls, indicating widespread host transcriptional remodeling.

**Figure 3**

`results/plots/Vene_diagram.jpg`

Distribution of significantly upregulated and downregulated genes.

**Figure 4**

`results/plots/heatmap_all_genes.png`

Hierarchical clustering heatmap illustrating global transcriptional differences between study groups.

---

## Functional Signaling & Pathway Enrichment

**Script**

- `scripts/enrichment_analysis.R`

### Approach

Functional interpretation of differentially expressed genes included:

- **Over-Representation Analysis (ORA)** using **gprofiler2**
- **Gene Set Enrichment Analysis (GSEA)** using **clusterProfiler**
- **MSigDB C2** curated pathway collections

### Results

Enrichment analyses identified extensive activation of innate immune responses together with strong B-cell–associated signaling pathways.

**Figure 5**

`results/plots/Upregulated_enrichment.jpg`

Manhattan plot of enriched biological pathways among upregulated genes.

**Figure 6**

`results/plots/Enrichment_bubble_plot.jpg`

Bubble plot summarizing significantly enriched Reactome pathways.

**Figure 7**

`results/plots/GSEA_result.png`

Representative GSEA enrichment plot highlighting activated follicular B-cell signaling pathways.

---

# 3. Core Biological Inferences

## Host Defense Mobilization

Cutaneous Leishmaniasis induces extensive host immune activation characterized by increased expression of genes involved in:

- Complement cascade
- Fc gamma receptor (FCGR) signaling
- Phagocytosis
- Innate immune defense

These findings support enhanced macrophage-mediated clearance of antibody-coated parasites.

---

## Humoral Immune Activation

Gene set enrichment revealed strong activation of pathways associated with:

- Follicular B-cell responses
- Extrafollicular B-cell activation
- Humoral immunity

High enrichment scores indicate active antibody-mediated immune responses within infected tissues.

---

## Immune Suppression Signature

In parallel with immune activation, several host pathways were significantly downregulated, suggesting localized immune suppression or parasite-driven immune evasion within lesion sites.

---

# Software & Dependencies

## Quantification

- FastQC
- MultiQC
- Kallisto

## R Packages

- edgeR
- limma
- ggplot2
- pheatmap
- EnhancedVolcano
- gprofiler2
- clusterProfiler
- msigdbr
- enrichplot
- fgsea

---

# Workflow Summary

```text
FASTQ Files
     │
     ▼
 FastQC
     │
     ▼
 MultiQC
     │
     ▼
 Kallisto Quantification
     │
     ▼
 Import into R
     │
     ▼
 Gene Filtering
     │
     ▼
 TMM Normalization
     │
     ▼
 PCA
     │
     ▼
 limma-voom Differential Expression
     │
     ▼
 ORA (gprofiler2)
     │
     ▼
 GSEA (clusterProfiler)
     │
     ▼
 Biological Interpretation
```

---

# Repository Structure

```text
results/
├── fastqc/
├── mappedReads/
├── multiqc/
├── plots/
│   ├── CPM_violin_TMM_filtering_comparison.jpg
│   ├── PCA_TMM_disease_VS_healthy.png
│   ├── Vene_diagram.jpg
│   ├── heatmap_all_genes.png
│   ├── Upregulated_enrichment.jpg
│   ├── Enrichment_bubble_plot.jpg
│   └── GSEA_result.png
└── tables/
```

---

# Citation

If you use this pipeline, please cite the corresponding software packages used for:

- Kallisto
- edgeR
- limma
- clusterProfiler
- gprofiler2
- MSigDB