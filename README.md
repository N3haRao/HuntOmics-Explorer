# 🚀 HuntOmics Explorer

**An interactive RNA-seq analysis dashboard for Huntington’s Disease transcriptomics.**

HuntOmics Explorer is a Shiny-based analytical platform for exploring transcriptomic alterations in Huntington’s Disease (HD). The application enables interactive investigation of RNA-seq data from human post-mortem BA9 brain tissue, integrating exploratory data analysis, differential expression, and functional enrichment within a unified interface.

This tool was developed to demonstrate reproducible RNA-seq workflows and interactive genomic data visualization in R.

---

## Scientific Context

Huntington’s Disease (HD) is a progressive neurodegenerative disorder caused by CAG repeat expansion in the *HTT* gene. Transcriptomic dysregulation in cortical regions, including Brodmann Area 9 (BA9), is a hallmark of disease pathology.

HuntOmics Explorer allows researchers to:

* Investigate global transcriptional structure (PCA, clustering)
* Explore gene-level expression patterns
* Identify differentially expressed genes (HD vs Control)
* Interpret functional enrichment signatures

The application is built around GEO dataset:

📁 **GSE64810**
*mRNA-Seq Expression profiling of human post-mortem BA9 brain tissue for Huntington’s Disease and neurologically normal individuals*
[https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64810](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64810)

---

## Key Features

### Sample Metadata Exploration

* Summary statistics of sample characteristics
* Distribution plots of numeric variables
* Interactive inspection of study design

### Gene Expression Matrix Exploration

* Summary of count distributions
* User-selected gene scatter plots
* Hierarchical clustering heatmaps
* Principal Component Analysis (PCA) visualization

### Differential Expression Analysis

* Interactive DE gene table
* Volcano plot visualization
* Significance filtering

### Gene Set Enrichment Analysis

* Gene Ontology (GO) enrichment
* Pathway-level summaries
* Bar plots and dot plots of enriched terms

---

## Tech Stack

* **R**
* **Shiny**
* **DESeq2**
* **ggplot2**
* **pheatmap**
* **clusterProfiler**
* **dplyr**
* **GEO-derived RNA-seq dataset**

---

## Repository Structure

```
HuntOmics-Explorer/
│
├── app.R
├── README.md
│
├── data/
│   ├── Counts.csv
│   ├── DE.csv
│   └── Metadata.csv
│
└── demo/
    └── demo.mp4
```

---

## Installation

### Clone the repository

```bash
git clone https://github.com/N3ha-Rao/HuntOmics-Explorer.git
cd HuntOmics-Explorer
```

### Install required R packages

```r
install.packages(c(
  "shiny",
  "ggplot2",
  "pheatmap",
  "clusterProfiler",
  "dplyr"
))

# If needed:
# BiocManager::install("DESeq2")
```

### Run the application

```r
library(shiny)
runApp()
```

---

## Data Preparation

The application expects:

* `Metadata.csv` — sample-level metadata
* `Counts.csv` — gene expression count matrix
* `DE.csv` — differential expression results

All input files correspond to GEO accession **GSE64810** and must follow the expected column formats defined within `app.R`.

---

## Reproducibility

All analyses are based on publicly available data (GSE64810).
The repository includes preprocessed count matrices and differential expression results derived using standard RNA-seq workflows.

This project emphasizes:

* Reproducible computational analysis
* Transparent data organization
* Modular visualization design
