# From Pipelines to Plots: Reproducible and Visual Data Analysis


Welcome to the GitHub repository for the **Australian and New Zealand Society for Immunology (ASI)** workshop:

> **From Pipelines to Plots: Reproducible and Visual Data Analysis**

This repository contains the slides, workshop code, datasets and supporting material used throughout the workshops delivered on the 9th and 10th of July 2026 by Dr Jennifer Currenti, Dr Aaron Beasley, and Dr Nicola Principe with assistance from Dr Nataliya Slater.

---

# Workshop Overview

Modern biological research increasingly relies on computational analyses that are reproducible, scalable and easy to interpret. This two-day workshop provides a practical introduction to building reproducible analysis workflows and producing publication-quality visualisations using R.

The workshop is designed for researchers who already have some basic experience with R and would like to develop more efficient, reproducible and automated analysis workflows.

---

# Workshop Structure

## Day 1 – Data Reproducibility and Automation

Run by:
- Dr Jennifer Currenti
- Dr Aaron Beazley

Topics include:

- Git and version control
- Managing computational environments with Conda
- Workflow automation with Nextflow
- Docker and Singularity
- Building reproducible analysis pipelines
- Best practices for reproducible computational research

---

## Day 2 – Data Fundamentals and Visualisation

Run by:
- Dr Jennifer Currenti
- Dr Aaron Beasley
- Dr Nicola Principe

Topics include:

- Writing reusable functions
- Automating analyses with loops
- Working with vectors, lists, data frames and tibbles
- Principal Component Analysis (PCA)
- Heatmaps with ComplexHeatmap
- Box plots, violin plots and ridge plots
- Dot plots
- UpSet plots
- Sankey diagrams
- Circos plots

---

# Repository Structure

```text
.
├── slides/                 # Workshop presentations
├── code/                   # R scripts used during the workshop
├── data/                   # Input datasets
├── results/                # Output generated during the workshop
├── environments/           # Conda environment files
├── README.md
└── LICENSE
```

---

# Workshop Files

The repository includes:

- Complete workshop slides
- Workshop scripts
- Fully annotated instructor solutions/html files
- Figures generated throughout the workshop
- Conda environment files (where applicable)

---

# Software Requirements

Participants should have installed prior to the workshop:

- R (≥ 4.4)
- RStudio
- Git
- Conda (Miniforge recommended)

The workshop has been tested on both Windows and macOS.

---

# Data

The visualisation exercises use publicly available single-cell RNA sequencing data from:

> Zhang Z. *et al.* **STARTRAC analyses of scRNA-seq data from tumor models reveal T cell dynamics and therapeutic targets.** *Journal of Experimental Medicine.* 2021.

Publication:
https://rupress.org/jem/article/218/6/e20201329/212026/STARTRAC-analyses-of-scRNAseq-data-from-tumor

Dataset (GEO):
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168944

The dataset includes:

- Gene expression (GEX)
- T-cell receptor (TCR) sequencing
- Cell metadata

---

# Learning Outcomes

By the end of the workshop you should be able to:

- Organise computational projects using reproducible workflows
- Manage software environments using Conda
- Use Git for version control
- Write reusable R functions
- Automate repetitive analyses with loops
- Manipulate data using the tidyverse
- Perform and interpret PCA
- Produce publication-quality figures in R
- Generate heatmaps with ComplexHeatmap
- Create a range of commonly used biological data visualisations

---

# Running the Workshop

Clone the repository:

```bash
git clone https://github.com/JenC9292/ASI-WA-branch-workshop-2026
```

Move into the repository:

```bash
cd ASI_Workshop
```

Open the RStudio Project or open the scripts located in the `code/` directory.

---

# Additional Resources

## R

- https://www.r-project.org/
- https://posit.co/download/rstudio-desktop/

## Git

- https://git-scm.com/

## Conda

- https://conda-forge.org/download/

## Nextflow

- https://www.nextflow.io/

## ComplexHeatmap

- https://jokergoo.github.io/ComplexHeatmap-reference/book/

---

# Citation

If you use these workshop materials in your teaching or training, please cite or acknowledge the workshop appropriately. Materials were created by Dr Jennifer Currenti, Dr Aaron Beasley, and Dr Nicola Principe.

---

# Contact

Email:
asiwa@immunology.org.au
jennifer.currenti@curtin.edu.au

---

These materials were developed for the **From Pipelines to Plots: Reproducible and Visual Data Analysis** workshop.