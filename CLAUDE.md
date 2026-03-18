# scRNA-seq Analysis Platform

Interactive Shiny-based platform for single-cell RNA sequencing analysis, progressing from individual-sample processing (Module I) through combined multi-sample analysis (future modules).

## Project Structure

- `module-I/` — Shiny app for individual-sample scRNA-seq analysis
- `test-data/` — Sample datasets (Cell Ranger MEX format and pre-processed Seurat .rds objects)

## Module I: Individual-Sample Analysis

A modular Shiny web application implementing a linear analysis pipeline for individual scRNA-seq samples:

**Import → QC Filtering → PCA → Doublet Removal → Cell Cycle Scoring → Biomarker Visualization**

### Modules (in `module-I/modules/`)

| Module | File | Purpose |
|--------|------|---------|
| Import | `mod_import.R` | Load Cell Ranger MEX output into Seurat object, compute QC metrics (mito/ribo/hemoglobin %) |
| QC | `mod_qc.R` | Interactive filtering with adjustable thresholds and visual pass/fail classification |
| PCA | `mod_pca.R` | Normalization (LogNormalize or SCTransform) and PCA with elbow/loading/heatmap plots |
| Doublet | `mod_doublet.R` | DoubletFinder-based doublet detection with UMAP visualization |
| Cell Cycle | `mod_cellcycle.R` | Cell cycle phase scoring (mouse/human) with UMAP display |
| Biomarker | `mod_biomarker.R` | Gene expression visualization from uploaded biomarker CSV lists |

**Disabled/incomplete modules:** `mod_singler_annotation.R`, `mod_annotation_manual.R`, `mod_cluster.R`

### Architecture

- Each module is a Shiny module (namespaced UI + server pair)
- Data flows as Seurat objects through reactive returns between modules
- Input: 10X Genomics Cell Ranger MEX format
- Output: Seurat .rds files downloadable at each step
- Entry point: `module-I/app.R`

## Tech Stack

- **Language:** R
- **Framework:** Shiny
- **Core packages:** Seurat, DoubletFinder, ggplot2, plotly, tidyverse, glmGamPoi, shinyFiles, shinyjs
- **Setup:** See `module-I/requirement.md` for dependencies and installation

## Planned Modules

- **Module I** — Individual-sample analysis (in development)
- **Module II** — TBD (future)
- **Module III** — TBD (future)
- **Module IV** — TBD (future)

The platform is designed to progress from individual-sample processing to combined multi-sample analysis across the four modules.
