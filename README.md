# scRNA-seq Analysis Platform

Interactive Shiny-based platform for single-cell RNA sequencing analysis.

## Overview

This platform provides a modular interface for analyzing scRNA-seq data, covering individual-sample processing, multi-sample integration, post-annotation exploration, and cell-cell communication.

## Modules

- **Module I** — Individual-sample scRNA-seq analysis (Import, QC, PCA, Doublet Removal, Cell Cycle Scoring, Annotation)
- **Module II** — Multi-sample integration and annotation (Upload & Merge, Integration, Benchmarking, Annotation)
- **Module III** — Post-annotation exploration (Upload, UMAP Cell-Type, Violin, Dot Plot, Heatmap, Differential Expression, Enrichment)
- **Module IV** — Cell-cell communication (Upload, CellChat compute + visualization, LIANA compute + visualization)

## Project Structure

- `R/` — Shared utilities used across modules (`utils.R`, `mod_save_config.R`)
- `module-I/` — Shiny app for individual-sample analysis
- `module-II/` — Shiny app for multi-sample integration and annotation
- `module-III/` — Shiny app for post-annotation exploration and visualization
- `module-IV/` — Shiny app for cell-cell communication analysis (CellChat + LIANA)
- `test-data/` — Sample datasets (Cell Ranger MEX format and pre-processed Seurat .rds objects)
- `manuscript_R/` — Standalone R scripts for manuscript figures and analysis

## Getting Started

See `module-I/requirement.md` for dependencies and setup instructions.

## License

MIT License — see [LICENSE](LICENSE) for details.
