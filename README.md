# scRNA-seq Analysis Platform

Interactive Shiny-based platform for single-cell RNA sequencing analysis.

## Overview

This platform provides a modular interface for analyzing scRNA-seq data, progressing from individual-sample processing through combined multi-sample analysis.

## Modules

- **Module I** — Individual-sample scRNA-seq analysis (Import, QC, PCA, Doublet Removal, Cell Cycle Scoring, Annotation)
- **Module II** — Multi-sample integration and annotation (Upload & Merge, Integration, Benchmarking, Annotation)
- **Module III–IV** — Planned

## Project Structure

- `R/` — Shared utilities used across modules (`utils.R`, `mod_save_config.R`)
- `module-I/` — Shiny app for individual-sample analysis
- `module-II/` — Shiny app for multi-sample integration and annotation
- `test-data/` — Sample datasets (Cell Ranger MEX format and pre-processed Seurat .rds objects)
- `manuscript_R/` — Standalone R scripts for manuscript figures and analysis

## Getting Started

See `module-I/requirement.md` for dependencies and setup instructions.

## License

MIT License — see [LICENSE](LICENSE) for details.
