# scRNA-seq Analysis Platform

Interactive Shiny-based platform for single-cell RNA sequencing analysis, progressing from individual-sample processing (Module I) through combined multi-sample analysis (Module II).

## Project Structure

- `R/` — Shared code used by all modules
  - `utils.R` — Package management (`load.or.install`) and common setup
  - `mod_save_config.R` — Reusable Shiny sub-module for configurable Seurat object export
- `module-I/` — Shiny app for individual-sample scRNA-seq analysis
  - `modules/` — Shiny module files (one per pipeline step)
  - `app.R` — Entry point: loads packages, defines UI layout, wires module servers
- `module-II/` — Shiny app for multi-sample integration and annotation
  - `modules/` — Shiny module files (one per pipeline step)
  - `app.R` — Entry point: loads packages, defines UI layout, wires module servers
- `test-data/` — Sample datasets (Cell Ranger MEX format and pre-processed Seurat .rds objects)
- `manuscript_R/` — Standalone R scripts for manuscript figures and analysis (not part of the Shiny app)

## Module I: Individual-Sample Analysis

A modular Shiny web application implementing a linear analysis pipeline for individual scRNA-seq samples:

**Import → QC Filtering → PCA → Doublet Removal → Cell Cycle Scoring → Annotation**

### Active Modules (in `module-I/modules/`)

| Module | File | Purpose |
|--------|------|---------|
| Import | `mod_import.R` | Load Cell Ranger MEX output into Seurat object, compute QC metrics (mito/ribo/hemoglobin %) |
| QC | `mod_qc.R` | Interactive filtering with adjustable thresholds and visual pass/fail classification |
| PCA | `mod_pca.R` | Normalization (LogNormalize or SCTransform) and PCA with elbow/loading/heatmap plots |
| Doublet | `mod_doublet.R` | DoubletFinder-based doublet detection with UMAP visualization |
| Cell Cycle | `mod_cellcycle.R` | Cell cycle phase scoring (mouse/human) with UMAP display |
| Annotation | `mod_annotation.R` | Orchestrates SingleR auto-annotation and manual annotation sub-modules |

### Architecture

- Each module is a Shiny module (namespaced UI + server pair)
- Data flows as Seurat objects through reactive returns between modules
- Tabs are progressively enabled as each pipeline step completes
- Input: 10X Genomics Cell Ranger MEX format
- Output: Seurat .rds files downloadable at each step (via shared `mod_save_config`)

## Module II: Multi-Sample Integration

A modular Shiny web application for integrating and annotating multiple scRNA-seq samples:

**Upload & Merge → Integration → Benchmarking → Annotation**

### Active Modules (in `module-II/modules/`)

| Module | File | Purpose |
|--------|------|---------|
| Upload & Merge | `mod_upload_merge.R` | Upload multiple Seurat .rds objects, assign metadata, merge |
| Integration | `mod_integrate.R` | Normalization (LogNormalize or SCTransform) and batch integration (CCA, RPCA, Harmony, FastMNN) |
| Benchmarking | `mod_benchmark.R` | Integration quality metrics (ASW, LISI) with ranked comparison |
| Annotation | `mod_annotation.R` | Interactive UMAP visualization, subset analysis, clustering, DE analysis, and manual annotation |

### Architecture

- Each module is a Shiny module (namespaced UI + server pair)
- Data flows as Seurat objects through reactive returns between modules
- Input: Pre-processed Seurat .rds objects (output of Module I)
- Output: Annotated Seurat .rds files (via shared `mod_save_config`)

## Tech Stack

- **Language:** R
- **Framework:** Shiny
- **Core packages:** Seurat, DoubletFinder, ggplot2, plotly, tidyverse, glmGamPoi, shinyFiles, shinyjs, ggpubr, shinycssloaders
- **Module II additional packages:** SeuratWrappers, SeuratObject, presto, cluster, lisi, Matrix

## Naming Conventions

All project code uses **dot-separated** naming following R convention. This applies to:

- **Function names:** `mod.import.ui`, `mod.qc.server`, `load.or.install`
- **Variable names:** `seurat.obj`, `folder.name`, `current.markers`
- **Function parameters:** `seurat.obj.qc`, `ui.testing`, `pca.dims`
- **Shiny input/output IDs:** `qc.metric.1`, `pca.run`, `marker.upload`, `plot.raw.vln`
- **Return list keys:** `seurat.obj`, `completed`
- **Constants:** `UI.TESTING`

**Exceptions** (do not rename):
- External package function names (`geom_point`, `theme_bw`, `case_when`, etc.)
- Seurat metadata column names (`nFeature_RNA`, `nCount_RNA`)
- Shiny element IDs used in JavaScript/CSS selectors (`main_tabs`)
- File names on disk (`mod_import.R`)

## Module Conventions

### Function naming
- UI function: `mod.<name>.ui(id)`
- Server function: `mod.<name>.server(id, ...)`

### Server return value
Every module server returns a named list:
```r
return(list(seurat.obj = <reactive>, completed = <reactiveVal>))
```
- `seurat.obj` — the processed Seurat object (reactive or eventReactive)
- `completed` — a `reactiveVal(FALSE)` set to `TRUE` when the step finishes

### Internal reactive naming
- Intermediate data reactives: `data.<module>` (e.g., `data.qc`, `data.pca`, `data.dbl`)
- Plot input reactives: `plot.input.<module>` (e.g., `plot.input.qc`, `plot.input.pca`)

### Plot rendering
- All `renderPlot` calls use `res = 96`

### UI patterns
- Action buttons use `class = "btn-success"`
- Each pipeline step is a `tabPanel` within the main `tabsetPanel`
- Sections within tabs are wrapped in `wellPanel`

## Planned Modules

- **Module I** — Individual-sample analysis (active)
- **Module II** — Multi-sample integration and annotation (active)
- **Module III** — TBD (future)
- **Module IV** — TBD (future)

The platform is designed to progress from individual-sample processing to combined multi-sample analysis across the four modules.

**Note:** Module II uses underscore-separated naming (`mod_upload_merge_ui`, `uploaded_seurat`) inherited from its original codebase, unlike Module I's dot-separated convention. New shared code in `R/` follows the dot-separated convention.
