# scNexus — scRNA-seq Analysis Platform

Interactive Shiny-based platform for single-cell RNA sequencing analysis, covering individual-sample processing (Module I), multi-sample integration (Module II), post-annotation exploration (Module III), and cell-cell communication (Module IV).

## Project Structure

- `R/` — Shared code used by all modules
  - `utils.R` — Package management (`load.or.install`), categorical-meta helper, mouse/human ortholog conversion helpers used by Module IV
  - `mod_save_config.R` — Reusable Shiny sub-module for configurable Seurat object export
- `module-I/` — Shiny app for individual-sample scRNA-seq analysis
  - `modules/` — Shiny module files (one per pipeline step)
  - `app.R` — Entry point: loads packages, defines UI layout, wires module servers
- `module-II/` — Shiny app for multi-sample integration and annotation
  - `modules/` — Shiny module files (one per pipeline step)
  - `app.R` — Entry point: loads packages, defines UI layout, wires module servers
- `module-III/` — Shiny app for post-annotation exploration and visualization
  - `modules/` — Shiny module files (one per exploration feature)
  - `app.R` — Entry point: loads packages, defines UI layout, wires module servers
- `module-IV/` — Shiny app for cell-cell communication analysis (CellChat + LIANA)
  - `modules/` — Shiny module files (upload, CellChat/LIANA compute + vis, shared vis helpers)
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
- All tabs are available at all times; each step's main compute action warns (with a proceed-anyway modal and a button hint) when its upstream step's `completed` signal is not yet set, rather than hard-blocking navigation
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

## Module III: Post-Annotation Exploration

A modular Shiny web application for exploring and visualizing annotated scRNA-seq data:

**Upload → UMAP Cell-Type → Violin Plot → Dot Plot → Heatmap → Differential Expression → Enrichment Analysis**

### Active Modules (in `module-III/modules/`)

| Module | File | Purpose |
|--------|------|---------|
| Upload | `mod_explore_upload.R` | Load slimmed-down Seurat .rds objects, validate structure, report summary |
| UMAP Cell-Type | `mod_explore_umap.R` | Interactive lasso selection with re-UMAP, cell type highlighting, gene (co-)expression |
| Violin Plot | `mod_explore_violin.R` | Gene expression by cell type with stacking, split.by, and figure download |
| Dot Plot | `mod_explore_dot.R` | Gene expression dot plots with per-identity color scaling when split.by is used |
| Heatmap | `mod_explore_heatmap.R` | Scaled expression heatmaps with HVG or custom genes, vars.to.regress support |
| DE | `mod_explore_de.R` | Differential expression between configurable cell populations with metadata subsetting |
| Enrichment | `mod_explore_enrich.R` | enrichR-based pathway analysis on DE results with bar plot visualization |

### Architecture

- Each module is a Shiny module (namespaced UI + server pair)
- Upload returns reactive Seurat object; all visualization modules consume it via `shared.data`
- DE module returns its result reactive, which feeds into the Enrichment module
- Input: Annotated Seurat .rds objects (output of Module II), slimmed down (data slot only)
- All tabs are available once data is uploaded (flat navigation, not a linear pipeline)

## Module IV: Cell-Cell Communication

A modular Shiny web application offering two parallel cell-cell communication engines — CellChat and LIANA — over a shared uploaded object:

**Upload → Run CellChat → Visualize CellChat → Run LIANA → Visualize LIANA**

### Active Modules (in `module-IV/modules/`)

| Module | File | Purpose |
|--------|------|---------|
| Upload | `mod_interact_upload.R` | Load annotated Seurat .rds, JoinLayers, validate data slot, assign user-picked metadata column as `group` |
| Vis helpers | `mod_interact_vis_helpers.R` | Shared UI + render helpers — `plot.grid`, `resolve.sel`, `cell.type.selector.*`, `common.controls.ui`, `make.render.and.download`, typed placeholders. Sourced before both vis modules |
| CellChat compute | `mod_interact_cellchat_comp.R` | Per-group CellChat pipeline on Secreted Signaling DB, followed by `mergeCellChat` for cross-group visualization |
| CellChat vis | `mod_interact_cellchat_vis.R` | Four subtabs: Global Network, Zoom-in (pathway / L-R pair), Signaling-Focused, Communication Patterns (incl. Manifold via `computeNetSimilarityPairwise` + `netEmbedding`) |
| LIANA compute | `mod_interact_liana_comp.R` | Per-group multi-method LIANA pipeline (`natmi`, `connectome`, `logfc`, `sca`, `cellphonedb`, `cytotalk`, `call_cellchat`) with lazy mouse ↔ human ortholog mapping when the chosen resource is not `MouseConsensus` |
| LIANA vis | `mod_interact_liana_vis.R` | Three subtabs: CCC Dot Plot (single-group), CCC Freq Heatmap (N-group grid via `liana::heat_freq`), CCC Freq Chord Diagram (N-group grid via `liana::chord_freq`) |

### Architecture

- Each module is a Shiny module (namespaced UI + server pair)
- A single upload is shared by all four downstream tabs via `shared.data`
- Compute modules return a reactive **result list** (not a Seurat object) holding per-group results + metadata, plus a `.rds` download
- Vis modules consume the upstream compute result first, falling back to a user-uploaded `.rds` for standalone use
- Every vis panel is gated behind a "Generate Plot" button (`shiny::bindEvent`), so parameter tweaks queue up without triggering expensive re-renders
- `plot.grid` routes heterogeneous panel items through the right render path: ggplots pass straight through, ComplexHeatmap via `grid::grid.grabExpr`, base graphics (CellChat circle/chord/hierarchy, LIANA chord) via `ggplotify::as.ggplot` so CellChat's internal `par(mfrow)` layout doesn't collide with the outer grid
- Input: annotated Seurat .rds files (output of Module II); the upload also tolerates Module III's slimmed `data`-only export
- Output: CellChat result list (`.rds`) or LIANA result list (`.rds`)
- Ensembl access (for ortholog conversion) uses `.with.ensembl.hosts` in `R/utils.R`, which tries `dec2021.archive` → `www.ensembl.org` → `useast` → `asia` end-to-end under a single tryCatch and only raises if every host fails

## Tech Stack

- **Language:** R
- **Framework:** Shiny
- **Core packages:** Seurat, DoubletFinder, ggplot2, plotly, tidyverse, glmGamPoi, shinyFiles, shinyjs, ggpubr, shinycssloaders
- **Module II additional packages:** SeuratWrappers, SeuratObject, presto, cluster, lisi, Matrix, batchelor (FastMNN backend)
- **Module III additional packages:** enrichR, openxlsx
- **Module IV additional packages:** CellChat (jinworks/CellChat), ComplexHeatmap, NMF, ggalluvial, ggplotify, future, liana (saezlab/liana), OmnipathR, biomaRt, entropy, digest

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

## Modules

- **Module I** — Individual-sample analysis (active)
- **Module II** — Multi-sample integration and annotation (active)
- **Module III** — Post-annotation exploration and visualization (active)
- **Module IV** — Cell-cell communication via CellChat + LIANA (active)

The platform progresses from individual-sample processing through integrated multi-sample analysis to downstream exploration and cell-cell communication.

## Branching Workflow

`main` is the single production branch — application code, Docker files (`docker/`, `docker-compose.yml`, `.dockerignore`), and docs all live there together. Work lands on `main` directly or via a short-lived feature branch off it.

The platform supports two install paths, both documented in `docs/`: Docker (`docs/docker.md`, the recommended path) and native R (`docs/native.md`, the fallback when a step needs more memory than the container has). A change that affects setup or dependencies needs both guides checked.

