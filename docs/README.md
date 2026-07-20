# scNexus — User Guide

A step-by-step walkthrough of the **happy path** through the platform: the
straight-line route from a raw Cell Ranger matrix to cell-cell communication
results, using the bundled sample data. Each module has its own guide with
screenshots of every step.

This guide documents the intended, working path — not every option or edge case.
Parameters that don't need changing for the demo are left at their defaults.

## The pipeline

The four modules chain together — each one's output is the next one's input:

```
Module I            Module II              Module III            Module IV
Individual sample   Multi-sample           Post-annotation       Cell-cell
processing          integration            exploration           communication
─────────────       ─────────────          ─────────────         ─────────────
Cell Ranger MEX  →  Processed .rds      →  Annotated .rds     →  Annotated .rds
      │             (one per sample)        │                     │
      ▼                    │                ▼                     ▼
Processed .rds  ──────────►│           Figures & tables      CellChat / LIANA
                    Annotated .rds                            result objects
```

You can also enter at any module if you already have data in the right shape —
each module's guide lists what it expects as input.

## Module guides

| Guide | App | Pipeline |
|-------|-----|----------|
| [Module I](module-I.md) | scNexus-Process | Import → QC → PCA → Doublet → Cell Cycle → Annotation |
| [Module II](module-II.md) | scNexus-Integrate | Upload & Merge → Integration → Benchmarking → Annotation |
| [Module III](module-III.md) | scNexus-Explore | Upload → UMAP → Violin → Dot → Heatmap → DE → Enrichment |
| [Module IV](module-IV.md) | scNexus-Interact | Upload → CellChat → LIANA |

## Data preparation

The demo datasets are not in this repository — download them before starting:

**[Download the demo datasets](https://drive.google.com/drive/folders/1ufBv0MfgJGrSATmonozJeY0CLJekHYnb)**

Put every file in the **`datasets/`** folder at the root of the repository. That
folder is where the platform reads and writes data; under Docker it is the only
host folder the apps can see.

| File | Used by |
|------|---------|
| `Demo_Module_I_21401X3.rds` | **Module I** — load with *Option 2 — Processed Seurat object (.rds)* on the Import tab |
| `Demo_Module_II_21401X3_control.rds` | **Module II** — one of three samples to merge and integrate |
| `Demo_Module_II_22713X2_disease.rds` | **Module II** |
| `Demo_Module_II_24143X4_control.rds` | **Module II** |

**Modules III and IV have no separate demo file** — both take the annotated
object that Module II produces. Run Module II first and save its output to
`datasets/`, then load that object in Module III or IV.

## Launching the platform

The platform runs as a single Docker image — see the [Docker setup
guide](docker.md) for build, run, and sample-data instructions. Once it's up,
each module opens in a browser (Module I at <http://localhost:3839>, and so on).

The Shiny apps can also be run directly from RStudio with a local R installation
— see the [RStudio setup guide](shiny.md). This is the fallback when a step needs
more memory than the container has; notably, Module IV's LIANA analysis may fail
under Docker.

## About the screenshots

Screenshots live in `docs/img/<module>/` and are named by step
(e.g. `01-import.png`). They are captured against the demo datasets, so anything
you see here is reproducible with them.
