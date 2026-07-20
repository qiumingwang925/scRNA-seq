# Setup — Running the Platform Natively

Docker is the recommended way to run the platform ([setup guide](docker.md)) —
it needs no R installation and no package builds. Run natively when you need a
module to use your machine's full memory: **Module IV's LIANA analysis can run
out of memory inside the container** and complete natively on the same machine.

## Prerequisites

- **R 4.4.2** — the version the Docker image is built against, so it is the
  version this stack is known to work on.
- **System libraries** for the single-cell stack (compilers, HDF5, libcurl,
  libxml2, and similar). These vary by operating system; on Linux they are the
  usual build dependencies for Seurat and Bioconductor packages.

## Install

There is no separate install step. Each module's `app.R` declares its own
packages and installs any that are missing on first launch (see
`load.or.install` in `R/utils.R`).

The **first launch is slow** — it compiles the whole Seurat + Bioconductor
stack, which can take considerably longer than a Docker build. Later launches
start immediately.

## Run

From the **repository root**, launch the module you want:

```bash
R -e 'shiny::runApp("module-I",   port = 3839)'
R -e 'shiny::runApp("module-II",  port = 3840)'
R -e 'shiny::runApp("module-III", port = 3841)'
R -e 'shiny::runApp("module-IV",  port = 3842)'
```

Each app must run with its own module folder as the working directory —
`shiny::runApp()` handles this, so launch by folder name rather than sourcing
`app.R` directly.

Then open the printed URL in a browser and follow the
[workflow guides](README.md).

## Sample data

Native runs read `test-data/` from the repository directly, so no extraction
step is needed. Module I's file browser is unrestricted here, unlike in the
container.
