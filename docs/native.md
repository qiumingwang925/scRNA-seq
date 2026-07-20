# Setup — Running the Platform Natively

Docker is the recommended way to run the platform ([setup guide](docker.md)) —
it needs no R installation and no package builds. Run natively when you need a
module to use your machine's full memory: **Module IV's LIANA analysis can run
out of memory inside the container** and complete natively on the same machine.

## Prerequisites

- **R 4.5.x**
- **RStudio** — recommended; the steps below use it.
- **System libraries** for the compiled packages in the Seurat and tidyverse
  stack.

On macOS, install them with [Homebrew](https://brew.sh):

```bash
brew install fribidi harfbuzz freetype libpng libtiff jpeg-turbo webp zlib gcc pkgconf
```

On Linux, install the equivalent development packages through your
distribution's package manager, along with a C/C++/Fortran toolchain.

## R packages

There is no separate install step. Each module's `app.R` declares its own
packages and installs any that are missing on first launch (see
`load.or.install` in `R/utils.R`).

The **first launch is slow** — it compiles the whole Seurat and Bioconductor
stack, which can take considerably longer than a Docker build. Later launches
start immediately.

## Run with RStudio

1. Open RStudio and set the working directory to the repository folder
   (*Session → Set Working Directory → Choose Directory*).
2. Open the `app.R` of the module you want — `module-I/app.R`, `module-II/app.R`,
   `module-III/app.R`, or `module-IV/app.R`.
3. Click **Run App** in the editor toolbar. RStudio runs the app from its own
   module folder, which is what the module's internal paths expect.
4. The app opens in a browser (or RStudio's viewer — use *Open in Browser* for a
   full window).

Run one module at a time. To switch modules, stop the running app first.

## Run from the command line

From the repository root:

```bash
R -e 'shiny::runApp("module-I",   port = 3839)'
R -e 'shiny::runApp("module-II",  port = 3840)'
R -e 'shiny::runApp("module-III", port = 3841)'
R -e 'shiny::runApp("module-IV",  port = 3842)'
```

Launch by folder name rather than sourcing `app.R` directly — `shiny::runApp()`
sets the working directory the module's internal paths depend on.

## There is no landing page

The Docker setup serves a landing page at <http://localhost:3838> that links to
all four modules. **That page does not exist when running natively.** Open the
URL the app itself prints on startup — RStudio assigns a port automatically, and
the command-line form above uses whichever port you pass it. Do not browse to
3838.

## Sample data

Native runs read `test-data/` from the repository directly, so no extraction
step is needed. Module I's file browser is unrestricted here, unlike in the
container.

## Next steps

With a module running, follow the [workflow guides](README.md).
