# Setup — Running the Shiny App from RStudio or a Headless Server

Docker is the recommended way to run the platform ([setup guide](docker.md)) —
it needs no R installation and no package builds. Run the Shiny apps directly
from RStudio when you need a module to use your machine's full memory:
**Module IV's LIANA analysis can run out of memory inside the container** and
complete from RStudio on the same machine.

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

These bind to `127.0.0.1`, reachable only from the machine running them. To
reach an app on a **remote headless server**, bind to all interfaces:

```bash
R -e 'shiny::runApp("module-I", host = "0.0.0.0", port = 3839)'
```

The apps have no authentication, so only do this on a trusted network. The safer
alternative is to leave the default binding and forward the port over SSH:

```bash
ssh -L 3839:localhost:3839 user@server
```

## There is no landing page

The Docker setup serves a landing page at <http://localhost:3838> that links to
all four modules. **That page does not exist outside Docker.** Open the URL the
app itself prints on startup — RStudio assigns a port automatically, and the
command-line form above uses whichever port you pass it. Do not browse to 3838.

## Sample data

Download the demo datasets into `datasets/` — see [data
preparation](README.md#data-preparation). Module I's file browser is
unrestricted here, unlike in the container, so it can reach any path on the
machine.

## Next steps

With a module running, follow the [workflow guides](README.md).
