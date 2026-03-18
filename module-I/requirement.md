# Requirements for scRNA-seq Shiny App

## Platform
- macOS (tested on Apple Silicon)
- Homebrew
- R 4.5.x

## R Package Dependencies
The app dynamically loads or installs these R packages in app.R:

- shiny
- shinyjs
- Seurat
- ggplot2
- tidyverse
- shinyFiles
- DoubletFinder
- plotly
- glmGamPoi
- ggpubr
- BiocManager

Install command (optional pre-install):

```bash
env PATH="<TOOLCHAIN_PATH>" \
R_LD_LIBRARY_PATH="<R_LD_LIBRARY_PATH>" \
R_MAKEVARS_USER="$PWD/.vscode/Makevars" \
R -q -e "options(repos=c(CRAN='https://cloud.r-project.org')); install.packages(c('shiny','shinyjs','Seurat','ggplot2','tidyverse','shinyFiles','DoubletFinder','plotly','glmGamPoi','ggpubr','BiocManager'))"
```

## System Dependencies (Homebrew)
These are required for compiled R package dependencies used by tidyverse/Seurat stack:

- fribidi
- harfbuzz
- freetype
- libpng
- libtiff
- jpeg-turbo
- webp
- zlib
- gcc
- pkgconf

Install command:

```bash
brew install fribidi harfbuzz freetype libpng libtiff jpeg-turbo webp zlib gcc pkgconf
```

## Workspace Build/Run Environment
The one-click VS Code run configuration expects:

- PATH=<TOOLCHAIN_PATH>
- R_LD_LIBRARY_PATH=<R_LD_LIBRARY_PATH>
- R_MAKEVARS_USER=${workspaceFolder}/.vscode/Makevars

And .vscode/Makevars contains FLIBS for Homebrew gcc runtime libs.

Placeholder examples:

- <TOOLCHAIN_PATH> = /opt/homebrew/bin:/opt/homebrew/sbin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin
- <R_LD_LIBRARY_PATH> = <MINIFORGE_PREFIX>/lib:<HOMEBREW_PREFIX>/lib/R/lib

## Run App
From workspace root:

```bash
env PATH="<TOOLCHAIN_PATH>" \
R_LD_LIBRARY_PATH="<R_LD_LIBRARY_PATH>" \
R_MAKEVARS_USER="$PWD/.vscode/Makevars" \
R -e "shiny::runApp('.', launch.browser = TRUE)"
```

Or use VS Code Run and Debug profile: Run Shiny App.
