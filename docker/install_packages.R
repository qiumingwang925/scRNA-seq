## ABOUTME: Installs the full package stack for all four scRNA-seq modules into the image.
## ABOUTME: CRAN (P3M binaries) -> Bioconductor -> GitHub; fails the build if anything is missing.

options(Ncpus = 4L, warn = 1)

# rocker's Rprofile.site points the CRAN repo at Posit Package Manager
# (precompiled binaries for this image's Ubuntu codename). Bootstrap BiocManager,
# then make CRAN (P3M binaries) AND Bioconductor visible for the whole run:
# several CRAN packages (e.g. NMF -> Biobase) depend on Bioconductor packages,
# so Bioc must be in the repo list before ANY install or those deps fail to
# resolve. BiocManager::repositories() keeps the existing P3M CRAN repo and only
# appends the Bioc repos, so CRAN binaries are still used.
install.packages("BiocManager")
options(repos = BiocManager::repositories())

cran.pkgs <- c(
  "shiny", "shinyjs", "shinycssloaders", "shinyFiles",
  "ggplot2", "tidyverse", "plotly", "DT", "patchwork", "scales",
  "ggnewscale", "ggpubr", "ggplotify", "ggalluvial", "future",
  "NMF", "cluster", "Matrix", "digest", "entropy", "harmony",
  "openxlsx", "enrichR", "remotes",
  "Seurat", "SeuratObject"
)
install.packages(cran.pkgs)

bioc.pkgs <- c(
  "glmGamPoi", "SingleR", "SingleCellExperiment",
  "ComplexHeatmap", "biomaRt", "OmnipathR"
)
BiocManager::install(bioc.pkgs, update = FALSE, ask = FALSE)

github.pkgs <- c(
  "satijalab/seurat-wrappers",
  "chris-mcginnis-ucsf/DoubletFinder",
  "immunogenomics/presto",
  "immunogenomics/lisi",
  "jinworks/CellChat",
  "saezlab/liana"
)
remotes::install_github(github.pkgs, upgrade = "never")

required <- c(
  "shiny", "shinyjs", "shinycssloaders", "shinyFiles", "ggplot2",
  "tidyverse", "plotly", "DT", "patchwork", "scales", "ggnewscale",
  "ggpubr", "ggplotify", "ggalluvial", "future", "NMF", "cluster",
  "Matrix", "digest", "entropy", "harmony", "openxlsx", "enrichR",
  "Seurat", "SeuratObject", "glmGamPoi", "SingleR", "SingleCellExperiment",
  "ComplexHeatmap", "biomaRt", "OmnipathR", "SeuratWrappers",
  "DoubletFinder", "presto", "lisi", "CellChat", "liana"
)
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Package installation failed for: ", paste(missing, collapse = ", "))
}

# Record exactly what got built (versions + GitHub commit SHAs) so these can be
# pinned later without guessing. Lives at /srv/BUILD_VERSIONS.txt in the image.
record <- vapply(required, function(pkg) {
  desc <- packageDescription(pkg)
  sha <- desc$RemoteSha
  sprintf("%-22s %-14s %s", pkg, desc$Version,
          if (!is.null(sha)) paste0("@", sha) else "")
}, character(1))
writeLines(c(paste("Built:", format(Sys.time())), "", record),
           "/srv/BUILD_VERSIONS.txt")
