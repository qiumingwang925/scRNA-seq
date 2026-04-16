## ABOUTME: Shared utility functions used across all modules in the scRNA-seq platform.
## ABOUTME: Provides package management (install-if-missing) and common setup.

cran.repo <- unname(getOption("repos")["CRAN"])
if (is.null(cran.repo) || length(cran.repo) == 0 || is.na(cran.repo) || cran.repo == "@CRAN@") {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

load.or.install <- function(pkg, github.url = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!is.null(github.url)) {
      if (!requireNamespace("remotes", quietly = TRUE))
        install.packages("remotes")
      remotes::install_github(github.url)
    } else {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install(pkg)
    }
  }
  library(pkg, character.only = TRUE)
}

# Returns non-numeric (categorical) metadata column names from a Seurat object
get.categorical.meta <- function(obj) {
  meta <- obj@meta.data
  cat.cols <- sapply(meta, function(x) !is.numeric(x))
  names(cat.cols)[cat.cols]
}
