## ABOUTME: Resumable installer for all four scRNA-seq modules' package stack.
## ABOUTME: Skips already-installed packages and continues past failures, so a re-run
## ABOUTME: against a BuildKit cache-mounted library only installs what is still missing.

options(Ncpus = 4L, warn = 1)

# Target library. The Dockerfile passes a BuildKit cache-mounted dir (persisted across
# builds, even failed ones) so a re-run resumes instead of reinstalling everything.
# Falls back to the default library when run without an argument.
args <- commandArgs(trailingOnly = TRUE)
lib <- if (length(args) >= 1L && nzchar(args[[1L]])) args[[1L]] else .libPaths()[1L]
dir.create(lib, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(lib, .libPaths()))

have <- function(pkg) requireNamespace(pkg, quietly = TRUE)
need <- function(pkgs) pkgs[!vapply(pkgs, have, logical(1L))]

# Bootstrap BiocManager, then expose CRAN (the frozen P3M snapshot the rocker base pins)
# AND Bioconductor 3.20 for the whole run: some CRAN packages (e.g. NMF -> Biobase)
# depend on Bioc, so Bioc must be in the repo list before ANY install. CRAN stays the
# P3M snapshot; BiocManager::repositories() only appends the Bioc repos.
if (!have("BiocManager")) install.packages("BiocManager", lib = lib)
options(repos = BiocManager::repositories())

cran.pkgs <- c(
  "shiny", "shinyjs", "shinycssloaders", "shinyFiles",
  "ggplot2", "tidyverse", "plotly", "DT", "patchwork", "scales",
  "ggnewscale", "ggpubr", "ggplotify", "ggalluvial", "future",
  "NMF", "cluster", "Matrix", "digest", "entropy", "harmony",
  "openxlsx", "enrichR", "remotes",
  "Seurat", "SeuratObject"
)

bioc.pkgs <- c(
  "glmGamPoi", "SingleR", "SingleCellExperiment",
  "ComplexHeatmap", "biomaRt", "OmnipathR",
  "batchelor"  # FastMNN backend for SeuratWrappers (Module II)
)

# GitHub deps pinned to the exact commits from a known-good build (BUILD_VERSIONS.txt,
# 2026-06-19). Default-branch HEAD drifts: an un-pinned build later pulls newer commits
# of the fast-moving packages (CellChat, liana) whose dependencies may no longer resolve
# against the frozen CRAN snapshot, breaking the build. Keyed by package namespace so an
# already-installed one can be skipped. CRAN (p3m.dev/cran/2025-02-27) and Bioc 3.20 (a
# closed release serving final versions) are already frozen, so only these need pinning.
github.pkgs <- c(
  SeuratWrappers = "satijalab/seurat-wrappers@ffaf74e306279b1ec16e31c9cb2142ebb2bc4bc1",
  DoubletFinder  = "chris-mcginnis-ucsf/DoubletFinder@1b244d8f0d54b4b1cb4365639931bbb16f01e1cd",
  presto         = "immunogenomics/presto@a24772a135c7895a8183b007376050556c60a05b",
  lisi           = "immunogenomics/lisi@a917556310d8d2c66833dcc35aa3d0f4d1b6e0f4",
  CellChat       = "jinworks/CellChat@75253cd0c9e68410e6e721a6d3a0419a1d7e358f",
  liana          = "saezlab/liana@6cab46c54234f861ea176c3de77c4b8aa45ecb3d"
)

# CRAN: install only what's missing. install.packages() warns (not errors) on a failed
# package and keeps going, so one bad build won't abort the rest.
cran.todo <- need(cran.pkgs)
if (length(cran.todo)) {
  message("[cran] installing: ", paste(cran.todo, collapse = ", "))
  tryCatch(install.packages(cran.todo, lib = lib),
           error = function(e) message("[cran] stage error: ", conditionMessage(e)))
} else message("[cran] all present")

# Bioconductor: same, restricted to missing packages.
bioc.todo <- need(bioc.pkgs)
if (length(bioc.todo)) {
  message("[bioc] installing: ", paste(bioc.todo, collapse = ", "))
  tryCatch(BiocManager::install(bioc.todo, lib = lib, update = FALSE, ask = FALSE),
           error = function(e) message("[bioc] stage error: ", conditionMessage(e)))
} else message("[bioc] all present")

# GitHub: per-package so one failure (install_github errors, unlike install.packages)
# doesn't abort the others; skip any already installed. No explicit lib: install_github
# can reject it; the prepended .libPaths() already makes /rlib the default location.
install.github.pinned <- function(ref) {
  # Try the API-based installer first (it resolves any Remotes: field), then fall back to
  # the commit archive tarball. The archive downloads from codeload.github.com, NOT the
  # rate-limited api.github.com. dependencies = FALSE is essential: without it remotes
  # still resolves the package's DESCRIPTION Remotes: field, which fires api.github.com
  # calls and re-triggers the 60 req/hr/IP rate limit even though the tarball itself came
  # from codeload. It is safe here because every dependency these packages need is already
  # installed (CRAN/Bioc, or an earlier entry in this ordered list).
  ok <- tryCatch({ remotes::install_github(ref, upgrade = "never"); TRUE },
                 error = function(e) { message("  install_github failed: ", conditionMessage(e)); FALSE })
  if (ok) return(TRUE)
  parts <- strsplit(ref, "@", fixed = TRUE)[[1L]]
  url <- sprintf("https://github.com/%s/archive/%s.tar.gz", parts[[1L]], parts[[2L]])
  for (attempt in seq_len(3L)) {
    message("  archive fallback ", attempt, "/3: ", url)
    ok <- tryCatch({ remotes::install_url(url, dependencies = FALSE, upgrade = "never"); TRUE },
                   error = function(e) { message("  install_url failed: ", conditionMessage(e)); FALSE })
    if (ok) return(TRUE)
    Sys.sleep(10)
  }
  FALSE
}
for (nm in names(github.pkgs)) {
  if (have(nm)) { message("[github] skip (installed): ", nm); next }
  message("[github] installing: ", github.pkgs[[nm]])
  if (!install.github.pinned(github.pkgs[[nm]])) message("[github] FAILED ", nm)
}

# Gate: everything was attempted above; now fail if anything is still missing so a broken
# image is never produced. The cache-mounted library keeps what succeeded, so re-running
# the build resumes and skips the already-installed packages.
required <- c(cran.pkgs, bioc.pkgs, names(github.pkgs))
missing <- need(required)
if (length(missing)) {
  stop("Still missing after this pass: ", paste(missing, collapse = ", "),
       "\nRe-run the build to resume — installed packages are cached and skipped.")
}

# Record exactly what got built (versions + GitHub commit SHAs) so these can be pinned
# later without guessing. Lives at /srv/BUILD_VERSIONS.txt in the image.
record <- vapply(required, function(pkg) {
  desc <- packageDescription(pkg)
  sha <- desc$RemoteSha
  sprintf("%-22s %-14s %s", pkg, desc$Version,
          if (!is.null(sha)) paste0("@", sha) else "")
}, character(1L))
writeLines(c(paste("Built:", format(Sys.time())), "", record),
           "/srv/BUILD_VERSIONS.txt")
