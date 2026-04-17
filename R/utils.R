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

# --- Mouse / human ortholog conversion (for LIANA non-MouseConsensus resources) ---
# biomaRt is namespace-qualified so sourcing utils.R does not require biomaRt.

.ensembl.mart <- function(dataset) {
  # dec2021 archive is the PI-pinned endpoint for reproducibility but it
  # occasionally fails with "Failed to perform HTTP request". Fall through
  # to the live site and the regional mirrors before giving up.
  hosts <- c("https://dec2021.archive.ensembl.org/",
             "https://www.ensembl.org/",
             "https://useast.ensembl.org/",
             "https://asia.ensembl.org/")
  last.err <- NULL
  for (host in hosts) {
    mart <- tryCatch(
      biomaRt::useMart("ensembl", dataset = dataset, host = host),
      error = function(e) { last.err <<- e; NULL }
    )
    if (!is.null(mart)) return(mart)
  }
  stop("All Ensembl hosts failed for dataset '", dataset,
       "'. Last error: ", conditionMessage(last.err))
}

# Deduplicate a getLDS result keeping one-to-one orthologs only
.dedup.ortho <- function(df, from, to) {
  df <- df[df[[to]] != "", ]
  df <- df[!duplicated(df), ]
  df <- df[!duplicated(df[[from]]), ]
  df <- df[!duplicated(df[[to]]), ]
  df
}

# Rewrite rownames of a Seurat object from mouse MGI symbols to human HGNC symbols.
# cache.env: optional environment() to avoid re-querying Ensembl for the same gene set.
# Returns list(obj = <converted Seurat>, map = <data.frame MGI.symbol / HGNC.symbol>).
convert.mouse.to.human.rownames <- function(obj, cache.env = NULL) {
  mouse.genes <- rownames(obj)
  cache.key <- paste0("m2h:", digest::digest(sort(mouse.genes)))

  gene.map <- if (!is.null(cache.env) && exists(cache.key, envir = cache.env)) {
    get(cache.key, envir = cache.env)
  } else {
    mouse <- .ensembl.mart("mmusculus_gene_ensembl")
    human <- .ensembl.mart("hsapiens_gene_ensembl")
    m2h <- biomaRt::getLDS(
      attributes = c("mgi_symbol"),
      filters = "mgi_symbol",
      values = mouse.genes,
      mart = mouse,
      attributesL = c("hgnc_symbol"),
      martL = human,
      uniqueRows = TRUE
    )
    m2h <- .dedup.ortho(m2h, "MGI.symbol", "HGNC.symbol")
    if (!is.null(cache.env)) assign(cache.key, m2h, envir = cache.env)
    m2h
  }

  aligned <- gene.map[match(rownames(obj), gene.map$MGI.symbol), ]
  keep <- !is.na(aligned$HGNC.symbol)
  obj <- obj[keep, ]
  rownames(obj) <- aligned$HGNC.symbol[keep]
  list(obj = obj, map = aligned[keep, ])
}

# Rewrite human HGNC gene symbols back to mouse MGI in a LIANA result tibble.
# Handles single-gene columns (ligand, receptor) and complex columns whose
# subunits are joined by "_" (ligand.complex, receptor.complex).
# Genes without a mouse ortholog fall back to their original human symbol.
convert.human.to.mouse.lr <- function(df.human, cache.env = NULL) {
  complex.cols <- c("ligand.complex", "receptor.complex")
  gene.cols <- c("ligand", "receptor")
  complex.cols <- intersect(complex.cols, colnames(df.human))
  gene.cols <- intersect(gene.cols, colnames(df.human))

  complex.parts <- unlist(lapply(complex.cols, function(cc) {
    strsplit(unique(df.human[[cc]]), "_")
  }))
  single.genes <- unlist(lapply(gene.cols, function(gc) df.human[[gc]]))
  human.genes <- unique(c(single.genes, complex.parts))
  human.genes <- human.genes[!is.na(human.genes) & nzchar(human.genes)]

  cache.key <- paste0("h2m:", digest::digest(sort(human.genes)))
  gene.map <- if (!is.null(cache.env) && exists(cache.key, envir = cache.env)) {
    get(cache.key, envir = cache.env)
  } else {
    mouse <- .ensembl.mart("mmusculus_gene_ensembl")
    human <- .ensembl.mart("hsapiens_gene_ensembl")
    h2m <- biomaRt::getLDS(
      attributes = c("hgnc_symbol"),
      filters = "hgnc_symbol",
      values = human.genes,
      mart = human,
      attributesL = c("mgi_symbol"),
      martL = mouse,
      uniqueRows = TRUE
    )
    h2m <- .dedup.ortho(h2m, "HGNC.symbol", "MGI.symbol")
    if (!is.null(cache.env)) assign(cache.key, h2m, envir = cache.env)
    h2m
  }

  lookup <- setNames(gene.map$MGI.symbol, gene.map$HGNC.symbol)
  map.gene <- function(g) {
    m <- lookup[g]
    ifelse(is.na(m), g, unname(m))
  }
  map.complex <- function(cx) {
    if (is.na(cx) || !nzchar(cx)) return(cx)
    paste(map.gene(strsplit(cx, "_")[[1]]), collapse = "_")
  }

  for (gc in gene.cols) df.human[[gc]] <- map.gene(df.human[[gc]])
  for (cc in complex.cols) {
    df.human[[cc]] <- vapply(df.human[[cc]], map.complex, character(1))
  }
  df.human
}
