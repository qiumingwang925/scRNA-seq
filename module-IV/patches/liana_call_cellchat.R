## ABOUTME: Runtime patch for liana::call_cellchat, which calls GetAssayData(slot = "data") --
## ABOUTME: defunct under SeuratObject >= 5, so `call_cellchat` errors out for every group.
## ABOUTME: patch.liana.call.cellchat() swaps our fixed copy into liana's namespace at startup.

# Against saezlab/liana@6cab46c54234f861ea176c3de77c4b8aa45ecb3d (R/cellchat_pipe.R), this
# copy changes exactly three things:
#   1. GetAssayData(slot = "data") -> layer = "data"   (the SeuratObject >= 5 fix)
#   2. future::plan("multiprocess") -> multisession    (defunct in future >= 1.32)
#   3. stop() on an unknown organism, instead of falling through with ccDB undefined
#
# Formal names below are an external contract: liana_wrap builds the argument list from
# liana_defaults()[["call_cellchat"]] and passes it by name, so none of them may be
# renamed to the project's dot-separated convention.
patched.call.cellchat <- function(sce,
                                  op_resource,
                                  .format = TRUE,
                                  exclude_anns = c(),
                                  nboot = 100,
                                  assay = "RNA",
                                  .seed = 1004,
                                  .normalize = FALSE,
                                  .do_parallel = FALSE,
                                  .raw_use = TRUE,
                                  expr_prop = 0,
                                  organism = "human",
                                  thresh = 1,
                                  de_thresh = 0.05,
                                  ...) {

  stringsAsFactors.old <- options("stringsAsFactors")[[1]]
  options(stringsAsFactors = FALSE)
  on.exit(options(stringsAsFactors = stringsAsFactors.old), add = TRUE)

  if (inherits(sce, "SingleCellExperiment")) {
    sce <- liana:::.liana_convert(sce, assay = assay)
  }

  labels <- Seurat::Idents(sce)
  meta <- data.frame(group = labels, row.names = names(labels))

  expr.data <- SeuratObject::GetAssayData(sce, assay = assay, layer = "data")
  if (.normalize) {
    expr.data <- CellChat::normalizeData(expr.data)
  }

  cellchat.omni <- CellChat::createCellChat(object = expr.data,
                                            meta = meta,
                                            group.by = "group")
  cellchat.omni <- CellChat::addMeta(cellchat.omni, meta = meta)
  cellchat.omni <- CellChat::setIdent(cellchat.omni, ident.use = "group")

  if (.do_parallel) {
    future::plan(future::multisession)
  }

  ccDB <- if (organism == "human") {
    CellChat::CellChatDB.human
  } else if (organism == "mouse") {
    CellChat::CellChatDB.mouse
  } else {
    stop("organism must be 'human' or 'mouse', got '", organism, "'")
  }

  if (!is.null(op_resource)) {
    ccDB <- liana::cellchat_formatDB(ccDB, op_resource, exclude_anns)
  } else {
    ccDB$interaction <- dplyr::filter(ccDB$interaction,
                                      !(annotation %in% exclude_anns))
  }
  cellchat.omni@DB <- ccDB

  cellchat.omni <- CellChat::subsetData(cellchat.omni)
  cellchat.omni <- CellChat::identifyOverExpressedGenes(cellchat.omni,
                                                        thresh.pc = expr_prop,
                                                        thresh.p = de_thresh)
  cellchat.omni <- CellChat::identifyOverExpressedInteractions(cellchat.omni)

  if (!.raw_use) {
    cellchat.omni <- CellChat::projectData(cellchat.omni, CellChat::PPI.human)
  }

  cellchat.omni <- CellChat::computeCommunProb(cellchat.omni,
                                               raw.use = .raw_use,
                                               seed.use = .seed,
                                               nboot = nboot)
  cellchat.omni <- CellChat::filterCommunication(cellchat.omni, min.cells = 1)

  df.omni <- CellChat::subsetCommunication(cellchat.omni, thresh = thresh, ...)

  if (.format) {
    df.omni <- tibble::as_tibble(
      dplyr::select(df.omni, source, target, ligand, receptor, prob, pval)
    )
  }

  df.omni
}

# liana_wrap dispatches through .select_method(), which yields the bare symbol
# `call_cellchat` for rlang::invoke to resolve inside liana's namespace at call time --
# so rebinding the name there is enough, no fork of the package required.
#
# assign() rebinds a name but does not change a function's closure environment, so
# without the environment() call below patched.call.cellchat would resolve its free
# variables in whatever environment sourced this file rather than in liana's namespace.
patch.liana.call.cellchat <- function() {
  tryCatch({
    patched <- patched.call.cellchat
    environment(patched) <- asNamespace("liana")
    unlockBinding("call_cellchat", asNamespace("liana"))
    assign("call_cellchat", patched, envir = asNamespace("liana"))
    lockBinding("call_cellchat", asNamespace("liana"))
    message("Patched liana::call_cellchat for SeuratObject >= 5.")
  }, error = function(e) {
    warning("Could not patch liana::call_cellchat: ", conditionMessage(e),
            "\nThe 'call_cellchat' method will fail; the other methods are unaffected.",
            call. = FALSE)
  })
  invisible(NULL)
}
