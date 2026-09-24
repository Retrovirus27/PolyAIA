# =============================================================================
#  Pseudobulk expression filter and gene-level differential expression (Libra)
# =============================================================================


# Internal: counts layer of an assay, joining Seurat v5 split layers if needed.
.PBCounts <- function(seu, assay) {
  a <- seu[[assay]]
  if (inherits(a, "Assay5")) a <- SeuratObject::JoinLayers(a)
  SeuratObject::LayerData(a, layer = "counts")
}


#' Pseudobulk counts and expression filter per cell type and comparison
#'
#' Aggregates raw counts into one pseudobulk sample per
#' (cell type, group, replicate) and applies an expression filter
#' (\code{edgeR::filterByExpr()} by default) to every cell type x comparison
#' block. Assay-agnostic: use the \code{RNA} assay for genes or the
#' \code{polyA} assay for peaks. The returned object is what
#' \code{\link{QCDensity}()} plots.
#'
#' Pseudobulk samples are identified by an internal id and described in
#' \code{sample_meta} built from \code{@meta.data}, rather than by parsing
#' \code{AggregateExpression()} column names (Seurat v5 rewrites \code{"_"}
#' to \code{"-"} and \code{make.names()} mangles names such as
#' \code{"L2/3 IT"}).
#'
#' @param seu A Seurat object.
#' @param comparisons \code{data.frame} with columns \code{treatment},
#'   \code{control} and \code{Condition}. \code{treatment}/\code{control} are
#'   values of \code{label_col} (e.g. \code{"B6_Alcohol"}, \code{"B6_Control"});
#'   \code{Condition} is the comparison label.
#' @param celltype_col Metadata column with the cell type (e.g.
#'   \code{"Subpopulation"}).
#' @param label_col Metadata column with the group compared (e.g.
#'   \code{"contrast"} = \code{paste0(Strain, "_", Treatment)}).
#' @param replicate_col Metadata column with the biological replicate (e.g.
#'   \code{"Number"}).
#' @param assay Assay to aggregate. \code{"RNA"} for genes, \code{"polyA"} for
#'   peaks. Default \code{"RNA"}.
#' @param features Optional subset of features (rows) to aggregate and test.
#' @param min_cells Minimum number of cells required in EACH group of a block
#'   for it to be tested. Default 10.
#' @param filter_method \code{"edgeR"} (default): \code{edgeR::filterByExpr()}
#'   with the block's design; \code{"manual"}: keep features with
#'   \code{rowSums(counts) >= mincounts}; \code{"none"}: keep all features.
#' @param mincounts Threshold for \code{filter_method = "manual"}. Default 10.
#' @param filterByExpr_args Named list of extra arguments for
#'   \code{edgeR::filterByExpr()} (e.g. \code{list(min.count = 5)}).
#' @param filter_scope Samples the filter is computed on.
#'   \code{"comparison"} (default): each cell type x comparison block on its
#'   own treatment + control samples. \code{"celltype"}: once per cell type on
#'   all its samples, and that feature set reused in every comparison (same
#'   universe across comparisons).
#' @param feature_type \code{"gene"} or \code{"peak"}; used in labels only.
#' @param verbose Print progress messages. Default \code{TRUE}.
#'
#' @return A list with
#'   \describe{
#'     \item{\code{pseudobulk}}{sparse features x pseudobulk-sample count matrix}
#'     \item{\code{sample_meta}}{\code{sample_id}, \code{Cells}, \code{label},
#'       \code{replicate}, \code{n_cells}}
#'     \item{\code{kept}}{features retained per block (\code{Cells},
#'       \code{Condition}, \code{feature})}
#'     \item{\code{summary}}{one row per block: cells and samples per group,
#'       features in / kept, status}
#'     \item{\code{comparisons}, \code{funnel}, \code{feature_type},
#'       \code{assay}, \code{filter_method}, \code{filter_scope}}{settings}
#'   }
#'
#' @examples
#' \dontrun{
#' seu$contrast <- paste0(seu$Strain, "_", seu$Treatment)
#' comparisons <- data.frame(
#'   treatment = c("B6_Alcohol", "3xTg_Control"),
#'   control   = c("B6_Control", "B6_Control"),
#'   Condition = c("B6_Alcohol_vs_B6_Control", "3xTg_Control_vs_B6_Control")
#' )
#' flt <- PseudobulkFilter(seu, comparisons,
#'                         celltype_col = "Subpopulation",
#'                         label_col = "contrast", replicate_col = "Number")
#' flt$summary
#' }
#'
#' @export
PseudobulkFilter <- function(seu,
                             comparisons,
                             celltype_col,
                             label_col,
                             replicate_col,
                             assay             = "RNA",
                             features          = NULL,
                             min_cells         = 10,
                             filter_method     = c("edgeR", "manual", "none"),
                             mincounts         = 10,
                             filterByExpr_args = list(),
                             filter_scope      = c("comparison", "celltype"),
                             feature_type      = c("gene", "peak"),
                             verbose           = TRUE) {
  feature_type  <- match.arg(feature_type)
  filter_method <- match.arg(filter_method)
  filter_scope  <- match.arg(filter_scope)
  if (filter_method == "manual" &&
      (!is.numeric(mincounts) || length(mincounts) != 1 || mincounts < 0)) {
    stop("`mincounts` must be a single non-negative number.")
  }

  req  <- c("treatment", "control", "Condition")
  miss <- setdiff(req, colnames(comparisons))
  if (length(miss)) {
    stop("`comparisons` is missing column(s): ", paste(miss, collapse = ", "))
  }
  comparisons <- as.data.frame(comparisons)[, req]
  comparisons[] <- lapply(comparisons, as.character)

  for (cl in c(celltype_col, label_col, replicate_col)) {
    if (!cl %in% colnames(seu@meta.data)) {
      stop("Column '", cl, "' not found in seu@meta.data.")
    }
  }
  if (!assay %in% SeuratObject::Assays(seu)) {
    stop("Assay '", assay, "' not found in seu.")
  }

  meta <- seu@meta.data
  used_labels <- unique(c(comparisons$treatment, comparisons$control))
  miss_lab <- setdiff(used_labels, unique(as.character(meta[[label_col]])))
  if (length(miss_lab)) {
    stop("Value(s) not found in '", label_col, "': ",
         paste(miss_lab, collapse = ", "),
         ". Present: ", paste(utils::head(unique(as.character(meta[[label_col]])), 10),
                              collapse = ", "))
  }
  meta <- meta[as.character(meta[[label_col]]) %in% used_labels, , drop = FALSE]
  cells_use <- rownames(meta)

  # ---- one pseudobulk sample per (cell type, label, replicate) --------------
  key <- data.frame(
    Cells     = as.character(meta[[celltype_col]]),
    label     = as.character(meta[[label_col]]),
    replicate = as.character(meta[[replicate_col]]),
    stringsAsFactors = FALSE
  )
  key_str <- paste(key$Cells, key$label, key$replicate, sep = "\r")
  first   <- !duplicated(key_str)

  sample_meta <- key[first, , drop = FALSE]
  sample_meta$key       <- key_str[first]
  sample_meta$sample_id <- paste0("pb", seq_len(nrow(sample_meta)))
  sample_meta$n_cells   <- as.integer(table(key_str)[sample_meta$key])
  rownames(sample_meta) <- NULL

  pb_id <- stats::setNames(sample_meta$sample_id[match(key_str, sample_meta$key)],
                           cells_use)

  # ---- aggregate: sum of raw counts through a sparse indicator matrix -------
  cnt <- .PBCounts(seu, assay)
  cnt <- cnt[, cells_use, drop = FALSE]
  if (!is.null(features)) {
    features <- intersect(features, rownames(cnt))
    if (!length(features)) stop("None of `features` is in assay '", assay, "'.")
    cnt <- cnt[features, , drop = FALSE]
  }
  ind <- Matrix::sparseMatrix(
    i    = seq_along(cells_use),
    j    = match(pb_id, sample_meta$sample_id),
    x    = 1,
    dims = c(length(cells_use), nrow(sample_meta)),
    dimnames = list(cells_use, sample_meta$sample_id)
  )
  pb <- cnt %*% ind

  # ---- filter ----------------------------------------------------------------
  cts <- sort(unique(sample_meta$Cells))
  kept_list <- list()
  summ_list <- list()

  .apply_filter <- function(x, grp) {
    keep <- switch(
      filter_method,
      edgeR  = do.call(edgeR::filterByExpr,
                       c(list(x, group = grp), filterByExpr_args)),
      manual = rowSums(x) >= mincounts,
      none   = rep(TRUE, nrow(x))
    )
    stats::setNames(as.logical(keep), rownames(x))
  }

  # filter_scope = "celltype": one filter per cell type on all its samples,
  # grouped by label (uses group sizes only, not differences between groups).
  keep_ct <- list()
  if (filter_scope == "celltype") {
    for (ct in cts) {
      s_ct <- sample_meta[sample_meta$Cells == ct, , drop = FALSE]
      if (!nrow(s_ct)) next
      x_ct <- as.matrix(pb[, s_ct$sample_id, drop = FALSE])
      keep_ct[[ct]] <- .apply_filter(x_ct, factor(s_ct$label))
      if (verbose) {
        message(sprintf("%s (all conditions): kept %d / %d %ss (%s)",
                        ct, sum(keep_ct[[ct]]), nrow(x_ct), feature_type, filter_method))
      }
    }
  }

  for (ci in seq_len(nrow(comparisons))) {
    trt  <- comparisons$treatment[ci]
    ctl  <- comparisons$control[ci]
    cond <- comparisons$Condition[ci]

    for (ct in cts) {
      s <- sample_meta[sample_meta$Cells == ct & sample_meta$label %in% c(trt, ctl), ,
                       drop = FALSE]
      n_cells_trt <- sum(s$n_cells[s$label == trt])
      n_cells_ctl <- sum(s$n_cells[s$label == ctl])

      row <- data.frame(
        Cells = ct, Condition = cond, treatment = trt, control = ctl,
        n_cells_trt = n_cells_trt, n_cells_ctl = n_cells_ctl,
        n_samples_trt = sum(s$label == trt), n_samples_ctl = sum(s$label == ctl),
        n_in = nrow(pb), n_kept = NA_integer_, status = "ok",
        stringsAsFactors = FALSE
      )

      if (n_cells_trt < min_cells || n_cells_ctl < min_cells) {
        row$status <- sprintf("skipped: < %d cells in a group", min_cells)
        summ_list[[length(summ_list) + 1]] <- row
        next
      }

      x    <- as.matrix(pb[, s$sample_id, drop = FALSE])
      grp  <- factor(s$label, levels = c(ctl, trt))
      keep <- if (filter_scope == "celltype") keep_ct[[ct]][rownames(x)] else .apply_filter(x, grp)

      row$n_kept <- sum(keep)
      summ_list[[length(summ_list) + 1]] <- row
      if (any(keep)) {
        kept_list[[length(kept_list) + 1]] <- data.frame(
          Cells = ct, Condition = cond, feature = rownames(x)[keep],
          stringsAsFactors = FALSE
        )
      }
      if (verbose && filter_scope == "comparison") {
        message(sprintf("%s | %s: kept %d / %d %ss (%s)",
                        cond, ct, sum(keep), length(keep), feature_type, filter_method))
      }
    }
  }

  summary_tbl <- dplyr::bind_rows(summ_list)
  if (verbose && any(summary_tbl$status != "ok")) {
    sk <- summary_tbl[summary_tbl$status != "ok", ]
    message("Skipped ", nrow(sk), " block(s): ",
            paste(paste0(sk$Cells, " (", sk$Condition, ")"), collapse = ", "))
  }

  list(
    pseudobulk    = pb,
    sample_meta   = sample_meta[, c("sample_id", "Cells", "label", "replicate", "n_cells")],
    kept          = dplyr::bind_rows(kept_list),
    summary       = summary_tbl,
    comparisons   = comparisons,
    funnel        = data.frame(
      stage = c(paste0(feature_type, "s in assay"), paste0(feature_type, "s tested")),
      n     = c(nrow(seu[[assay]]), nrow(pb))
    ),
    feature_type  = feature_type,
    assay         = assay,
    filter_method = filter_method,
    filter_scope  = filter_scope
  )
}


#' Gene-level differential expression per cell type and comparison (Libra)
#'
#' For every cell type x comparison block: builds pseudobulk counts, filters
#' lowly expressed genes (\code{\link{PseudobulkFilter}()}), and runs
#' \code{Libra::run_de()} on the retained genes only. Library sizes are then
#' computed on the retained genes (equivalent to \code{keep.lib.sizes = FALSE}
#' after \code{filterByExpr}). The group factor is ordered
#' \code{c(control, treatment)}, so \code{avg_logFC} is treatment vs control.
#'
#' Libra names its per-group columns after the group values (e.g.
#' \code{"B6_Control.pct"}); these are renamed by name to \code{group1.*}
#' (control) and \code{group2.*} (treatment), so every comparison returns the
#' same columns.
#'
#' @inheritParams PseudobulkFilter
#' @param de_family,de_method,de_type Passed to \code{Libra::run_de()}.
#'   Defaults \code{"pseudobulk"} / \code{"edgeR"} / \code{"LRT"}. Examples:
#'   \code{de_method = "DESeq2", de_type = "Wald"};
#'   \code{de_method = "limma", de_type = "voom"}.
#' @param libra_min_cells Libra's own \code{min_cells} (a filter inside
#'   \code{run_de()}, different from this function's \code{min_cells}).
#'   \code{NULL} (default) keeps Libra's default.
#' @param return_filter If \code{TRUE}, also return the
#'   \code{\link{PseudobulkFilter}()} object (for \code{\link{QCDensity}()}).
#'   Default \code{FALSE}.
#' @param ... Other arguments for \code{Libra::run_de()} (e.g.
#'   \code{min_reps}, \code{n_threads}).
#'
#' @return A list with \code{de}, a tibble with columns \code{cell_type},
#'   \code{gene}, \code{avg_logFC}, \code{group1.pct}, \code{group2.pct},
#'   \code{group1.exp}, \code{group2.exp}, \code{p_val}, \code{p_val_adj},
#'   \code{de_family}, \code{de_method}, \code{de_type}, \code{group1}
#'   (control), \code{group2} (treatment), \code{Condition},
#'   \code{regulation} and \code{n_tested} (genes tested in that block); plus
#'   \code{filter} when \code{return_filter = TRUE}.
#'
#' @examples
#' \dontrun{
#' seu$contrast <- paste0(seu$Strain, "_", seu$Treatment)
#' degs <- DEGsMatrix(seu, comparisons,
#'                    celltype_col = "Subpopulation",
#'                    label_col = "contrast", replicate_col = "Number",
#'                    return_filter = TRUE)
#' head(degs$de)
#' QCDensity(degs$filter, stage = "after", rows_by_condition = TRUE)
#' }
#'
#' @export
DEGsMatrix <- function(seu,
                       comparisons,
                       celltype_col,
                       label_col,
                       replicate_col,
                       assay             = "RNA",
                       features          = NULL,
                       min_cells         = 10,
                       filter_method     = c("edgeR", "manual", "none"),
                       mincounts         = 10,
                       filterByExpr_args = list(),
                       filter_scope      = c("comparison", "celltype"),
                       de_family         = "pseudobulk",
                       de_method         = "edgeR",
                       de_type           = "LRT",
                       libra_min_cells   = NULL,
                       return_filter     = FALSE,
                       verbose           = TRUE,
                       ...) {
  if (!requireNamespace("Libra", quietly = TRUE)) {
    stop("Package 'Libra' is required: remotes::install_github('neurorestore/Libra')")
  }

  flt <- PseudobulkFilter(
    seu, comparisons, celltype_col, label_col, replicate_col,
    assay = assay, features = features, min_cells = min_cells,
    filter_method = match.arg(filter_method), mincounts = mincounts,
    filterByExpr_args = filterByExpr_args,
    filter_scope = match.arg(filter_scope),
    feature_type = "gene", verbose = verbose
  )

  cnt    <- .PBCounts(seu, assay)
  blocks <- flt$summary[flt$summary$status == "ok" & flt$summary$n_kept > 0, , drop = FALSE]

  if (verbose) {
    message("Running Libra (", de_family, " / ", de_method, " / ", de_type,
            ") on ", nrow(blocks), " block(s)...")
  }

  libra_args <- c(
    list(de_family = de_family, de_method = de_method, de_type = de_type),
    if (!is.null(libra_min_cells)) list(min_cells = libra_min_cells),
    list(...)
  )

  .rename_groups <- function(out, ctl, trt) {
    for (suffix in c("pct", "exp")) {
      for (g in list(c(ctl, "group1"), c(trt, "group2"))) {
        cand <- c(paste0(g[1], ".", suffix), paste0(make.names(g[1]), ".", suffix))
        hit  <- intersect(cand, colnames(out))
        if (length(hit)) colnames(out)[colnames(out) == hit[1]] <- paste0(g[2], ".", suffix)
      }
    }
    out
  }

  md_all <- seu@meta.data
  res <- lapply(seq_len(nrow(blocks)), function(i) {
    b <- blocks[i, ]
    genes <- flt$kept$feature[flt$kept$Cells == b$Cells &
                                flt$kept$Condition == b$Condition]
    cells <- rownames(md_all)[as.character(md_all[[celltype_col]]) == b$Cells &
                              as.character(md_all[[label_col]]) %in% c(b$treatment, b$control)]

    md <- md_all[cells, c(celltype_col, label_col, replicate_col), drop = FALSE]
    md[[label_col]] <- factor(as.character(md[[label_col]]),
                              levels = c(b$control, b$treatment))

    obj <- suppressWarnings(Seurat::CreateSeuratObject(
      counts = cnt[genes, cells, drop = FALSE], meta.data = md
    ))

    out <- tryCatch(
      do.call(Libra::run_de, c(
        list(obj,
             cell_type_col = celltype_col,
             replicate_col = replicate_col,
             label_col     = label_col),
        libra_args
      )),
      error = function(e) {
        warning(sprintf("Libra failed for %s | %s: %s",
                        b$Condition, b$Cells, conditionMessage(e)), call. = FALSE)
        NULL
      }
    )
    if (is.null(out) || !nrow(out)) return(NULL)

    out <- .rename_groups(as.data.frame(out), b$control, b$treatment)
    out$group1    <- b$control
    out$group2    <- b$treatment
    out$Condition <- b$Condition
    out$n_tested  <- length(genes)
    out
  })

  de <- dplyr::bind_rows(res)
  if (nrow(de) && "avg_logFC" %in% colnames(de)) {
    de$regulation <- ifelse(de$avg_logFC > 0, "Up", "Down")
    front <- c("cell_type", "gene", "avg_logFC",
               "group1.pct", "group2.pct", "group1.exp", "group2.exp",
               "p_val", "p_val_adj", "de_family", "de_method", "de_type",
               "group1", "group2", "Condition", "regulation", "n_tested")
    de <- de[, c(intersect(front, colnames(de)), setdiff(colnames(de), front)), drop = FALSE]
  }
  de <- tibble::as_tibble(de)

  if (isTRUE(return_filter)) list(de = de, filter = flt) else list(de = de)
}
