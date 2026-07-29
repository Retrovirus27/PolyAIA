#' Filter cells with per-group (e.g. per-cell-type) QC thresholds
#'
#' Applies QC cutoffs that vary by a grouping variable (typically
#' \code{CellType}) instead of one global threshold for the whole object. Each
#' cell is matched to its group's row in a \code{thresholds} table, and cells
#' failing any supplied bound (min/max total counts, min/max detected features,
#' max mitochondrial percentage) for their own group are removed.
#'
#' Any threshold column absent from \code{thresholds} is simply not applied, so
#' you can filter on only the metrics you care about.
#'
#' @param seurat_obj A Seurat object carrying the QC metric columns below.
#' @param thresholds A \code{data.frame} of per-group cutoffs: one row per group,
#'   a column named \code{group_col} holding the group labels, and any of the
#'   optional numeric columns \code{min_counts}, \code{max_counts},
#'   \code{min_features}, \code{max_features}, \code{mito_cutoff}. Only the
#'   columns present are enforced.
#' @param group_col Metadata/thresholds column used to match each cell to its
#'   group. Default \code{"CellType"}.
#' @param count_col,feature_col,mito_col Metadata columns holding total counts,
#'   detected features, and mitochondrial percentage. Defaults
#'   \code{"nCount_RNA"}, \code{"nFeature_RNA"}, \code{"percent_mito"}.
#' @param unmatched What to do with cells whose group has no row in
#'   \code{thresholds}: \code{"drop"} (default -- removes them, matching a plain
#'   \code{subset()} against \code{NA} thresholds) or \code{"keep"} (retains
#'   them unfiltered).
#' @param store_thresholds Logical. If \code{TRUE}, keep the per-cell threshold
#'   columns (\code{min_counts}, ...) on the returned object's metadata. Default
#'   \code{FALSE} (they're used internally, then dropped).
#' @param verbose Logical. Print a kept/removed summary (overall and per group).
#'   Default \code{TRUE}.
#'
#' @return The filtered Seurat object.
#'
#' @examples
#' \dontrun{
#' thresholds <- data.frame(
#'   CellType     = c("Astro", "End", "Exc", "Inh", "Mic", "Oligo", "OPC"),
#'   min_counts   = c(700, 300, 2200, 1500, 250, 250, 250),
#'   max_counts   = c(5000, 5000, 9000, 9000, 5000, 5000, 5000),
#'   min_features = c(700, 300, 2200, 1500, 250, 250, 250),
#'   max_features = c(3000, 2000, 4000, 4000, 2000, 2000, 3000),
#'   mito_cutoff  = c(3, 3, 3, 3, 3, 3, 3)
#' )
#' seurat_obj <- QCFilter(seurat_obj, thresholds, group_col = "CellType")
#' }
#'
#' @export
QCFilter <- function(
    seurat_obj,
    thresholds,
    group_col        = "CellType",
    count_col        = "nCount_RNA",
    feature_col      = "nFeature_RNA",
    mito_col         = "percent_mito",
    unmatched        = c("drop", "keep"),
    store_thresholds = FALSE,
    verbose          = TRUE
) {
  unmatched <- match.arg(unmatched)

  if (!methods::is(seurat_obj, "Seurat")) stop("`seurat_obj` must be a Seurat object.")
  thresholds <- as.data.frame(thresholds, stringsAsFactors = FALSE)

  if (!group_col %in% colnames(thresholds)) {
    stop("`thresholds` must contain the group column '", group_col, "'.")
  }
  if (!group_col %in% colnames(seurat_obj@meta.data)) {
    stop("`group_col` '", group_col, "' not found in seurat_obj metadata.")
  }
  if (anyDuplicated(thresholds[[group_col]])) {
    stop("`thresholds` has duplicate '", group_col, "' values -- one row per group is expected.")
  }

  # Which bounds were supplied, the metric column each needs, and the comparison.
  bound_specs <- list(
    min_counts   = list(metric = count_col,   op = `>`),
    max_counts   = list(metric = count_col,   op = `<`),
    min_features = list(metric = feature_col, op = `>`),
    max_features = list(metric = feature_col, op = `<`),
    mito_cutoff  = list(metric = mito_col,    op = `<`)
  )
  active <- intersect(names(bound_specs), colnames(thresholds))
  if (length(active) == 0) {
    stop("`thresholds` supplies none of: ",
         paste(names(bound_specs), collapse = ", "), ".")
  }
  needed_metrics <- unique(vapply(bound_specs[active], function(x) x$metric, character(1)))
  miss_metric <- setdiff(needed_metrics, colnames(seurat_obj@meta.data))
  if (length(miss_metric) > 0) {
    stop("Metric column(s) required by the supplied thresholds not found in ",
         "metadata: ", paste(miss_metric, collapse = ", "), ".")
  }

  meta <- seurat_obj@meta.data
  n0   <- ncol(seurat_obj)
  idx  <- match(as.character(meta[[group_col]]), as.character(thresholds[[group_col]]))
  matched <- !is.na(idx)

  # Apply each supplied bound using the per-cell threshold for its group. NA
  # (unmatched cells / any NA threshold) propagates and is resolved below.
  keep <- rep(TRUE, nrow(meta))
  for (b in active) {
    spec      <- bound_specs[[b]]
    bound_vec <- thresholds[[b]][idx]
    metric    <- meta[[spec$metric]]
    keep      <- keep & spec$op(metric, bound_vec)
  }

  # Resolve unmatched cells per policy, then treat any residual NA as fail.
  keep[!matched] <- (unmatched == "keep")
  keep[is.na(keep)] <- FALSE

  if (verbose) {
    n_unmatched <- sum(!matched)
    if (n_unmatched > 0) {
      message("QCFilter: ", n_unmatched, " cell(s) had no matching '", group_col,
              "' row in thresholds -> ", unmatched, ".")
    }
    kept_tab   <- table(as.character(meta[[group_col]])[keep])
    total_tab  <- table(as.character(meta[[group_col]]))
    for (g in names(total_tab)) {
      k <- kept_tab[g]
      if (is.na(k)) k <- 0L
      message(sprintf("  %-14s kept %d / %d", g, as.integer(k), as.integer(total_tab[g])))
    }
  }

  # Optionally stash the per-cell thresholds before subsetting.
  if (store_thresholds) {
    for (b in active) seurat_obj@meta.data[[b]] <- thresholds[[b]][idx]
  }

  seurat_obj <- seurat_obj[, keep]

  if (verbose) {
    message(sprintf("QCFilter: kept %d / %d cells (%.1f%%).",
                    ncol(seurat_obj), n0, 100 * ncol(seurat_obj) / n0))
  }

  seurat_obj
}

#' Create comprehensive QC plots for Seurat objects
#'
#' @param seurat_obj Seurat object with QC metrics
#' @param count_col Column name for total counts (default "nCount_RNA")
#' @param feature_col Column name for features (default "nFeature_RNA")
#' @param mito_col Column name for mitochondrial percentage (default "percent_mito")
#' @param plot_types Vector of plot types to generate. Options:
#'   "histogram", "violin", "scatter", or "all" (default)
#' @param split_by split "histogram", "violin", "scatter" by a variable
#' @param count_max Maximum x-axis value for count histogram (auto if NULL)
#' @param count_breaks Break interval for count histogram (auto-calculated if NULL)
#' @param feature_max Maximum x-axis value for features histogram (auto if NULL)
#' @param feature_breaks Break interval for features histogram (auto-calculated if NULL)
#' @param mt_max Maximum x-axis value for mt histogram (auto if NULL)
#' @param mt_breaks Break interval for mt histogram (auto-calculated if NULL)
#' @param min_counts Optional lower filtering cutoff for total counts. If supplied,
#'   drawn as a dashed reference line on the count histogram/scatter/violin plots.
#'   Default NULL (no line drawn)
#' @param max_counts Optional upper filtering cutoff for total counts, drawn the
#'   same way as \code{min_counts}. Default NULL (no line drawn)
#' @param min_features Optional lower filtering cutoff for detected features,
#'   drawn as a dashed reference line on the feature histogram/scatter/violin
#'   plots. Default NULL (no line drawn)
#' @param max_features Optional upper filtering cutoff for detected features,
#'   drawn the same way as \code{min_features}. Default NULL (no line drawn)
#' @param mito_cutoff Optional filtering cutoff for mitochondrial percentage,
#'   drawn as a dashed reference line on the mt histogram/violin plots. Default
#'   NULL (no line drawn)
#' @param bins Number of bins for histograms (default 50)
#' @param base_size Base font size (default 9)
#' @param point_alpha Transparency for scatter/violin points (default 0.3)
#' @param n_col Number of columns for layout (default NULL, auto-arranged)
#' @param n_row Number of rows for layout (default NULL, auto-arranged)
#'
#'
#' @return Combined patchwork plot object
#' @export
QCPlot <- function(
    seurat_obj,
    count_col = "nCount_RNA",
    feature_col = "nFeature_RNA",
    mito_col = "percent_mito",
    plot_types = "all",
    split_by = NULL,
    count_max = NULL,
    count_breaks = NULL,
    feature_max = NULL,
    feature_breaks = NULL,
    mt_max = NULL,
    mt_breaks = NULL,
    min_counts = NULL,
    max_counts = NULL,
    min_features = NULL,
    max_features = NULL,
    mito_cutoff = NULL,
    bins = 50,
    base_size = 9,
    point_alpha = 0.3,
    n_col = NULL,
    n_row = NULL
) {

  # Extract metadata
  meta_data <- seurat_obj@meta.data

  # Validate required columns exist up front, with a clear message, instead of
  # silently propagating -Inf/NA through the auto-calculated axis ranges.
  required_cols <- c(count_col, feature_col, mito_col)
  missing_cols <- required_cols[!required_cols %in% names(meta_data)]
  if (length(missing_cols) > 0) {
    stop("Column(s) not found in seurat_obj@meta.data: ", paste(missing_cols, collapse = ", "))
  }

  # Determine plot types
  valid_types <- c("histogram", "violin", "scatter", "all")
  invalid_types <- setdiff(plot_types, valid_types)
  if (length(invalid_types) > 0) {
    warning("Ignoring unrecognized plot_types: ", paste(invalid_types, collapse = ", "),
            ". Valid options: ", paste(valid_types, collapse = ", "))
  }
  if ("all" %in% plot_types) {
    plot_types <- c("histogram", "violin", "scatter")
  }

  if (!is.null(split_by)) {
    if (length(split_by) > 2) {
      stop("split_by supports at most 2 variables (got ", length(split_by), ").")
    }
    missing_split <- split_by[!split_by %in% names(meta_data)]
    if (length(missing_split) > 0) {
      stop("Column(s) ", paste(missing_split, collapse = ", "), " (split_by) not found in data.")
    }
  }

  # Auto-calculate x axis if not provided
  if (is.null(count_max))   { count_max   <- max(meta_data[[count_col]], na.rm = TRUE)   }
  if (is.null(feature_max)) { feature_max <- max(meta_data[[feature_col]], na.rm = TRUE) }
  if (is.null(mt_max))      { mt_max      <- max(meta_data[[mito_col]], na.rm = TRUE)    }

  if (is.null(count_breaks))   count_breaks   <- diff(pretty(c(0, count_max),   n = 6))[1]
  if (is.null(feature_breaks)) feature_breaks <- diff(pretty(c(0, feature_max), n = 6))[1]
  if (is.null(mt_breaks))      mt_breaks      <- diff(pretty(c(0, mt_max),      n = 5))[1]

  # Draw the actual filtering cutoffs (min_counts/max_counts/min_features/
  # max_features/mito_cutoff), when supplied, as dashed reference lines --
  # distinct from count_max/feature_max/mt_max above, which only control the
  # axis range. `values` already drops NULLs via c(), so passing e.g.
  # c(min_counts, max_counts) when only one is set still works.
  .cutoff_lines <- function(values, orientation = c("v", "h")) {
    orientation <- match.arg(orientation)
    if (length(values) == 0) return(list())
    if (orientation == "v") {
      list(ggplot2::geom_vline(xintercept = values, linetype = "dashed", color = "red3", linewidth = 0.6))
    } else {
      list(ggplot2::geom_hline(yintercept = values, linetype = "dashed", color = "red3", linewidth = 0.6))
    }
  }
  count_cutoffs   <- c(min_counts, max_counts)
  feature_cutoffs <- c(min_features, max_features)

  # Initialize plot list
  plot_list <- list()

  # === HISTOGRAM PLOTS ===
  if ("histogram" %in% plot_types) {

    # Count histogram
    hist_count <- meta_data %>%
      ggplot2::ggplot(ggplot2::aes(x = .data[[count_col]])) +
      ggplot2::geom_histogram(bins = bins, fill = "steelblue", color = "black") +
      .cutoff_lines(count_cutoffs, "v") +
      ggplot2::scale_x_continuous(
        limits = c(0, count_max),
        breaks = seq(0, count_max, count_breaks)
      ) +
      ggplot2::labs(x = NULL, y = NULL, title = "nCount") +
      th +
      ggplot2::theme( axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1) )

    # Feature histogram
    hist_feature <- meta_data %>%
      ggplot2::ggplot(ggplot2::aes(x = .data[[feature_col]])) +
      ggplot2::geom_histogram(bins = bins, fill = "skyblue", color = "black") +
      .cutoff_lines(feature_cutoffs, "v") +
      ggplot2::scale_x_continuous(
        limits = c(0, feature_max),
        breaks = seq(0, feature_max, feature_breaks)
      ) +
      ggplot2::labs(x = NULL, y = NULL, title = "nFeature") +
      th +
      ggplot2::theme( axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1) )

    # Mitochondrial percentage histogram
    hist_mito <- meta_data %>%
      ggplot2::ggplot(ggplot2::aes(x = .data[[mito_col]])) +
      ggplot2::geom_histogram(bins = bins, fill = "salmon", color = "black") +
      .cutoff_lines(mito_cutoff, "v") +
      ggplot2::scale_x_continuous(
        limits = c(0, mt_max),
        breaks = seq(0, mt_max, mt_breaks)
      ) +
      ggplot2::labs(x = NULL, y = NULL, title = "% Mt") +
      th +
      ggplot2::theme( axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1) )

    # Add facet if split_by is set
    if (!is.null(split_by)) {
      # Ensure the variable(s) exist
      if (any(!split_by %in% names(meta_data))) {
        stop(paste("Column(s)", paste(split_by[!split_by %in% names(meta_data)], collapse = ", "), "not found in data."))
      }
      # If 2 variables -> facet_grid() dynamically
      if (length(split_by) == 2) {
        facet_formula <- stats::as.formula(
          paste(split_by[1], "~", split_by[2])
        )
        hist_count = hist_count + ggh4x::facet_nested_wrap(facet_formula, nrow = n_row, ncol = n_col, scales = "free", labeller = ggplot2::label_both) +
          ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5),
                         strip.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = 0.5))
        hist_feature = hist_feature + ggh4x::facet_nested_wrap(facet_formula, nrow = n_row, ncol = n_col, scales = "free", labeller = ggplot2::label_both) +
          ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5),
                         strip.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = 0.5))
        hist_mito = hist_mito + ggh4x::facet_nested_wrap(facet_formula, nrow = n_row, ncol = n_col, scales = "free", labeller = ggplot2::label_both) +
          ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5),
                         strip.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = 0.5))
      } else {
        # Otherwise facet_wrap as before
        hist_count = hist_count + ggh4x::facet_nested_wrap(ggplot2::vars(.data[[split_by[1]]]), nrow = n_row, ncol = n_col, scales = "free", labeller = ggplot2::label_both) +
          ggplot2::theme(strip.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = 0.5))
        hist_feature = hist_feature + ggh4x::facet_nested_wrap(ggplot2::vars(.data[[split_by[1]]]), nrow = n_row, ncol = n_col, scales = "free", labeller = ggplot2::label_both) +
          ggplot2::theme(strip.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = 0.5))
        hist_mito = hist_mito + ggh4x::facet_nested_wrap(ggplot2::vars(.data[[split_by[1]]]), nrow = n_row, ncol = n_col, scales = "free", labeller = ggplot2::label_both) +
          ggplot2::theme(strip.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = 0.5))
      }
    }

    plot_list <- c(plot_list, list(
      hist_count = hist_count,
      hist_feature = hist_feature,
      hist_mito = hist_mito
    ))
  }

  # === SCATTER PLOT ===
  if ("scatter" %in% plot_types) {

    scatter_plot <- meta_data %>%
      ggplot2::ggplot(ggplot2::aes(x = .data[[count_col]],
                 y = .data[[feature_col]],
                 color = .data[[mito_col]])) +
      ggplot2::geom_point(alpha = point_alpha + 0.4, size = 1) +
      .cutoff_lines(count_cutoffs, "v") +
      .cutoff_lines(feature_cutoffs, "h") +
      ggplot2::scale_color_gradient(
        low = "darkblue",
        high = "gold",
        limits = c(0, mt_max),
        breaks = seq(0, mt_max, mt_breaks)
      ) +
      ggplot2::scale_x_continuous(
        limits = c(0, count_max),
        breaks = seq(0, count_max, count_breaks)
      ) +
      ggplot2::scale_y_continuous(
        limits = c(0, feature_max),
        breaks = seq(0, feature_max, feature_breaks)
      ) +
      ggplot2::labs(x = NULL, y = NULL, title = "Linear relation", color = "% Mt") +
      th +
      ggplot2::theme( axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1) )

    # Add facet if split_by is set
    if (!is.null(split_by)) {
      # Ensure the variable(s) exist
      if (any(!split_by %in% names(meta_data))) {
        stop(paste("Column(s)", paste(split_by[!split_by %in% names(meta_data)], collapse = ", "), "not found in data."))
      }
      # If 2 variables -> facet_grid() dynamically
      if (length(split_by) == 2) {
        facet_formula <- stats::as.formula(
          paste(split_by[1], "~", split_by[2])
        )
        scatter_plot = scatter_plot + ggh4x::facet_nested_wrap(facet_formula, nrow = n_row, ncol = n_col, labeller = ggplot2::label_both) +
          ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5),
                         strip.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = 0.5))
      } else {
        # Otherwise facet_wrap as before
        scatter_plot = scatter_plot + ggh4x::facet_nested_wrap(ggplot2::vars(.data[[split_by[1]]]), nrow = n_row, ncol = n_col, labeller = ggplot2::label_both) +
          ggplot2::theme(strip.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = 0.5))
      }
    }

    plot_list <- c(plot_list, list(scatter = scatter_plot))
  }

  # === VIOLIN PLOTS ===
  if ("violin" %in% plot_types) {

    # Count violin
    violin_count <- meta_data %>%
      ggplot2::ggplot(ggplot2::aes(x = "Counts", y = .data[[count_col]])) +
      ggplot2::geom_violin(fill = "lightblue", alpha = 0.7) +
      ggplot2::geom_jitter(alpha = point_alpha, width = 0.1, size = 0.5) +
      .cutoff_lines(count_cutoffs, "h") +
      ggplot2::scale_y_continuous(
        limits = c(0, count_max),
        breaks = seq(0, count_max, count_breaks)
      ) +
      ggplot2::labs(x = NULL, y = NULL, title = NULL) +
      th +
      ggplot2::theme( axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1) )

    # Feature violin
    violin_feature <- meta_data %>%
      ggplot2::ggplot(ggplot2::aes(x = "Features", y = .data[[feature_col]])) +
      ggplot2::geom_violin(fill = "lightblue", alpha = 0.7) +
      ggplot2::geom_jitter(alpha = point_alpha, width = 0.1, size = 0.5) +
      .cutoff_lines(feature_cutoffs, "h") +
      ggplot2::scale_y_continuous(
        limits = c(0, feature_max),
        breaks = seq(0, feature_max, feature_breaks)
      ) +
      ggplot2::labs(x = NULL, y = NULL, title = NULL) +
      th +
      ggplot2::theme( axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1) )

    # Mitochondrial violin
    violin_mito <- meta_data %>%
      ggplot2::ggplot(ggplot2::aes(x = "% Mt", y = .data[[mito_col]])) +
      ggplot2::geom_violin(fill = "lightblue", alpha = 0.7) +
      ggplot2::geom_jitter(alpha = point_alpha, width = 0.1, size = 0.5) +
      .cutoff_lines(mito_cutoff, "h") +
      ggplot2::scale_y_continuous(
        limits = c(0, mt_max),
        breaks = seq(0, mt_max, mt_breaks)
      ) +
      ggplot2::labs(x = NULL, y = NULL, title = NULL) +
      th +
      ggplot2::theme( axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1) )

    # Cell count barplot
    if (!is.null(split_by)) {
      bar_cells <- meta_data %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(split_by))) %>%
        dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
        ggplot2::ggplot(ggplot2::aes(x = if (length(split_by) == 2) .data[[split_by[2]]] else .data[[split_by[1]]], y = n)) +
        ggplot2::geom_col(fill = "lightblue", alpha = 0.7) +
        ggplot2::geom_text(ggplot2::aes(label = scales::comma(n)), vjust = -0.5, size = 3) +
        ggplot2::labs(x = NULL, y = NULL, title = NULL) +
        th +
        ggplot2::theme( axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1) )

      if (length(split_by) == 2) {
        bar_cells <- bar_cells +
          ggplot2::facet_wrap(stats::as.formula(paste("~", split_by[1])), nrow = n_row, ncol = n_col, labeller = ggplot2::label_both) +
          ggplot2::theme(strip.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = 0.5))
      } else {
        bar_cells <- bar_cells +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
      }
    } else {
      bar_cells <- meta_data %>%
        dplyr::summarise(n = dplyr::n()) %>%
        ggplot2::ggplot(ggplot2::aes(x = "Number of cells", y = n)) +
        ggplot2::geom_col(fill = "lightblue", alpha = 0.7) +
        ggplot2::geom_text(ggplot2::aes(label = scales::comma(n)), vjust = -0.5, size = 3) +
        ggplot2::labs(x = NULL, y = NULL, title = NULL) +
        th +
        ggplot2::theme( axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1) )
    }

    # Add facet if split_by is set
    if (!is.null(split_by)) {
      # Ensure the variable(s) exist
      if (any(!split_by %in% names(meta_data))) {
        stop(paste("Column(s)", paste(split_by[!split_by %in% names(meta_data)], collapse = ", "), "not found in data."))
      }
      # If 2 variables -> facet_grid() dynamically
      if (length(split_by) == 2) {
        facet_formula <- stats::as.formula(
          paste(split_by[1], "~", split_by[2])
        )
        violin_count = violin_count + ggh4x::facet_nested_wrap(facet_formula, nrow = n_row, ncol = n_col, labeller = ggplot2::label_both) +
          ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5),
                         strip.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = 0.5))
        violin_feature = violin_feature + ggh4x::facet_nested_wrap(facet_formula, nrow = n_row, ncol = n_col, labeller = ggplot2::label_both) +
          ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5),
                         strip.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = 0.5))
        violin_mito = violin_mito + ggh4x::facet_nested_wrap(facet_formula, nrow = n_row, ncol = n_col, labeller = ggplot2::label_both) +
          ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5),
                         strip.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = 0.5))
      } else {
        # Otherwise facet_wrap as before
        violin_count = violin_count + ggh4x::facet_nested_wrap(ggplot2::vars(.data[[split_by[1]]]), nrow = n_row, ncol = n_col, labeller = ggplot2::label_both) +
          ggplot2::theme(strip.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = 0.5))
        violin_feature = violin_feature + ggh4x::facet_nested_wrap(ggplot2::vars(.data[[split_by[1]]]), nrow = n_row, ncol = n_col, labeller = ggplot2::label_both) +
          ggplot2::theme(strip.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = 0.5))
        violin_mito = violin_mito + ggh4x::facet_nested_wrap(ggplot2::vars(.data[[split_by[1]]]), nrow = n_row, ncol = n_col, labeller = ggplot2::label_both) +
          ggplot2::theme(strip.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = 0.5))
      }
    }

    plot_list <- c(plot_list, list(
      violin_count = violin_count,
      violin_feature = violin_feature,
      violin_mito = violin_mito,
      bar_cells = bar_cells
    ))
  }


  # === COMBINE PLOTS ===
  if (length(plot_list) == 0) {
    stop("No valid plot types selected")
  }

  # Arrange plots based on what was generated
  if (all(c("histogram", "violin", "scatter") %in% plot_types)) {
    # Full layout: histograms + scatter on top, violins on bottom
    combined_plot <- (
      plot_list$hist_count |
        plot_list$hist_feature |
        plot_list$hist_mito |
        plot_list$scatter
    ) / (
      plot_list$violin_count |
        plot_list$violin_feature |
        plot_list$violin_mito |
        plot_list$bar_cells
    ) + patchwork::plot_annotation(
      title = "Quality Control Metrics Overview",
      tag_levels = "A"
    )
  } else {
    # Auto-arrange based on selected plots
    if (is.null(n_col)) {
      n_col <- min(4, length(plot_list))
    }
    combined_plot <- patchwork::wrap_plots(plot_list, ncol = n_col, nrow = n_row) +
      patchwork::plot_annotation(
        title = "Quality Control Metrics",
        tag_levels = "A"
      )
  }

  return(combined_plot)
}

#' Create UMAP/dimensionality reduction plots
#'
#' Generate customizable UMAP or other dimensionality reduction plots with
#' optional faceting, labels, and color schemes for continuous or discrete variables.
#'
#' @param seurat_obj Seurat object
#' @param reduction.method Name of the reduction method (e.g., "Harmony", "UMAP", "PCA").
#'   Looked up in \code{seurat_obj} as \code{"umap_<reduction.method>"} (lowercased,
#'   e.g. \code{"umap_harmony"}) since that's this pipeline's convention for the 2D
#'   UMAP embedding, separate from the higher-dimensional integration/PCA space
#'   itself. Falls back to the bare lowercased name if no \code{"umap_"}-prefixed
#'   reduction is found. Coordinate column names are read directly from the
#'   resulting embedding, not guessed from this argument.
#' @param color Column name to use for coloring points. Can be continuous or discrete
#' @param label_by Column name to use for labels (independent of color). Required if show_labels = TRUE
#' @param split_by Optional column name(s) for faceting. Can be a single variable or
#'   a vector of two variables for grid faceting
#' @param palette_d Viridis palette option for discrete variables. Default is "turbo"
#' @param saturation Desaturation amount for discrete palettes. Default is 0.2
#' @param alpha_point Transparency of the points. Default is 0.5
#' @param size_point Size of the points. Default is 0.1
#' @param palette_c Color palette for continuous variables. A vector of three colors
#'   for low, mid, and high values. Default is c("blue","lightgray","red")
#' @param limits Limits for continuous color scale. Default is c(-0.1, 1.1)
#' @param breaks Breaks for continuous color scale. Default is seq(0, 1, 0.2)
#' @param title.prefix Optional prefix for plot title
#' @param legend.position Position of legend. Default is "top"
#' @param legend_point_size Size of points in legend for discrete variables. Default is 4
#' @param show_labels Logical. If TRUE, adds text labels at the centroid of each group. Default is FALSE
#' @param show_numbers Logical. If TRUE, includes cell counts in labels. Default is FALSE
#' @param n_col Number of columns for facet layout. Default is NULL (automatic)
#' @param n_row Number of rows for facet layout. Default is NULL (automatic)
#' @param sizes Numeric value for base font size. Default is c(9,0.5,11,10,10) for the following: base_size, base_line_size, legend.title, legend.text, strip.text, strip.text.x and strip.text.y
#'
#' @return A ggplot2 object
#'
#' @export
#' @concept visualization
#'
UmapPlot <-
  function(seurat_obj,
           reduction.method = "Harmony",
           color = "Population",
           label_by = NULL,
           split_by = NULL,
           palette_d = "turbo",
           saturation = 0.2,
           alpha_point = 0.5,
           size_point = 0.1,
           palette_c = c("blue","lightgray","red"),
           limits = c(-0.1, 1.1),
           breaks = seq(0, 1, 0.2),
           title.prefix = NULL,
           legend.position = "top",
           legend_point_size = 4,
           show_labels = FALSE,
           show_numbers = FALSE,
           n_col = NULL,
           n_row = NULL,
           sizes = c(9,0.5,11,10,10)) {

    # UmapPlot() is for 2D UMAP-style visualizations, and this pipeline's
    # convention is to name that embedding "umap_<method>" (e.g. "umap_harmony"),
    # keeping it separate from the higher-dimensional integration/PCA space
    # itself (e.g. "harmony", 50 dims -- plotting that directly would just show
    # its first two dimensions, which looks like a PCA plot, not a UMAP).
    # Try the "umap_" prefix first; fall back to the bare name in case
    # reduction.method already *is* the exact reduction slot name.
    available_reductions <- SeuratObject::Reductions(seurat_obj)
    reduction_key <- paste0("umap_", tolower(reduction.method))
    if (!reduction_key %in% available_reductions) {
      fallback_key <- tolower(reduction.method)
      if (fallback_key %in% available_reductions) {
        reduction_key <- fallback_key
      } else {
        stop("Reduction '", reduction_key, "' (or '", fallback_key, "') not found in seurat_obj. ",
             "Available reductions: ", paste(available_reductions, collapse = ", "))
      }
    }

    embeddings <- as.data.frame(SeuratObject::Embeddings(seurat_obj, reduction = reduction_key))

    # Read the x/y coordinate column names directly off the embeddings instead
    # of guessing them from reduction.method (e.g. "Harmony_1") -- the actual
    # key set on the reduction when it was created (often lowercase, e.g.
    # "harmony_1") doesn't necessarily match the capitalization of the display
    # name passed in here, and guessing wrong throws a cryptic dplyr error.
    x_col <- colnames(embeddings)[1]
    y_col <- colnames(embeddings)[2]

    data <- embeddings %>%
      dplyr::bind_cols(seurat_obj@meta.data)

    # Validate split_by / label_by up front -- these are used later to build
    # label_df and facets, and a clear error here beats a cryptic dplyr/facet
    # error deep in the function.
    if (!is.null(split_by)) {
      if (length(split_by) > 2) {
        stop("split_by supports at most 2 variables (got ", length(split_by), ").")
      }
      missing_split <- split_by[!split_by %in% names(data)]
      if (length(missing_split) > 0) {
        stop("Column(s) ", paste(missing_split, collapse = ", "), " (split_by) not found in data.")
      }
    }
    if (show_labels && !is.null(label_by) && !label_by %in% names(data)) {
      stop("Column '", label_by, "' (label_by) not found in data.")
    }

    th <- ggprism::theme_prism(
      base_size = sizes[1],
      base_line_size = sizes[2]
    )
    th$legend.title <- ggplot2::element_text(size = sizes[3], face = "bold")
    th$legend.text  <- ggplot2::element_text(size = sizes[4])
    th$strip.text   <- ggplot2::element_text(size = sizes[5], face = "bold")  # Tamaño de facet
    th$strip.text.x <- ggplot2::element_text(size = sizes[5], face = "bold")  # Facet horizontal
    th$strip.text.y <- ggplot2::element_text(size = sizes[5], face = "bold", angle = 90)  # Facet vertical

    # Create arrows for axes
    arrow_x <- grid::segmentsGrob(
      x0 = grid::unit(0.05, "npc"),
      y0 = grid::unit(0.05, "npc"),
      x1 = grid::unit(0.05, "npc") + grid::unit(1, "cm"),  # longitud fija
      y1 = grid::unit(0.05, "npc"),
      arrow = grid::arrow(length = grid::unit(0.15, "cm"), type = "closed"),
      gp = grid::gpar(col = "black", fill = "black", lwd = 1.5)
    )

    arrow_y <- grid::segmentsGrob(
      x0 = grid::unit(0.05, "npc"),
      y0 = grid::unit(0.05, "npc"),
      x1 = grid::unit(0.05, "npc"),
      y1 = grid::unit(0.05, "npc") + grid::unit(1, "cm"),  # longitud fija
      arrow = grid::arrow(length = grid::unit(0.15, "cm"), type = "closed"),
      gp = grid::gpar(col = "black", fill = "black", lwd = 1.5)
    )

    p <-
      (if (is.null(color)) {
        ggplot2::ggplot(data, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]]))
      } else {
        ggplot2::ggplot(data, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]], color = .data[[color]]))
      }) +
      ggplot2::geom_point(alpha = alpha_point, size = size_point) +
      ggplot2::labs(title = if (!is.null(title.prefix)) title.prefix else NULL,x = NULL, y = NULL) +
      th +
      ggplot2::theme(legend.position = legend.position,
                     axis.line = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     axis.text = ggplot2::element_blank(),
                     plot.title = ggplot2::element_blank(),
                     plot.margin = ggplot2::margin(0, 0, 0, 0)  ) +
      ggplot2::annotation_custom(arrow_x, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
      ggplot2::annotation_custom(arrow_y, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)

    # auto-detect discrete vs continuous (only relevant if a color variable was given)
    if (!is.null(color)) {
      var <- if (is.numeric(data[[color]])) "continuous" else "discrete"

      p = p +
        # Color scales based on variable type
        switch(
          var,
          "continuous" = ggplot2::scale_color_gradient2(
            low = palette_c[1], mid = palette_c[2], high = palette_c[3], midpoint = 0,
            limits = limits, breaks = breaks, oob = scales::squish
          ),
          "discrete" = ggplot2::scale_color_manual(values = colorspace::desaturate(viridisLite::viridis(length(unique(data[[color]])) + 1, option = palette_d)[2:(length(unique(data[[color]])) + 1)], amount = saturation))
        ) +
        switch(
          var,
          "continuous" = ggplot2::guides(color = ggplot2::guide_colorbar(barwidth = 0.6, barheight = 5)),
          "discrete"   = ggplot2::guides(color = ggplot2::guide_legend(
            override.aes = list(size = legend_point_size, alpha = 1)
          ))
        )
    }

    # Add labels if requested
    if (show_labels) {
      # Check that label_by is specified
      if (is.null(label_by)) {
        warning("show_labels is TRUE but label_by is NULL. No labels will be added.")
      } else {
        # If split_by is defined, group by both label_by and split variables
        if (!is.null(split_by)) {
          label_df <- data %>%
            dplyr::group_by(dplyr::across(tidyselect::all_of(c(split_by, label_by)))) %>%
            dplyr::summarise(
              x_mean = mean(.data[[x_col]], na.rm = TRUE),
              y_mean = mean(.data[[y_col]], na.rm = TRUE),
              n_cells = dplyr::n(),
              .groups = "drop"
            )
        } else {
          label_df <- data %>%
            dplyr::group_by(.data[[label_by]]) %>%
            dplyr::summarise(
              x_mean = mean(.data[[x_col]], na.rm = TRUE),
              y_mean = mean(.data[[y_col]], na.rm = TRUE),
              n_cells = dplyr::n(),
              .groups = "drop"
            )
        }

        # Build label text (with or without counts)
        label_df$label_text <- if (show_numbers) {
          paste0(label_df[[label_by]], " (", label_df$n_cells, ")")
        } else {
          as.character(label_df[[label_by]])
        }

        # Add labels, respecting facet mapping
        p <- p + ggrepel::geom_text_repel(
          data = label_df,
          ggplot2::aes(x = x_mean, y = y_mean, label = label_text),
          size = 7, color = "black", fontface = "bold",
          box.padding = 0.3, max.overlaps = Inf,
          inherit.aes = FALSE
        )
      }
    }

    # Add facet if split_by is set
    if (!is.null(split_by)) {
      # Ensure the variable(s) exist
      if (any(!split_by %in% names(data))) {
        stop(paste("Column(s)", paste(split_by[!split_by %in% names(data)], collapse = ", "), "not found in data."))
      }
      # If 2 variables -> facet_grid() dynamically
      if (length(split_by) == 2) {
        facet_formula <- stats::as.formula(
          paste(split_by[1], "~", split_by[2])
        )
        p = p + ggh4x::facet_nested_wrap(facet_formula, nrow = n_row, ncol = n_col) +
          ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5),
                strip.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = 0.5))
      } else {
        # Otherwise facet_wrap as before
        p = p + ggh4x::facet_nested_wrap(ggplot2::vars(.data[[split_by[1]]]), nrow = n_row, ncol = n_col) +
          ggplot2::theme(strip.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = 0.5))
      }
    }

    return(p)
  }


#' Create a Dot Plot for Gene Expression Visualization
#'
#' Generates a dot plot showing average gene expression and percentage of cells
#' expressing each gene across different groups. Supports optional conditioning
#' by a secondary variable and grouping of genes into categories.
#'
#' @param seu A Seurat object containing single-cell RNA-seq data. If
#'   \code{rownames(seu)} aren't already gene symbols, they are automatically
#'   swapped for the assay's symbol feature-metadata column (matched
#'   case-insensitively against \code{"SYMBOL"}, as set by
#'   \code{AnnotateGenes()}/\code{UploadSce()}), when available. Requested
#'   \code{features} are also matched case-insensitively as a fallback (e.g. a
#'   human-style \code{"GFAP"} will still match a mouse object's \code{"Gfap"}).
#' @param group_var Character string specifying the column name in the Seurat
#'   object metadata to use for grouping (e.g., cell types, clusters).
#' @param features A named list of character vectors, where each element contains
#'   gene names for a feature group. Names will be used as facet labels.
#'   Example: \code{list(Markers1 = c("Gene1", "Gene2"), Markers2 = c("Gene3", "Gene4"))}.
#' @param condition_var Character string specifying an optional column name for
#'   a secondary grouping variable (e.g., treatment condition). Default is \code{NULL}.
#' @param palette Character vector of length 3 specifying colors for the gradient
#'   (low, mid, high). Default is \code{c("blue", "white", "red")}.
#' @param limit_fc Numeric vector of length 2 specifying the limits for the
#'   color scale (z-score range). Default is \code{c(-2, 2)}.
#' @param dot_scale Numeric value controlling the maximum size of dots.
#'   Default is \code{5}.
#' @param rotate_x Logical indicating whether to rotate x-axis labels.
#'   Default is \code{TRUE}.
#' @param title Character string for the plot title. Default is \code{NULL}.
#' @param poslegend Character string specifying the position of the legend
#'   (e.g., \code{"right"}, \code{"bottom"}, \code{"left"}, \code{"top"}).
#'   Default is \code{"right"}.
#' @param sizes Numeric vector of length 4 controlling font sizes for various
#'   plot elements: axis text, axis line size, title text, and legend text,
#'   respectively. Default is \code{c(9, 0.5, 11, 10)}.
#'
#' @return A \code{ggplot2} object showing a dot plot where:
#'   \itemize{
#'     \item Dot size represents the percentage of cells expressing each gene.
#'     \item Dot color represents scaled average expression (z-score).
#'     \item Genes are grouped by feature categories displayed as facets.
#'     \item If \code{condition_var} is provided, the plot is additionally
#'           faceted by the condition variable.
#'   }
#'
#' @examples
#' \dontrun{
#' # Basic usage with cell type markers
#' markers <- list(
#'   Tcells = c("CD3D", "CD3E"),
#'   Bcells = c("CD19", "MS4A1")
#' )
#' DotPlot(seurat_obj, group_var = "celltype", features = markers)
#'
#' # With condition variable and custom legend position
#' DotPlot(
#'   seurat_obj,
#'   group_var  = "celltype",
#'   features   = markers,
#'   condition_var = "treatment",
#'   poslegend  = "bottom",
#'   sizes      = c(9, 0.5, 11, 10)
#' )
#' }
#'
#' @export
DotPlot <- function(
    seu,
    group_var,
    features,
    condition_var = NULL,
    palette = c("blue","white","red"),
    limit_fc = c(-2,2),
    dot_scale = 5,
    rotate_x = TRUE,
    title = NULL,
    poslegend = "right",
    sizes = c(9,0.5,11,10)
){
  th <- ggprism::theme_prism(
    base_size = sizes[1],
    base_line_size = sizes[2]
  )
  th$legend.title <- ggplot2::element_text(size = sizes[3], face = "bold")
  th$legend.text  <- ggplot2::element_text(size = sizes[4])

  # Validate grouping columns up front instead of failing deep inside FetchData
  meta_cols <- colnames(seu@meta.data)
  if (!group_var %in% meta_cols) {
    stop("group_var '", group_var, "' not found in seu metadata.")
  }
  if (!is.null(condition_var) && !condition_var %in% meta_cols) {
    stop("condition_var '", condition_var, "' not found in seu metadata.")
  }

  # `features` is normally supplied as gene symbols (e.g. "Gfap"), but rownames(seu)
  # may still be whatever ID type UploadSce()/AnnotateGenes() used (e.g. Ensembl,
  # via IDtype), or the feature-metadata symbol column may have been dropped by an
  # intervening split()/JoinLayers() (Seurat v5 doesn't always carry per-feature
  # metadata through that round-trip). Look for a "SYMBOL"-like feature-metadata
  # column (matched case-insensitively, in case it was renamed) and switch to it
  # as the source of truth when found -- this replaces having to manually build a
  # symbol-renamed copy of `seu` (e.g.
  # `rownames(seurat_obj2) <- seurat_obj@assays[["RNA"]]@meta.data$SYMBOL`)
  # before calling DotPlot().
  assay_meta <- seu[[SeuratObject::DefaultAssay(seu)]]@meta.data
  symbol_col <- grep("^symbol$", colnames(assay_meta), ignore.case = TRUE, value = TRUE)
  symbol_col <- if (length(symbol_col) > 0) symbol_col[1] else NA_character_
  used_symbol_col <- FALSE
  if (!is.na(symbol_col) && !identical(rownames(seu), assay_meta[[symbol_col]])) {
    symbols <- assay_meta[[symbol_col]]
    # Guard against NA/blank/duplicated symbols, which would otherwise produce
    # NA or non-unique rownames
    missing_symbol <- is.na(symbols) | symbols == ""
    symbols[missing_symbol] <- rownames(seu)[missing_symbol]
    if (anyDuplicated(symbols)) symbols <- make.unique(symbols)
    rownames(seu) <- symbols
    used_symbol_col <- TRUE
  }

  # Genes that exist in the object
  genes_in_seu <- rownames(seu)

  # Resolve requested gene names against genes_in_seu: exact match first, then
  # a case-insensitive fallback (handles e.g. a human-style "GFAP" supplied
  # against a mouse object's "Gfap"), returning names in genes_in_seu's actual
  # case so downstream FetchData()/pivot_longer() always get a real match.
  resolve_genes <- function(x) {
    exact <- x[x %in% genes_in_seu]
    unmatched <- setdiff(x, exact)
    if (length(unmatched) > 0) {
      lut <- stats::setNames(genes_in_seu, toupper(genes_in_seu))
      ci_matches <- lut[toupper(unmatched)]
      exact <- c(exact, unname(ci_matches[!is.na(ci_matches)]))
    }
    unique(exact)
  }

  # Flatten the marker list
  feature_vec <- resolve_genes(unique(unlist(features)))
  if (length(feature_vec) == 0) {
    stop(
      "None of the requested `features` were found in `seu` (checked against gene ",
      "symbols, case-insensitively). ",
      if (used_symbol_col) {
        paste0("Renamed rownames(seu) using its '", symbol_col, "' feature-metadata column. ")
      } else if (!is.na(symbol_col)) {
        paste0("A '", symbol_col, "' feature-metadata column exists but already matched rownames(seu). ")
      } else {
        paste0("No SYMBOL-like feature-metadata column was found (it may have been dropped ",
               "by an intervening split()/JoinLayers()), so rownames(seu) were used as-is. ")
      },
      "First few rownames(seu): ", paste(utils::head(genes_in_seu, 5), collapse = ", "), ". ",
      "Check for typos or a species/ID-type mismatch."
    )
  }
  # Fix the feature list to remove invalid genes (same exact/case-insensitive resolution)
  features <- lapply(features, resolve_genes)
  # Now safely build vars
  vars <- c(group_var, feature_vec)
  if (!is.null(condition_var)) vars <- c(vars, condition_var)
  # Now FetchData will never fail
  df <- SeuratObject::FetchData(seu, vars = vars)
  # Pivot safer
  df_long <- df %>%
    tidyr::pivot_longer(
      cols = tidyselect::all_of(feature_vec),
      names_to = "gene",
      values_to = "expr"
    )
  # Agrupar dependiendo de condición
  if (is.null(condition_var)) {
    stats_df <- df_long %>%
      dplyr::group_by(!!dplyr::sym(group_var), gene) %>%
      dplyr::summarise(
        avg_exp = mean(expr, na.rm = TRUE),
        pct_exp = mean(expr > 0, na.rm = TRUE) * 100,
        .groups = "drop"
      )
  } else {
    stats_df <- df_long %>%
      dplyr::group_by(!!dplyr::sym(group_var), !!dplyr::sym(condition_var), gene) %>%
      dplyr::summarise(
        avg_exp = mean(expr, na.rm = TRUE),
        pct_exp = mean(expr > 0, na.rm = TRUE) * 100,
        .groups = "drop"
      )
  }
  # Convierte tu lista de marcadores en un dataframe
  marker_df <- utils::stack(features)
  colnames(marker_df) <- c("gene", "GeneGroup")
  # Extraer datos
  stats_df <- stats_df %>%
    dplyr::left_join(marker_df, by = "gene") %>%
    dplyr::group_by(gene) %>%
    dplyr::mutate(avg_xp_scaled = as.numeric(scale(avg_exp)))

  # ---- Plot ----
  p <- ggplot2::ggplot(stats_df, ggplot2::aes(x = gene, y = !!dplyr::sym(group_var),
                         size = pct_exp, color = avg_xp_scaled)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_gradient2(low = palette[1], mid = palette[2], high = palette[3],
                          midpoint = 0, limits = limit_fc, oob = scales::squish) +
    ggplot2::scale_size(range = c(1, dot_scale), name = "% Cells") +
    ggplot2::labs(x = NULL, y = NULL, color = "Z score", title = title) +
    th +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      panel.spacing = grid::unit(0.5, "lines"),
      legend.position = poslegend
    )
  if (is.null(condition_var)) {
    p <- p + ggplot2::facet_wrap("GeneGroup", scales = "free_x", nrow = 1)
  }
  if (!is.null(condition_var)) {
    p <- p + ggplot2::facet_grid( rows  = ggplot2::vars(!!dplyr::sym(condition_var)), cols  = ggplot2::vars(GeneGroup),
                         scales = "free_x", space  = "free_x" )
  }
  p <- p + ggplot2::theme(strip.text = ggplot2::element_blank())
  return(p)
}


# Suggest a number of PCs/dimensions to carry into downstream steps (e.g.
# IntegrateLayers()/RunUMAP() in SeuratPipeline(), see 1_AnnotationSceData.R)
# by locating the ACTUAL elbow ("codo") of the PCA standard-deviation curve --
# i.e. the point where the ElbowPlot bends most sharply.
#
# Uses the standard geometric "distance to the chord" method (a.k.a. the
# Kneedle idea): draw the straight line from the first PC's stdev to the last
# PC's stdev, then pick the PC whose stdev point lies farthest from that line.
# That maximal-curvature point is the elbow. Unlike the previous
# percent-change-below-a-cutoff heuristic, this always returns a real bend
# within 1:n.pcs (no cutoff to tune, no fall-back-to-n.pcs warning).
# @noRd
.suggest_n_dims <- function(object, reduction.method = "pca", n.pcs = 30) {
  stdev_vec <- object@reductions[[reduction.method]]@stdev
  if (length(stdev_vec) < n.pcs) {
    n.pcs <- length(stdev_vec)
  }
  y <- stdev_vec[seq_len(n.pcs)]

  # Degenerate cases: with fewer than 3 points there's no interior bend.
  if (n.pcs < 3) return(n.pcs)

  x <- seq_len(n.pcs)

  # Perpendicular distance from each (x, y) point to the chord joining the
  # first point (x1, y1) and the last point (xn, yn).
  x1 <- x[1]; y1 <- y[1]
  xn <- x[n.pcs]; yn <- y[n.pcs]
  # Line through (x1,y1)-(xn,yn) as a*x + b*y + c = 0.
  a <- yn - y1
  b <- x1 - xn
  cc <- xn * y1 - x1 * yn
  denom <- sqrt(a^2 + b^2)
  if (denom == 0) return(n.pcs)
  dist_to_line <- abs(a * x + b * y + cc) / denom

  which.max(dist_to_line)
}

#' Elbow and JackStraw plots for PCA dimensionality selection
#'
#' Generates diagnostic plots to assist in selecting the optimal number of
#' principal components for downstream clustering and embedding. The number
#' of PCs is suggested based on the point at which the percent change in
#' standard deviation between consecutive PCs falls below a minimal
#' threshold (indicating diminishing returns). The cumulative standard
#' deviation cutoff is shown for visual reference only and does not affect
#' the suggested PC.
#'
#' @param object A Seurat object with PCA computed.
#' @param reduction.method Reduction method to evaluate. Default is \code{"pca"}.
#' @param cumulative.stdev.cutoff Cumulative standard deviation threshold (as a
#'   percentage of total), shown on the plot for visual reference only.
#'   Default is \code{51}.
#' @param pct.change.cutoff Minimum percent change in standard deviation
#'   between consecutive PCs. The first PC at which the change drops below
#'   this threshold is used to suggest the number of PCs. Default is
#'   \code{0.1}.
#' @param n.pcs Number of leading PCs to display/evaluate. Default is \code{30}.
#' @param JackStraw Logical. Whether to compute and plot JackStraw results.
#'   Reserved for future implementation; currently has no effect. Default is
#'   \code{FALSE}.
#'
#' @return A patchwork of Elbow and (optionally) JackStraw plots. The
#'   suggested PC (based on percent change) is shown in the plot subtitle.
#'
#' @examples
#' # ClusterPlots(seurat_obj)
#' # ClusterPlots(seurat_obj, pct.change.cutoff = 0.1)
ClusterPlots <- function(object,
                         reduction.method = "pca",
                         cumulative.stdev.cutoff = 51,
                         pct.change.cutoff = 0.1,
                         n.pcs = 30,
                         JackStraw = FALSE) {

  stdev_vec <- object@reductions[[reduction.method]]@stdev

  if (length(stdev_vec) < n.pcs) {
    n.pcs <- length(stdev_vec)
  }

  pc_data <- stdev_vec %>%
    as.data.frame() %>%
    stats::setNames("Stdev") %>%
    tibble::rownames_to_column(var = "PC") %>%
    dplyr::slice(1:n.pcs) %>%
    dplyr::mutate(
      PC_num          = as.numeric(PC),
      CumulativeStdev = cumsum(Stdev),
      PctChange       = c(NA, abs(diff(Stdev) / Stdev[-length(Stdev)]) * 100)
    )

  # Suggested PC: first point where consecutive % change drops below threshold
  suggested_pc <- which(pc_data$PctChange < pct.change.cutoff)[1]

  if (is.na(suggested_pc)) {
    warning(
      "No PC met the pct.change.cutoff of ", pct.change.cutoff,
      "% within the first ", n.pcs, " PCs. Consider increasing n.pcs ",
      "or relaxing pct.change.cutoff."
    )
  }

  # cumulative.stdev.cutoff used only for visual reference (coloring), not for the decision
  pc_data <- pc_data %>%
    dplyr::mutate(color = ifelse(CumulativeStdev < cumulative.stdev.cutoff, "TRUE", "FALSE"))

  subtitle_text <- paste0(
    "Suggested PCs (\u0394stdev < ", pct.change.cutoff, "%): ",
    ifelse(is.na(suggested_pc), "not reached", suggested_pc)
  )

  ClusterElbow <- pc_data %>%
    ggplot2::ggplot(ggplot2::aes(x = CumulativeStdev, y = Stdev, colour = color)) +
    ggplot2::geom_point() +
    ggplot2::geom_text(ggplot2::aes(label = PC), vjust = -0.5, size = 3) +
    ggplot2::scale_x_continuous(breaks = seq(0, 100, 10)) +
    ggplot2::scale_color_manual(values = c("TRUE" = "black", "FALSE" = "lightgray")) +
    ggplot2::labs(
      title = paste0("Elbow Plot - ", toupper(reduction.method)),
      subtitle = subtitle_text,
      x = "Cumulative SD (color reference only, cutoff = cumulative.stdev.cutoff)",
      y = "Standard Deviation (SD)"
    ) +
    th +
    ggplot2::theme(legend.position = "none")

  if (!is.na(suggested_pc)) {
    ClusterElbow <- ClusterElbow +
      ggplot2::geom_vline(
        xintercept = pc_data$CumulativeStdev[suggested_pc],
        linetype = "dashed", colour = "gray40"
      )
  }

  # JackStraw plotting reserved for future implementation
  # if (JackStraw) {
  #   ClusterJack <- object@reductions[[reduction.method]]@jackstraw$overall.p.values %>%
  #     as.data.frame() %>%
  #     ggplot2::ggplot(ggplot2::aes(x = PC, y = log10(Score))) +
  #     ggplot2::geom_point() +
  #     ggplot2::labs(title = paste0("JackStraw Plot - ", toupper(reduction.method))) +
  #     Universal_theme
  #   return(ClusterElbow + ClusterJack)
  # }

  return(ClusterElbow)
}

