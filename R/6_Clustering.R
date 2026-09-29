# =============================================================================
#  Cluster exploration and decontamination: ClusterTree(), ClusterComposition(),
#  FlagContaminants(). Clustering itself runs in SeuratPipeline(cluster = TRUE).
# =============================================================================


#' Clustering tree across resolutions
#'
#' Draws a \code{clustree} tree of the clusterings produced by
#' \code{SeuratPipeline(cluster = TRUE)} (one metadata column per resolution).
#' By default each node is coloured by the mean \code{Subpopulation_Score}
#' (Azimuth confidence) and labelled with its most frequent
#' \code{Subpopulation}.
#'
#' @param seurat_obj Seurat object with cluster columns (\code{prefix}*).
#' @param prefix Prefix of the cluster columns. Default \code{"RNA_snn_res."}.
#' @param node_colour Metadata column for the node colour, or
#'   \code{"sc3_stability"} (then set \code{node_colour_aggr = NULL}). Default
#'   \code{"Subpopulation_Score"}.
#' @param node_colour_aggr Aggregation of \code{node_colour} per node. Default
#'   \code{"mean"}.
#' @param node_label Metadata column for the node label. Default
#'   \code{"Subpopulation"}.
#' @param node_label_aggr Aggregation of \code{node_label} per node. Default
#'   \code{"getMode"}: the most frequent value (provided by this function if it
#'   is not already defined in the global environment).
#' @param node_text_colour Colour of the node text. Default \code{"white"}.
#' @param colour_title Title of the colour bar. Default \code{"SC3 Stability"}.
#' @param low,high Colours of the node colour gradient. Default \code{"blue"},
#'   \code{"red"}.
#' @param ... Passed to \code{clustree::clustree()}.
#'
#' @return A \code{ggplot}/\code{ggraph} object.
#'
#' @details Requires the \code{clustree} package (Suggests). \code{clustree}
#'   looks up the aggregation functions by name, so \code{getMode} is created
#'   in the global environment for the duration of the call when it does not
#'   exist there already. \code{ggraph} is attached if needed (ggplot2 finds
#'   ggraph's legends by name).
#'
#' @examples
#' \dontrun{
#' seurat_obj <- SeuratPipeline(seurat_obj)          # cluster = TRUE by default
#' ClusterTree(seurat_obj)
#' ClusterTree(seurat_obj, node_label = "CellType")
#' ClusterTree(seurat_obj, node_colour = "sc3_stability", node_colour_aggr = NULL)
#' }
#'
#' @export
ClusterTree <- function(seurat_obj,
                        prefix           = "RNA_snn_res.",
                        node_colour      = "Subpopulation_Score",
                        node_colour_aggr = "mean",
                        node_label       = "Subpopulation",
                        node_label_aggr  = "getMode",
                        node_text_colour = "white",
                        colour_title     = "SC3 Stability",
                        low              = "blue",
                        high             = "red",
                        ...) {
  if (!requireNamespace("clustree", quietly = TRUE)) {
    stop("Package 'clustree' is required. Install with: install.packages('clustree')")
  }
  if (!any(startsWith(colnames(seurat_obj@meta.data), prefix))) {
    stop("No metadata columns start with '", prefix, "'. Run SeuratPipeline(cluster = TRUE) ",
         "(or FindClusters) first.")
  }
  # library(clustree) attaches ggraph; ggplot2 needs it attached to find
  # ggraph's legends ("Unknown guide: edge_colourbar" otherwise).
  if (!"package:ggraph" %in% search()) {
    suppressPackageStartupMessages(attachNamespace("ggraph"))
  }
  # clustree looks up the aggregation function by name from the global
  # environment; provide getMode there for this call if it does not exist.
  genv <- globalenv()
  if (identical(node_label_aggr, "getMode") &&
      !exists("getMode", envir = genv, inherits = FALSE)) {
    assign("getMode", function(x) {
      ux <- unique(x)
      ux[which.max(tabulate(match(x, ux)))]
    }, envir = genv)
    on.exit(rm("getMode", envir = genv), add = TRUE)
  }

  clustree::clustree(seurat_obj,
                     prefix           = prefix,
                     node_colour      = node_colour,
                     node_colour_aggr = node_colour_aggr,
                     node_label       = node_label,
                     node_label_aggr  = node_label_aggr,
                     node_text_colour = node_text_colour,
                     ...) +
    ggplot2::scale_colour_gradient(low = low, high = high) +
    ggplot2::guides(colour = ggplot2::guide_colorbar(title = colour_title))
}


#' Cell type composition of clusters
#'
#' Counts cells per cluster x cell type and draws a heatmap of proportions.
#'
#' @param seurat_obj Seurat object.
#' @param cluster_col Metadata column with the clusters (e.g.
#'   \code{"RNA_snn_res.0.1"}).
#' @param celltype_col Metadata column with the cell type. Default
#'   \code{"CellType"}.
#' @param normalize What the heatmap proportions add up to:
#'   \code{"celltype"} (default) -- each cell type's cells across clusters
#'   (each row sums to 100\%; where each cell type goes);
#'   \code{"cluster"} -- each cluster's cells across cell types (each column
#'   sums to 100\%; what each cluster is made of).
#'
#' @return A list with \code{table} (columns \code{cluster}, \code{celltype},
#'   \code{n}, \code{prop_in_celltype}, \code{prop_in_cluster}) and \code{plot}.
#'
#' @examples
#' \dontrun{
#' comp <- ClusterComposition(seurat_obj, "RNA_snn_res.0.1")
#' comp$plot
#' ClusterComposition(seurat_obj, "RNA_snn_res.0.1", normalize = "cluster")$plot
#' }
#'
#' @export
ClusterComposition <- function(seurat_obj,
                               cluster_col,
                               celltype_col = "CellType",
                               normalize    = c("celltype", "cluster")) {
  normalize <- match.arg(normalize)
  md <- seurat_obj@meta.data
  for (cl in c(cluster_col, celltype_col)) {
    if (!cl %in% colnames(md)) stop("Column '", cl, "' not found in seurat_obj metadata.")
  }
  tab <- data.frame(cluster  = as.character(md[[cluster_col]]),
                    celltype = as.character(md[[celltype_col]]),
                    stringsAsFactors = FALSE) %>%
    dplyr::count(cluster, celltype, name = "n") %>%
    dplyr::group_by(celltype) %>%
    dplyr::mutate(prop_in_celltype = n / sum(n)) %>%
    dplyr::group_by(cluster) %>%
    dplyr::mutate(prop_in_cluster = n / sum(n)) %>%
    dplyr::ungroup()

  cl_levels <- unique(tab$cluster)
  num <- suppressWarnings(as.numeric(cl_levels))
  cl_levels <- if (!anyNA(num)) cl_levels[order(num)] else sort(cl_levels)
  tab$cluster <- factor(tab$cluster, levels = cl_levels)

  fill_col <- if (normalize == "celltype") "prop_in_celltype" else "prop_in_cluster"
  plot <- ggplot2::ggplot(tab, ggplot2::aes(x = .data$cluster, y = .data$celltype,
                                            fill = .data[[fill_col]])) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = scales::percent(.data[[fill_col]], accuracy = 0.1)),
                       size = 3) +
    ggplot2::scale_fill_gradient2(low = "white", mid = "orange", high = "darkred",
                                  midpoint = 0.5, labels = scales::percent_format()) +
    ggplot2::labs(
      x = paste0("Cluster (", cluster_col, ")"), y = celltype_col, fill = "Proportion",
      title = if (normalize == "celltype") "Cluster distribution per cell type"
              else "Cell type composition per cluster",
      subtitle = if (normalize == "celltype") "Each row (cell type) sums to 100%"
                 else "Each column (cluster) sums to 100%"
    ) +
    th

  list(table = tab, plot = plot)
}


#' Flag cells whose cell type disagrees with their cluster
#'
#' Each cluster gets an expected cell type by majority vote; cells of any other
#' cell type in that cluster are flagged as contaminants. The object is only
#' annotated -- remove the flagged cells yourself, e.g.
#' \code{subset(seurat_obj, Is_Contaminant == FALSE)}.
#'
#' @param seurat_obj Seurat object.
#' @param cluster_col Metadata column with the clusters (e.g.
#'   \code{"RNA_snn_res.0.1"}).
#' @param celltype_col Metadata column with the cell type. Default
#'   \code{"CellType"}.
#' @param protected Character vector of cell types (values of
#'   \code{celltype_col}) that are never flagged, whatever cluster they fall in
#'   -- e.g. \code{c("End")} for a rare population without a cluster of its
#'   own. Defined by cell type, so it does not depend on cluster numbers.
#'   Default \code{NULL}.
#' @param assign How the expected cell type of a cluster is chosen:
#'   \code{"cluster"} (default) -- the most abundant cell type within the
#'   cluster; \code{"celltype"} -- the cell type with the largest share of its
#'   own cells in that cluster (proportions normalised per cell type). The
#'   latter can let a rare cell type that sits entirely in a large cluster
#'   "win" it, flagging the cluster's majority; kept for comparison.
#' @param na_contaminant Logical. Flag cells with \code{NA} cell type. Default
#'   \code{TRUE}.
#' @param verbose Logical. Print the per-cluster summary. Default \code{TRUE}.
#'
#' @return \code{seurat_obj} with metadata columns \code{Expected_CellType} and
#'   \code{Is_Contaminant}. A per-cluster summary (\code{cells},
#'   \code{contaminants}, \code{pct_contaminant}, \code{protected_kept} = cells
#'   of a protected type kept although they differ from the expected cell type)
#'   is stored in \code{seurat_obj@misc$contaminant_summary}.
#'
#' @examples
#' \dontrun{
#' seurat_obj <- FlagContaminants(seurat_obj, "RNA_snn_res.0.1", protected = c("End"))
#' seurat_obj_clean <- subset(seurat_obj, Is_Contaminant == FALSE)
#' }
#'
#' @export
FlagContaminants <- function(seurat_obj,
                             cluster_col,
                             celltype_col   = "CellType",
                             protected      = NULL,
                             assign         = c("cluster", "celltype"),
                             na_contaminant = TRUE,
                             verbose        = TRUE) {
  assign <- match.arg(assign)
  md <- seurat_obj@meta.data
  for (cl in c(cluster_col, celltype_col)) {
    if (!cl %in% colnames(md)) stop("Column '", cl, "' not found in seurat_obj metadata.")
  }
  clu <- as.character(md[[cluster_col]])
  ct  <- as.character(md[[celltype_col]])

  if (is.list(protected)) {
    stop("`protected` must be a character vector of cell types, e.g. c(\"End\").")
  }
  protected <- as.character(protected)
  bad_ct <- setdiff(protected, unique(ct))
  if (length(bad_ct)) {
    warning("Protected cell type(s) not present in '", celltype_col, "': ",
            paste(bad_ct, collapse = ", "), call. = FALSE)
  }

  # expected cell type per cluster (majority vote)
  comp <- ClusterComposition(seurat_obj, cluster_col, celltype_col)$table %>%
    dplyr::filter(!is.na(celltype))
  score_col <- if (assign == "cluster") "prop_in_cluster" else "prop_in_celltype"
  expected <- comp %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(.data[[score_col]], n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  exp_map <- stats::setNames(expected$celltype, as.character(expected$cluster))

  expected_ct <- unname(exp_map[clu])
  is_cont   <- ct != expected_ct
  is_prot   <- ct %in% protected
  kept_prot <- is_prot & (is_cont %in% TRUE)   # would have been flagged
  is_cont[is_prot]   <- FALSE                  # protected cell types: never flagged
  is_cont[is.na(ct)] <- isTRUE(na_contaminant)

  seurat_obj$Expected_CellType <- expected_ct
  seurat_obj$Is_Contaminant    <- is_cont

  summ <- data.frame(cluster = clu, expected = expected_ct, contaminant = is_cont,
                     kept_prot = kept_prot, stringsAsFactors = FALSE) %>%
    dplyr::group_by(cluster, expected) %>%
    dplyr::summarise(cells           = dplyr::n(),
                     contaminants    = sum(contaminant, na.rm = TRUE),
                     pct_contaminant = round(100 * contaminants / cells, 1),
                     protected_kept  = sum(kept_prot),
                     .groups = "drop")
  num <- suppressWarnings(as.numeric(summ$cluster))
  if (!anyNA(num)) summ <- summ[order(num), ]
  seurat_obj@misc$contaminant_summary <- summ

  if (verbose) {
    message("FlagContaminants (", cluster_col, ", assign = \"", assign, "\"): ",
            sum(is_cont, na.rm = TRUE), " of ", length(is_cont), " cells flagged (",
            round(100 * mean(is_cont, na.rm = TRUE), 1), "%)",
            if (length(protected)) paste0("; ", sum(kept_prot), " cell(s) of protected type(s) {",
                                          paste(protected, collapse = ", "), "} kept") else "",
            ".")
    print(as.data.frame(summ), row.names = FALSE)
  }
  seurat_obj
}
