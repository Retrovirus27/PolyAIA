# =============================================================================
#  Internal helpers for PolyAPlot(): coverage tracks built from the fragment
#  files (UMI deduplication, per-gene normalisation, overlay, replicates), the
#  site-usage panel and the gene-expression panel. Used only when PolyAPlot()
#  is called with any of dedup_umi / normalize != "signac" / overlay /
#  show_replicates / usage_panel / expression_panel; otherwise PolyAPlot()
#  keeps using Signac::CoveragePlot().
# =============================================================================


#' Read the fragments of one region for a set of cells
#'
#' Reads every fragment file of \code{assay} over \code{chr:start-end} with
#' tabix and maps file barcodes to object cell names.
#'
#' @return \code{data.frame(cell, start, end, umi)} (1-based), or \code{NULL}.
#' @noRd
.PolyAFragmentReads <- function(seu, assay, chr, start, end, cells) {
  gr  <- GenomicRanges::GRanges(chr, IRanges::IRanges(start, end))
  frs <- Signac::Fragments(seu[[assay]])
  if (!length(frs)) stop("Assay '", assay, "' has no fragment files attached.")
  out <- lapply(frs, function(fr) {
    path <- Signac::GetFragmentData(fr, "path")
    cmap <- Signac::GetFragmentData(fr, "cells")  # names = object cells, values = file barcodes
    l <- tryCatch(Rsamtools::scanTabix(Rsamtools::TabixFile(path), param = gr)[[1]],
                  error = function(e) character(0))
    if (!length(l)) return(NULL)
    f  <- strsplit(l, "\t", fixed = TRUE)
    bc <- vapply(f, `[`, character(1), 4)
    cell <- if (!is.null(cmap) && length(cmap)) names(cmap)[match(bc, cmap)] else bc
    ok <- !is.na(cell) & cell %in% cells
    if (!any(ok)) return(NULL)
    data.frame(
      cell  = cell[ok],
      start = as.integer(vapply(f[ok], `[`, character(1), 2)) + 1L,  # BED 0-based -> 1-based
      end   = as.integer(vapply(f[ok], `[`, character(1), 3)),
      umi   = vapply(f[ok], function(x) if (length(x) >= 5) x[5] else NA_character_,
                     character(1)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}


#' Build coverage tracks and site usage for PolyAPlot()
#'
#' @param tracks Named list: track name -> character vector of cells.
#' @param group_of Named character vector: cell -> group (tracks are split /
#'   overlaid by group).
#' @param group_levels Order of the groups.
#' @param rep_of Named character vector: cell -> replicate, or \code{NULL}.
#' @param sites \code{data.frame(site, peak, start, end)} of the gene's sites,
#'   labelled S1..Sn from 5' to 3'.
#' @param highlight \code{GRanges} of regions to shade, or \code{NULL}.
#' @return list(plots, usage_plot, usage, usage_rep, summary)
#' @noRd
.PolyACoverageTracks <- function(seu, assay, tracks, group_of, group_levels,
                                 rep_of, show_replicates, sites, chr,
                                 roi_start, roi_end, dedup_umi, normalize,
                                 overlay, bin_size, highlight, colors,
                                 usage_panel, text_size, axis_text_size,
                                 axis_title_size, strip_text_size,
                                 legend_text_size, verbose,
                                 legend_position = "top", usage_type = "bar",
                                 gene = NULL, expression_panel = FALSE,
                                 expression_source = "RNA", usage_scale = "percent") {

  all_cells <- unique(unlist(tracks, use.names = FALSE))
  rd_start  <- min(roi_start, sites$start)
  rd_end    <- max(roi_end,   sites$end)

  reads <- .PolyAFragmentReads(seu, assay, chr, rd_start, rd_end, all_cells)
  if (is.null(reads) || !nrow(reads)) {
    stop("No fragments from the selected cells in ", chr, ":", rd_start, "-", rd_end, ".")
  }
  n_reads <- nrow(reads)
  if (isTRUE(dedup_umi)) {
    if (all(is.na(reads$umi))) {
      warning("Fragment files have no UMI column; dedup_umi ignored.", call. = FALSE)
    } else {
      reads <- reads[!duplicated(paste(reads$cell, reads$umi)), , drop = FALSE]
    }
  }
  if (verbose) message("Coverage: ", n_reads, " reads -> ", nrow(reads),
                       if (isTRUE(dedup_umi)) " molecules (UMI-deduplicated)" else " reads")
  reads$group     <- unname(group_of[reads$cell])
  reads$replicate <- if (!is.null(rep_of)) unname(rep_of[reads$cell]) else NA_character_

  # molecules overlapping each site window / any site of the gene
  ov_site <- lapply(seq_len(nrow(sites)), function(i) reads$start <= sites$end[i] & reads$end >= sites$start[i])
  in_gene <- if (nrow(sites)) Reduce(`|`, ov_site) else rep(TRUE, nrow(reads))

  tot_counts <- if (normalize == "counts" || usage_scale == "cpm") {
    Matrix::colSums(SeuratObject::LayerData(seu, assay = assay, layer = "counts")[, all_cells, drop = FALSE])
  } else NULL
  # per-site log2(CPM + 1): site molecules / polyA counts of the unit's cells
  site_cpm <- function(n, unit_cells) {
    if (is.null(tot_counts)) return(rep(NA_real_, length(n)))
    log2(1e6 * n / max(sum(tot_counts[unit_cells]), 1) + 1)
  }

  # binned coverage over the plotted region only
  W     <- roi_end - roi_start + 1
  edges <- seq(roi_start, roi_end + bin_size, by = bin_size)
  bins  <- factor(cut(seq_len(W) + roi_start - 1, edges, right = FALSE, labels = FALSE),
                  levels = seq_len(length(edges) - 1))
  mids  <- edges[-length(edges)] + bin_size / 2

  cov_bins <- function(sel) {
    r <- reads[sel & reads$end >= roi_start & reads$start <= roi_end, , drop = FALSE]
    if (!nrow(r)) return(numeric(length(mids)))
    s0 <- pmax(r$start, roi_start) - roi_start + 1
    e0 <- pmax(pmin(r$end, roi_end) - roi_start + 1, s0)
    cv <- as.numeric(IRanges::coverage(IRanges::IRanges(s0, e0), width = W))
    v  <- as.numeric(tapply(cv, bins, mean))
    v[is.na(v)] <- 0
    v
  }
  scale_of <- function(sel, unit_cells) {
    switch(normalize,
           gene   = 1 / max(sum(sel & in_gene), 1),
           counts = 1e6 / max(sum(tot_counts[unit_cells]), 1))
  }
  site_n <- function(sel) vapply(ov_site, function(o) sum(o & sel), numeric(1))

  sig <- list(); sig_rep <- list(); usage <- list(); usage_rep <- list(); summ <- list()
  for (tr in names(tracks)) {
    tc <- tracks[[tr]]
    in_tr <- reads$cell %in% tc
    for (g in group_levels) {
      gc <- tc[group_of[tc] == g]
      if (!length(gc)) next
      sel <- in_tr & reads$group == g
      sig[[length(sig) + 1]] <- data.frame(track = tr, group = g, pos = mids,
                                           value = cov_bins(sel) * scale_of(sel, gc))
      n <- site_n(sel)
      usage[[length(usage) + 1]] <- data.frame(track = tr, group = g, site = sites$site,
                                               peak = sites$peak, n = n,
                                               pct = 100 * n / max(sum(n), 1),
                                               cpm = site_cpm(n, gc))
      summ[[length(summ) + 1]] <- data.frame(track = tr, group = g, cells = length(gc),
                                             molecules_in_gene = sum(sel & in_gene))
      if (!is.null(rep_of)) {
        for (rp in unique(rep_of[gc])) {
          rc   <- gc[rep_of[gc] == rp]
          selr <- sel & reads$replicate == rp
          if (isTRUE(show_replicates)) {
            sig_rep[[length(sig_rep) + 1]] <- data.frame(
              track = tr, group = g, replicate = rp, pos = mids,
              value = cov_bins(selr) * scale_of(selr, rc))
          }
          nr <- site_n(selr)
          usage_rep[[length(usage_rep) + 1]] <- data.frame(
            track = tr, group = g, replicate = rp, site = sites$site, peak = sites$peak,
            n = nr, pct = 100 * nr / max(sum(nr), 1), cpm = site_cpm(nr, rc))
        }
      }
    }
  }
  fix <- function(d) {
    if (is.null(d) || !length(d)) return(NULL)
    d <- do.call(rbind, d)
    d$group <- factor(d$group, levels = group_levels)
    d$track <- factor(d$track, levels = names(tracks))
    if ("site" %in% names(d)) d$site <- factor(d$site, levels = sites$site)
    if ("replicate" %in% names(d)) d$line <- paste(d$group, d$replicate, sep = "|")
    d
  }
  sig <- fix(sig); sig_rep <- fix(sig_rep); usage <- fix(usage); usage_rep <- fix(usage_rep)
  summ <- do.call(rbind, summ)

  # common y-axis across tracks
  y_max <- max(c(sig$value, if (!is.null(sig_rep)) sig_rep$value), na.rm = TRUE)
  if (!is.finite(y_max) || y_max <= 0) y_max <- 1
  # Same title as Signac for every normalisation
  y_name <- paste0("Normalized signal\n(range 0 - ", signif(y_max, 2), ")")

  hl_df <- if (!is.null(highlight) && length(highlight)) {
    data.frame(xmin = IRanges::start(highlight), xmax = IRanges::end(highlight))
  } else NULL

  plots <- lapply(names(tracks), function(tr) {
    d  <- sig[sig$track == tr, , drop = FALSE]
    dr <- if (!is.null(sig_rep)) sig_rep[sig_rep$track == tr, , drop = FALSE] else NULL
    has_rep <- !is.null(dr) && nrow(dr) > 0

    p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$pos, y = .data$value,
                                         fill = .data$group, colour = .data$group))
    if (!is.null(hl_df)) {
      p <- p + ggplot2::geom_rect(data = hl_df,
                                  ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax),
                                  ymin = -Inf, ymax = Inf, inherit.aes = FALSE,
                                  fill = "grey", alpha = 0.4)
    }
    if (isTRUE(overlay)) {
      p <- p +
        ggplot2::geom_area(position = "identity", alpha = if (has_rep) 0.15 else 0.25,
                           colour = NA)
      if (has_rep) {
        p <- p + ggplot2::geom_line(data = dr, ggplot2::aes(group = .data$line),
                                    linewidth = 0.3, alpha = 0.6)
      }
      p <- p + ggplot2::geom_line(linewidth = if (has_rep) 0.8 else 0.5)
    } else {
      p <- p +
        ggplot2::geom_area(stat = "identity", alpha = if (has_rep) 0.5 else 1, colour = NA)
      if (has_rep) {
        p <- p + ggplot2::geom_line(data = dr, ggplot2::aes(group = .data$line),
                                    linewidth = 0.3, alpha = 0.8)
      }
      p <- p +
        ggplot2::geom_hline(yintercept = 0, linewidth = 0.1) +   # Signac baseline
        ggplot2::facet_wrap(ggplot2::vars(.data$group), strip.position = "left", ncol = 1)
    }
    p <- p +
      ggplot2::coord_cartesian(xlim = c(roi_start, roi_end)) +
      ggplot2::scale_y_continuous(limits = c(0, y_max), name = y_name) +
      ggplot2::labs(x = NULL, fill = NULL, colour = NULL) +
      ggplot2::theme_classic(base_size = text_size) +
      ggplot2::theme(
        axis.text.x       = ggplot2::element_blank(),
        axis.ticks.x      = ggplot2::element_blank(),
        axis.line.x       = ggplot2::element_blank(),
        axis.text.y       = ggplot2::element_text(size = axis_text_size),
        axis.title.y      = ggplot2::element_text(size = axis_title_size),
        strip.text.y.left = ggplot2::element_text(angle = 0, size = strip_text_size),
        strip.background  = ggplot2::element_blank(),
        strip.placement   = "outside",
        panel.spacing.y   = ggplot2::unit(0, "lines"),
        # legend shown when groups share a panel (overlay) or when the usage
        # panel needs it; the usage panel itself has no legend
        legend.position   = if (isTRUE(overlay) || isTRUE(usage_panel) || isTRUE(expression_panel))
                              legend_position else "none",
        legend.text       = ggplot2::element_text(size = legend_text_size)
      )
    if (!isTRUE(overlay)) {
      # Exactly Signac's CoverageTrack look (theme_browser): theme_classic at
      # its default size, group names as left strips inside the axis (fixed
      # size in pt, they don't move when the figure is resized), no y tick
      # labels, the shared scale given only in the title "(range 0 - max)".
      p <- p +
        ggplot2::theme_classic() +
        ggplot2::theme(
          axis.text.x       = ggplot2::element_blank(),
          axis.ticks.x      = ggplot2::element_blank(),
          axis.line.x       = ggplot2::element_blank(),
          axis.text.y       = ggplot2::element_blank(),
          axis.title.y      = ggplot2::element_text(size = axis_title_size),
          strip.background  = ggplot2::element_blank(),
          strip.text.y.left = ggplot2::element_text(angle = 0),
          panel.spacing.y   = ggplot2::unit(0, "lines"),
          legend.position   = "none"
        )
    }
    if (!is.null(colors)) {
      p <- p + ggplot2::scale_fill_manual(values = colors) +
        ggplot2::scale_colour_manual(values = colors)
    }
    p
  })
  names(plots) <- names(tracks)

  # site-usage panel: bars = pooled group, points = replicates
  usage_plot <- NULL
  if (isTRUE(usage_panel) && !is.null(usage)) {
    dodge <- ggplot2::position_dodge(width = 0.8)
    use_y <- if (usage_scale == "cpm") "cpm" else "pct"
    if (usage_type == "box") {
      # one box per group x site over the replicates
      usage_plot <- ggplot2::ggplot(usage_rep, ggplot2::aes(x = .data$site, y = .data[[use_y]],
                                                            fill = .data$group)) +
        ggplot2::geom_boxplot(position = dodge, width = 0.7, alpha = 0.7,
                              outlier.shape = NA, colour = "grey20", linewidth = 0.3)
    } else {
      # bars = pooled group
      usage_plot <- ggplot2::ggplot(usage, ggplot2::aes(x = .data$site, y = .data[[use_y]],
                                                        fill = .data$group)) +
        ggplot2::geom_col(position = dodge, width = 0.75, alpha = 0.7)
    }
    if (isTRUE(show_replicates) && !is.null(usage_rep)) {
      usage_plot <- usage_plot +
        ggplot2::geom_point(
          data = usage_rep,
          ggplot2::aes(x = .data$site, y = .data[[use_y]], colour = .data$group),
          position = ggplot2::position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8,
                                                   seed = 1),
          size = 1.5, inherit.aes = FALSE, show.legend = FALSE)
    }
    if (length(tracks) > 1) {
      usage_plot <- usage_plot + ggplot2::facet_wrap(ggplot2::vars(.data$track), ncol = 1)
    }
    usage_plot <- usage_plot +
      ggplot2::labs(x = NULL, fill = NULL, colour = NULL,
                    y = if (usage_scale == "cpm") "Site expression\nlog2(CPM + 1)"
                        else "% of gene molecules") +
      ggplot2::theme_classic(base_size = text_size) +
      ggplot2::theme(
        axis.text        = ggplot2::element_text(size = axis_text_size),
        axis.title       = ggplot2::element_text(size = axis_title_size),
        strip.text       = ggplot2::element_text(size = strip_text_size, face = "bold"),
        strip.background = ggplot2::element_blank(),
        legend.position  = "none"
      )
    if (!is.null(colors)) {
      usage_plot <- usage_plot + ggplot2::scale_fill_manual(values = colors) +
        ggplot2::scale_colour_manual(values = colors)
    }
  }

  # ---- gene-expression panel: log2(CPM + 1) per group (and replicate) --------
  #   "RNA":   counts of the gene in the RNA assay / RNA library size of the cells
  #   "polyA": molecules of the gene (reads over its sites, after dedup) / polyA
  #            counts of the cells -- same data as the coverage
  expr <- NULL; expr_rep <- NULL; expr_plot <- NULL
  if (isTRUE(expression_panel)) {
    if (expression_source == "RNA") {
      # Switch to the RNA assay to read it (local copy of the object), joining
      # split Seurat v5 count layers (e.g. counts.1, counts.2) if needed.
      SeuratObject::DefaultAssay(seu) <- "RNA"
      rna <- seu[["RNA"]]
      if (inherits(rna, "Assay5") &&
          length(SeuratObject::Layers(rna, search = "counts")) > 1) {
        rna <- SeuratObject::JoinLayers(rna)
      }
      rna_cm <- SeuratObject::LayerData(rna, layer = "counts")
      SeuratObject::DefaultAssay(seu) <- assay   # back to the polyA assay
      feat <- gene
      if (!feat %in% rownames(rna_cm)) {
        rmeta <- seu[["RNA"]][[]]
        sc <- grep("^symbol$", colnames(rmeta), ignore.case = TRUE, value = TRUE)
        hit <- if (length(sc)) which(toupper(rmeta[[sc[1]]]) == toupper(gene)) else integer(0)
        if (!length(hit)) stop("Gene '", gene, "' not found in the RNA assay (rownames or SYMBOL column).")
        feat <- rownames(rmeta)[hit[1]]
      }
      g_cnt <- rna_cm[feat, all_cells]
      lib   <- Matrix::colSums(rna_cm[, all_cells, drop = FALSE])
    } else {
      mol_cell <- table(factor(reads$cell[in_gene], levels = all_cells))
      g_cnt <- stats::setNames(as.numeric(mol_cell), all_cells)
      lib   <- Matrix::colSums(SeuratObject::LayerData(seu, assay = assay, layer = "counts")[, all_cells, drop = FALSE])
    }
    cpm <- function(cells) log2(1e6 * sum(g_cnt[cells]) / max(sum(lib[cells]), 1) + 1)
    e_list <- list(); er_list <- list()
    for (tr in names(tracks)) {
      tc <- tracks[[tr]]
      for (g in group_levels) {
        gc <- tc[group_of[tc] == g]
        if (!length(gc)) next
        e_list[[length(e_list) + 1]] <- data.frame(track = tr, group = g, expr = cpm(gc))
        if (!is.null(rep_of)) {
          for (rp in unique(rep_of[gc])) {
            er_list[[length(er_list) + 1]] <- data.frame(track = tr, group = g, replicate = rp,
                                                         expr = cpm(gc[rep_of[gc] == rp]))
          }
        }
      }
    }
    expr <- fix(e_list); expr_rep <- fix(er_list)

    dodge <- ggplot2::position_dodge(width = 0.8)
    if (usage_type == "box" && !is.null(expr_rep)) {
      expr_plot <- ggplot2::ggplot(expr_rep, ggplot2::aes(x = .data$group, y = .data$expr,
                                                          fill = .data$group)) +
        ggplot2::geom_boxplot(width = 0.7, alpha = 0.7, outlier.shape = NA,
                              colour = "grey20", linewidth = 0.3)
    } else {
      expr_plot <- ggplot2::ggplot(expr, ggplot2::aes(x = .data$group, y = .data$expr,
                                                      fill = .data$group)) +
        ggplot2::geom_col(width = 0.75, alpha = 0.7)
    }
    if (isTRUE(show_replicates) && !is.null(expr_rep)) {
      expr_plot <- expr_plot +
        ggplot2::geom_point(data = expr_rep,
                            ggplot2::aes(x = .data$group, y = .data$expr, colour = .data$group),
                            position = ggplot2::position_jitter(width = 0.12, height = 0, seed = 1),
                            size = 1.5, inherit.aes = FALSE, show.legend = FALSE)
    }
    if (length(tracks) > 1) {
      expr_plot <- expr_plot + ggplot2::facet_wrap(ggplot2::vars(.data$track), ncol = 1)
    }
    expr_plot <- expr_plot +
      ggplot2::labs(x = NULL, fill = NULL, colour = NULL,
                    y = paste0(gene, " expression\nlog2(CPM + 1)",
                               if (expression_source == "polyA") " [polyA]" else "")) +
      ggplot2::theme_classic(base_size = text_size) +
      ggplot2::theme(
        axis.text        = ggplot2::element_text(size = axis_text_size),
        axis.text.x      = ggplot2::element_blank(),   # groups identified by colour
        axis.ticks.x     = ggplot2::element_blank(),
        axis.title       = ggplot2::element_text(size = axis_title_size),
        strip.text       = ggplot2::element_text(size = strip_text_size, face = "bold"),
        strip.background = ggplot2::element_blank(),
        legend.position  = "none"
      )
    if (!is.null(colors)) {
      expr_plot <- expr_plot + ggplot2::scale_fill_manual(values = colors) +
        ggplot2::scale_colour_manual(values = colors)
    }
  }

  list(plots = plots, usage_plot = usage_plot, usage = usage,
       usage_rep = usage_rep, expr_plot = expr_plot, expr = expr,
       expr_rep = expr_rep, summary = summ)
}
