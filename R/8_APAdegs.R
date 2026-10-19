# =============================================================================
#  Gene-level differential expression from polyA-site counts (APAdeg)
#  Complements DEGsMatrix() (RNA) and DEPsMatrix() (site usage / RED).
#  APAdeg: https://github.com/XuLabSJTU/APAdeg (GPL-3, called, not copied).
# =============================================================================


# Size factors as APAdeg computes them: DESeq-style median of ratios on the
# gene sums (genes x samples matrix).
.APAdegSizeFactors <- function(g) {
  g  <- g[rowSums(g) > 0, , drop = FALSE]
  gm <- apply(g, 1, function(v) exp(mean(log(v[v > 0]))))
  apply(g / gm, 2, function(r) stats::median(r[r > 0 & is.finite(r)]))
}


# pmax(<number>, log(x)) -> log(x) anywhere inside a function body.
.UnclampLog <- function(e) {
  if (!is.call(e)) return(e)
  if (identical(e[[1]], as.name("pmax")) && length(e) == 3 &&
      is.call(e[[3]]) && identical(e[[3]][[1]], as.name("log"))) {
    return(e[[3]])
  }
  for (i in seq_along(e)[-1]) {
    if (is.call(e[[i]])) e[[i]] <- .UnclampLog(e[[i]])
  }
  e
}


# APAdeg builds the GLMM offset as pmax(1e-8, log(size_factor)), so samples
# with a size factor < 1 get offset 0 (not normalised). This swaps it for
# log(size_factor) in APAdeg:::reshape_long_format for the current call.
# Returns the original function (to restore with on.exit) or NULL when the
# clamp is not there (e.g. fixed upstream).
.PatchAPAdegOffset <- function(verbose = TRUE) {
  orig  <- utils::getFromNamespace("reshape_long_format", "APAdeg")
  fixed <- orig
  body(fixed) <- .UnclampLog(body(orig))
  if (identical(body(fixed), body(orig))) {
    if (verbose) message("fix_offset: APAdeg's offset is not clamped in this version; nothing to fix.")
    return(NULL)
  }
  utils::assignInNamespace("reshape_long_format", fixed, ns = "APAdeg")
  orig
}


# Evaluate quietly (APAdegFit() always prints a progress bar).
.QuietOutput <- function(expr) {
  utils::capture.output(res <- expr)
  res
}


#' Gene-level differential expression from polyA-site counts (APAdeg)
#'
#' Runs APAdeg (Briefings in Bioinformatics, 2026;
#' \url{https://github.com/XuLabSJTU/APAdeg}) on every cell type x comparison
#' block of the polyA assay, as a complement to \code{\link{DEGsMatrix}()}
#' (RNA) and \code{\link{DEPsMatrix}()} (site usage).
#'
#' For each block:
#' \enumerate{
#'   \item pseudobulk polyA-site counts (cell type x group x replicate) and the
#'     same expression filter as \code{DEGsMatrix()} / \code{DEPsMatrix()}
#'     (\code{\link{PseudobulkFilter}()}), applied per site;
#'   \item sites are grouped by gene (\code{gene_col});
#'   \item \code{APAdeg::MetaData()} and \code{APAdeg::APAdegFit()} with the
#'     control as reference.
#' }
#'
#' How APAdeg works: per gene it fits a negative binomial GLMM
#' \code{count ~ group * is_major + (1 | apa_id)} with offset
#' \eqn{\log}(size factor) on the site counts plus an added "total gene" row.
#' \code{is_major} is the row with the most reads, in practice the total row,
#' so \code{log2FC} is the treatment-vs-control change in total gene counts;
#' the site rows mainly help estimate the gene's dispersion. Dispersions are
#' shrunk towards a trend across genes (needs hundreds of genes per block).
#' Sites with zero counts in a whole group (and >= 10 reads in the other) are
#' imputed with small values, so on/off sites give large but finite log2FC.
#'
#' Known issue in APAdeg: the offset is \code{pmax(1e-8, log(size factor))},
#' so samples with size factor < 1 are not normalised. \code{fix_offset = TRUE}
#' (default) corrects this during the call (APAdeg's internal function is
#' restored afterwards). \code{$blocks$offset_bias_log2} estimates the log2FC
#' bias the clamp would cause, and \code{$samples} lists the size factors.
#'
#' APAdeg fits two GLMMs per gene: about 0.4 s per gene on one core (about
#' 80 min for 10,000 genes). Use \code{cells} / \code{condition} to run only
#' some blocks, \code{parallel = TRUE}, and \code{min_sites = 2} (genes with
#' one site reduce to a standard NB test on the total).
#'
#' Requires APAdeg: \code{remotes::install_github("XuLabSJTU/APAdeg")}.
#'
#' @inheritParams PseudobulkFilter
#' @param comparisons \code{data.frame} with \code{treatment}, \code{control}
#'   and \code{Condition} (values of \code{label_col}), as in
#'   \code{DEGsMatrix()}. Extra columns are ignored, so the
#'   \code{DEPsMatrix()} table works too.
#' @param cells,condition Optional subset of blocks: cell types
#'   (\code{celltype_col} values) and \code{Condition} values. \code{NULL}
#'   (default) = all.
#' @param assay polyA assay. Default \code{"polyA"}.
#' @param genomic_positions \code{Intron.exon.location} values used as sites.
#'   Default \code{c("3' most exon", "Intron")}.
#' @param gene_col,symbol_col Gene id / symbol columns of the polyA assay
#'   \code{meta.features}. Defaults \code{"ensembl_gene_id"},
#'   \code{"Gene_Symbol"}.
#' @param filterByExpr_args Named list for \code{edgeR::filterByExpr()}
#'   (applied per site). Default \code{list(min.count = 1,
#'   min.total.count = 30)}, the adjusted filter used with
#'   \code{DEPsMatrix()} / PASTA, so all methods start from the same sites.
#'   \code{list()} = edgeR defaults (10 / 15).
#' @param min_sites Minimum number of kept sites per gene. Default 1 (every
#'   gene); 2 keeps the genes where the site structure adds information and
#'   roughly halves the time.
#' @param max_genes Only the \code{max_genes} most expressed genes per block
#'   (quick timing tests). \code{NULL} (default) = all. Keep >= 200: the
#'   dispersion trend needs many genes.
#' @param method GLMM engine for the final fit: \code{"lme4"} (default) or
#'   \code{"glmmTMB"}. With \code{"glmmTMB"} APAdeg reports the group
#'   coefficient (site rows) instead of the total-row contrast.
#' @param fitMethod Dispersion trend: \code{"NLS"} (default) or \code{"GAM"}.
#' @param Design Optional extra terms for the model, using the per-sample
#'   columns \code{replicate}, \code{n_cells} and \code{log_cells}, e.g.
#'   \code{"log_cells"} or \code{"1|replicate"}. Default \code{NULL}.
#' @param shrink log2FC shrinkage: \code{"none"} (default), \code{"apeGLMM"}
#'   or \code{"ash"}.
#' @param fix_offset Correct APAdeg's clamped offset (see Details). Default
#'   \code{TRUE}.
#' @param parallel Logical. Fit genes in parallel (BiocParallel; forking on
#'   macOS/Linux). Default \code{FALSE}.
#' @param ncores Workers when \code{parallel = TRUE}. \code{NULL} (default) =
#'   all cores but one.
#' @param BPPARAM Optional BiocParallel back-end (overrides
#'   \code{parallel}/\code{ncores}).
#' @param fit If \code{FALSE}, only build the pseudobulk and apply the filter
#'   (no APAdeg, seconds): \code{$blocks} gives the sites and genes each block
#'   would test, to compare \code{filterByExpr_args} settings before a long
#'   run. Default \code{TRUE}.
#' @param return_filter If \code{TRUE}, also return the
#'   \code{\link{PseudobulkFilter}()} object (for \code{\link{QCDensity}()}).
#'   Default \code{FALSE}.
#' @param verbose Print progress. Default \code{TRUE}.
#'
#' @return A list with
#'   \describe{
#'     \item{\code{de}}{one row per block x gene: \code{Cells},
#'       \code{Condition}, \code{group1} (control), \code{group2} (treatment),
#'       \code{gene}, \code{ensembl_gene_id}, \code{n_sites}, \code{baseMean},
#'       \code{log2FC} (treatment vs control), \code{pvalue}, \code{padj},
#'       \code{regulation}, \code{n_tested}}
#'     \item{\code{blocks}}{one row per block: samples, sites, genes, genes
#'       with padj < 0.05, run time, \code{offset_bias_log2},
#'       \code{offset_fixed}}
#'     \item{\code{samples}}{pseudobulk samples with their size factor and
#'       \code{offset_clamped} (not normalised by unpatched APAdeg)}
#'     \item{\code{disp}}{APAdeg dispersion table per block}
#'     \item{\code{filter}}{the \code{PseudobulkFilter()} object, when
#'       \code{return_filter = TRUE} or \code{fit = FALSE}}
#'   }
#'
#' @seealso \code{\link{CompareDEGs}()}, \code{\link{CompareDEGsPlot}()}
#'
#' @examples
#' \dontrun{
#' seurat_polyA$contrast <- paste0(seurat_polyA$Strain, "_", seurat_polyA$Treatment)
#' comparisons <- data.frame(treatment = "3xTg_Alcohol", control = "B6_Control",
#'                           Condition = "3xTg_Alcohol_vs_B6_Control")
#'
#' apadegs <- APAdegsMatrix(seurat_polyA, comparisons,
#'                          celltype_col      = "Subpopulation",
#'                          label_col         = "contrast",
#'                          replicate_col     = "Number",
#'                          cells             = "L2/3 IT",
#'                          parallel          = TRUE,
#'                          ncores            = 8)
#' apadegs$blocks
#' head(apadegs$de)
#'
#' # Compare filters first (no APAdeg, seconds): sites / genes per block
#' APAdegsMatrix(seurat_polyA, comparisons, "Subpopulation", "contrast", "Number",
#'               fit = FALSE)$blocks                                   # default (1 / 30)
#' APAdegsMatrix(seurat_polyA, comparisons, "Subpopulation", "contrast", "Number",
#'               filterByExpr_args = list(min.count = 5, min.total.count = 30),
#'               fit = FALSE)$blocks
#' APAdegsMatrix(seurat_polyA, comparisons, "Subpopulation", "contrast", "Number",
#'               filterByExpr_args = list(), fit = FALSE)$blocks       # edgeR defaults
#' }
#'
#' @export
APAdegsMatrix <- function(seu,
                          comparisons,
                          celltype_col,
                          label_col,
                          replicate_col,
                          cells             = NULL,
                          condition         = NULL,
                          assay             = "polyA",
                          genomic_positions = c("3' most exon", "Intron"),
                          gene_col          = "ensembl_gene_id",
                          symbol_col        = "Gene_Symbol",
                          min_cells         = 10,
                          filter_method     = c("edgeR", "manual", "none"),
                          mincounts         = 10,
                          filterByExpr_args = list(min.count = 1, min.total.count = 30),
                          filter_scope      = c("comparison", "celltype"),
                          min_sites         = 1,
                          max_genes         = NULL,
                          method            = c("lme4", "glmmTMB"),
                          fitMethod         = c("NLS", "GAM"),
                          Design            = NULL,
                          shrink            = c("none", "apeGLMM", "ash"),
                          fix_offset        = TRUE,
                          parallel          = FALSE,
                          ncores            = NULL,
                          BPPARAM           = NULL,
                          fit               = TRUE,
                          return_filter     = FALSE,
                          verbose           = TRUE) {
  if (isTRUE(fit) && !requireNamespace("APAdeg", quietly = TRUE)) {
    stop("Package 'APAdeg' is required: remotes::install_github('XuLabSJTU/APAdeg')")
  }
  filter_method <- match.arg(filter_method)
  filter_scope  <- match.arg(filter_scope)
  method        <- match.arg(method)
  fitMethod     <- match.arg(fitMethod)
  shrink        <- match.arg(shrink)
  if (!assay %in% SeuratObject::Assays(seu)) stop("Assay '", assay, "' not found in seu.")

  # ---- sites with a gene, in genomic_positions --------------------------------
  mf <- seu[[assay]]@meta.features
  miss <- setdiff(c(gene_col, symbol_col, "Intron.exon.location"), colnames(mf))
  if (length(miss)) {
    stop("Column(s) not in the ", assay, " meta.features: ", paste(miss, collapse = ", "),
         ". Annotate the assay first (AnnotatePolyADb / PolyAPipeline).")
  }
  ok <- mf$Intron.exon.location %in% genomic_positions &
    !is.na(mf[[gene_col]]) & mf[[gene_col]] != ""
  sites     <- rownames(mf)[ok]
  site_gene <- stats::setNames(as.character(mf[[gene_col]]), rownames(mf))
  gene_sym  <- stats::setNames(as.character(mf[[symbol_col]]), as.character(mf[[gene_col]]))
  if (!length(sites)) {
    stop("No site with a gene in Intron.exon.location {",
         paste(genomic_positions, collapse = ", "), "}.")
  }
  if (verbose) {
    message(length(sites), " polyA sites with a gene in {",
            paste(genomic_positions, collapse = ", "), "}")
  }

  # ---- pseudobulk + expression filter (same as DEGsMatrix / DEPsMatrix) ------
  flt <- PseudobulkFilter(
    seu, comparisons, celltype_col, label_col, replicate_col,
    assay             = assay,
    features          = sites,
    min_cells         = min_cells,
    filter_method     = filter_method,
    mincounts         = mincounts,
    filterByExpr_args = filterByExpr_args,
    filter_scope      = filter_scope,
    feature_type      = "peak",
    verbose           = FALSE
  )

  blocks <- flt$summary
  if (!is.null(cells))     blocks <- blocks[blocks$Cells %in% cells, , drop = FALSE]
  if (!is.null(condition)) blocks <- blocks[blocks$Condition %in% condition, , drop = FALSE]
  if (!nrow(blocks)) {
    stop("No block left. Cells: ", paste(unique(flt$summary$Cells), collapse = ", "),
         " | Conditions: ", paste(unique(flt$summary$Condition), collapse = ", "))
  }
  skipped <- blocks[blocks$status != "ok" | is.na(blocks$n_kept) | blocks$n_kept == 0, , drop = FALSE]
  if (verbose && nrow(skipped)) {
    message("Skipped ", nrow(skipped), " block(s): ",
            paste0(skipped$Cells, " (", skipped$Condition, ", ", skipped$status, ")",
                   collapse = "; "))
  }
  blocks <- blocks[blocks$status == "ok" & !is.na(blocks$n_kept) & blocks$n_kept > 0, , drop = FALSE]
  if (!nrow(blocks)) stop("No block passed min_cells / the expression filter.")

  # ---- parallel back-end (passed to APAdegFit) --------------------------------
  if (isTRUE(parallel) && is.null(BPPARAM)) {
    if (is.null(ncores)) ncores <- max(1L, parallel::detectCores() - 1L)
    BPPARAM <- if (.Platform$OS.type == "windows") {
      BiocParallel::SnowParam(workers = ncores)
    } else {
      BiocParallel::MulticoreParam(workers = ncores)
    }
  }
  use_bp <- !is.null(BPPARAM)
  if (!use_bp) BPPARAM <- BiocParallel::SerialParam()

  # ---- APAdeg offset fix (restored on exit) -----------------------------------
  offset_fixed <- FALSE
  if (isTRUE(fit) && isTRUE(fix_offset)) {
    orig <- .PatchAPAdegOffset(verbose = verbose)
    if (!is.null(orig)) {
      offset_fixed <- TRUE
      on.exit(utils::assignInNamespace("reshape_long_format", orig, ns = "APAdeg"), add = TRUE)
    }
  }

  if (verbose && isTRUE(fit)) {
    message("Running APAdeg on ", nrow(blocks), " block(s)",
            if (use_bp) paste0(" (", BiocParallel::bpnworkers(BPPARAM), " workers)"), "...")
  }

  # ---- one APAdeg fit per block -----------------------------------------------
  out <- lapply(seq_len(nrow(blocks)), function(i) {
    b <- blocks[i, ]

    s <- flt$sample_meta %>%
      dplyr::filter(Cells == b$Cells, label %in% c(b$control, b$treatment)) %>%
      dplyr::mutate(group = ifelse(label == b$control, "control", "treatment")) %>%
      dplyr::arrange(group, replicate) %>%
      dplyr::mutate(sample = make.unique(make.names(
        paste0(ifelse(group == "control", "ctl_", "trt_"), replicate))))
    is_c <- s$group == "control"

    peaks <- flt$kept$feature[flt$kept$Cells == b$Cells & flt$kept$Condition == b$Condition]
    x <- as.matrix(flt$pseudobulk[peaks, s$sample_id, drop = FALSE])
    colnames(x) <- s$sample

    site_tbl <- data.frame(peak = peaks, gene_id = unname(site_gene[peaks]),
                           stringsAsFactors = FALSE) %>%
      dplyr::add_count(gene_id, name = "n_sites") %>%
      dplyr::filter(n_sites >= min_sites)
    if (!nrow(site_tbl)) {
      warning(sprintf("%s | %s: no gene with >= %d kept sites.",
                      b$Condition, b$Cells, min_sites), call. = FALSE)
      return(NULL)
    }

    gsum  <- rowsum(x[site_tbl$peak, , drop = FALSE], site_tbl$gene_id)
    n_all <- nrow(gsum)
    if (!is.null(max_genes) && n_all > max_genes) {
      top      <- names(sort(rowSums(gsum), decreasing = TRUE))[seq_len(max_genes)]
      site_tbl <- site_tbl[site_tbl$gene_id %in% top, , drop = FALSE]
      gsum     <- gsum[top, , drop = FALSE]
    }

    # fit = FALSE: filter only, report what the block would test
    if (!isTRUE(fit)) {
      return(list(block = data.frame(
        Cells = b$Cells, Condition = b$Condition,
        n_ctl = sum(is_c), n_trt = sum(!is_c),
        n_sites_in = b$n_in, n_sites_kept = length(peaks),
        n_sites = nrow(site_tbl), n_genes = nrow(gsum), n_genes_block = n_all,
        n_genes_1site = sum(table(site_tbl$gene_id) == 1),
        median_gene_reads = stats::median(rowSums(gsum)),
        status = "filter only",
        stringsAsFactors = FALSE
      )))
    }

    # size factors APAdeg will use, and the log2FC bias of the clamped offset
    sf   <- .APAdegSizeFactors(gsum)[s$sample]
    bias <- (mean(pmin(0, log(sf[!is_c]))) - mean(pmin(0, log(sf[is_c])))) / log(2)

    apa_in <- data.frame(gene_id = site_tbl$gene_id,
                         x[site_tbl$peak, , drop = FALSE],
                         check.names = FALSE, stringsAsFactors = FALSE, row.names = NULL)
    md <- APAdeg::MetaData(APAseq = apa_in, group1 = "control", group2 = "treatment",
                           n_group1 = sum(is_c), n_group2 = sum(!is_c))
    pheno <- data.frame(sample = s$sample, replicate = as.character(s$replicate),
                        n_cells = s$n_cells, log_cells = log(s$n_cells),
                        stringsAsFactors = FALSE)

    if (verbose) {
      message(sprintf("[%d/%d] %s | %s: %d genes, %d sites, %d vs %d samples",
                      i, nrow(blocks), b$Condition, b$Cells, nrow(gsum), nrow(site_tbl),
                      sum(is_c), sum(!is_c)))
    }

    t0 <- Sys.time()
    run_fit <- function() {
      APAdeg::APAdegFit(md,
                        Reference  = "control",
                        phenoData  = pheno,
                        Design     = Design,
                        method     = method,
                        fitMethod  = fitMethod,
                        Parallel   = use_bp,
                        BPPARAM    = BPPARAM,
                        Shrinklfc  = shrink != "none",
                        ShrinkType = shrink)
    }
    apa_fit <- tryCatch(
      if (verbose) run_fit() else .QuietOutput(run_fit()),
      error = function(e) {
        warning(sprintf("APAdeg failed for %s | %s: %s",
                        b$Condition, b$Cells, conditionMessage(e)), call. = FALSE)
        NULL
      }
    )
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    de <- disp <- NULL
    if (!is.null(apa_fit)) {
      de <- as.data.frame(apa_fit@APAdegResult) %>%
        dplyr::left_join(dplyr::distinct(site_tbl, gene_id, n_sites), by = "gene_id") %>%
        dplyr::transmute(Cells = b$Cells, Condition = b$Condition,
                         group1 = b$control, group2 = b$treatment,
                         gene = unname(gene_sym[gene_id]), ensembl_gene_id = gene_id,
                         n_sites, baseMean, log2FC, pvalue, padj,
                         regulation = ifelse(log2FC > 0, "Up", "Down"),
                         n_tested = nrow(gsum))
      disp <- as.data.frame(apa_fit@APAdegDispData) %>%
        dplyr::mutate(Cells = b$Cells, Condition = b$Condition, .before = 1)
    }

    list(
      de   = de,
      disp = disp,
      samples = data.frame(
        Cells = b$Cells, Condition = b$Condition, sample = s$sample, group = s$group,
        label = s$label, replicate = s$replicate, n_cells = s$n_cells,
        size_factor = round(unname(sf), 3),
        offset_clamped = !offset_fixed & unname(sf) < 1,
        stringsAsFactors = FALSE
      ),
      block = data.frame(
        Cells = b$Cells, Condition = b$Condition,
        n_ctl = sum(is_c), n_trt = sum(!is_c),
        n_sites = nrow(site_tbl), n_genes = nrow(gsum), n_genes_block = n_all,
        n_fitted = if (is.null(de)) 0L else sum(!is.na(de$log2FC)),
        n_padj05 = if (is.null(de)) 0L else sum(de$padj < 0.05, na.rm = TRUE),
        minutes  = round(secs / 60, 1),
        offset_bias_log2 = round(bias, 3),
        offset_fixed     = offset_fixed,
        status   = if (is.null(apa_fit)) "APAdeg failed" else "ok",
        stringsAsFactors = FALSE
      )
    )
  })

  if (!isTRUE(fit)) {
    return(list(blocks = dplyr::bind_rows(lapply(out, `[[`, "block")), filter = flt))
  }

  res <- list(
    de      = tibble::as_tibble(dplyr::bind_rows(lapply(out, `[[`, "de"))),
    blocks  = dplyr::bind_rows(lapply(out, `[[`, "block")),
    samples = dplyr::bind_rows(lapply(out, `[[`, "samples")),
    disp    = dplyr::bind_rows(lapply(out, `[[`, "disp"))
  )
  if (isTRUE(return_filter)) res$filter <- flt
  res
}


#' Compare APAdeg with DEGsMatrix (and DEPsMatrix) per block and gene
#'
#' Joins \code{\link{APAdegsMatrix}()} with \code{\link{DEGsMatrix}()} on the
#' blocks APAdeg was run on, classifies every gene, and optionally flags genes
#' with significant APA from \code{\link{DEPsMatrix}()}.
#'
#' Cell types are matched with \code{make.names()} (\code{"L2/3 IT"} =
#' \code{"L2.3.IT"}); genes by symbol, or by Ensembl id when the DEGsMatrix
#' genes are Ensembl ids.
#'
#' @param apadegs \code{APAdegsMatrix()} output (or its \code{$de}).
#' @param degs \code{DEGsMatrix()} output (or its \code{$de}).
#' @param deps Optional \code{DEPsMatrix()} output (or its
#'   \code{$red_scores}): adds \code{RED_sig} and \code{RED_min_padj} per gene.
#' @param alpha Significance threshold. Default 0.05.
#' @param sig_by \code{"padj"} (default, BH-adjusted) or \code{"pvalue"} (raw
#'   p-value), used for both methods.
#' @param deps_padj RED p-value column used for \code{RED_sig}. Default
#'   \code{"RED_LRT_padj"}.
#'
#' @return A list with \code{table} (one row per block x gene: log2FC,
#'   p-values and padj of both methods, \code{tested} = both / DEGsMatrix only
#'   / APAdeg only, \code{sig_DEG}, \code{sig_APAdeg}, \code{class} = Both,
#'   same direction / Both, opposite direction / DEGsMatrix only / APAdeg only
#'   / Neither, and the RED columns), \code{summary} (genes per class and the
#'   Spearman correlation of log2FC per block), \code{sig_by}, \code{alpha}.
#'
#' @examples
#' \dontrun{
#' cmp <- CompareDEGs(apadegs, degs, deps = deps)
#' cmp$summary
#' CompareDEGsPlot(cmp, sig_only = TRUE)
#' }
#'
#' @export
CompareDEGs <- function(apadegs, degs, deps = NULL, alpha = 0.05,
                        sig_by = c("padj", "pvalue"),
                        deps_padj = "RED_LRT_padj") {
  sig_by <- match.arg(sig_by)
  ap <- if (is.data.frame(apadegs)) apadegs else apadegs$de
  dg <- if (is.data.frame(degs)) degs else degs$de
  if (!nrow(ap)) stop("`apadegs` has no rows.")
  if (!nrow(dg)) stop("`degs` has no rows.")

  by_ens <- mean(grepl("^ENS", dg$gene)) > 0.5      # DEGsMatrix genes = Ensembl ids?

  ap <- ap %>%
    dplyr::mutate(key_ct = make.names(Cells),
                  key_gene = if (by_ens) ensembl_gene_id else gene) %>%
    dplyr::filter(!is.na(key_gene), key_gene != "") %>%
    dplyr::arrange(padj) %>%
    dplyr::distinct(key_ct, Condition, key_gene, .keep_all = TRUE) %>%
    dplyr::transmute(key_ct, Condition, key_gene, Cells, gene, ensembl_gene_id, n_sites,
                     logFC_APAdeg = log2FC, pval_APAdeg = pvalue, padj_APAdeg = padj)

  dg_all <- dg %>%
    dplyr::transmute(key_ct = make.names(as.character(cell_type)),
                     Condition = as.character(Condition),
                     key_gene = gene, Cells_DEG = as.character(cell_type),
                     logFC_DEG = avg_logFC, pval_DEG = p_val, padj_DEG = p_val_adj)
  dg <- dg_all %>%
    dplyr::semi_join(dplyr::distinct(ap, key_ct, Condition), by = c("key_ct", "Condition"))

  if (!nrow(dg)) {
    stop("No DEGsMatrix rows for the APAdeg blocks.\n  APAdeg:     ",
         paste(unique(paste(ap$Cells, "|", ap$Condition)), collapse = "; "),
         "\n  DEGsMatrix: ",
         paste(utils::head(unique(paste(dg_all$Cells_DEG, "|", dg_all$Condition)), 6),
               collapse = "; "), call. = FALSE)
  }
  if (!nrow(dplyr::semi_join(dg, ap, by = c("key_ct", "Condition", "key_gene")))) {
    stop("Same blocks but no shared genes.\n  DEGsMatrix genes: ",
         paste(utils::head(dg$key_gene, 4), collapse = ", "),
         "\n  APAdeg genes:     ", paste(utils::head(ap$key_gene, 4), collapse = ", "),
         call. = FALSE)
  }

  tbl <- dplyr::full_join(dg, ap, by = c("key_ct", "Condition", "key_gene")) %>%
    dplyr::mutate(
      Cells = dplyr::coalesce(Cells, Cells_DEG),
      gene  = dplyr::coalesce(gene, key_gene),
      tested = dplyr::case_when(
        is.na(logFC_APAdeg) ~ "DEGsMatrix only",
        is.na(logFC_DEG)    ~ "APAdeg only",
        TRUE                ~ "both"),
      sig_DEG    = dplyr::coalesce(if (sig_by == "padj") padj_DEG < alpha else pval_DEG < alpha, FALSE),
      sig_APAdeg = dplyr::coalesce(if (sig_by == "padj") padj_APAdeg < alpha else pval_APAdeg < alpha, FALSE),
      class = dplyr::case_when(
        sig_DEG & sig_APAdeg & sign(logFC_DEG) == sign(logFC_APAdeg) ~ "Both, same direction",
        sig_DEG & sig_APAdeg ~ "Both, opposite direction",
        sig_DEG              ~ "DEGsMatrix only",
        sig_APAdeg           ~ "APAdeg only",
        TRUE                 ~ "Neither")
    )

  if (!is.null(deps)) {
    red <- if (is.data.frame(deps)) deps else deps$red_scores
    if (!deps_padj %in% colnames(red)) stop("`", deps_padj, "` not in the RED table.")
    red_gene <- red %>%
      dplyr::filter(!is.na(.data[[deps_padj]])) %>%
      dplyr::group_by(key_ct = make.names(as.character(Cells)), Condition,
                      key_gene = if (by_ens) Ensembl_ID else Gene_symbol) %>%
      dplyr::summarise(RED_min_padj = min(.data[[deps_padj]]),
                       RED_sig      = any(.data[[deps_padj]] < alpha),
                       .groups = "drop")
    tbl <- tbl %>%
      dplyr::left_join(red_gene, by = c("key_ct", "Condition", "key_gene")) %>%
      dplyr::mutate(RED_sig = dplyr::coalesce(RED_sig, FALSE))
  }

  tbl <- tbl %>%
    dplyr::select(-key_ct, -key_gene, -Cells_DEG) %>%
    dplyr::relocate(Cells, Condition, gene, ensembl_gene_id, tested, class)

  .spearman <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3) NA_real_ else round(stats::cor(x[ok], y[ok], method = "spearman"), 3)
  }

  summ <- tbl %>%
    dplyr::count(Cells, Condition, class) %>%
    tidyr::pivot_wider(names_from = class, values_from = n, values_fill = 0) %>%
    dplyr::left_join(
      tbl %>%
        dplyr::filter(tested == "both") %>%
        dplyr::group_by(Cells, Condition) %>%
        dplyr::summarise(n_tested_both   = dplyr::n(),
                         spearman_log2FC = .spearman(logFC_DEG, logFC_APAdeg),
                         .groups = "drop"),
      by = c("Cells", "Condition"))

  list(table = tbl, summary = summ, sig_by = sig_by, alpha = alpha)
}


#' log2FC of DEGsMatrix vs APAdeg
#'
#' Scatter of the log2FC of both methods for the genes tested by both,
#' coloured by \code{class} (\code{\link{CompareDEGs}()}), one panel per cell
#' type x comparison.
#'
#' @param cmp \code{CompareDEGs()} output (or its \code{$table}).
#' @param sig_only If \code{TRUE}, only genes significant in at least one
#'   method, coloured Both / DEGsMatrix only / APAdeg only; genes significant
#'   but tested by one method only are counted in the subtitle. Default
#'   \code{FALSE} (all genes, five classes).
#' @param cells,condition Optional subset of panels.
#' @param label_top Label the \code{label_top} genes significant in both with
#'   the largest |log2FC| per panel (0 = no labels). Default 10.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' CompareDEGsPlot(cmp)
#' CompareDEGsPlot(cmp, sig_only = TRUE)
#' }
#'
#' @export
CompareDEGsPlot <- function(cmp, sig_only = FALSE, cells = NULL, condition = NULL,
                            label_top = 10) {
  d <- if (is.data.frame(cmp)) cmp else cmp$table
  rule <- if (is.list(cmp) && !is.data.frame(cmp) && !is.null(cmp$sig_by)) {
    paste0(" (", cmp$sig_by, " < ", cmp$alpha, ")")
  } else ""
  if (!is.null(cells))     d <- d[d$Cells %in% cells, , drop = FALSE]
  if (!is.null(condition)) d <- d[d$Condition %in% condition, , drop = FALSE]

  sub <- NULL
  if (isTRUE(sig_only)) {
    n_one <- d %>% dplyr::filter(tested != "both", sig_DEG | sig_APAdeg) %>% dplyr::count(tested)
    if (nrow(n_one)) {
      sub <- paste0("Significant but tested by one method only (not shown): ",
                    paste0(n_one$tested, " = ", n_one$n, collapse = ", "))
    }
    d <- d %>%
      dplyr::filter(tested == "both", sig_DEG | sig_APAdeg) %>%
      dplyr::mutate(group = dplyr::case_when(sig_DEG & sig_APAdeg ~ "Both",
                                             sig_DEG              ~ "DEGsMatrix only",
                                             TRUE                 ~ "APAdeg only"))
    lev <- c("Both", "DEGsMatrix only", "APAdeg only")
    pal <- c("#7c3aed", "#2563eb", "#ea580c")
  } else {
    d <- d %>% dplyr::filter(tested == "both") %>% dplyr::mutate(group = class)
    lev <- c("Both, same direction", "Both, opposite direction",
             "DEGsMatrix only", "APAdeg only", "Neither")
    pal <- c("#7c3aed", "#dc2626", "#2563eb", "#ea580c", "grey80")
  }
  if (!nrow(d)) stop("No genes to plot.")

  n   <- table(factor(d$group, levels = lev))
  lab <- stats::setNames(paste0(lev, " (", n, ")"), lev)
  d <- d %>%
    dplyr::mutate(group = factor(lab[group], levels = lab)) %>%
    dplyr::arrange(dplyr::desc(group))                  # significant drawn on top

  top <- d %>%
    dplyr::filter(sig_DEG, sig_APAdeg) %>%
    dplyr::group_by(Cells, Condition) %>%
    dplyr::slice_max(order_by = abs(logFC_DEG) + abs(logFC_APAdeg), n = label_top,
                     with_ties = FALSE) %>%
    dplyr::ungroup()

  ggplot2::ggplot(d, ggplot2::aes(logFC_DEG, logFC_APAdeg, colour = group)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey80") +
    ggplot2::geom_vline(xintercept = 0, colour = "grey80") +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
    ggplot2::geom_point(size = if (isTRUE(sig_only)) 1.5 else 0.9, alpha = 0.8) +
    ggrepel::geom_text_repel(data = top, ggplot2::aes(label = gene), size = 3,
                             colour = "black", max.overlaps = Inf, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = stats::setNames(pal, lab), drop = FALSE) +
    ggplot2::facet_grid(Condition ~ Cells) +
    ggplot2::labs(x = "log2FC  DEGsMatrix (RNA)", y = "log2FC  APAdeg (polyA sites)",
                  colour = paste0("Significant in", rule), subtitle = sub) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")
}
