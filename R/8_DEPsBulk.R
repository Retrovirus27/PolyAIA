# =============================================================================
#  Differential polyA site usage (RED) for bulk 3' RNA-seq
#  Same RED core as DEPsMatrix(); the input is a sites x samples count matrix
#  (e.g. polyApipe_bulk.py) instead of a single-cell Seurat object.
# =============================================================================


#' Differential polyA site usage (RED) for bulk 3' RNA-seq
#'
#' Bulk counterpart of \code{\link{DEPsMatrix}()}: same RED scores, pair
#' choice, DPAU/gDPAU and tests (\code{RED_LRT}, \code{RED_QL},
#' \code{RED_fisher}, \code{RED_ttest}, with covariates, paired designs and
#' interactions), from a polyA sites x samples count matrix (e.g. the
#' \code{*_counts_matrix.tsv.gz} of \code{polyApipe_bulk.py}). Each sample is
#' one replicate; there are no cells, cell types or PASTA residuals.
#'
#' Steps:
#' \enumerate{
#'   \item sites annotated with polyA_DB (same rule as
#'     \code{\link{AnnotatePolyADb}()}: the site's position, nearest polyA_DB
#'     site on the same strand within \code{max.dist}); sites kept if their
#'     \code{Intron.exon.location} is in \code{genomic_positions};
#'   \item optional removal of mispriming / internal-priming sites and
#'     RepeatMasker (TE) handling, as in \code{DEPsMatrix()};
#'   \item per comparison: expression filter, proximal/distal pairs, RED and
#'     tests (the \code{DEPsMatrix()} RED core).
#' }
#'
#' The function does not interpret sample names: \code{samples} says which
#' column of \code{counts} is which sample and gives its design (condition,
#' line/donor, covariates). Build it in your own script from the names of your
#' dataset.
#'
#' @param counts Sites x samples counts: a \code{data.frame}/matrix with a
#'   \code{peak} column (or rownames) and one column per sample, or the path
#'   to such a tab-separated file (\code{.tsv} or \code{.tsv.gz}).
#' @param sites Site coordinates: a \code{data.frame} or path with columns
#'   \code{peak}, \code{chrom}, \code{site} (1-based position of the polyA
#'   site) and \code{strand}, e.g. the \code{*_polyA_peaks_info.tsv.gz} of
#'   \code{polyApipe_bulk.py}. Optional \code{misprime} /
#'   \code{internal_priming} columns are used by \code{exclude_misprime} /
#'   \code{exclude_internal_priming}.
#' @param samples \code{data.frame} with a \code{sample} column equal to the
#'   column names of \code{counts}, the \code{condition_col} column and any
#'   design column used in \code{covariate_formula} / \code{interaction}
#'   (e.g. \code{line}, \code{genotype}). Samples of \code{counts} not listed
#'   are not used.
#' @param comparisons \code{data.frame} with \code{treatment}, \code{control}
#'   (values of \code{condition_col}; several joined by \code{"+"} to pool
#'   them, e.g. for an \code{interaction}) and \code{Condition} (label). Extra
#'   columns are stamped onto the output rows.
#' @param condition_col Column of \code{samples} holding the groups compared.
#' @param polyAdb.file,polyAdb polyA_DB table (path, or the table itself), same
#'   genome build as the counts.
#' @param max.dist Maximum distance (nt) between a site and its polyA_DB site.
#'   Default 50.
#' @param genome_build \code{"hg38"}, \code{"hg19"} or \code{"mm10"}: gene
#'   start/end from the EnsDb of that build (hg38 = v86, hg19 = v75, mm10 =
#'   v79), needed for \code{RLDi} (intronic RLD). \code{NULL} (default): no
#'   gene coordinates, \code{RLDi} is \code{NA}; everything else is computed.
#' @param genomic_positions \code{Intron.exon.location} values kept. Default
#'   \code{c("3' most exon", "Intron")}.
#' @param exclude_misprime,exclude_internal_priming Drop sites flagged in the
#'   \code{misprime} / \code{internal_priming} column of \code{sites}. Default
#'   \code{FALSE}.
#' @param repeat_masker,rmsk,repeat_masker_positions TE handling, as in
#'   \code{DEPsMatrix()}: \code{"none"} (default), \code{"Remove"} or
#'   \code{"Keep"} (flag), with a UCSC RepeatMasker table \code{rmsk}, for the
#'   sites whose location is in \code{repeat_masker_positions} (default
#'   \code{"Intron"}).
#' @param filter_method,mincounts,filterByExpr_args Expression filter of the
#'   sites, as in \code{DEPsMatrix()}.
#' @param filter_scope \code{"comparison"} (default): filter on each
#'   comparison's own samples. \code{"all"}: one filter on the samples of all
#'   groups in \code{comparisons}, the same sites for every comparison.
#' @param covariate_formula,interaction Design of \code{RED_LRT} /
#'   \code{RED_QL}, as in \code{DEPsMatrix()}, using columns of
#'   \code{samples} (e.g. \code{covariate_formula = "line"} for a paired
#'   design, \code{interaction = "genotype"}).
#' @param pair_rule,red_tests,gdpau_min_reads,min_site_reads,min_site_usage,min_rep_pairs,drop_na_cols
#'   As in \code{DEPsMatrix()}.
#' @param return_filter If \code{TRUE}, also return the per-comparison filter
#'   (for \code{\link{QCDensity}()}). Default \code{FALSE}.
#' @param verbose Print progress. Default \code{TRUE}.
#'
#' @return A list with \code{red_scores} (same columns as
#'   \code{DEPsMatrix()}; \code{Cells = "bulk"}), \code{polyAdb} (the analyzed
#'   sites with their polyA_DB annotation, \code{match_distance} and
#'   \code{Repeated_Masker}), \code{sites} (all input sites with their match,
#'   or \code{NA}), \code{samples} (sample, the key used inside, group) and,
#'   with \code{return_filter = TRUE}, \code{filter}.
#'
#' @seealso \code{\link{DEPsMatrix}()}
#'
#' @examples
#' \dontrun{
#' mat <- read.delim("PolyA_counts_matrix.tsv", check.names = FALSE)
#'
#' # sample design, built from this dataset's names (GSM..._<line>_<infection>)
#' samples <- data.frame(sample = colnames(mat)[-1]) |>
#'   tidyr::separate(sample, into = c("gsm", "line", "infection"), sep = "_",
#'                   remove = FALSE, extra = "merge") |>
#'   dplyr::mutate(genotype = sub("^f(WT|AD)[0-9]+$", "\\1", line),
#'                 group    = paste(genotype, infection, sep = "_"))
#'
#' comparisons <- data.frame(
#'   treatment = c("WT_TBEV", "WT_HSV-1", "AD_TBEV", "AD_HSV-1"),
#'   control   = c("WT_mock", "WT_mock",  "AD_mock", "AD_mock"),
#'   Condition = c("WT TBEV", "WT HSV-1", "AD TBEV", "AD HSV-1")
#' )
#'
#' deps_bulk <- DEPsMatrixBulk(
#'   counts            = mat,
#'   sites             = "PolyA_polyA_peaks_info.tsv",
#'   samples           = samples,
#'   comparisons       = comparisons,
#'   condition_col     = "group",
#'   polyAdb.file      = "polyAdb_hg38.txt",
#'   genome_build      = "hg38",
#'   covariate_formula = "line"          # paired: each line mock vs infected
#' )
#'
#' # does the infection effect differ between WT and AD?
#' inter <- data.frame(treatment = "WT_TBEV+AD_TBEV", control = "WT_mock+AD_mock",
#'                     Condition = "TBEV: AD vs WT")
#' deps_int <- DEPsMatrixBulk(mat, "PolyA_polyA_peaks_info.tsv", samples, inter,
#'                            condition_col = "group", polyAdb.file = "polyAdb_hg38.txt",
#'                            covariate_formula = "line", interaction = "genotype")
#' }
#'
#' @export
DEPsMatrixBulk <- function(
    counts,
    sites,
    samples,
    comparisons,
    condition_col,
    polyAdb.file      = NULL,
    polyAdb           = NULL,
    max.dist          = 50,
    genome_build      = NULL,
    genomic_positions = c("3' most exon", "Intron"),
    exclude_misprime  = FALSE,
    exclude_internal_priming = FALSE,
    repeat_masker     = c("none", "Remove", "Keep"),
    rmsk              = NULL,
    repeat_masker_positions = "Intron",
    filter_method     = c("edgeR", "manual", "none"),
    filter_scope      = c("comparison", "all"),
    mincounts         = 10,
    filterByExpr_args = list(),
    covariate_formula = NULL,
    interaction       = NULL,
    pair_rule         = c("max", "top"),
    red_tests         = c("LRT", "QL", "fisher", "ttest"),
    gdpau_min_reads   = 5,
    min_site_reads    = 10,
    min_site_usage    = 0.05,
    min_rep_pairs     = 0.9,
    drop_na_cols      = TRUE,
    return_filter     = FALSE,
    verbose           = TRUE
) {
  repeat_masker <- match.arg(repeat_masker)
  filter_method <- match.arg(filter_method)
  filter_scope  <- match.arg(filter_scope)
  pair_rule     <- match.arg(pair_rule)
  red_tests     <- match.arg(red_tests, several.ok = TRUE)
  ct <- "bulk"

  .read_tab <- function(x, what) {
    if (is.character(x) && length(x) == 1) {
      if (!file.exists(x)) stop("`", what, "` file not found: ", x)
      x <- utils::read.delim(x, check.names = FALSE, stringsAsFactors = FALSE)
    }
    as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
  }

  # ---- counts --------------------------------------------------------------
  counts <- .read_tab(counts, "counts")
  if ("peak" %in% colnames(counts)) {
    rownames(counts) <- as.character(counts$peak)
    counts$peak <- NULL
  }
  if (is.null(rownames(counts)) || all(rownames(counts) == as.character(seq_len(nrow(counts))))) {
    stop("`counts` needs a `peak` column or site names as rownames.")
  }

  # ---- samples -------------------------------------------------------------
  samples <- as.data.frame(samples, stringsAsFactors = FALSE)
  if (!"sample" %in% colnames(samples)) stop("`samples` needs a `sample` column.")
  if (!condition_col %in% colnames(samples)) {
    stop("`condition_col` '", condition_col, "' is not a column of `samples`.")
  }
  samples$sample <- as.character(samples$sample)
  if (anyDuplicated(samples$sample)) {
    stop("Duplicated names in samples$sample: ",
         paste(unique(samples$sample[duplicated(samples$sample)]), collapse = ", "))
  }
  miss_s <- setdiff(samples$sample, colnames(counts))
  if (length(miss_s)) {
    stop("samples$sample not found in the columns of `counts`: ", paste(miss_s, collapse = ", "),
         ". Columns of counts: ", paste(utils::head(colnames(counts), 5), collapse = ", "), " ...")
  }
  unused <- setdiff(colnames(counts), samples$sample)
  if (verbose && length(unused)) {
    message("Not in `samples`, not used: ", paste(unused, collapse = ", "))
  }
  # Key used inside (last "_" field of the internal column names): no "_".
  samples$key <- make.unique(gsub("_", "-", samples$sample, fixed = TRUE))
  groups <- as.character(samples[[condition_col]])
  if (anyNA(groups) || any(!nzchar(groups))) stop("`", condition_col, "` has empty/NA values.")

  # ---- comparisons ---------------------------------------------------------
  comparisons <- as.data.frame(comparisons, stringsAsFactors = FALSE)
  req_comp <- c("treatment", "control", "Condition")
  miss_c <- setdiff(req_comp, colnames(comparisons))
  if (length(miss_c)) stop("`comparisons` is missing column(s): ", paste(miss_c, collapse = ", "))
  wanted <- unique(.SplitLabels(c(comparisons$treatment, comparisons$control)))
  bad_g  <- setdiff(wanted, unique(groups))
  if (length(bad_g)) {
    stop("treatment/control value(s) not found in samples$", condition_col, ": ",
         paste(bad_g, collapse = ", "), ". Present: ", paste(unique(groups), collapse = ", "))
  }
  extra_cols <- setdiff(colnames(comparisons), req_comp)

  # ---- covariates / interaction (columns of `samples`) ---------------------
  if (!is.null(interaction) &&
      (!is.character(interaction) || length(interaction) != 1 || !interaction %in% colnames(samples))) {
    stop("`interaction` must be the name of one column of `samples`.")
  }
  if (!is.null(covariate_formula)) {
    miss_v <- setdiff(all.vars(stats::as.formula(paste("~", covariate_formula))), colnames(samples))
    if (length(miss_v)) stop("Not columns of `samples`: ", paste(miss_v, collapse = ", "))
  }
  sample_metadata <- data.frame(Sample_key = samples$key,
                                samples[, setdiff(colnames(samples), c("key")), drop = FALSE],
                                stringsAsFactors = FALSE, check.names = FALSE)

  # ---- sites + polyA_DB annotation -----------------------------------------
  sites <- .read_tab(sites, "sites")
  miss_site <- setdiff(c("peak", "chrom", "site", "strand"), colnames(sites))
  if (length(miss_site)) stop("`sites` is missing column(s): ", paste(miss_site, collapse = ", "))
  sites$peak <- as.character(sites$peak)
  sites <- sites[sites$peak %in% rownames(counts) & sites$strand %in% c("+", "-"), , drop = FALSE]
  if (!nrow(sites)) stop("No site of `sites` is in `counts` (check the `peak` names).")

  if (is.null(polyAdb)) {
    if (is.null(polyAdb.file) || !file.exists(polyAdb.file)) {
      stop("Give `polyAdb.file` (path to the polyA_DB table) or `polyAdb`.")
    }
    if (verbose) message("Reading polyA_DB: ", polyAdb.file)
    polyAdb <- utils::read.table(file = polyAdb.file, header = TRUE)
  }
  polyAdb <- as.data.frame(polyAdb)

  rg <- GenomicRanges::GRanges(
    seqnames = as.character(sites$chrom),
    ranges   = IRanges::IRanges(start = as.integer(sites$site), width = 1L),
    strand   = as.character(sites$strand))
  m <- .MatchPolyADb(rg, polyAdb, max.dist = max.dist)

  ann <- data.frame(
    seqnames = as.character(sites$chrom),
    start    = as.integer(sites$site),
    end      = as.integer(sites$site),
    width    = 1L,
    strand   = as.character(sites$strand),
    peak     = sites$peak,
    match_distance = m$distance,
    stringsAsFactors = FALSE
  )
  for (cl in setdiff(colnames(polyAdb), c("hg38_Chromosome_format", "Strand"))) {
    ann[[cl]] <- polyAdb[[cl]][m$row]
  }
  for (cl in intersect(c("misprime", "internal_priming"), colnames(sites))) ann[[cl]] <- sites[[cl]]

  .is_true <- function(x) toupper(as.character(x)) %in% c("TRUE", "T", "1")
  keep <- !is.na(m$row) & ann$Intron.exon.location %in% genomic_positions &
    !is.na(ann$ensembl_gene_id) & nzchar(as.character(ann$ensembl_gene_id))
  if (verbose) {
    message("Sites: ", nrow(ann), " in counts; ", sum(!is.na(m$row)), " matched polyA_DB within ",
            max.dist, " nt; ", sum(keep), " in {", paste(genomic_positions, collapse = ", "),
            "} with a gene.")
  }
  if (isTRUE(exclude_misprime)) {
    if (!"misprime" %in% colnames(ann)) stop("`sites` has no `misprime` column.")
    drop <- keep & .is_true(ann$misprime)
    if (verbose) message("exclude_misprime: removed ", sum(drop), " site(s).")
    keep <- keep & !drop
  }
  if (isTRUE(exclude_internal_priming)) {
    if (!"internal_priming" %in% colnames(ann)) {
      stop("`sites` has no `internal_priming` column (run polyApipe_bulk.py with --genome_fasta).")
    }
    drop <- keep & .is_true(ann$internal_priming)
    if (verbose) message("exclude_internal_priming: removed ", sum(drop), " site(s).")
    keep <- keep & !drop
  }
  if (!any(keep)) stop("No site left after annotation/filters.")

  Peaks_all <- ann[keep, , drop = FALSE]
  if (!"symbol" %in% colnames(Peaks_all)) {
    sym_src <- intersect(c("Gene_Symbol", "SYMBOL", "Symbol", "gene_name"), colnames(Peaks_all))
    if (length(sym_src)) Peaks_all$symbol <- Peaks_all[[sym_src[1]]]
  }
  # Gene start/end (used by the RED core for RLDi). In DEPsMatrix they come
  # from the RNA assay's gene annotation; here from the EnsDb of genome_build.
  Peaks_all$gene_start <- NA_real_
  Peaks_all$gene_end   <- NA_real_
  if (!is.null(genome_build)) {
    ensdb <- switch(genome_build,
      hg38 = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
      hg19 = EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,
      mm10 = EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79,
      stop("`genome_build` must be \"hg38\", \"hg19\" or \"mm10\"."))
    gid <- sub("\\..*$", "", as.character(Peaks_all$ensembl_gene_id))
    gc  <- tryCatch(
      suppressMessages(AnnotationDbi::select(ensdb, keys = unique(gid), keytype = "GENEID",
                                             columns = c("GENEID", "GENESEQSTART", "GENESEQEND"))),
      error = function(e) NULL)
    if (!is.null(gc) && nrow(gc)) {
      gc <- gc[!duplicated(gc$GENEID), , drop = FALSE]
      mi <- match(gid, gc$GENEID)
      Peaks_all$gene_start <- as.numeric(gc$GENESEQSTART[mi])
      Peaks_all$gene_end   <- as.numeric(gc$GENESEQEND[mi])
    }
    if (verbose) {
      message("Gene coordinates (", genome_build, ") for ",
              length(unique(gid[!is.na(Peaks_all$gene_start)])), " of ", length(unique(gid)),
              " genes (RLDi is NA for the rest).")
    }
  } else if (verbose) {
    message("genome_build = NULL: no gene coordinates, RLDi will be NA.")
  }

  # PASTA statistic columns the RED core expects (no PASTA in bulk)
  for (cl in c("Estimate", "p.value", "p_val_adj", "percent.1", "percent.2")) Peaks_all[[cl]] <- NA_real_
  Peaks_all$CellType <- ct
  red_required <- c("symbol", "ensembl_gene_id", "seqnames", "hg38_Position", "strand",
                    "Intron.exon.location", "label", "peak", "PAS_hexamer",
                    "PolyaStrength_percentile")
  miss_red <- setdiff(red_required, colnames(Peaks_all))
  if (length(miss_red)) {
    stop("The polyA_DB table lacks column(s) the RED core needs: ", paste(miss_red, collapse = ", "))
  }

  # ---- RepeatMasker (TE) ---------------------------------------------------
  te_peaks <- NULL
  if (repeat_masker != "none") {
    if (is.null(rmsk)) stop("`rmsk` (UCSC RepeatMasker table) is required when repeat_masker = '",
                            repeat_masker, "'.")
    rmsk <- as.data.frame(rmsk, stringsAsFactors = FALSE)
    miss_rm <- setdiff(c("genoName", "genoStart", "genoEnd", "strand"), colnames(rmsk))
    if (length(miss_rm)) stop("`rmsk` is missing column(s): ", paste(miss_rm, collapse = ", "))
    qdf <- Peaks_all[Peaks_all$Intron.exon.location %in% repeat_masker_positions,
                     c("seqnames", "start", "end", "strand", "peak"), drop = FALSE]
    te_peaks <- character(0)
    if (nrow(qdf)) {
      q_gr <- GenomicRanges::makeGRangesFromDataFrame(qdf, keep.extra.columns = TRUE)
      r_gr <- GenomicRanges::makeGRangesFromDataFrame(
        rmsk, seqnames.field = "genoName", start.field = "genoStart", end.field = "genoEnd",
        strand.field = "strand", starts.in.df.are.0based = TRUE)
      hits <- suppressWarnings(GenomicRanges::findOverlaps(q_gr, r_gr, ignore.strand = FALSE))
      te_peaks <- unique(qdf$peak[S4Vectors::queryHits(hits)])
    }
    if (repeat_masker == "Remove") {
      Peaks_all <- Peaks_all[!Peaks_all$peak %in% te_peaks, , drop = FALSE]
      if (verbose) message("repeat_masker = 'Remove': dropped ", length(te_peaks), " TE site(s).")
      if (!nrow(Peaks_all)) stop("All sites were removed by the repeat_masker filter.")
    } else if (verbose) {
      message("repeat_masker = 'Keep': flagged ", length(te_peaks), " TE site(s).")
    }
  }
  .te_flag <- function(pk) {
    if (is.null(te_peaks)) return(rep(NA, length(pk)))
    ifelse(is.na(pk), NA, pk %in% te_peaks)
  }

  # ---- counts in the RED core's layout: "<group>_bulk_<key>" ---------------
  cnt <- counts[Peaks_all$peak, samples$sample, drop = FALSE]
  cnt[] <- lapply(cnt, function(v) { v <- as.numeric(v); v[is.na(v)] <- 0; v })
  colnames(cnt) <- paste0(groups, "_", ct, "_", samples$key)

  # filter_scope = "all": one expression filter on the samples of every group
  keep_all <- NULL
  if (filter_scope == "all") {
    cols <- lapply(wanted, function(l) .PseudobulkCols(colnames(cnt), l, ct))
    grp  <- rep(wanted, lengths(cols))
    k <- .FilterFeatures(as.matrix(cnt[, unlist(cols), drop = FALSE]), grp,
                         filter_method = filter_method, mincounts = mincounts,
                         filterByExpr_args = filterByExpr_args)
    keep_all <- names(k)[k]
    if (verbose) message("Common filter (all groups): ", length(keep_all), " / ", length(k), " sites.")
  }

  # ---- RED per comparison --------------------------------------------------
  if (verbose) message("Calculating RED scores for ", nrow(comparisons), " comparison(s)...")
  fails <- character(0)
  block_filter <- list()
  res <- lapply(seq_len(nrow(comparisons)), function(ci) {
    cmp <- comparisons[ci, , drop = FALSE]
    err <- NULL
    red <- tryCatch(
      .RedScoresCellType(
        Peaks             = Peaks_all,
        counts_polyA      = cnt,
        CellType          = ct,
        treatment         = cmp$treatment,
        control           = cmp$control,
        Condition         = NULL,
        filter_method     = filter_method,
        mincounts         = mincounts,
        filterByExpr_args = filterByExpr_args,
        sample_metadata   = sample_metadata,
        covariate_formula = covariate_formula,
        interaction       = interaction,
        gdpau_min_reads   = gdpau_min_reads,
        pair_rule         = pair_rule,
        min_site_reads    = min_site_reads,
        min_site_usage    = min_site_usage,
        min_rep_pairs     = min_rep_pairs,
        keep_counts       = isTRUE(return_filter),
        keep_peaks        = keep_all,
        red_tests         = red_tests
      ),
      error = function(e) { err <<- paste0("[", cmp$Condition, "] ", conditionMessage(e)); NULL }
    )
    if (!is.null(err)) fails <<- c(fails, err)
    fi <- if (!is.null(red)) attr(red, "red_filter") else NULL
    block_filter[[ci]] <<- list(
      Cells = ct, Condition = cmp$Condition, treatment = cmp$treatment, control = cmp$control,
      red_status = if (is.null(red)) "error" else if (nrow(red) == 0) "empty" else "ok",
      n_peaks = if (is.null(fi)) NA_integer_ else fi$n_peaks,
      n_kept  = if (is.null(fi)) NA_integer_ else fi$n_kept,
      kept    = if (is.null(fi)) character(0) else fi$kept,
      counts  = if (is.null(fi)) NULL else fi$counts,
      samples = if (is.null(fi)) NULL else fi$samples
    )
    if (is.null(red) || !nrow(red)) return(NULL)
    red$Cells     <- ct
    red$Condition <- cmp$Condition
    for (col in extra_cols) red[[col]] <- cmp[[col]]
    red
  })
  if (length(fails)) {
    warning("DEPsMatrixBulk: ", length(fails), " comparison(s) failed:\n  ",
            paste(fails, collapse = "\n  "), call. = FALSE)
  }

  out <- dplyr::bind_rows(res)
  if (nrow(out)) {
    out$p_Repeated_Masker <- .te_flag(out$p_peak)
    out$d_Repeated_Masker <- .te_flag(out$d_peak)
    if (isTRUE(drop_na_cols)) {
      all_na <- vapply(out, function(col) all(is.na(col)), logical(1))
      out <- out[, !all_na, drop = FALSE]
    }
    out <- dplyr::relocate(out, dplyr::any_of(c("Cells", extra_cols, "Condition")))
  }
  if (verbose) {
    message("DEPsMatrixBulk complete: ", nrow(out), " RED rows across ",
            nrow(comparisons) - length(fails), " of ", nrow(comparisons), " comparison(s).")
  }

  polyA_sites <- Peaks_all[, setdiff(colnames(Peaks_all),
                                     c("CellType", "Estimate", "p.value", "p_val_adj",
                                       "percent.1", "percent.2")), drop = FALSE]
  polyA_sites$Repeated_Masker <- .te_flag(polyA_sites$peak)
  rownames(polyA_sites) <- NULL

  result <- list(
    red_scores = out,
    polyAdb    = polyA_sites,
    sites      = ann,
    samples    = data.frame(sample = samples$sample, key = samples$key, group = groups,
                            stringsAsFactors = FALSE)
  )
  if (isTRUE(return_filter)) {
    pseudo <- cnt
    result$filter <- .BuildDEPsFilter(
      block_filter, comparisons, assay = "bulk", filter_method,
      filter_scope   = filter_scope,
      pseudo         = pseudo,
      input_features = Peaks_all$peak,
      sample_info    = data.frame(sample_id = colnames(cnt), n_cells = NA_integer_,
                                  replicate = samples$sample, stringsAsFactors = FALSE)
    )
  }
  result
}
