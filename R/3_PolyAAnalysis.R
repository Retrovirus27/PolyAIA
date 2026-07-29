#' Read raw polyA counts and build per-sample PolyA assays
#'
#' Reads each sample's raw polyA site counts (\code{PASTA::ReadPolyApipe()})
#' and immediately builds its \code{PolyAAssay} (\code{PASTA::CreatePolyAAssay()})
#' one sample at a time, discarding that sample's raw counts before moving to
#' the next. This keeps peak memory to roughly one sample's raw counts at a
#' time, instead of reading every sample's raw counts into a list first (often
#' hundreds of MB to several GB each) and only building assays afterwards in a
#' second pass over that list -- which needs all of them in memory
#' simultaneously, on top of the assays being built from them.
#'
#' @param samples Character vector of sample names. Used to build each
#'   sample's count/fragment file paths (see \code{counts_dir}/\code{fragments_dir})
#'   and as the barcode-tagging prefix.
#' @param counts_dir Directory containing each sample's raw polyA count file,
#'   named \code{"<sample>.tab.gz"}.
#' @param peaks_file Path to the shared peak annotation (\code{.gff}) file,
#'   passed to \code{ReadPolyApipe()} for every sample.
#' @param fragments_dir Directory containing each sample's fragment/barcode
#'   file, named \code{"<sample>.blocks.sort.bed.barcode.gz"}.
#' @param genome Genome build passed to \code{CreatePolyAAssay()}. Default \code{"hg38"}.
#' @param filter.chromosomes,min.features,min.cells Passed through to
#'   \code{ReadPolyApipe()}. Defaults \code{TRUE}, \code{10}, \code{25}.
#' @param gc_each_sample Logical. Whether to force garbage collection after
#'   each sample's raw counts are discarded -- a little slower, but keeps peak
#'   memory down when samples are large or numerous. Default \code{TRUE}.
#' @param verbose Logical. Whether to print progress messages. Default \code{TRUE}.
#'
#' @return A named list of \code{PolyAAssay} objects, one per sample, ready
#'   for \code{FilterPeaks()}/\code{RequantifyPolyA()}.
#'
#' @examples
#' \dontrun{
#' polyA.assays <- UploadPolyAAssays(
#'   samples       = samples,
#'   counts_dir    = "3_polyA_counts/",
#'   peaks_file    = "4_MetaData/polyA_Bam_Files_polyA_peaks.gff",
#'   fragments_dir = "1_Bed_files/",
#'   genome        = "mm10"
#' )
#' }
#'
#' @export
UploadPolyAAssays <- function(samples,
                               counts_dir,
                               peaks_file,
                               fragments_dir,
                               genome = "hg38",
                               filter.chromosomes = TRUE,
                               min.features = 10,
                               min.cells = 25,
                               gc_each_sample = TRUE,
                               verbose = TRUE) {

  polyA.assays <- stats::setNames(vector("list", length(samples)), samples)

  for (x in samples) {
    if (verbose) message("Reading sample: ", x)
    counts.file <- file.path(counts_dir, paste0(x, ".tab.gz"))
    counts <- PASTA::ReadPolyApipe(
      counts.file = counts.file,
      peaks.file = peaks_file,
      filter.chromosomes = filter.chromosomes,
      min.features = min.features,
      min.cells = min.cells
    )
    colnames(counts) <- paste(x, colnames(counts), sep = "_")
    colnames(counts) <- paste0(colnames(counts), "-1")

    fragment.file <- file.path(fragments_dir, paste0(x, ".blocks.sort.bed.barcode.gz"))

    if (verbose) message("Building PolyA assay: ", x)
    polyA.assays[[x]] <- PASTA::CreatePolyAAssay(
      counts = counts,
      genome = genome,
      fragments = fragment.file
    )

    # Discard this sample's raw counts before moving to the next -- they're
    # no longer needed once the assay is built, and are usually the largest
    # object in memory at this point.
    rm(counts)
    if (gc_each_sample) invisible(gc(verbose = FALSE))
  }

  polyA.assays
}

#' Build a peak-level annotation table from a polyA assay's own metadata
#'
#' Filters a Seurat object's polyA assay to sites already annotated as 3'-most
#' exon or intronic, drops intronic sites that overlap a transposable-element
#' (repeat-masker) region (these are frequently priming/mapping artifacts
#' rather than genuine polyA sites), and joins in gene-level metadata from the
#' object's \code{RNA} assay via \code{ensembl_gene_id}. This produces the
#' \code{Peaks} table consumed by \code{DEPsMatrix()}.
#'
#' Note this does not look anything up in an external reference polyA
#' database -- it relies on the polyA assay's own \code{meta.features}
#' already carrying \code{Intron.exon.location}/\code{ensembl_gene_id}/etc.
#' (e.g. from whatever annotated the assay upstream).
#'
#' @param seu A Seurat object containing the polyA assay (\code{"polyA"}) to
#'   filter/annotate, and an \code{RNA} assay to join gene metadata from.
#' @param rmsk A repeat-masker table (as read in via e.g. \code{read.delim()})
#'   for the SAME genome build as \code{seu}'s polyA assay coordinates --
#'   mismatched builds/species will silently produce meaningless overlaps.
#'   Expected to have \code{genoName}/\code{genoStart}/\code{genoEnd}/\code{strand}
#'   columns (UCSC RepeatMasker track format).
#' @param genomic_positions Character vector of \code{Intron.exon.location}
#'   values to keep (e.g. \code{c("3' most exon", "Intron")}).
#' @param remove_te Logical. Whether to drop intronic sites overlapping a
#'   transposable element (repeat masker). Default \code{TRUE}. Set to
#'   \code{FALSE} to keep every site and let \code{DEPsMatrix(repeat_masker =
#'   ...)} handle TE removal/flagging downstream instead (\code{rmsk} is then
#'   not needed here).
#'
#' @return A \code{data.frame} of peaks restricted to \code{genomic_positions},
#'   with transposable-element-overlapping intronic sites removed (when
#'   \code{remove_te = TRUE}), joined with \code{RNA} assay gene metadata.
#'
#' @examples
#' \dontrun{
#' rmsk <- read.delim("4_MetaData/mm10.rmsk.txt")
#' polyAdb <- Peaksdb(
#'   seu = seurat_obj_filtered,
#'   rmsk = rmsk,
#'   genomic_positions = c("3' most exon", "Intron")
#' )
#' }
#'
#' @export
Peaksdb <- function(seu, rmsk = NULL, genomic_positions, remove_te = TRUE) {

  if (!"polyA" %in% SeuratObject::Assays(seu)) {
    stop("`seu` has no 'polyA' assay.")
  }
  if (!"RNA" %in% SeuratObject::Assays(seu)) {
    stop("`seu` has no 'RNA' assay.")
  }

  # Load polyAdb (the polyA assay's own feature metadata)
  polyAdb <- data.frame(seu[["polyA"]]@meta.features)
  polyAdb <- dplyr::filter(polyAdb, !is.na(Intron.exon.location))
  polyAdb <- polyAdb[polyAdb$Intron.exon.location %in% genomic_positions, ]

  if (isTRUE(remove_te)) {
    if (is.null(rmsk)) {
      stop("`rmsk` is required when remove_te = TRUE. ",
           "Pass a RepeatMasker table, or set remove_te = FALSE to defer TE ",
           "handling to DEPsMatrix(repeat_masker = ...).")
    }

    te_query <- polyAdb %>%
      dplyr::filter(Intron.exon.location == "Intron") %>%
      GenomicRanges::makeGRangesFromDataFrame(
        seqnames.field     = "seqnames",
        start.field        = "start",
        end.field          = "end",
        strand.field       = "strand",
        keep.extra.columns = TRUE
      )

    # Repeat masker
    gr <- GenomicRanges::makeGRangesFromDataFrame(
      rmsk,
      seqnames.field = "genoName",
      start.field = "genoStart",
      end.field = "genoEnd",
      strand.field = "strand",
      starts.in.df.are.0based = TRUE,
      keep.extra.columns = TRUE
    )

    # Overlap: drop intronic sites that overlap a repeat/TE region
    hits     <- GenomicRanges::findOverlaps(te_query, gr, ignore.strand = FALSE)
    te_peaks <- unique(te_query$peak[S4Vectors::queryHits(hits)])

    polyAdb <- polyAdb %>%
      dplyr::filter(!(Intron.exon.location == "Intron" & peak %in% te_peaks))
  }

  # Merge polyAdb with RNA assay gene metadata via ensembl_gene_id
  meta_peak <- polyAdb %>%
    dplyr::rename(ID = ensembl_gene_id)

  meta_gene <- seu[["RNA"]]@meta.data[1:13] %>%
    dplyr::filter(!is.na(ID), ID != "") %>%
    dplyr::distinct(ID, .keep_all = TRUE)

  Peaks <- meta_peak %>%
    dplyr::left_join(meta_gene, by = "ID") %>%
    dplyr::distinct(peak, .keep_all = TRUE)
  Peaks$strand <- Peaks$strand.x

  Peaks
}

#' Clean PolyA Metadata by Removing Duplicate and Positional Columns
#'
#' @param obj Seurat object with polyA assay
#' @param remove_coords Remove coordinate columns? Default TRUE
#' @param coord_cols Vector of coordinate column names to remove
#' @param verbose Print messages? Default TRUE
#'
#' @return Seurat object with cleaned polyA metadata
#'
#' @export
CleanPolyAMetadata <- function(
    obj,
    remove_coords = TRUE,
    coord_cols = c("start", "end", "width", "strand", "seqnames"),
    verbose = TRUE
) {

  if (verbose) message("Cleaning polyA metadata...")

  meta <- obj[["polyA"]]@meta.features

  # Remove coordinate columns if requested
  if (remove_coords) {
    cols_to_remove <- coord_cols[coord_cols %in% colnames(meta)]
    if (length(cols_to_remove) > 0) {
      if (verbose) message("Removing ", length(cols_to_remove), " coordinate columns")
      meta <- meta[, !colnames(meta) %in% cols_to_remove, drop = FALSE]
    }
  }

  # Remove duplicate column names
  dup_cols <- duplicated(colnames(meta))
  if (any(dup_cols)) {
    if (verbose) message("Removing ", sum(dup_cols), " duplicate columns")
    meta <- meta[, !dup_cols, drop = FALSE]
  }

  obj[["polyA"]]@meta.features <- meta

  if (verbose) message("Metadata cleaned. Remaining columns: ", ncol(meta))

  return(obj)
}

#' Filter and Standardize Genomic Peaks Across Samples
#'
#' This function standardizes peaks across all samples by creating a unified peak set,
#' reducing overlaps, and counting sample representation. It also generates visualization
#' plots showing peak distribution across samples.
#'
#' @param polyA.assays A list or collection of assay objects containing genomic ranges.
#'   Each element should be compatible with \code{GenomicRanges::granges()}.
#' @param multiple_samples set TRUE if multiple studies are analyzed together.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{Ranges}{A \code{GRanges} object containing the unified peak set after
#'       reducing overlaps across all samples (strand-aware).}
#'     \item{CommonPeaks}{An integer vector indicating the number of samples that
#'       overlap with each peak in the unified set.}
#'   }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Extracts all GRanges from individual assays
#'   \item Creates a unified peak set by reducing overlapping ranges (strand-aware)
#'   \item Counts how many samples each peak appears in
#'   \item Generates two diagnostic plots:
#'     \itemize{
#'       \item Number of peaks per sample (bar plot)
#'       \item Distribution of peaks by number of samples (bar plot)
#'     }
#' }
#'
#' @note The function prints combined plots to the current graphics device.
#'   Requires the \code{patchwork} package for plot combination (using \code{|} operator).
#'   A theme object \code{th} must be defined in the global environment.
#'
#' @examples
#' \dontrun{
#' # Assuming you have a list of polyA assays
#' results <- FilterPeaks(my_polyA_assays)
#' unified_peaks <- results$Ranges
#' sample_counts <- results$CommonPeaks
#' }
#'
#' @importFrom GenomicRanges granges reduce countOverlaps GRangesList
#' @importFrom ggplot2 ggplot aes geom_bar labs theme element_text scale_y_continuous scale_x_continuous
#' @importFrom dplyr count
#' @importFrom scales comma
#'
#' @export
FilterPeaks <- function(polyA.assays, multiple_samples = NULL){
  # Before merging, standardize peaks across all samples
  # 1. Extract all GRanges from individual assays
  all.ranges <- lapply(polyA.assays, GenomicRanges::granges)
  unified.ranges <- unlist(GenomicRanges::GRangesList(all.ranges))
  # 2. Create a unified peak set by reducing overlaps
  unified.peaks <- GenomicRanges::reduce(unified.ranges, ignore.strand = FALSE)

  if (!is.null(multiple_samples)) {
    # Count how many samples have at least one range overlapping each unified peak
    peak_sample_counts <- sapply(all.ranges, function(sample_ranges) {
      countOverlaps(unified.peaks, sample_ranges, ignore.strand = FALSE) > 0
    })
    # Sum across samples (each column is a sample, each row is a unified peak)
    peak_sample_counts <- rowSums(peak_sample_counts)
  } else {
    # For single sample, count total overlapping ranges per peak
    peak_sample_counts <- GenomicRanges::countOverlaps(unified.peaks, unified.ranges, ignore.strand = FALSE)
  }

  # Number of peaks per sample
  samples <- names(polyA.assays)
  Sample_peaks <-
    data.frame(sample = samples, n_peaks = sapply(polyA.assays, function(x) length(granges(x)))) %>%
    ggplot2::ggplot(ggplot2::aes(x = sample, y = n_peaks)) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
    ggplot2::labs(
      title = "Peaks Across Samples",
      x = NULL,
      y = "Number of Peaks" ) +
    th +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45 , hjust = 1, vjust = 1)) +
    ggplot2::scale_y_continuous(
      labels = scales::comma,
      n.breaks = 6 )
  # Number of peaks shared by sample
  Sample_shared <-
    data.frame(peak_sample_counts) %>%
    dplyr::count(peak_sample_counts) %>%
    ggplot2::ggplot(ggplot2::aes(x = peak_sample_counts, y = n)) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
    ggplot2::labs(
      title = "Common peaks Across Samples",
      x = "Number of Samples",
      y = "Number of Peaks" ) +
    ggplot2::scale_x_continuous(breaks = 1:max(peak_sample_counts)) +
    th +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45 , hjust = 1, vjust = 1)) +
    ggplot2::scale_y_continuous(
      labels = scales::comma,
      n.breaks = 6  )

  # Display plots side by side - FOR R MARKDOWN
  combined_plot <- Sample_peaks | Sample_shared
  print(combined_plot)

  return(list(Ranges=unified.peaks, CommonPeaks=peak_sample_counts))
}

#' Requantify PolyA Assays Using Unified Peaks
#'
#' This function requantifies multiple PolyA assays using a unified set of peaks.
#' It maps original peak counts to unified peaks by finding overlaps and summing
#' counts from overlapping regions.
#'
#' @param polyA.assays A named list of PolyA assay objects to be requantified.
#'   Each assay should contain counts data, genomic ranges, fragments, and annotations.
#' @param unified.peaks The \code{list} returned by \code{\link{FilterPeaks}()},
#'   with elements \code{Ranges} (a \code{GRanges} of the unified peak set) and
#'   \code{CommonPeaks} (an integer vector of per-peak sample counts). Strand
#'   information is preserved during overlap detection.
#' @param samples A character vector of sample names corresponding to names in
#'   \code{polyA.assays}. Defaults to all names in \code{polyA.assays}.
#' @param threshold Threshold of the minimum peaks taken for requantification.
#' @param genome Character. Genome build to tag the requantified assays with
#'   (passed to \code{PASTA::CreatePolyAAssay(genome = ...)}). Default \code{"hg19"}.
#'
#' @return A named list of requantified PolyA assay objects with the same names as
#'   the input samples. Each assay contains:
#'   \itemize{
#'     \item counts matrix with rows corresponding to unified peaks
#'     \item unified genomic ranges
#'     \item original fragments
#'     \item original annotations
#'   }
#'
#' @examples
#' \dontrun{
#' # Requantify all samples
#' requantified <- RequantifyPolyA(
#'   polyA.assays = my_assays,
#'   unified.peaks.filtered = unified_peaks
#' )
#'
#' # Requantify specific samples
#' requantified <- RequantifyPolyA(
#'   polyA.assays = my_assays,
#'   unified.peaks.filtered = unified_peaks,
#'   samples = c("sample1", "sample2")
#' )
#' }
#'
#' @export
RequantifyPolyA <- function(polyA.assays,
                            unified.peaks,
                            threshold,
                            samples = names(polyA.assays),
                            genome = "hg19") {

  # Validate inputs
  if (!is.list(polyA.assays) || is.null(names(polyA.assays))) {
    stop("polyA.assays must be a named list")
  }

  unified.peaks.filtered <- unified.peaks$Ranges[unified.peaks$CommonPeaks >= threshold]

  message("A total of ", length(unified.peaks.filtered), " common peaks")

  # Process each sample
  polyA.assays.unified <- lapply(samples, function(x) {
    message("Requantifying sample: ", x)

    original.counts <- SeuratObject::GetAssayData(polyA.assays[[x]], slot = "counts")
    original.ranges <- GenomicRanges::granges(polyA.assays[[x]])

    # Find overlaps (strand-specific)
    ov <- GenomicRanges::findOverlaps(unified.peaks.filtered, original.ranges, ignore.strand = FALSE)

    # Build a mapping matrix: unified_peaks x original_peaks
    # This sums multiple original peaks into each unified peak
    n_unified <- length(unified.peaks.filtered)
    n_original <- length(original.ranges)

    # Create sparse mapping matrix
    mapping <- Matrix::sparseMatrix(
      i = S4Vectors::queryHits(ov),   # unified peak index
      j = S4Vectors::subjectHits(ov), # original peak index
      x = 1,
      dims = c(n_unified, n_original)
    )

    # Matrix multiplication: mapping %*% original.counts
    # This sums counts from overlapping original peaks
    new.counts <- mapping %*% original.counts

    rownames(new.counts) <- names(unified.peaks.filtered)
    colnames(new.counts) <- colnames(original.counts)

    message("  ", nrow(new.counts), " peaks x ", ncol(new.counts), " cells")

    PASTA::CreatePolyAAssay(
      counts = new.counts,
      ranges = unified.peaks.filtered,
      genome = genome,
      fragments = Signac::Fragments(polyA.assays[[x]]),
      annotation = Signac::Annotation(polyA.assays[[x]])
    )
  })

  # Preserve names
  names(polyA.assays.unified) <- samples

  return(polyA.assays.unified)
}

#' Plot PolyA Coverage for Differentially Expressed Sites
#'
#' Creates a combined track plot showing coverage, polyA sites, peaks, and gene
#' annotation for a specified gene. Allows customization of region boundaries
#' and text sizes for all plot elements.
#'
#' @param seu Seurat object containing a polyA assay. If its \code{RNA} assay's
#'   \code{meta.data} already has \code{chromosome}/\code{strand}/
#'   \code{gene_start}/\code{gene_end}/\code{transcript_info} for \code{gene}
#'   (as set by \code{UploadSce()}/\code{AnnotateGenes()} when the object was
#'   built), that's reused directly instead of re-querying \code{AnnotateGenes()}.
#' @param gene Gene symbol to plot
#' @param deg_list A single \code{data.frame}/tibble of \code{DEPsMatrix()}
#'   results across all cell types (e.g. \code{rbind}/\code{bind_rows} of
#'   per-comparison runs) -- NOT a named list. Must contain \code{Gene_symbol},
#'   \code{Cells} (cell type label; add this yourself if \code{DEPsMatrix()}
#'   was run per cell type, e.g. \code{red$Cells <- celltype}), \code{strand},
#'   \code{p_pos}/\code{d_pos}, \code{p_peak}/\code{d_peak}, and \code{Condition}
#'   (see \code{DEPsMatrix()}'s \code{Condition} argument/output column).
#' @param polyAdb PolyA database object with gene annotations
#' @param cell.type Character. Restrict BOTH the plotted region (via
#'   \code{deg_list$Cells}) and the coverage cells to this single cell
#'   type/subpopulation. Its metadata column is auto-resolved (it need not be
#'   \code{"CellType"} -- e.g. \code{"L2/3 IT"} is found in a subpopulation
#'   column). If \code{NULL} (default), all cell types are used.
#' @param group.by Metadata column name used to group cells
#' @param split.by Metadata column name used to split the plot
#' @param colors Named color palette for groups. If \code{NULL}, uses default colors.
#'   Default \code{NULL}
#' @param expand_bp Integer. Number of base pairs to expand the plotted region beyond
#'   gene boundaries. Default \code{0}
#' @param polyARegion Logical. Whether to zoom the coverage plot over changing polyA
#'   sites only. Default \code{FALSE}
#' @param highlight_region Logical. Whether to highlight the significant proximal/distal
#'   polyA site region(s) on the coverage tracks. Default \code{FALSE}
#' @param region_start Integer. Custom start position for the plotted region. If
#'   \code{NULL}, uses the gene start coordinate. Default \code{NULL}
#' @param region_end Integer. Custom end position for the plotted region. If
#'   \code{NULL}, uses the gene end coordinate. Default \code{NULL}
#' @param species Species for genome annotation.
#'   Options: \code{"human"} or \code{"mouse"}
#' @param version Genome build version.
#'   Options: \code{"hg38"}, \code{"hg19"}, \code{"mm10"}, or \code{"mm39"}
#' @param cell_filter Character vector of values used to filter which cells are
#'   shown in the coverage tracks (matched against \code{filter_col}). Accepts a
#'   broad class (e.g. \code{"Exc"}) or a list of subpopulations. Only the
#'   coverage tracks are affected -- the plotted region is defined by the gene's
#'   polyA sites regardless. \code{NULL} (default) keeps all cells.
#' @param filter_col Metadata column that \code{cell_filter} is matched against
#'   (and, when \code{bulk = FALSE}, the column tracks are split by). Default
#'   \code{"CellType"}. Set to your subpopulation column to filter/split by
#'   subpopulation.
#' @param condition Optional character vector of \code{deg_list$Condition}
#'   values (e.g. \code{"B6_Alcohol_vs_B6_Control"}) to restrict the plotted
#'   region, highlighted sites, and polyA-site condition labelling to those
#'   comparison(s). \code{NULL} (default) uses all conditions.
#' @param condition_group_cols Metadata columns whose values (joined by
#'   \code{condition_group_sep}) reproduce the group tokens on either side of a
#'   \code{condition} label (e.g. \code{c("Strain", "Treatment")} ->
#'   \code{"B6_Alcohol"}). When supplied together with \code{condition}, the
#'   COVERAGE cells are also restricted to just the compared groups. If
#'   \code{NULL} (default), only the region/labels are condition-filtered, not
#'   the coverage tracks.
#' @param condition_group_sep Separator joining \code{condition_group_cols} into
#'   a per-cell token. Default \code{"_"}.
#' @param condition_vs Delimiter separating the two compared groups inside a
#'   \code{condition} label. Default \code{"_vs_"}.
#' @param cell_order Deprecated alias of \code{cell_filter} (it always acted as
#'   a filter, never an ordering). Default \code{NULL}.
#' @param plot_layers Character. Controls which plot layers are rendered.
#'   Options: \code{"All"} or \code{"Tracks"}. Default \code{"All"}
#' @param bulk Logical. If \code{TRUE} (default), the selected cells are shown as
#'   a single combined "Bulk" coverage track (e.g. all subpopulations of
#'   \code{"Exc"} together). If \code{FALSE}, one track per unique
#'   \code{filter_col} value among the selected cells. Combine with
#'   \code{cell_filter} to bulk a chosen subset.
#' @param remove_na Logical. If \code{TRUE} (default), drop sites for the gene
#'   that have \code{NA} in \code{p_pos}/\code{d_pos}/\code{p_peak}/\code{d_peak}
#'   (e.g. single-site RED entries) before building the plotted region and
#'   highlights. This avoids the "'start' or 'end' cannot contain NAs" error
#'   and removes the need to pre-filter \code{deg_list} with
#'   \code{complete.cases()}.
#' @param verbose Logical. Whether to print progress messages. Default \code{TRUE}
#' @param text_size Numeric. General text size for plot elements. Default \code{14}
#' @param axis_text_size Numeric. Text size for axis tick labels. Default \code{12}
#' @param axis_title_size Numeric. Text size for axis titles. Default \code{14}
#' @param strip_text_size Numeric. Text size for facet strip labels. Default \code{12}
#' @param legend_text_size Numeric. Text size for legend labels. Default \code{12}
#' @param legend_title_size Numeric. Text size for legend title. Default \code{14}
#'
#' @return A combined track plot (\code{patchwork} object) with coverage, polyA
#'   sites, peaks, and gene annotation tracks
#'
#' @export
#'
PolyAPlot <- function(
    seu,
    gene,
    deg_list,
    polyAdb,
    cell.type       = NULL,
    group.by        = NULL,
    split.by        = NULL,
    colors          = NULL,
    expand_bp       = 0,
    polyARegion     = FALSE,
    highlight_region = FALSE,
    region_start    = NULL,
    region_end      = NULL,
    species         = c("human", "mouse"),
    version         = c("hg38", "hg19", "mm10", "mm39"),
    cell_filter     = NULL,
    filter_col      = "CellType",
    condition       = NULL,
    condition_group_cols = NULL,
    condition_group_sep  = "_",
    condition_vs         = "_vs_",
    cell_order      = NULL,
    plot_layers     = c("All", "Tracks"),
    bulk            = TRUE,
    remove_na       = TRUE,
    verbose         = TRUE,
    text_size       = 12,
    axis_text_size  = 10,
    axis_title_size = 12,
    strip_text_size = 10,
    legend_text_size  = 8,
    legend_title_size = 10
) {

  # Resolve multiple-choice arguments up front so mistakes are caught
  # immediately instead of erroring later with "condition has length > 1".
  species     <- match.arg(species)
  version     <- match.arg(version)
  plot_layers <- match.arg(plot_layers)

  # `cell_order` was renamed to `cell_filter` (it always acted as a filter, not
  # an ordering). Keep the old name working as a deprecated alias.
  if (!is.null(cell_order)) {
    if (is.null(cell_filter)) cell_filter <- cell_order
    warning("`cell_order` is deprecated; use `cell_filter` instead.", call. = FALSE)
  }

  # Resolve which metadata column actually holds a set of filter values. The
  # same cell label can live in different columns depending on annotation --
  # e.g. a broad class "Exc" in "CellType" but a subpopulation "L2/3 IT" in a
  # separate column -- so `filter_col`/`cell.type` shouldn't have to assume
  # "CellType". Prefer `preferred` when it contains all values, otherwise search
  # the object's other character/factor columns for one that does. Returns
  # NA_character_ if no column contains all the values.
  resolve_meta_col <- function(values, preferred) {
    md <- seu@meta.data
    values <- as.character(values)
    if (!is.null(preferred) && preferred %in% colnames(md) &&
        all(values %in% as.character(md[[preferred]]))) {
      return(preferred)
    }
    for (cn in colnames(md)) {
      col <- md[[cn]]
      if ((is.character(col) || is.factor(col)) &&
          all(values %in% as.character(col))) {
        return(cn)
      }
    }
    NA_character_
  }

  if (is.null(cell.type)) {
    # Search across all cell types
    if (verbose) message("Searching for gene '", gene, "' across all cell types...")

    # Filter the dataframe for the gene
    results <- deg_list %>%
      dplyr::filter(Gene_symbol == gene)

    if (nrow(results) == 0) {
      stop("Gene '", gene, "' not found in any cell type")
    }

    if (verbose) message("Found ", nrow(results), " significant sites in ",
                         length(unique(results$Cells)), " cell type(s)")
  } else {
    # Single cell type
    if (verbose) message("Subsetting to cell type: ", cell.type)

    # Subset Seurat object. Resolve which column `cell.type` lives in instead of
    # assuming "CellType" -- e.g. "L2/3 IT" is typically a Subpopulation value,
    # so `seu$CellType == "L2/3 IT"` used to match nothing and error with
    # "No cells found".
    ct_col <- resolve_meta_col(cell.type, filter_col)
    if (is.na(ct_col)) {
      stop("`cell.type` value '", cell.type, "' was not found in any metadata ",
           "column of `seu` (checked '", filter_col, "' and the rest). Check the ",
           "spelling/casing -- names come from the Seurat object's metadata.")
    }
    if (verbose && ct_col != filter_col) {
      message("  (matched cell.type '", cell.type, "' in column '", ct_col, "')")
    }
    seu <- seu[, as.character(seu@meta.data[[ct_col]]) == cell.type]

    # Filter deg_list for specific cell type and gene
    results <- deg_list %>%
      dplyr::filter(Cells == cell.type, Gene_symbol == gene)
    if (nrow(results) == 0) {
      stop("Gene '", gene, "' not found in cell type '", cell.type, "'")
    }

    if (verbose) message("Found ", nrow(results), " significant sites")
  }

  # Optional Condition filter: restrict the plotted region / highlighted sites
  # (and, below, the polyA-site condition labelling) to one or more comparisons
  # from deg_list$Condition (e.g. "B6_Alcohol_vs_B6_Control").
  if (!is.null(condition)) {
    if (!"Condition" %in% colnames(deg_list)) {
      stop("`condition` was supplied but `deg_list` has no `Condition` column.")
    }
    bad_cond <- setdiff(condition, unique(as.character(deg_list$Condition)))
    if (length(bad_cond) > 0) {
      stop("`condition` value(s) not found in deg_list$Condition: ",
           paste(bad_cond, collapse = ", "), ". Available: ",
           paste(unique(as.character(deg_list$Condition)), collapse = ", "))
    }
    results  <- results  %>% dplyr::filter(Condition %in% condition)
    deg_list <- deg_list %>% dplyr::filter(Condition %in% condition)
    if (nrow(results) == 0) {
      stop("Gene '", gene, "' has no sites for condition(s): ",
           paste(condition, collapse = ", "))
    }
    if (verbose) message("Filtered to condition(s) {", paste(condition, collapse = ", "),
                         "}: ", nrow(results), " site(s) remaining.")
  }

  # Drop rows lacking the coordinate columns used to build the plotted region
  # and highlights. Single-site RED entries carry NA proximal fields
  # (p_pos/p_peak) and NA RED, which otherwise propagate into
  # IRanges()/min()/max() below and error with
  # "'start' or 'end' cannot contain NAs". With remove_na = TRUE (default)
  # these are filtered out here so callers no longer need to pre-filter
  # deg_list with complete.cases().
  if (isTRUE(remove_na)) {
    n_before_na <- nrow(results)
    results <- results %>%
      dplyr::filter(
        !is.na(p_pos), !is.na(d_pos),
        !is.na(p_peak), !is.na(d_peak)
      )
    if (verbose && nrow(results) < n_before_na) {
      message("Removed ", n_before_na - nrow(results),
              " site(s) with NA coordinates; ", nrow(results), " remaining.")
    }
    if (nrow(results) == 0) {
      stop("All sites for gene '", gene, "' had NA coordinates after removing NAs. ",
           "Nothing left to plot (set remove_na = FALSE to inspect, or check that ",
           "this gene has a valid proximal/distal polyA pair).")
    }
  }

  # Extract coordinates and create region.highlight
  region.highlight <- NULL
  if (isTRUE(highlight_region)) {
    coordinates <- c(results$p_peak, results$d_peak)
    coordinates <- do.call(rbind, strsplit(coordinates, "-"))

    # Bug fix: start/end used to be de-duplicated INDEPENDENTLY
    # (`unique(coordinates[, 2])`, `unique(coordinates[, 3])`), which breaks
    # the row-wise pairing between a peak's start and its own end whenever
    # two different peaks happen to share a start (or an end) but not both
    # -- the two unique() calls can come back different lengths, silently
    # mismatching/recycling start against the wrong end. Deduplicate whole
    # (chr, start, end) rows instead, so each kept row is still a real peak.
    coord_df <- unique(as.data.frame(coordinates, stringsAsFactors = FALSE))

    region.highlight <- GenomicRanges::GRanges(
      seqnames = sub(".*\\.", "", coord_df[, 1]),
      ranges = IRanges::IRanges(
        start = as.numeric(coord_df[, 2]),
        end = as.numeric(coord_df[, 3])
      )
    )

    if (verbose) message("Highlighting ", length(region.highlight), " region(s)")
  }

  # Create ROI from gene coordinates.
  # Note: `species`/`version` default to "human"/"hg38" (the first option in
  # each match.arg()), NOT the mouse data most of this package's other
  # scripts (G1-G3) work with -- if you're plotting a mouse gene and forget
  # to pass species = "mouse", version = "mm10", AnnotationDbi will fail
  # with an unhelpful "None of the keys entered are valid keys for 'SYMBOL'"
  # (it queried the human OrgDb for a mouse-cased symbol). Wrapping the call
  # to surface a clearer, actionable message for exactly that case.
  # Reuse cached gene-level annotation from the RNA assay's own meta.data if
  # it's already there, instead of always hitting AnnotateGenes() (which
  # re-queries AnnotationDbi/an EnsDb from scratch every call -- the slow,
  # repeated "downloading coordinates every time" step). UploadSce()/
  # AnnotateGenes() already stash exactly these columns (chromosome, strand,
  # gene_start, gene_end, transcript_info) onto seu[["RNA"]]@meta.data,
  # keyed by gene symbol, the first time the Seurat object is built --
  # there's no need to requery for a gene that's already annotated there.
  required_meta_cols <- c("chromosome", "strand", "gene_start", "gene_end", "transcript_info")
  rna_meta <- tryCatch(seu[["RNA"]]@meta.data, error = function(e) NULL)

  gene_row <- NULL
  if (!is.null(rna_meta) && all(required_meta_cols %in% colnames(rna_meta))) {
    if (gene %in% rownames(rna_meta)) {
      # rownames(seu[["RNA"]]) already are gene symbols
      gene_row <- rna_meta[gene, , drop = FALSE]
    } else {
      # Bug fix: rownames(seu[["RNA"]]) are frequently NOT gene symbols (e.g.
      # Ensembl IDs) -- same underlying issue DotPlot() had to handle. The
      # first version of this cache check only tried `gene %in%
      # rownames(rna_meta)`, which silently never matched here, so it fell
      # through to AnnotateGenes() (and its slow Signac::GetGRangesFromEnsDb()
      # call) on every single run. Fall back to matching the "SYMBOL"
      # feature-metadata column (case-insensitive column name, exact-then-
      # case-insensitive value match), same convention as DotPlot().
      symbol_col <- grep("^symbol$", colnames(rna_meta), ignore.case = TRUE, value = TRUE)
      if (length(symbol_col) > 0) {
        symbol_col <- symbol_col[1]
        match_idx <- which(rna_meta[[symbol_col]] == gene)
        if (length(match_idx) == 0) {
          match_idx <- which(toupper(rna_meta[[symbol_col]]) == toupper(gene))
        }
        if (length(match_idx) >= 1) {
          gene_row <- rna_meta[match_idx[1], , drop = FALSE]
        }
      }
    }
  }

  if (!is.null(gene_row) &&
      !anyNA(gene_row[1, c("chromosome", "strand", "gene_start", "gene_end")])) {

    if (verbose) message("Using cached annotation for '", gene, "' from seu[[\"RNA\"]]@meta.data (skipping AnnotateGenes())")
    gene_meta <- gene_row

  } else {
    gene_meta <- tryCatch(
      suppressWarnings(suppressMessages(
        AnnotateGenes(data = gene, IDtype = "SYMBOL", species = species, genome_build = version)
      )),
      error = function(e) {
        stop(
          "AnnotateGenes() failed to find '", gene, "' as a SYMBOL for species = '",
          species, "' (genome_build = '", version, "'). If this is a mouse gene ",
          "(e.g. 'Camk2a'), pass species = \"mouse\", version = \"mm10\" (or \"mm39\") ",
          "explicitly -- PolyAPlot()'s species/version default to human/hg38. ",
          "Original error: ", conditionMessage(e),
          call. = FALSE
        )
      }
    )
  }

  # Use custom region coordinates if provided, otherwise use gene coordinates
  if (!is.null(region_start) && !is.null(region_end) && isTRUE(polyARegion)) {
    stop("region_start/region_end and polyARegion cannot be used together")

  } else if (!is.null(region_start) && !is.null(region_end)) {
    if (verbose) message("Using custom region coordinates: ", region_start, "-", region_end)
    roi <- GenomicRanges::GRanges(
      seqnames = gene_meta$chromosome,
      ranges = IRanges::IRanges(start = region_start, end = region_end),
      strand = gene_meta$strand
    )

  } else if (!is.null(region_start) || !is.null(region_end)) {
    stop("Both region_start and region_end must be provided together")

  } else if (isTRUE(polyARegion)) {
    # `unique(results$strand == "+")` used to be handed straight to `if()` --
    # fine as long as every row for this gene agrees on strand, but if
    # `results` ever contained a mixed/inconsistent strand annotation for the
    # same gene (a data problem, but one that shouldn't be swallowed
    # silently), that returns a length > 1 logical and `if()` errors with an
    # opaque "the condition has length > 1". Check explicitly instead.
    gene_strands <- unique(results$strand)
    if (length(gene_strands) > 1) {
      stop("Gene '", gene, "' has inconsistent strand values in `deg_list` (",
           paste(gene_strands, collapse = ", "), ") -- expected a single strand.")
    }
    is_plus_strand <- gene_strands == "+"
    roi <- GenomicRanges::GRanges(
      seqnames = gene_meta$chromosome,
      ranges = IRanges::IRanges(
        start = if (is_plus_strand) min(results$p_pos) else max(results$d_pos),
        end   = if (is_plus_strand) max(results$d_pos) else min(results$p_pos)
      ),
      strand = gene_meta$strand
    )

  } else {
    roi <- GenomicRanges::GRanges(
      seqnames = gene_meta$chromosome,
      ranges = IRanges::IRanges(start = gene_meta$gene_start, end = gene_meta$gene_end),
      strand = gene_meta$strand
    )
  }

  # Expand roi by specified bp
  roi_expanded <- roi + expand_bp
  roi_start <- IRanges::start(roi_expanded)
  roi_end <- IRanges::end(roi_expanded)

  if (verbose) message("Final region: ", as.character(GenomicRanges::seqnames(roi))[1], ":", roi_start, "-", roi_end)


  gene_annot_all<- GenomicRanges::GRanges(
    tx_id = as.data.frame(gene_meta$transcript_info)$tx_id,
    seqnames = as.data.frame(gene_meta$transcript_info)$tx_chromosome,
    ranges = IRanges::IRanges(start = as.data.frame(gene_meta$transcript_info)$tx_start, end = as.data.frame(gene_meta$transcript_info)$tx_end),
    strand = as.data.frame(gene_meta$transcript_info)$tx_strand,
    type = as.data.frame(gene_meta$transcript_info)$tx_type
  )

  if (verbose) {
    message("Found ", length(gene_annot_all), " annotation entries for gene ", gene)
    if (length(gene_annot_all) > 0) {
      message("Types present: ", paste(unique(gene_annot_all$type), collapse = ", "))
      # Check for transcript IDs
      if ("tx_id" %in% names(S4Vectors::mcols(gene_annot_all))) {
        message("Unique transcripts: ", length(unique(gene_annot_all$tx_id)))
      }
    }
  }

  # Extract polyA sites for this gene
  polyA_gene <- polyAdb[polyAdb$Gene_Symbol %in% gene, ]

  # Bug fix: "Shared" used to be hardcoded to `length(conds) == 3` with a
  # comment claiming that means "all 4 conditions" (a leftover/wrong comment
  # -- 3 happens to be right for this specific 3-comparison dataset, but it
  # silently stops meaning "all conditions" the moment deg_list is built from
  # a different number of comparisons). Derive it from how many distinct
  # Condition values actually exist in deg_list instead of hardcoding it.
  n_conditions <- dplyr::n_distinct(deg_list$Condition, na.rm = TRUE)

  peak_conditions <- deg_list %>%
    dplyr::select(p_peak, d_peak, Condition) %>%
    tidyr::pivot_longer(
      cols = c(p_peak, d_peak),
      names_to = "site",
      values_to = "peak"
    ) %>%
    dplyr::filter(!is.na(peak)) %>%
    dplyr::group_by(peak) %>%
    dplyr::summarize(
      conditions = {
        conds <- unique(Condition)
        if (length(conds) == n_conditions) "Shared"  # found in every comparison
        else if (length(conds) == 1) conds            # unique to one condition
        else paste(sort(conds), collapse = " & ")     # shared between some but not all
      },
      .groups = "drop"
    )

  polyA_gene <- polyA_gene %>%
    dplyr::left_join(peak_conditions, by = "peak")

  if (nrow(polyA_gene) == 0) {
    warning("No polyA sites found in polyAdb for gene '", gene, "'")
    polyA_gr <- GenomicRanges::GRanges()
    df <- data.frame()
  } else {
    polyA_gr <- GenomicRanges::GRanges(
      seqnames = polyA_gene$seqnames,
      ranges = IRanges::IRanges(
        start = polyA_gene$hg38_Position,
        end = polyA_gene$hg38_Position
      ),
      strand = polyA_gene$strand,
      Condition = polyA_gene$conditions
    )

    #peak.intersect <- subsetByOverlaps(x = polyA_gr, ranges = roi_expanded)
    df <- as.data.frame(x = polyA_gr)

    if (verbose && nrow(df) > 0) {
      message("Found ", nrow(df), " polyA sites in the selected region")
    }
  }

  # Create plots
  gene_plot <- plot_annotation_from_granges(
    gene = gene,
    gene_annot_all  = gene_annot_all,
    roi_expanded    = roi_expanded
  )

  # Check version
  if (version %in% c("hg38","mm10")) {
    roi_expanded <- GenomeInfoDb::renameSeqlevels(roi_expanded, stats::setNames(paste0("chr", GenomeInfoDb::seqlevels(roi_expanded)), GenomeInfoDb::seqlevels(roi_expanded)))
  }

  # Coverage plot per cell
  seu <- CleanPolyAMetadata(seu)

  # --- cell_filter: keep only the requested cells for the coverage tracks -----
  # This does NOT change the plotted region (that's defined by the gene's polyA
  # sites in `results`) -- it only controls which cells' coverage is shown.
  # cell_filter values are matched against `filter_col` (default the broad
  # "CellType", e.g. "Exc"); to select subpopulations instead, pass
  # filter_col = "<your subpopulation column>" and a vector of subpopulations.
  if (!is.null(cell_filter)) {
    # Resolve the column cell_filter's values live in (same logic as cell.type),
    # so callers don't have to hand-match filter_col to broad-class vs
    # subpopulation names. Falls back to filter_col if nothing contains all
    # values (then the "no cells matched" error below fires with context).
    resolved_col <- resolve_meta_col(cell_filter, filter_col)
    if (is.na(resolved_col)) resolved_col <- filter_col
    if (!resolved_col %in% colnames(seu@meta.data)) {
      stop("`filter_col` '", resolved_col, "' not found in seu metadata.")
    }
    if (verbose && !identical(resolved_col, filter_col)) {
      message("cell_filter: matched values in column '", resolved_col,
              "' (not '", filter_col, "').")
    }
    filter_col <- resolved_col  # use downstream (also for bulk = FALSE track split)

    keep_cells <- as.character(seu@meta.data[[filter_col]]) %in% cell_filter
    if (!any(keep_cells)) {
      stop("No cells matched cell_filter = {", paste(cell_filter, collapse = ", "),
           "} in column '", filter_col, "'. Present values include: ",
           paste(utils::head(unique(as.character(seu@meta.data[[filter_col]])), 10),
                 collapse = ", "), " ...")
    }
    if (verbose) message("cell_filter kept ", sum(keep_cells), " / ", ncol(seu),
                         " cells (", filter_col, " in {",
                         paste(cell_filter, collapse = ", "), "}).")
    seu <- seu[, keep_cells]
  }

  # --- Restrict coverage cells to the group(s) being compared in `condition` --
  # The Condition label (e.g. "B6_Alcohol_vs_B6_Control") is split on
  # `condition_vs` into its two group tokens ("B6_Alcohol", "B6_Control"); a
  # per-cell token built from `condition_group_cols` (joined by
  # `condition_group_sep`) is matched against them, so the coverage tracks show
  # only the cells actually being compared. Requires `condition_group_cols` so
  # PolyAPlot knows how the label maps to metadata (e.g. c("Strain","Treatment")
  # -> "B6_Alcohol"); without it the region/labels are still condition-filtered
  # but the coverage cells are not.
  if (!is.null(condition)) {
    if (is.null(condition_group_cols)) {
      if (verbose) {
        message("Note: coverage cells NOT filtered by condition -- pass ",
                "`condition_group_cols` (e.g. c(\"Strain\",\"Treatment\")) to ",
                "also restrict the coverage tracks to the compared groups.")
      }
    } else {
      if (!all(condition_group_cols %in% colnames(seu@meta.data))) {
        stop("Not all `condition_group_cols` found in seu metadata: ",
             paste(setdiff(condition_group_cols, colnames(seu@meta.data)), collapse = ", "))
      }
      group_tokens <- unique(unlist(
        strsplit(as.character(condition), condition_vs, fixed = TRUE)
      ))
      cell_tokens <- do.call(paste, c(seu@meta.data[condition_group_cols],
                                      sep = condition_group_sep))
      keep_cond <- cell_tokens %in% group_tokens
      if (!any(keep_cond)) {
        stop("No cells matched the condition group tokens {",
             paste(group_tokens, collapse = ", "), "} built from columns {",
             paste(condition_group_cols, collapse = ", "), "}. Example cell token: '",
             cell_tokens[1], "'. Check `condition_group_cols`/`condition_group_sep`/",
             "`condition_vs` -- the token order/separator must reproduce your ",
             "Condition label halves.")
      }
      if (verbose) {
        message("Condition coverage filter kept ", sum(keep_cond), " / ", ncol(seu),
                " cells (", paste(condition_group_cols, collapse = condition_group_sep),
                " in {", paste(group_tokens, collapse = ", "), "}).")
      }
      seu <- seu[, keep_cond]
    }
  }

  if (isTRUE(bulk)) {
    # One combined "Bulk" coverage track over all (filtered) cells. When the
    # caller doesn't supply group.by we force a single constant grouping column
    # -- otherwise Signac::CoveragePlot() groups by Idents(seu) (here the
    # fine-grained celltype/condition identity) and emits one track per level,
    # i.e. the "one track per subpopulation of Exc" behavior. group.by/split.by
    # are still honored when explicitly provided.
    if (is.null(group.by)) {
      seu@meta.data[[".polyaplot_bulk"]] <- "Bulk"
      bulk_group <- ".polyaplot_bulk"
    } else {
      bulk_group <- group.by
    }
    cov_plot <- list(
      Bulk = suppressWarnings(
        Signac::CoveragePlot(
          object           = seu,
          region           = roi_expanded,
          region.highlight = region.highlight,
          annotation       = FALSE,
          peaks            = FALSE,
          group.by         = bulk_group,
          split.by         = split.by,
          expression.assay = "polyA",
          expression.slot  = "scale.data"
        )
      )
    )
    cell_types <- "Bulk"

  } else {
    # One coverage track per unique `filter_col` value among the (filtered)
    # cells -- per broad CellType by default, or per subpopulation when
    # filter_col is a subpopulation column.
    if (!filter_col %in% colnames(seu@meta.data)) {
      stop("`filter_col` '", filter_col, "' not found in seu metadata.")
    }
    track_vals <- if (!is.null(cell_filter)) {
      intersect(cell_filter, unique(as.character(seu@meta.data[[filter_col]])))
    } else {
      unique(as.character(seu@meta.data[[filter_col]]))
    }
    seu@meta.data[[filter_col]] <- factor(as.character(seu@meta.data[[filter_col]]),
                                          levels = track_vals)
    cell_types <- track_vals

    cov_plot <- sapply(cell_types, function(x) {
      if (verbose) message("Generating coverage plot for: ", x)
      keep_x <- !is.na(seu@meta.data[[filter_col]]) & seu@meta.data[[filter_col]] == x
      suppressWarnings(
        Signac::CoveragePlot(
          object           = seu[, keep_x],
          region           = roi_expanded,
          region.highlight = region.highlight,
          annotation       = FALSE,
          peaks            = FALSE,
          group.by         = group.by,
          split.by         = split.by,
          expression.assay = "polyA",
          expression.slot  = "scale.data"
        )
      )
    }, simplify = FALSE)
  }

  # Bug fix: this used to run unconditionally with `values = colors`, so
  # when `colors` is left at its default NULL (as in most calls), it applied
  # a manual scale with ZERO color values on top of Signac::CoveragePlot()'s
  # own default fill (which is typically Idents(seu) -- e.g. 4 treatment
  # groups here), immediately erroring with "Insufficient values in manual
  # scale. N needed but only 0 provided." Only override the fill scale when
  # the caller actually supplied colors.
  if (!is.null(colors)) {
    for (i in seq_along(cov_plot)) {
      cov_plot[[i]] <- cov_plot[[i]] & ggplot2::scale_fill_manual(values = colors)
    }
  }

  # Normalize y-axis across all plots using the reference from the max plot
  y_maxes <- sapply(cov_plot, function(p) {
    built <- ggplot2::ggplot_build(p[[1]])
    max(built$layout$panel_params[[1]]$y.range, na.rm = TRUE)
  })

  max_idx <- which.max(y_maxes)
  ref_scale <- ggplot2::ggplot_build(cov_plot[[max_idx]][[1]])$layout$panel_params[[1]]
  ref_breaks <- ref_scale$y.breaks
  ref_limits <- ref_scale$y.range

  # Extract and normalize the inner ggplot from each patchwork
  cov_plot_norm <- lapply(seq_along(cov_plot), function(i) {
    x     <- cell_types[i]
    inner <- cov_plot[[i]][[1]]

    suppressMessages(
      inner & ggplot2::scale_y_continuous(
        limits = ref_limits,
        breaks = ref_breaks,
        name   = paste0("Normalized signal\n(range 0 - ", round(ref_limits[2]), ")")
      )
    ) &
      ggplot2::labs(tag = x) &
      ggplot2::theme(
        plot.tag          = ggplot2::element_text(
          angle = 0,
          size  = text_size,
          hjust = 0.5,
          vjust = 0.5,
          face  = "bold"
        ),
        plot.tag.position = "top"
      )
  })
  names(cov_plot_norm) <- cell_types

  peak_plot <- suppressWarnings(
    Signac::PeakPlot(
      object = seu,
      region = roi_expanded
    ) +
      ggplot2::labs(y = "\nPeaks") +
      ggplot2::theme(
        text = ggplot2::element_text(size = text_size),
        axis.text = ggplot2::element_text(size = axis_text_size),
        axis.title = ggplot2::element_text(size = axis_title_size)
      )
  )

  df<- df[stats::complete.cases(df),]

  polyA_plot <- ggplot2::ggplot(df) +
    ggplot2::geom_point(
      ggplot2::aes(x = start, y = 0.5, fill = Condition, color = Condition),
      size = 2.5,
      shape = 25,
    ) +
    ggplot2::coord_cartesian(xlim = c(roi_start, roi_end), ylim = c(0.25, 0.75)) +
    ggplot2::labs(y = "\nPAS") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "none",
      axis.ticks.y  = ggplot2::element_blank(),
      axis.text.y   = ggplot2::element_blank(),
      text          = ggplot2::element_text(size = text_size),
      axis.text.x   = ggplot2::element_text(size = axis_text_size),
      axis.title    = ggplot2::element_text(size = axis_title_size),
      legend.text   = ggplot2::element_text(size = legend_text_size),
      legend.title  = ggplot2::element_text(size = legend_title_size)
    )

  if (plot_layers == "All") {
    Signac::CombineTracks(
      plotlist = c(
        cov_plot_norm,
        list(
          polyA_plot,
          peak_plot,
          gene_plot
        )
      ),
      heights = c(rep(1.6, length(cov_plot_norm)), 0.5, 0.5, 1.5)
    )
  } else {
    Signac::CombineTracks(
      plotlist = c(
        cov_plot_norm,
        list(
          gene_plot
        )
      ),
      heights = c(rep(1.6, length(cov_plot_norm)),1)
    )
  }

}


#' Plot Gene Annotation from GRanges Object
#'
#' Generates a gene annotation track plot from a \code{GRanges} object for a
#' specified region of interest. Intended as an internal track component within
#' \code{PolyAPlot()}.
#'
#' @param gene Character. Gene symbol, used only for axis/plot labeling.
#' @param gene_annot_all \code{GRanges} object containing gene annotation features
#'   (e.g., exons, transcripts) for the region of interest
#' @param roi_expanded \code{GRanges} object defining the expanded region of
#'   interest to plot
#' @param text_size Numeric. General text size for plot elements. Default \code{14}
#' @param axis_text_size Numeric. Text size for axis tick labels. Default \code{12}
#' @param axis_name_size Numeric. Text size for gene/feature name labels. Default \code{6}
#' @param axis_title_size Numeric. Text size for axis titles. Default \code{14}
#'
#' @return A \code{ggplot2} object with the gene annotation track
#'
#' @keywords internal
#'
plot_annotation_from_granges <- function(
    gene,
    gene_annot_all,
    roi_expanded,
    text_size       = 12,
    axis_text_size  = 10,
    axis_name_size  = 6,
    axis_title_size = 12
) {

  # Convert to df and filter to roi
  annot_df <- as.data.frame(gene_annot_all, row.names = NULL) %>%
    dplyr::rename(chromosome = seqnames) %>%
    dplyr::filter(
      start <= IRanges::end(IRanges::ranges(roi_expanded)),
      end   >= IRanges::start(IRanges::ranges(roi_expanded))
    )

  if (nrow(annot_df) == 0) {
    message("No annotation entries overlap the roi")
    return(ggplot2::ggplot() + ggplot2::theme_void())
  }

  roi_start <- IRanges::start(IRanges::ranges(roi_expanded))
  roi_end   <- IRanges::end(IRanges::ranges(roi_expanded))

  # Clip coordinates to roi
  annot_df <- annot_df %>%
    dplyr::mutate(
      start = pmax(start, roi_start),
      end   = pmin(end,   roi_end)
    )

  # Assign y position per transcript (stagger them)
  tx_levels <- unique(annot_df$tx_id)
  tx_y      <- stats::setNames(seq_along(tx_levels), tx_levels)
  annot_df  <- annot_df %>%
    dplyr::mutate(y = tx_y[tx_id])

  n_tx <- length(tx_levels)

  # Backbone per transcript (thin line spanning min-max)
  backbone_df <- annot_df %>%
    dplyr::group_by(tx_id) %>%
    dplyr::summarise(
      x    = min(start),
      xend = max(end),
      y    = unique(y),
      .groups = "drop"
    )

  # Exons / UTRs / CDS heights
  height_map <- c(exon = 0.35, utr = 0.15, cds = 0.35, gap = 0.05)

  annot_df <- annot_df %>%
    dplyr::mutate(
      block_h = dplyr::recode(as.character(type), !!!height_map, .default = 0.2) * 8  # scale factor, tune to taste
    )

  # Gene name label position
  label_df <- backbone_df %>%
    dplyr::slice(1) %>%
    dplyr::mutate(
      x     = (x + xend) / 2,
      label = unique(annot_df$gene_name)
    )

  plot_df <- annot_df %>%
    dplyr::filter(type != "gap") %>%
    dplyr::group_by(tx_id) %>%
    dplyr::mutate(has_cds = any(type == "cds")) %>%
    dplyr::ungroup() %>%
    dplyr::filter(
      (has_cds & type %in% c("utr", "cds")) |
        (!has_cds & type == "exon")
    )


  gene_plot <- ggplot2::ggplot() +

    ggplot2::geom_segment(
      data = backbone_df,
      ggplot2::aes(x = x, xend = xend, y = y, yend = y),
      linewidth   = 0.1,
      show.legend = FALSE
    ) +
    ggplot2::geom_segment(
      data = plot_df %>% dplyr::filter(type != "gap"),
      ggplot2::aes(x = start, xend = end, y = y, yend = y, linewidth = block_h),
      show.legend = FALSE
    ) +
    ggplot2::scale_linewidth_identity() +
    ggplot2::coord_cartesian(xlim = c(roi_start, roi_end), ylim = c(0, max(backbone_df$y) + 1)) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.05)) +
    #ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = backbone_df$y, labels = backbone_df$tx_id,) +
    #ggplot2::ylim(c(0.9, max(tx_y) + 0.4)) +
    ggplot2::labs(x = paste0(unique(annot_df$chromosome), " position (bp)"), y = paste0(gene,"\nTranscripts")) +
    ggplot2::theme_classic(base_size = text_size) +
    ggplot2::theme(
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.y  = ggplot2::element_text(size = axis_name_size, hjust = 1, vjust = 0.25, margin = ggplot2::margin(r = -80)),
      axis.text.x  = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.title   = ggplot2::element_text(size = axis_title_size),
      axis.line.x.bottom = ggplot2::element_line(size = 0.1),
      axis.line.y.left   = ggplot2::element_line(size = 0.5)
    )

  return(gene_plot)
}



#' Get Cell Type-Specific PolyA Differential Expression
#'
#' Performs differential expression analysis on polyA peaks for a specified
#' cell type, comparing two identity groups using a linear model. Optionally
#' restricts the analysis to a subset of features and supports covariate
#' correction.
#'
#' @param seu_obj A Seurat object containing a polyA assay.
#' @param celltype Character string specifying the cell type to subset and test.
#' @param ident1 Character string specifying the first identity group
#'   (e.g., \code{"Treatment"}). Default is \code{"Treatment"}.
#' @param ident2 Character string specifying the second identity group
#'   (e.g., \code{"Control"}). Default is \code{"Control"}.
#' @param features Character vector of polyA peak features to test. If
#'   \code{NULL}, all features in the polyA assay are used. Default is
#'   \code{NULL}.
#' @param covariates Character vector of metadata column names to include as
#'   covariates in the linear model. Default is \code{NULL}.
#' @param verbose Logical indicating whether to print progress messages.
#'   Default is \code{TRUE}.
#'
#' @return A data frame with differential expression results for polyA peaks,
#'   including statistical estimates comparing \code{ident1} vs \code{ident2}
#'   within the specified cell type.
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' DEPolyAPeaks(seu_obj, celltype = "Astrocytes")
#'
#' # With custom identities and covariates
#' DEPolyAPeaks(
#'   seu_obj,
#'   celltype   = "Microglia",
#'   ident1     = "Alcohol",
#'   ident2     = "Control",
#'   covariates = c("sex", "age"),
#'   verbose    = FALSE
#' )
#' }
#'
#' @export
DEPolyAPeaks <- function(
    seu_obj,
    celltype,
    ident1 = "Treatment",
    ident2 = "Control",
    features = NULL,
    covariates = NULL,
    verbose = TRUE
) {

  # Build identity names
  ident_1 <- paste0(celltype, "_", ident1)
  ident_2 <- paste0(celltype, "_", ident2)

  if (verbose) message("Running DEGs for: ", ident_1, " vs ", ident_2)

  # Run differential expression test
  if (!requireNamespace("PASTA", quietly = TRUE)) {
    stop("Package 'PASTA' is required for FindDifferentialPolyA")
  }
  res <- PASTA::FindDifferentialPolyA(
    seu_obj,
    ident.1 = ident_1,
    ident.2 = ident_2,
    features = features,
    covariates = covariates
  )

  # Filter complete cases
  res <- res[stats::complete.cases(res), ]

  if (verbose) message("Found ", nrow(res), " complete results")

  # Add metadata
  res$CellType <- celltype
  res$peak <- rownames(res)

  return(res)
}

#' Merge differential-polyA results with peak- and gene-level metadata
#'
#' Joins a \code{DEPolyAPeaks()} result table (per cell type) with the polyA
#' assay's peak metadata (\code{meta.features}) and the RNA assay's gene
#' metadata, then stamps a \code{Condition} label and de-duplicates by peak.
#' Internal helper for \code{DEPsMatrix()}; mirrors the "Merge DEGs with
#' metadata" step previously written inline in the analysis scripts.
#'
#' @param deg_df A per-cell-type \code{DEPolyAPeaks()} result (must contain a
#'   \code{peak} column).
#' @param meta_peak Peak-level metadata (\code{seu[[assay]]@meta.features});
#'   must contain \code{peak} and \code{peak_gene_key}.
#' @param meta_gene Gene-level metadata slice from the RNA assay (columns 1
#'   through the last one needed, e.g. \code{transcript_info}); must contain
#'   \code{rna_gene_key}.
#' @param condition Character label written to a \code{Condition} column
#'   (\code{NULL} leaves it unset).
#' @param peak_gene_key Column in \code{meta_peak}/\code{deg_df} holding the
#'   Ensembl gene id. Default \code{"ensembl_gene_id"}.
#' @param rna_gene_key Column in \code{meta_gene} holding the Ensembl gene id.
#'   Default \code{"ID"}.
#'
#' @return A merged, peak-deduplicated \code{data.frame}.
#'
#' @noRd
.MergePeakGeneMetadata <- function(deg_df, meta_peak, meta_gene, condition = NULL,
                                   peak_gene_key = "ensembl_gene_id",
                                   rna_gene_key = "ID") {
  if (!"peak" %in% colnames(meta_peak)) meta_peak$peak <- rownames(meta_peak)

  df <- dplyr::left_join(deg_df, meta_peak, by = "peak")
  join_spec <- stats::setNames(rna_gene_key, peak_gene_key)
  df <- dplyr::left_join(df, meta_gene, by = join_spec, relationship = "many-to-many")

  if (!is.null(condition)) df$Condition <- condition

  df <- dplyr::distinct(df, peak, .keep_all = TRUE)
  rownames(df) <- df$peak
  df
}

#' Build a peak table without PASTA differential testing
#'
#' For \code{Pasta = FALSE}: constructs the peak table \code{DEPsMatrix()}'s
#' RED core consumes directly from peak- and gene-level metadata, replicated
#' across cell types, with the PASTA-derived statistic columns filled with
#' \code{NA} (no residuals / no \code{FindDifferentialPolyA()} are run). RED
#' scores themselves are still computed downstream from pseudobulk counts.
#'
#' @param meta_peak Peak-level metadata (\code{seu[[assay]]@meta.features}).
#' @param meta_gene Gene-level metadata slice from the RNA assay.
#' @param celltypes Character vector of cell types to replicate rows across.
#' @param features Optional character vector restricting which peaks to keep.
#' @param peak_gene_key,rna_gene_key Join keys, see \code{.MergePeakGeneMetadata}.
#'
#' @return A \code{data.frame} of peaks x cell types with \code{NA} PASTA stats.
#'
#' @noRd
.BuildPeaksNoPASTA <- function(meta_peak, meta_gene, celltypes, features = NULL,
                               peak_gene_key = "ensembl_gene_id",
                               rna_gene_key = "ID") {
  if (!"peak" %in% colnames(meta_peak)) meta_peak$peak <- rownames(meta_peak)

  base <- meta_peak
  if (!is.null(features)) base <- base[base$peak %in% features, , drop = FALSE]

  # PASTA / FindDifferentialPolyA output columns the RED core expects -- filled
  # with NA since PASTA isn't being run in this branch.
  base$Estimate  <- NA_real_
  base$p.value   <- NA_real_
  base$p_val_adj <- NA_real_
  base$percent.1 <- NA_real_
  base$percent.2 <- NA_real_

  # The RED core needs a lowercase `symbol` column. In the PASTA path this comes
  # from the differential-test result; here it isn't produced, so map it from
  # the polyA assay's gene-symbol column (Gene_Symbol / SYMBOL / Symbol).
  if (!"symbol" %in% colnames(base)) {
    sym_src <- intersect(c("Gene_Symbol", "SYMBOL", "Symbol", "gene_name"), colnames(base))
    if (length(sym_src) > 0) {
      base$symbol <- base[[sym_src[1]]]
    }
  }

  join_spec <- stats::setNames(rna_gene_key, peak_gene_key)
  base <- dplyr::left_join(base, meta_gene, by = join_spec, relationship = "many-to-many")
  base <- dplyr::distinct(base, peak, .keep_all = TRUE)

  do.call(rbind, lapply(celltypes, function(ct) {
    d <- base
    d$CellType <- ct
    d
  }))
}

#' Differential polyA usage + RED scores, end to end
#'
#' Integrated pipeline that takes an annotated Seurat object (with a polyA
#' assay and an RNA assay) and returns a combined RED-score table across all
#' requested cell types and comparisons. It folds together what were
#' previously separate manual steps:
#' \enumerate{
#'   \item (\code{Pasta = TRUE}) polyA residual calculation via
#'     \code{PASTA::CalcPolyAResiduals()} and per-comparison differential
#'     testing via \code{DEPolyAPeaks()} (\code{PASTA::FindDifferentialPolyA});
#'   \item merging the differential results with peak- and gene-level metadata
#'     (\code{.MergePeakGeneMetadata});
#'   \item pseudobulk aggregation via \code{Seurat::AggregateExpression()};
#'   \item per-(cell type x comparison) RED-score computation via the internal
#'     RED core (formerly the exported \code{DEPsMatrix()}).
#' }
#' When \code{Pasta = FALSE}, residual calculation and differential testing are
#' skipped: the peak table is still built (from metadata, with \code{NA} PASTA
#' statistics) and RED scores are computed from the pseudobulk counts alone.
#'
#' @param seu A Seurat object containing a polyA assay (\code{assay}) and an
#'   \code{RNA} assay carrying gene-level metadata.
#' @param comparisons A \code{data.frame}/list of comparisons, one row each,
#'   with required columns \code{ident1}, \code{ident2} (identity tokens for
#'   \code{DEPolyAPeaks()}, matched against \code{condition_cols}),
#'   \code{treatment}, \code{control} (pseudobulk count-column prefixes for the
#'   RED core), and \code{Condition} (a label). Any additional columns (e.g.
#'   \code{Strain}, \code{Treatment}) are stamped onto every output row for
#'   that comparison.
#' @param celltype_col Metadata column defining cell types (e.g.
#'   \code{"Subpopulation"}).
#' @param condition_cols Metadata column(s) whose values, pasted with
#'   \code{condition_sep}, reproduce the \code{ident1}/\code{ident2} tokens
#'   (e.g. \code{c("Treatment", "Strain")} -> \code{"Control_B6"}). Used to
#'   build the differential-testing identity as
#'   \code{paste0(<celltype>, "_", <condition token>)}.
#' @param group.by Grouping passed to \code{Seurat::AggregateExpression()} for
#'   pseudobulk counts (e.g.
#'   \code{c("Strain", "Treatment", "Subpopulation", "Number")}).
#' @param features Optional polyA features (peaks) to analyze. If \code{NULL}
#'   (default), they are derived from the polyA assay's own \code{meta.features}
#'   as the sites whose \code{Intron.exon.location} is in
#'   \code{genomic_positions} -- i.e. the peak-list step formerly done with a
#'   separate \code{Peaksdb()} call is folded in here. Pass a vector to override.
#' @param genomic_positions \code{Intron.exon.location} values kept when
#'   deriving \code{features} automatically (ignored if \code{features} is
#'   supplied). Default \code{c("3' most exon", "Intron")}.
#' @param assay Name of the polyA assay. Default \code{"polyA"}.
#' @param gene.names Gene-annotation column in the polyA assay's meta.features,
#'   passed to \code{CalcPolyAResiduals()}. Default \code{"Gene_Symbol"}.
#' @param background Identity value used as the background distribution for
#'   \code{CalcPolyAResiduals()} (only when \code{Pasta = TRUE}).
#' @param background_cols Optional metadata column(s) concatenated (with
#'   \code{background_sep}) to set \code{Idents(seu)} for the residual/background
#'   step. If \code{NULL}, the object's current \code{Idents} are used as-is.
#' @param background_sep Separator for \code{background_cols}. Default \code{""}
#'   (e.g. \code{Treatment}+\code{Strain} -> \code{"ControlB6"}).
#' @param min.counts.background Passed to \code{CalcPolyAResiduals()}. Default 5.
#' @param Pasta Logical. If \code{TRUE} (default), run PASTA residuals +
#'   differential testing; if \code{FALSE}, skip both and compute RED from
#'   counts only (PASTA stat columns are \code{NA}).
#' @param de_features Optional features passed to \code{DEPolyAPeaks()} for the
#'   differential test. Default \code{NULL} = test every peak that has a
#'   residual. Do NOT set this to \code{features}: residual calculation drops
#'   peaks lacking a gene annotation or below \code{min.counts.background}, so
#'   \code{features} is a superset of the residual-bearing peaks and passing it
#'   here yields all-NA rows that get filtered out (collapsing results to 0).
#' @param min_cells Integer. Native minimum-cells threshold: cell types (in
#'   \code{celltype_col}) with fewer than this many cells are dropped up front
#'   (too few for stable pseudobulk / RED estimates), replacing the previously
#'   manual "drop tiny subpopulations" step. Those cells still contribute to
#'   the residual background; they are only excluded from differential testing
#'   and RED. Additionally, a cell type that clears this threshold overall but
#'   lacks the required pseudobulk columns for a given comparison (e.g. no cells
#'   in one condition) is skipped for that comparison with a warning rather than
#'   erroring. Default 10.
#' @param filter_method,mincounts Passed to the RED core. Default
#'   \code{"edgeR"}, 10.
#' @param gdpau_min_reads Minimum total supporting reads a group must have (for
#'   the gene, for gDPAU; for the site pair, for DPAU) to compute the index
#'   there; below this that group's value is \code{NA}, making the
#'   \code{DPAU}/\code{gDPAU} difference \code{NA}. Default 5 (Wang et al. 2022).
#'   See \strong{Value} for the \code{DPAU}/\code{gDPAU} columns.
#' @param covariates Optional covariate columns passed to \code{DEPolyAPeaks()}.
#' @param sample_metadata,covariate_formula Optional covariate handling passed
#'   to the RED core.
#' @param gene_meta_last_col Name of the last RNA-metadata column to include in
#'   the gene-metadata join; columns 1 through this one are kept. Default
#'   \code{"transcript_info"} (the last column carrying transcript ranges).
#' @param peak_gene_key,rna_gene_key Join keys between peak metadata and RNA
#'   gene metadata. Defaults \code{"ensembl_gene_id"} / \code{"ID"}.
#' @param condition_sep Separator for \code{condition_cols}. Default \code{"_"}.
#' @param aggregate_slot Slot/layer passed to \code{AggregateExpression()}.
#'   Default \code{"counts"}.
#' @param repeat_masker RepeatMasker (transposable-element) handling for the
#'   peak table, applied after residuals/differential testing but before RED
#'   scoring (so the \code{return_seurat} object keeps every site):
#'   \code{"none"} (default, no TE handling; \code{rmsk} not needed),
#'   \code{"Remove"} (drop peaks overlapping a TE), or \code{"Keep"} (retain all
#'   but flag them). The result gains \code{p_Repeated_Masker} /
#'   \code{d_Repeated_Masker} columns (\code{TRUE} = the site overlaps a TE;
#'   \code{NA} when \code{"none"}). Requires \code{rmsk} for \code{"Remove"}/
#'   \code{"Keep"}; the peak coordinates come from the object's polyA
#'   \code{meta.features}, so no external polyA database is needed.
#' @param rmsk A UCSC RepeatMasker table (data.frame with
#'   \code{genoName}/\code{genoStart}/\code{genoEnd}/\code{strand} columns) for
#'   the SAME genome build as the polyA coordinates. Required when
#'   \code{repeat_masker} is \code{"Remove"} or \code{"Keep"}.
#' @param repeat_masker_positions Character vector of \code{Intron.exon.location}
#'   values checked for TE overlap. Default \code{"Intron"} (matches
#'   \code{Peaksdb()}); include \code{"3' most exon"} to also flag/remove exonic
#'   TE sites. Sites outside this set are never flagged (\code{FALSE}).
#' @param drop_na_cols Logical. Drop columns of \code{red_scores} that are
#'   entirely \code{NA} -- e.g. with \code{Pasta = FALSE} all PASTA statistic
#'   columns are \code{NA} by construction, and \code{RED_LRT_cov_*} is
#'   \code{NA} without a \code{covariate_formula}. Default \code{TRUE}; set
#'   \code{FALSE} to keep a fixed column layout.
#' @param return_seurat Logical. The result is always a list with
#'   \code{$red_scores} and \code{$polyAdb} (see \strong{Value}). If
#'   \code{TRUE}, the list additionally includes \code{$seu} -- the
#'   (residual-bearing, when \code{Pasta = TRUE}) Seurat object for
#'   \code{PolyAPlot()} coverage tracks. Default \code{FALSE}.
#' @param verbose Logical. Print progress. Default \code{TRUE}.
#'
#' @return A combined \code{data.frame} of RED scores across all cell types and
#'   comparisons, with \code{Cells} (cell type), \code{Condition}, and any
#'   extra \code{comparisons} label columns stamped on. Two distal-usage
#'   difference columns (treatment - control; Wang et al. 2022,
#'   doi:10.1073/pnas.2113504118) are included: \code{DPAU} -- the classic
#'   2-site index for the row's proximal/distal PAIR
#'   (\eqn{distal / (distal + proximal)} from each group's stacked reads),
#'   defined per pair; and \code{gDPAU} -- the gene-level location-index-weighted
#'   generalization (\eqn{\frac{1}{n-1}\sum_{i=1}^{n}(i-1)p_i} over all of the
#'   gene's sites; equals DPAU when \eqn{n=2}). A group with fewer than
#'   \code{gdpau_min_reads} supporting reads gives \code{NA}.
#'
#'   Always returns a \strong{list}: \code{$red_scores} (the RED table above,
#'   suitable as \code{PolyAPlot()}'s \code{deg_list}) and \code{$polyAdb} (the
#'   analyzed peak-level table -- one row per peak, TE handling applied, with a
#'   \code{Repeated_Masker} flag -- suitable as \code{PolyAPlot()}'s
#'   \code{polyAdb}). With \code{return_seurat = TRUE} the list also has
#'   \code{$seu}, the (residual-bearing, when \code{Pasta = TRUE}) object.
#'
#' @examples
#' \dontrun{
#' comparisons <- data.frame(
#'   ident1    = c("Alcohol_B6", "Control_3xTg", "Alcohol_3xTg"),
#'   ident2    = c("Control_B6", "Control_B6",   "Control_B6"),
#'   treatment = c("B6_Alcohol", "3xTg_Control", "3xTg_Alcohol"),
#'   control   = c("B6_Control", "B6_Control",   "B6_Control"),
#'   Condition = c("B6_Alcohol_vs_B6_Control",
#'                 "3xTg_Control_vs_B6_Control",
#'                 "3xTg_Alcohol_vs_B6_Control"),
#'   Strain    = c("B6", "3xTg", "3xTg"),
#'   Treatment = c("Alcohol", "Control", "Alcohol"),
#'   stringsAsFactors = FALSE
#' )
#'
#' red_scores <- DEPsMatrix(
#'   seu             = seurat_polyA,
#'   comparisons     = comparisons,
#'   celltype_col    = "Subpopulation",
#'   condition_cols  = c("Treatment", "Strain"),
#'   group.by        = c("Strain", "Treatment", "Subpopulation", "Number"),
#'   features        = features.last.exon,
#'   background      = "ControlB6",
#'   background_cols = c("Treatment", "Strain"),
#'   Pasta           = TRUE
#' )
#' }
#'
#' @export
DEPsMatrix <- function(
    seu,
    comparisons,
    celltype_col,
    condition_cols,
    group.by,
    features            = NULL,
    genomic_positions   = c("3' most exon", "Intron"),
    assay               = "polyA",
    gene.names          = "Gene_Symbol",
    background          = NULL,
    background_cols     = NULL,
    background_sep      = "",
    min.counts.background = 5,
    Pasta               = TRUE,
    de_features         = NULL,
    min_cells           = 10,
    filter_method       = "edgeR",
    mincounts           = 10,
    gdpau_min_reads     = 5,
    covariates          = NULL,
    sample_metadata     = NULL,
    covariate_formula   = NULL,
    gene_meta_last_col  = "transcript_info",
    peak_gene_key       = "ensembl_gene_id",
    rna_gene_key        = "ID",
    condition_sep       = "_",
    aggregate_slot      = "counts",
    repeat_masker       = c("none", "Remove", "Keep"),
    rmsk                = NULL,
    repeat_masker_positions = "Intron",
    drop_na_cols        = TRUE,
    return_seurat       = FALSE,
    verbose             = TRUE
) {

  # ---- Validation -----------------------------------------------------------
  repeat_masker <- match.arg(repeat_masker)
  if (!methods::is(seu, "Seurat")) stop("`seu` must be a Seurat object.")
  if (repeat_masker != "none" && is.null(rmsk)) {
    stop("`rmsk` (a RepeatMasker table with genoName/genoStart/genoEnd/strand ",
         "columns) is required when repeat_masker = '", repeat_masker, "'.")
  }
  if (!assay %in% SeuratObject::Assays(seu)) {
    stop("Assay '", assay, "' not found in `seu`.")
  }
  if (!"RNA" %in% SeuratObject::Assays(seu)) {
    stop("`seu` needs an 'RNA' assay carrying gene-level metadata.")
  }

  comparisons <- as.data.frame(comparisons, stringsAsFactors = FALSE)
  req_comp <- c("ident1", "ident2", "treatment", "control", "Condition")
  miss_comp <- setdiff(req_comp, colnames(comparisons))
  if (length(miss_comp) > 0) {
    stop("`comparisons` is missing required columns: ", paste(miss_comp, collapse = ", "))
  }

  if (!celltype_col %in% colnames(seu@meta.data)) {
    stop("`celltype_col` '", celltype_col, "' not found in seu metadata.")
  }
  if (!all(condition_cols %in% colnames(seu@meta.data))) {
    stop("Not all `condition_cols` found in seu metadata: ",
         paste(setdiff(condition_cols, colnames(seu@meta.data)), collapse = ", "))
  }

  celltypes_raw <- unique(as.character(seu@meta.data[[celltype_col]]))

  # ---- Native minimum-cells threshold --------------------------------------
  # Replaces the previously-manual "drop tiny subpopulations" step: cell types
  # with fewer than `min_cells` cells are dropped up front (too few cells for
  # stable pseudobulk / RED estimates). Nothing is subset out of the object
  # itself -- those cells still contribute to the residual background -- they
  # are just not iterated over for differential testing / RED.
  ct_counts <- table(as.character(seu@meta.data[[celltype_col]]))
  keep_ct   <- names(ct_counts)[ct_counts >= min_cells]
  dropped_ct <- setdiff(celltypes_raw, keep_ct)
  if (length(dropped_ct) > 0 && verbose) {
    message("Dropping ", length(dropped_ct), " cell type(s) with < ", min_cells,
            " cells: ",
            paste0(dropped_ct, " (", ct_counts[dropped_ct], ")", collapse = ", "))
  }
  celltypes_raw <- intersect(celltypes_raw, keep_ct)
  if (length(celltypes_raw) == 0) {
    stop("No cell type in '", celltype_col, "' has >= min_cells (", min_cells,
         ") cells. Lower `min_cells` or check `celltype_col`.")
  }

  # Canonical (Seurat-object) cell-type names are the source of truth. Cell-type
  # strings picked up elsewhere (pseudobulk column names, differential results)
  # get sanitized to dotted form (spaces/slashes/hyphens -> "."), so we keep a
  # sanitized <-> canonical lookup and always report the canonical name in the
  # output `Cells` column -- callers no longer need to remap dotted names back.
  ct_sanitize   <- function(x) make.names(x)
  canonical_by_san <- stats::setNames(celltypes_raw, ct_sanitize(celltypes_raw))

  # ---- Peak- and gene-level metadata ---------------------------------------
  meta_peak <- seu[[assay]]@meta.features
  if (!"peak" %in% colnames(meta_peak)) meta_peak$peak <- rownames(meta_peak)

  rna_meta <- seu[["RNA"]]@meta.data
  # Keep columns 1 through `gene_meta_last_col` (e.g. transcript_info, which is
  # the last column carrying the per-gene transcript ranges needed downstream).
  up_to <- match(gene_meta_last_col, colnames(rna_meta))
  if (is.na(up_to)) {
    up_to <- min(12L, ncol(rna_meta))
    warning("`gene_meta_last_col` = '", gene_meta_last_col,
            "' not found in RNA meta.data; falling back to the first ", up_to,
            " columns.")
  }
  meta_gene <- rna_meta[, seq_len(up_to), drop = FALSE]
  if (!rna_gene_key %in% colnames(meta_gene)) {
    stop("`rna_gene_key` '", rna_gene_key, "' not among the first ", up_to,
         " RNA meta.data columns (", paste(colnames(meta_gene), collapse = ", "), ").")
  }
  if (!peak_gene_key %in% colnames(meta_peak)) {
    stop("`peak_gene_key` '", peak_gene_key, "' not found in ", assay,
         " meta.features.")
  }

  # ---- Derive features if not supplied -------------------------------------
  # Fold the peak-list building (formerly a separate Peaksdb() call) into
  # DEPsMatrix: when `features` is NULL, take the sites whose
  # `Intron.exon.location` is in `genomic_positions` (3'-most exon / intron by
  # default -- the same position filter Peaksdb() uses) straight from the polyA
  # assay's own meta.features. TE removal still happens later via
  # `repeat_masker`, so this returns every position-matched site.
  if (is.null(features)) {
    if (!"Intron.exon.location" %in% colnames(meta_peak)) {
      stop("Cannot derive `features` automatically: the ", assay,
           " meta.features has no 'Intron.exon.location' column. Pass `features` ",
           "explicitly, or annotate the assay first.")
    }
    keep_pos <- !is.na(meta_peak$Intron.exon.location) &
      meta_peak$Intron.exon.location %in% genomic_positions
    features <- unique(meta_peak$peak[keep_pos])
    if (length(features) == 0) {
      stop("No peaks matched genomic_positions {",
           paste(genomic_positions, collapse = ", "), "} in the ", assay,
           " meta.features.")
    }
    if (verbose) {
      message("Derived ", length(features), " features from ", assay,
              " meta.features (Intron.exon.location in {",
              paste(genomic_positions, collapse = ", "), "}).")
    }
  }

  # ---- Build the peak table -------------------------------------------------
  if (isTRUE(Pasta)) {

    # 1. Residuals -- run under the background grouping identity.
    if (!is.null(background_cols)) {
      if (!all(background_cols %in% colnames(seu@meta.data))) {
        stop("Not all `background_cols` found in seu metadata: ",
             paste(setdiff(background_cols, colnames(seu@meta.data)), collapse = ", "))
      }
      seu@meta.data[[".deps_bg_ident"]] <-
        do.call(paste, c(seu@meta.data[background_cols], sep = background_sep))
      Seurat::Idents(seu) <- seu@meta.data[[".deps_bg_ident"]]
    }

    if (!requireNamespace("PASTA", quietly = TRUE)) {
      stop("Package 'PASTA' is required when Pasta = TRUE.")
    }
    if (verbose) message("Calculating polyA residuals (PASTA::CalcPolyAResiduals)...")
    seu <- PASTA::CalcPolyAResiduals(
      seu,
      assay                 = assay,
      features              = features,
      gene.names            = gene.names,
      background            = background,
      min.counts.background = min.counts.background,
      verbose               = verbose
    )

    # 2. Differential-testing identity: paste0(<celltype>, "_", <condition token>)
    cond_token <- do.call(paste, c(seu@meta.data[condition_cols], sep = condition_sep))
    seu@meta.data[[".deps_de_ident"]] <-
      paste0(seu@meta.data[[celltype_col]], "_", cond_token)
    Seurat::Idents(seu) <- seu@meta.data[[".deps_de_ident"]]

    # 3. Per comparison x cell type: DEPolyAPeaks + metadata merge.
    # Track outcomes so an all-empty result yields an actionable error rather
    # than a silent 0x0 data.frame.
    de_errors    <- character(0)
    n_de_ok      <- 0L
    n_de_empty   <- 0L

    # Sanity check: the idents DEPolyAPeaks builds (paste0(celltype, "_",
    # ident1/2)) must actually exist in Idents(seu). Warn early with a concrete
    # example if none of the requested comparison idents are present.
    present_idents <- unique(as.character(Seurat::Idents(seu)))
    wanted_idents  <- unique(unlist(lapply(seq_len(nrow(comparisons)), function(ci) {
      c(paste0(celltypes_raw, "_", comparisons$ident1[ci]),
        paste0(celltypes_raw, "_", comparisons$ident2[ci]))
    })))
    if (length(intersect(wanted_idents, present_idents)) == 0) {
      stop("None of the comparison identities exist in Idents(seu). ",
           "DEPolyAPeaks builds idents as paste0(<celltype>, \"_\", <ident1/ident2>). ",
           "Example wanted: '", wanted_idents[1], "'. ",
           "Example present: '", paste(utils::head(present_idents, 3), collapse = "', '"), "'. ",
           "Check `condition_cols` (currently c(", paste(condition_cols, collapse = ", "),
           ")) and `condition_sep` ('", condition_sep, "') reproduce your ident1/ident2 tokens.")
    }

    peaks_list <- lapply(seq_len(nrow(comparisons)), function(ci) {
      cmp <- comparisons[ci, ]
      per_ct <- lapply(celltypes_raw, function(ct) {
        deg <- tryCatch(
          # NOTE: DE features default to NULL (test every peak that has a
          # residual), NOT `features`. CalcPolyAResiduals() above already
          # restricted residuals to `features` minus those dropped for lacking
          # a gene annotation or failing min.counts.background -- so
          # `features` is a SUPERSET of the residual-bearing peaks. Passing it
          # here made FindDifferentialPolyA emit all-NA rows for the missing
          # peaks, which complete.cases() then removed, collapsing results to
          # 0. Mirror the standalone workflow (which passed no features here).
          DEPolyAPeaks(seu, ct, ident1 = cmp$ident1, ident2 = cmp$ident2,
                       features = de_features, covariates = covariates, verbose = verbose),
          error = function(e) {
            de_errors <<- c(de_errors, paste0("[", ct, " / ", cmp$Condition, "] ",
                                              conditionMessage(e)))
            NULL
          }
        )
        if (is.null(deg)) return(NULL)
        if (nrow(deg) == 0) { n_de_empty <<- n_de_empty + 1L; return(NULL) }
        n_de_ok <<- n_de_ok + 1L
        .MergePeakGeneMetadata(deg, meta_peak, meta_gene, condition = cmp$Condition,
                               peak_gene_key = peak_gene_key, rna_gene_key = rna_gene_key)
      })
      dplyr::bind_rows(per_ct[!vapply(per_ct, is.null, logical(1))])
    })
    Peaks_all <- dplyr::bind_rows(peaks_list)

    if (verbose) {
      message("Differential testing: ", n_de_ok, " (celltype x comparison) blocks returned peaks, ",
              n_de_empty, " returned 0 rows, ", length(de_errors), " errored.")
    }

    # If nothing came back, surface the actual reason instead of an empty frame.
    if (nrow(Peaks_all) == 0) {
      if (length(de_errors) > 0) {
        stop("DEPsMatrix: every DEPolyAPeaks() call failed. First error:\n  ",
             de_errors[1],
             if (length(de_errors) > 1) paste0("\n  (", length(de_errors) - 1,
                                               " more similar errors)") else "",
             "\nCommon causes: residuals not present in the polyA assay's scale.data, ",
             "or `features` not matching the residual matrix rownames.")
      }
      stop("DEPsMatrix: differential testing ran but every block returned 0 rows ",
           "after complete.cases filtering (n empty = ", n_de_empty, "). This usually ",
           "means FindDifferentialPolyA produced all-NA results -- check that polyA ",
           "residuals were computed (scale.data) and that `gene.names`/`background` are correct.")
    }

  } else {
    if (verbose) {
      message("Pasta = FALSE: skipping residual calculation and differential ",
              "testing; building peak table from metadata (PASTA stats = NA).")
    }
    Peaks_all <- .BuildPeaksNoPASTA(meta_peak, meta_gene, celltypes_raw,
                                    features = features,
                                    peak_gene_key = peak_gene_key,
                                    rna_gene_key = rna_gene_key)
  }

  if (nrow(Peaks_all) == 0) {
    stop("No peaks produced. See messages above for where the pipeline collapsed.")
  }

  # Peak-level strand: the peak/gene joins produce strand.x (peak) / strand.y
  # (gene); the RED core needs a single `strand` column set to the peak strand.
  if (!"strand" %in% colnames(Peaks_all) && "strand.x" %in% colnames(Peaks_all)) {
    Peaks_all$strand <- Peaks_all$strand.x
  }

  # Up-front check that the peak table has every column the RED core needs.
  # Without this the RED core stops per (cell type x comparison), those stops
  # get swallowed by the per-block tryCatch into warnings, and the run ends with
  # a silent 0-row result. Fail loudly here naming the missing columns instead.
  red_required <- c(
    "symbol", "ensembl_gene_id", "seqnames", "hg38_Position", "strand",
    "Intron.exon.location", "label", "peak", "PAS_hexamer",
    "PolyaStrength_percentile", "Estimate", "p.value", "p_val_adj",
    "percent.1", "percent.2", "CellType"
  )
  miss_red <- setdiff(red_required, colnames(Peaks_all))
  if (length(miss_red) > 0) {
    stop("The peak table is missing column(s) the RED core needs: ",
         paste(miss_red, collapse = ", "),
         ". These come from the polyA assay's meta.features (peak/gene ",
         "annotation) plus the differential-test result; in the Pasta = FALSE ",
         "path the PASTA stat columns are set to NA and `symbol` is mapped from ",
         "a gene-symbol column. If one is still absent, annotate the polyA ",
         "assay (e.g. via GetPolyADbAnnotation/Peaksdb) so meta.features carries it.")
  }

  # ---- RepeatMasker (TE) overlap on the peak table -------------------------
  # Folds Peaksdb()'s intronic-TE handling into DEPsMatrix, applied AFTER
  # residuals/differential testing (so the returned residual object keeps every
  # site) but BEFORE the RED pairing/scoring. Sites whose `Intron.exon.location`
  # is in `repeat_masker_positions` (default "Intron") are overlapped against
  # `rmsk`; a per-peak flag is built and, in "Remove" mode, TE peaks are dropped
  # here so they never form a proximal/distal pair.
  #   - "Remove": drop TE peaks; retained peaks flagged FALSE.
  #   - "Keep":   keep all peaks; TE peaks flagged TRUE, the rest FALSE.
  #   - "none":   flag is NA (no rmsk needed).
  # The flag surfaces on the output as p_Repeated_Masker / d_Repeated_Masker
  # (the result rows are proximal/distal pairs). With the default intron-only
  # positions, distal sites (never intronic) are always FALSE.
  te_peaks_all <- NULL
  if (repeat_masker != "none") {
    need_cols <- c("seqnames", "start", "end", "strand", "Intron.exon.location", "peak")
    miss_pk <- setdiff(need_cols, colnames(Peaks_all))
    if (length(miss_pk) > 0) {
      stop("repeat_masker overlap needs column(s) not present in the peak table: ",
           paste(miss_pk, collapse = ", "),
           ". These come from the polyA assay's meta.features.")
    }
    rmsk <- as.data.frame(rmsk, stringsAsFactors = FALSE)
    miss_rm <- setdiff(c("genoName", "genoStart", "genoEnd", "strand"), colnames(rmsk))
    if (length(miss_rm) > 0) {
      stop("`rmsk` is missing UCSC RepeatMasker column(s): ",
           paste(miss_rm, collapse = ", "), ".")
    }

    in_pos <- Peaks_all$Intron.exon.location %in% repeat_masker_positions
    qdf    <- Peaks_all[in_pos, c("seqnames", "start", "end", "strand", "peak"), drop = FALSE]
    qdf    <- qdf[!duplicated(qdf$peak), , drop = FALSE]

    te_peaks <- character(0)
    if (nrow(qdf) > 0) {
      q_gr <- GenomicRanges::makeGRangesFromDataFrame(
        qdf, seqnames.field = "seqnames", start.field = "start", end.field = "end",
        strand.field = "strand", keep.extra.columns = TRUE
      )
      r_gr <- GenomicRanges::makeGRangesFromDataFrame(
        rmsk, seqnames.field = "genoName", start.field = "genoStart",
        end.field = "genoEnd", strand.field = "strand",
        starts.in.df.are.0based = TRUE, keep.extra.columns = TRUE
      )
      hits <- GenomicRanges::findOverlaps(q_gr, r_gr, ignore.strand = FALSE)
      te_peaks <- unique(qdf$peak[S4Vectors::queryHits(hits)])
    }

    if (repeat_masker == "Remove") {
      n_pk_before <- length(unique(Peaks_all$peak))
      Peaks_all <- Peaks_all[!(Peaks_all$peak %in% te_peaks), , drop = FALSE]
      if (verbose) {
        message("repeat_masker = 'Remove': dropped ", length(te_peaks),
                " TE-overlapping peak(s) (of ", n_pk_before, " unique).")
      }
      if (nrow(Peaks_all) == 0) {
        stop("All peaks were removed by the repeat_masker filter. Check `rmsk`/",
             "`repeat_masker_positions` (and that rmsk's genome build matches).")
      }
    } else if (verbose) {
      message("repeat_masker = 'Keep': flagged ", length(te_peaks),
              " TE-overlapping peak(s) in Repeated_Masker.")
    }

    # Keep the TE peak IDs themselves and test membership with %in% downstream.
    # (A named-lookup vector keyed on the analyzed peaks returns NA for any ID
    # not found, which silently turned most rows' flags into NA; %in% correctly
    # yields FALSE for "peak exists and does not overlap a TE".)
    te_peaks_all <- te_peaks
  }

  # Sanitize cell-type labels with the SAME transform R applies to the
  # pseudobulk column names (as.data.frame() -> make.names()), so treatment/
  # control prefixes match the pseudobulk columns regardless of whether the
  # object's cell-type names contain dots, spaces, slashes or hyphens.
  Peaks_all$CellType <- ct_sanitize(Peaks_all$CellType)
  celltypes_red <- unique(Peaks_all$CellType)

  # ---- Pseudobulk aggregation ----------------------------------------------
  # suppressMessages hides Seurat's benign "group.by variable starts with a
  # number, appending 'g'" notice (from sanitizing e.g. Strain "3xTg").
  if (verbose) message("Aggregating pseudobulk polyA counts...")
  pseudo <- suppressMessages(Seurat::AggregateExpression(
    seu,
    assays               = assay,
    group.by             = group.by,
    slot                 = aggregate_slot,
    normalization.method = "none",
    scale.factor         = 0
  ))
  pseudo <- as.data.frame(pseudo)
  pseudo$peak <- rownames(pseudo)

  # ---- RED scores per (cell type x comparison) -----------------------------
  extra_cols <- setdiff(colnames(comparisons), req_comp)

  # Track per-block outcomes so failures are reported loudly at the end instead
  # of only as deferred warnings (which are easy to miss, and made whole cell
  # types silently disappear from the result).
  block_fail  <- character(0)
  block_empty <- character(0)

  results <- lapply(seq_len(nrow(comparisons)), function(ci) {
    cmp <- comparisons[ci, ]
    red_ct <- lapply(celltypes_red, function(ct) {
      # A cell type may pass the overall min_cells threshold yet still lack
      # enough cells in one condition of THIS comparison, so its pseudobulk
      # treatment/control columns can be missing -- the RED core stops in that
      # case. Skip (recording why) rather than aborting the whole run.
      red <- tryCatch(
        .RedScoresCellType(
          Peaks             = Peaks_all,
          counts_polyA      = pseudo,
          CellType          = ct,
          treatment         = cmp$treatment,
          control           = cmp$control,
          # Only filter Peaks by Condition when PASTA populated it per-comparison;
          # in the no-PASTA branch the same peak rows serve every comparison.
          Condition         = if (isTRUE(Pasta)) cmp$Condition else NULL,
          filter_method     = filter_method,
          mincounts         = mincounts,
          sample_metadata   = sample_metadata,
          covariate_formula = covariate_formula,
          gdpau_min_reads   = gdpau_min_reads
        ),
        error = function(e) {
          block_fail <<- c(block_fail,
                           paste0("[", canonical_by_san[[ct]], " / ", cmp$Condition, "] ",
                                  conditionMessage(e)))
          NULL
        }
      )
      if (is.null(red)) return(NULL)
      if (nrow(red) == 0) {
        block_empty <<- c(block_empty,
                          paste0(canonical_by_san[[ct]], " / ", cmp$Condition))
        return(NULL)
      }
      # Report the canonical (Seurat-object) cell-type name, not the dotted one.
      red$Cells     <- canonical_by_san[[ct]]
      red$Condition <- cmp$Condition
      for (col in extra_cols) red[[col]] <- cmp[[col]]
      red
    })
    dplyr::bind_rows(red_ct[!vapply(red_ct, is.null, logical(1))])
  })

  out <- dplyr::bind_rows(results)

  # Report skipped blocks up front -- a cell type missing from the result is
  # almost always one of these, not a real biological absence.
  if (length(block_fail) > 0) {
    n_show <- min(5L, length(block_fail))
    warning("DEPsMatrix: ", length(block_fail),
            " (cell type x comparison) block(s) errored and were skipped. First ",
            n_show, ":\n  ", paste(block_fail[seq_len(n_show)], collapse = "\n  "),
            call. = FALSE)
    if (verbose) {
      message("Skipped ", length(block_fail), " block(s) due to errors; first ",
              n_show, ":\n  ", paste(block_fail[seq_len(n_show)], collapse = "\n  "))
    }
  }
  if (length(block_empty) > 0 && verbose) {
    message("Note: ", length(block_empty),
            " block(s) ran but produced no gene passed filtering/pairing: ",
            paste(utils::head(block_empty, 10), collapse = "; "),
            if (length(block_empty) > 10) " ..." else "")
  }

  # Attach the per-peak RepeatMasker flag to each pair's proximal/distal site.
  # NA when repeat_masker = "none"; otherwise a clean logical -- TRUE if that
  # site overlaps a TE, FALSE if it doesn't (so in "Remove" mode every retained
  # row is FALSE). Only a genuinely absent site (p_peak NA on single-site rows)
  # stays NA.
  .te_flag <- function(pk) {
    if (is.null(te_peaks_all)) return(rep(NA, length(pk)))
    ifelse(is.na(pk), NA, pk %in% te_peaks_all)
  }
  if (nrow(out) > 0) {
    out$p_Repeated_Masker <- .te_flag(out$p_peak)
    out$d_Repeated_Masker <- .te_flag(out$d_peak)
  }

  # Drop columns that are entirely NA -- e.g. with Pasta = FALSE every PASTA
  # statistic column (Estimate/p.value/percent.*/...) is NA by construction, and
  # RED_LRT_cov_* is NA unless a covariate_formula was supplied. Keeps the table
  # readable; set drop_na_cols = FALSE to retain a fixed column layout.
  if (isTRUE(drop_na_cols) && nrow(out) > 0) {
    all_na <- vapply(out, function(col) all(is.na(col)), logical(1))
    if (any(all_na)) {
      if (verbose) {
        message("Dropping ", sum(all_na), " all-NA column(s): ",
                paste(names(out)[all_na], collapse = ", "))
      }
      out <- out[, !all_na, drop = FALSE]
    }
  }

  # Bring the identifying/label columns to the front: cell type, then the extra
  # comparison label columns (e.g. Strain, Treatment), then Condition.
  # any_of() keeps this safe if some are absent.
  out <- out %>%
    dplyr::relocate(dplyr::any_of(c("Cells", extra_cols, "Condition")))

  if (verbose) {
    n_out_ct <- if (nrow(out) > 0) length(unique(out$Cells)) else 0L
    message("DEPsMatrix complete: ", nrow(out), " RED rows from ", n_out_ct,
            " of ", length(celltypes_red), " cell type(s) across ",
            nrow(comparisons), " comparison(s).")
    if (n_out_ct < length(celltypes_red)) {
      absent <- setdiff(unname(canonical_by_san[celltypes_red]), unique(out$Cells))
      message("  Cell type(s) with no rows in the result: ",
              paste(absent, collapse = ", "))
    }
  }

  # A peak-level table (one row per peak) of exactly the sites that were
  # analyzed, with TE handling already applied (TE peaks removed in "Remove"
  # mode) and a per-peak `Repeated_Masker` flag. This is the `polyAdb` argument
  # PolyAPlot() expects -- returning it means callers don't have to keep a
  # separate Peaksdb() output around, and the plotted sites match the analyzed
  # set.
  polyAdb <- dplyr::distinct(Peaks_all, peak, .keep_all = TRUE)
  polyAdb$Repeated_Masker <- .te_flag(polyAdb$peak)

  # Always return a list so `polyAdb` (the analyzed peak-level table, one row per
  # peak, TE handling applied) is available regardless of `return_seurat`.
  # `return_seurat = TRUE` additionally includes the (residual-bearing, when
  # Pasta = TRUE) Seurat object for PolyAPlot() coverage tracks.
  if (isTRUE(return_seurat)) {
    return(list(red_scores = out, seu = seu, polyAdb = polyAdb))
  }
  list(red_scores = out, polyAdb = polyAdb)
}

#' Compute RED scores for a single cell type / comparison (internal RED core)
#'
#' The per-cell-type, per-comparison RED-matrix engine. Processes polyA site
#' data for one cell type to create pairwise comparisons between proximal (p)
#' and distal (d) polyA sites within each gene, restricted to sites in introns
#' or the 3' most exon, and calculates Relative Expression Difference (RED)
#' scores to quantify alternative polyadenylation changes.
#'
#' This used to be the exported \code{DEPsMatrix()} function. It is now an
#' internal helper called once per (cell type x comparison) by the top-level
#' \code{DEPsMatrix()} orchestrator, which additionally handles residual
#' calculation, differential testing, metadata merging, and pseudobulk
#' aggregation. The body is unchanged from the previous \code{DEPsMatrix()}.
#'
#' @param Peaks A \code{data.frame} containing raw polyA peak data to be processed.
#'   Expected columns include \code{symbol}, \code{ensembl_gene_id}, \code{strand},
#'   \code{hg38_Position}, \code{Intron.exon.location}, \code{PAS_hexamer},
#'   \code{PolyaStrength}, count columns for treatment and control conditions,
#'   and PASTA output columns (\code{estimate}, \code{percent.1}, \code{padj}).
#'   Genes are grouped by \code{ensembl_gene_id} (not \code{symbol}) since gene
#'   symbols can be duplicated across distinct Ensembl genes -- this is the
#'   column produced by \code{Peaksdb()}/the RNA assay's Ensembl join, NOT the
#'   \code{Ensembl_ID} column in this function's own output (see below).
#' @param counts_polyA Numeric. Pseudocount value added to polyA site counts
#'   before proportion and RED score calculations to avoid division by zero
#' @param CellType Character. Cell type to analyze. Must match a level present
#'   in the data. Default \code{NULL}
#' @param treatment Character. Prefix identifying the treatment condition's
#'   columns in \code{counts_polyA} (matched as \code{"<treatment>_<CellType>"}).
#'   Default \code{NULL}
#' @param control Character. Prefix identifying the control condition's
#'   columns in \code{counts_polyA} (matched as \code{"<control>_<CellType>"}).
#'   Default \code{NULL}
#' @param Condition Character. A value to filter \code{Peaks$Condition} by
#'   (\code{Peaks} is expected to already have a literal \code{Condition}
#'   column, e.g. set via \code{dplyr::mutate(Condition = "TreatmentA_vs_Control")}
#'   before calling this function) -- NOT the name of a column to look up.
#'   If \code{NULL}, no \code{Condition} filtering is applied. Default \code{NULL}
#'
#' @return A \code{data.frame} with pairwise proximal–distal comparisons between
#'   polyA sites within each gene. Each row represents one proximal (p) / distal (d)
#'   site pair, where the distal site is never intronic and has a higher positional
#'   index than the proximal site. For each RED type (\code{REDi}/\code{REDu}),
#'   only the pair with the greatest combined proportion change is retained.
#'   Returned columns:
#'   \describe{
#'     \item{\code{Gene_symbol}}{Gene symbol}
#'     \item{\code{Ensembl_ID}}{Ensembl gene ID -- the key genes are grouped by
#'       (copied from the input's \code{ensembl_gene_id} column)}
#'     \item{\code{strand}}{Genomic strand (\code{"+"} or \code{"-"})}
#'     \item{\code{p_site_id}, \code{d_site_id}}{Site indices within the gene,
#'       ordered 5' to 3'}
#'     \item{\code{p_peak}, \code{d_peak}}{Peak identifiers}
#'     \item{\code{p_pos}, \code{d_pos}}{Genomic positions (\code{hg38_Position})}
#'     \item{\code{p_region}, \code{d_region}}{Intron/exon location
#'       (from \code{Intron.exon.location})}
#'     \item{\code{p_region_label}, \code{d_region_label}}{Region labels}
#'     \item{\code{p_PAS_hexamer}, \code{d_PAS_hexamer}}{PAS hexamer sequences}
#'     \item{\code{p_PolyaStrength}, \code{d_PolyaStrength}}{PolyA signal
#'       strength scores}
#'     \item{\code{p_treatment_count}, \code{d_treatment_count}}{Raw polyA
#'       counts in the treatment condition}
#'     \item{\code{p_control_count}, \code{d_control_count}}{Raw polyA counts
#'       in the control condition}
#'     \item{\code{p_treatment_prop}, \code{d_treatment_prop}}{Proportions in
#'       the treatment condition}
#'     \item{\code{p_control_prop}, \code{d_control_prop}}{Proportions in the
#'       control condition}
#'     \item{\code{p_percent_change}, \code{d_percent_change}}{Absolute
#'       percentage change in proportions between conditions}
#'     \item{\code{p_estimate_PASTA}, \code{d_estimate_PASTA}}{PASTA log-fold
#'       change estimates}
#'     \item{\code{p_percent1_PASTA}, \code{d_percent1_PASTA}}{PASTA
#'       \code{percent.1} values}
#'     \item{\code{p_padj_PASTA}, \code{d_padj_PASTA}}{PASTA adjusted
#'       p-values}
#'     \item{\code{RED_type}}{\code{"REDi"} if the proximal site is intronic;
#'       \code{"REDu"} if the proximal site is exonic}
#'     \item{\code{RED}}{Relative Expression Difference score:
#'       \eqn{\log_2(d\_treatment / p\_treatment) - \log_2(d\_control / p\_control)}}
#'     \item{\code{Direction}}{\code{"Lengthened"} if \code{RED > 0};
#'       \code{"Shortened"} otherwise}
#'   }
#'   Genes with fewer than 2 polyA sites in introns or the 3' most exon, or
#'   without at least one valid proximal–distal pair, are excluded.
#'
#' @noRd
.RedScoresCellType <- function(
    Peaks,
    counts_polyA,
    CellType    = NULL,
    treatment   = NULL,
    control     = NULL,
    Condition   = NULL,
    filter_method = "edgeR",
    mincounts     = 10,
    sample_metadata   = NULL,
    covariate_formula = NULL,
    gdpau_min_reads   = 5,
    verbose           = FALSE
) {

  # --- Input validation ---
  # Note: the input column is "ensembl_gene_id" (lowercase, as produced by
  # Peaksdb()/the RNA assay's ENSEMBL join) -- not to be confused with this
  # function's own OUTPUT column "Ensembl_ID", which is just a
  # nicer-cased copy of the same value (set further down via
  # `Ensembl_ID = x` once split by ensembl_gene_id).
  required_peaks_cols <- c(
    "symbol", "ensembl_gene_id", "seqnames", "hg38_Position", "strand",
    "Intron.exon.location", "label", "peak",
    "PAS_hexamer", "PolyaStrength_percentile",
    "Estimate", "p.value", "p_val_adj", "percent.1", "percent.2", "CellType"
  )

  missing_cols <- setdiff(required_peaks_cols, colnames(Peaks))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in Peaks: ", paste(missing_cols, collapse = ", "))
  }

  if (is.null(CellType))   stop("`CellType` must be provided.")
  if (is.null(treatment))  stop("`treatment` must be provided.")
  if (is.null(control))    stop("`control` must be provided.")
  #if (is.null(Condition))  stop("`Condition` must be provided.")

  if (!CellType %in% unique(Peaks$CellType)) {
    stop("`CellType` '", CellType, "' not found in Peaks$CellType.")
  }

  if (!is.null(covariate_formula) && is.null(sample_metadata)) {
    stop("`sample_metadata` must be provided when `covariate_formula` is specified.")
  }
  if (!is.null(sample_metadata) && !all(c("Sample_key") %in% colnames(sample_metadata))) {
    stop("`sample_metadata` must contain a `Sample_key` column matching sanitized sample IDs.")
  }

  #if (!Condition %in% unique(Peaks$Condition)) {
  #  stop("`Condition` '", Condition, "' not found in Peaks$Condition.")
  #}

  # Validate counts_polyA has matching columns for treatment/control
  counts_obj <- if (CellType %in% names(counts_polyA)) {
    counts_polyA[[CellType]]
  } else {
    counts_polyA
  }

  has_treatment <- any(grepl(paste0(treatment, "_", CellType), colnames(counts_obj)))
  has_control   <- any(grepl(paste0(control,   "_", CellType), colnames(counts_obj)))

  if (!has_treatment) {
    stop("No column matching '", treatment, "_", CellType, "' found in counts_polyA.")
  }
  if (!has_control) {
    stop("No column matching '", control, "_", CellType, "' found in counts_polyA.")
  }

  if (!is.null(Condition)) {
    Peaks_filtered <- Peaks %>%
      dplyr::filter(
        .data$Condition == .env$Condition,
        .data$CellType  == .env$CellType )
  } else {
    Peaks_filtered <- Peaks %>%
      dplyr::filter(
        .data$CellType  == .env$CellType )
  }

  if (!is.null(covariate_formula)) {
    trt_orig_cols  <- grep(paste0(treatment, "_", CellType), colnames(counts_obj), value = TRUE)
    ctrl_orig_cols <- grep(paste0(control,   "_", CellType), colnames(counts_obj), value = TRUE)

    # Strip assay prefix (e.g. "polyA.") before extracting the sample token
    trt_orig_cols_clean  <- sub("^polyA\\.", "", trt_orig_cols)
    ctrl_orig_cols_clean <- sub("^polyA\\.", "", ctrl_orig_cols)

    trt_sample_ids  <- sub(paste0("^", treatment, "_", CellType, "_"), "", trt_orig_cols_clean)
    ctrl_sample_ids <- sub(paste0("^", control,   "_", CellType, "_"), "", ctrl_orig_cols_clean)

    trt_covariates  <- sample_metadata[match(trt_sample_ids,  sample_metadata$Sample_key), all.vars(stats::as.formula(paste("~", covariate_formula)))]
    ctrl_covariates <- sample_metadata[match(ctrl_sample_ids, sample_metadata$Sample_key), all.vars(stats::as.formula(paste("~", covariate_formula)))]

    # sanity check — should be no NAs after matching sample_metadata
    stopifnot(!anyNA(trt_covariates), !anyNA(ctrl_covariates))
  }

  counts_polyA <- counts_obj %>%
    dplyr::rename_with(
      ~ paste0("polyA_treatment_", seq_along(.x)),
      dplyr::matches(paste0(treatment, "_", CellType))
    ) %>%
    dplyr::rename_with(
      ~ paste0("polyA_control_", seq_along(.x)),
      dplyr::matches(paste0(control, "_", CellType))
    ) %>%
    dplyr::mutate(peak = rownames(.))

  n_trt  <- sum(grepl("^polyA_treatment_", colnames(counts_polyA)))
  n_ctrl <- sum(grepl("^polyA_control_",   colnames(counts_polyA)))

  PolyA_DEP <- Peaks_filtered %>%
    dplyr::left_join(counts_polyA, by = "peak") %>%
    dplyr::select(
      dplyr::all_of(colnames(Peaks_filtered)),
      dplyr::matches("^polyA_treatment_"),
      dplyr::matches("^polyA_control_")
    ) %>%
    dplyr::mutate(peak_id = peak) %>%
    tibble::column_to_rownames("peak_id")

  # Extract count matrix
  count_cols <- grep("^polyA_treatment_|^polyA_control_", colnames(PolyA_DEP), value = TRUE)
  count_mat  <- as.matrix(PolyA_DEP[, count_cols])

  # Filter
  if (filter_method == "edgeR") {
    group <- dplyr::case_when(
      grepl("^polyA_treatment_", count_cols) ~ "treatment",
      grepl("^polyA_control_",   count_cols) ~ "control"
    )
    keep <- edgeR::filterByExpr(count_mat, group = group)

  } else if (filter_method == "manual") {
    keep <- rowSums(count_mat) >= mincounts

  } else {
    stop("filter_method must be 'edgeR' or 'manual'")
  }

  if (verbose) message(sprintf("Keeping %d / %d peaks after '%s' filtering", sum(keep), nrow(count_mat), filter_method))

  PolyA_DEP <- PolyA_DEP[keep, ]

  # Drop peaks with no usable gene id BEFORE splitting. split() would otherwise
  # create a group named "" (or drop NA rows), and R's `[[` never matches an
  # empty or NA name -- so data_list[[""]] returns NULL and the per-gene
  # `NULL %>% filter(...)` errors with "no applicable method for 'filter'
  # applied to an object of class NULL", killing the whole block. Unannotated
  # peaks can't be assigned to a proximal/distal pair anyway.
  n_pre_gene <- nrow(PolyA_DEP)
  PolyA_DEP <- PolyA_DEP[!is.na(PolyA_DEP$ensembl_gene_id) &
                           trimws(as.character(PolyA_DEP$ensembl_gene_id)) != "", , drop = FALSE]
  if (verbose && nrow(PolyA_DEP) < n_pre_gene) {
    message(sprintf("Dropping %d peak(s) with no ensembl_gene_id",
                    n_pre_gene - nrow(PolyA_DEP)))
  }
  if (nrow(PolyA_DEP) == 0) {
    return(tibble::tibble())
  }

  # Split by ensembl_gene_id rather than gene symbol: symbols can be
  # duplicated across distinct Ensembl genes, which would otherwise
  # incorrectly merge unrelated genes into a single proximal/distal analysis
  # unit. Each split group's name (`x` below) becomes the output's
  # "Ensembl_ID" column.
  data_list <- PolyA_DEP %>%
    split(.$ensembl_gene_id)
  # Guard against empty groups (unused factor levels, etc.)
  data_list <- data_list[vapply(data_list, function(d) !is.null(d) && nrow(d) > 0, logical(1))]

  # --- Dynamic column name builders ---
  trt_count_cols  <- paste0("polyA_treatment_", 1:n_trt)
  ctrl_count_cols <- paste0("polyA_control_",   1:n_ctrl)

  p_nread_trt_cols  <- paste0("p_nread_treatment_", 1:n_trt)
  p_nread_ctrl_cols <- paste0("p_nread_control_",   1:n_ctrl)
  d_nread_trt_cols  <- paste0("d_nread_treatment_", 1:n_trt)
  d_nread_ctrl_cols <- paste0("d_nread_control_",   1:n_ctrl)

  # Use a progress-bar lapply only when verbose; otherwise plain lapply so the
  # per-(cell type x comparison) "0% ~calculating" bars don't flood the console.
  .red_lapply <- if (verbose) pbapply::pblapply else lapply
  # Iterate by POSITION (not by name): looking elements up by name breaks for
  # any group whose name R's `[[` can't match (empty string / NA). `x` is still
  # the gene id, used below for the Ensembl_ID column.
  gene_ids <- names(data_list)
  Matrix <- .red_lapply(seq_along(data_list), function(gi){
    x  <- gene_ids[gi]
    dl <- data_list[[gi]]
    if (is.null(dl) || nrow(dl) == 0) return(NULL)

    polyA_sites <- dl %>%
      dplyr::filter(Intron.exon.location %in% c("3' most exon","Intron"))

    # Single-site gene: return as distal-only row with NA for proximal fields
    if(nrow(polyA_sites) == 1) {
      site <- polyA_sites %>%
        dplyr::mutate(
          dplyr::across(dplyr::matches("^polyA_treatment_"),
                        ~ .x / sum(.x), .names = "prop_{.col}"),
          dplyr::across(dplyr::matches("^polyA_control_"),
                        ~ .x / sum(.x), .names = "prop_{.col}"),
          prop_treatment  = rowMeans(dplyr::pick(dplyr::matches("^prop_polyA_treatment_")), na.rm = TRUE),
          prop_control    = rowMeans(dplyr::pick(dplyr::matches("^prop_polyA_control_")), na.rm = TRUE),
          prop_change     = abs(prop_treatment - prop_control),
          prop_pct_change = prop_change * 100,
          site_id = 1L, Ensembl_ID = x
        )
      return(tibble::tibble(
        Gene_symbol       = site$symbol,
        Ensembl_ID        = x,
        chr               = site$seqnames,
        strand            = site$strand,
        npas              = 1L,
        RED_type          = dplyr::if_else(site$Intron.exon.location == "Intron", "REDi", "REDu"),
        RLD               = NA_real_,
        RLD_direction     = NA_character_,
        RED               = NA_real_,
        RED_Direction     = NA_character_,
        RED_fisher_pval   = NA_real_,
        RED_fisher_padj   = NA_real_,
        RED_LRT_pval      = NA_real_,
        RED_LRT_cov_pval  = NA_real_,
        PASTA_pval        = NA_real_,
        PASTA_padj        = NA_real_,
        p_site_id         = NA_integer_,
        d_site_id         = 1L,
        p_peak            = NA_character_,
        d_peak            = site$peak,
        p_pos             = NA_real_,
        d_pos             = site$hg38_Position,
        p_region          = NA_character_,
        d_region          = site$Intron.exon.location,
        p_PAS_hexamer     = NA_character_,
        d_PAS_hexamer     = site$PAS_hexamer,
        p_PolyaStrength   = NA_real_,
        d_PolyaStrength   = site$PolyaStrength_percentile,
        p_estimate_PASTA  = NA_real_,
        d_estimate_PASTA  = site$Estimate,
        p_pval_PASTA      = NA_real_,
        d_pval_PASTA      = site$p.value,
        p_padj_PASTA      = NA_real_,
        d_padj_PASTA      = site$p_val_adj,
        p_percent1_PASTA  = NA_real_,
        d_percent1_PASTA  = site$percent.1,
        p_percent2_PASTA  = NA_real_,
        d_percent2_PASTA  = site$percent.2,
        p_treatment_prop  = NA_real_,
        d_treatment_prop  = site$prop_treatment,
        p_control_prop    = NA_real_,
        d_control_prop    = site$prop_control,
        p_percent_change  = NA_real_,
        d_percent_change  = site$prop_pct_change,
        # gDPAU/DPAU are undefined for a single-site gene (need a pair / n >= 2).
        gDPAU             = NA_real_,
        DPAU              = NA_real_
      ))
    }

    # Calculate proportions
    gene_processed <- polyA_sites %>%
      dplyr::arrange(dplyr::if_else(strand == "+", hg38_Position, dplyr::desc(hg38_Position))) %>%
      dplyr::mutate(
        # Per-replicate proportions
        dplyr::across(dplyr::matches("^polyA_treatment_"),
               ~ .x / sum(.x), .names = "prop_{.col}"),
        dplyr::across(dplyr::matches("^polyA_control_"),
               ~ .x / sum(.x), .names = "prop_{.col}"),
        # Average the per-replicate proportions
        prop_treatment  = rowMeans(dplyr::pick(dplyr::matches("^prop_polyA_treatment_")), na.rm = TRUE),
        prop_control    = rowMeans(dplyr::pick(dplyr::matches("^prop_polyA_control_")), na.rm = TRUE),
        prop_change     = abs(prop_treatment - prop_control),
        prop_pct_change = prop_change * 100
      ) %>%
      dplyr::mutate(site_id = dplyr::row_number(), Ensembl_ID = x)

    # --- gDPAU (Wang et al. 2022) ------------------------------------------
    # Location-index-weighted distal-usage index over ALL of this gene's polyA
    # sites (ordered 5'->3' by site_id above):
    #   gDPAU = 1/(n-1) * sum_{i=1..n} (i-1) * p_i,  p_i = site's read fraction.
    # p_i uses STACKED counts within the group (summed across replicates, the
    # paper's convention), NOT the mean of per-replicate proportions. When n = 2
    # this reduces to DPAU (distal-site fraction). If a group's total supporting
    # reads across the gene's sites are < gdpau_min_reads, its gDPAU is NA
    # (noise), and so is gDPAU_diff.
    .gdpau <- function(count_cols) {
      site_reads <- rowSums(as.matrix(gene_processed[, count_cols, drop = FALSE]), na.rm = TRUE)
      tot <- sum(site_reads)
      n_s <- length(site_reads)
      if (n_s < 2 || tot < gdpau_min_reads) return(NA_real_)
      p_i <- site_reads / tot
      sum((gene_processed$site_id - 1) * p_i) / (n_s - 1)
    }
    trt_count_cols_g  <- grep("^polyA_treatment_", names(gene_processed), value = TRUE)
    ctrl_count_cols_g <- grep("^polyA_control_",   names(gene_processed), value = TRUE)
    gDPAU_treatment_g <- .gdpau(trt_count_cols_g)
    gDPAU_control_g   <- .gdpau(ctrl_count_cols_g)
    gDPAU_diff_g      <- gDPAU_treatment_g - gDPAU_control_g

    # --- DPAU (classic 2-site index) ---------------------------------------
    # The "normal" distal-usage index for the specific proximal/distal PAIR in
    # a row: distal / (distal + proximal) reads (stacked within the group). This
    # is per-tested-pair, so it's defined for every row; for a gene with exactly
    # 2 sites the pair spans the whole gene and DPAU == gDPAU. NA if the pair's
    # total reads in a group are < gdpau_min_reads.
    .dpau <- function(count_cols, p_sid, d_sid) {
      if (is.na(p_sid) || is.na(d_sid)) return(NA_real_)
      p_row <- which(gene_processed$site_id == p_sid)
      d_row <- which(gene_processed$site_id == d_sid)
      if (length(p_row) != 1 || length(d_row) != 1) return(NA_real_)
      # Use as.matrix() (like .gdpau above), NOT as.numeric(): gene_processed is
      # a data.frame, so gene_processed[row, <several cols>] returns a
      # data.frame and as.numeric() on it errors with "'list' object cannot be
      # coerced to type 'double'" -- an error that escapes the per-gene loop and
      # kills the whole (cell type x comparison) block.
      cm <- as.matrix(gene_processed[, count_cols, drop = FALSE])
      p_reads <- sum(cm[p_row, ], na.rm = TRUE)
      d_reads <- sum(cm[d_row, ], na.rm = TRUE)
      tot <- p_reads + d_reads
      if (!is.finite(tot) || tot < gdpau_min_reads) return(NA_real_)
      d_reads / tot
    }

    # RLDu: from gene_processed, UTR sites only
    utr_sites <- gene_processed %>%
      dplyr::filter(Intron.exon.location == "3' most exon") %>%
      dplyr::arrange(dplyr::if_else(strand == "+", hg38_Position, -hg38_Position)) %>%
      dplyr::mutate(l_utr = ifelse(dplyr::row_number() == 1, 1, abs(hg38_Position - hg38_Position[1])))

    RLDu <- if(nrow(utr_sites) >= 2){
      log2( sum(utr_sites$prop_treatment * utr_sites$l_utr) / sum(utr_sites$prop_control   * utr_sites$l_utr) )
    } else { NA_real_ }

    # RLDi: intron sites vs exon sites, global proportions
    intron_sites <- gene_processed %>%
      dplyr::filter(Intron.exon.location == "Intron") %>%
      dplyr::arrange(dplyr::if_else(strand == "+", hg38_Position, -hg38_Position)) %>%
      dplyr::mutate(l = abs(hg38_Position - dplyr::if_else(strand == "+", gene_start, gene_end)))

    exon_sites <- gene_processed %>%
      dplyr::filter(Intron.exon.location == "3' most exon") %>%
      dplyr::arrange(dplyr::if_else(strand == "+", hg38_Position, -hg38_Position)) %>%
      dplyr::mutate(l = abs(hg38_Position - dplyr::if_else(strand == "+", gene_start, gene_end)))

    RLDi <- if(nrow(intron_sites) >= 1 && nrow(exon_sites) >= 1){
      log2(sum(exon_sites$prop_treatment * exon_sites$l) / sum(intron_sites$prop_treatment * intron_sites$l)) -
        log2(sum(exon_sites$prop_control * exon_sites$l) / sum(intron_sites$prop_control   * intron_sites$l))
    } else { NA_real_ }

    # Need at least one non-intron site to be distal
    if(sum(gene_processed$label != "Intron") < 1) return(NULL)

    # Get all ordered pairs where d_idx > p_idx and d is not an intron
    pairs <- expand.grid(p = gene_processed$site_id, d = gene_processed$site_id) %>%
      dplyr::filter(
        d > p,
        gene_processed$label[d] != "Intron"   # introns are never distal
      ) %>%
      dplyr::arrange(dplyr::desc(d), dplyr::desc(p))

    # Need at least one valid pair after filtering
    if(nrow(pairs) < 1) return(NULL)

    pairwise_comparisons <- tibble::tibble(
      Gene_symbol = gene_processed$symbol[pairs$d],
      Ensembl_ID  = gene_processed$Ensembl_ID[pairs$d],
      chr         = gene_processed$seqnames[pairs$d],
      p_pos       = gene_processed$hg38_Position[pairs$p],
      d_pos       = gene_processed$hg38_Position[pairs$d],
      strand      = gene_processed$strand[pairs$d],
      npas        = nrow(gene_processed),

      p_site_id = gene_processed$site_id[pairs$p],
      d_site_id = gene_processed$site_id[pairs$d],

      p_peak = gene_processed$peak[pairs$p],
      d_peak = gene_processed$peak[pairs$d],

      p_region = gene_processed$Intron.exon.location[pairs$p],
      d_region = gene_processed$Intron.exon.location[pairs$d],

      p_PAS_hexamer = gene_processed$PAS_hexamer[pairs$p],
      d_PAS_hexamer = gene_processed$PAS_hexamer[pairs$d],

      p_PolyaStrength = gene_processed$PolyaStrength_percentile[pairs$p],
      d_PolyaStrength = gene_processed$PolyaStrength_percentile[pairs$d],

      p_treatment_count = gene_processed$polyA_treatment[pairs$p],
      d_treatment_count = gene_processed$polyA_treatment[pairs$d],

      p_control_count = gene_processed$polyA_control[pairs$p],
      d_control_count = gene_processed$polyA_control[pairs$d],

      p_treatment_prop = gene_processed$prop_treatment[pairs$p],
      d_treatment_prop = gene_processed$prop_treatment[pairs$d],

      p_control_prop = gene_processed$prop_control[pairs$p],
      d_control_prop = gene_processed$prop_control[pairs$d],

      p_percent_change = gene_processed$prop_pct_change[pairs$p],
      d_percent_change = gene_processed$prop_pct_change[pairs$d],

      p_estimate_PASTA = gene_processed$Estimate[pairs$p],
      d_estimate_PASTA = gene_processed$Estimate[pairs$d],

      p_pval_PASTA = gene_processed$p.value[pairs$p],
      d_pval_PASTA = gene_processed$p.value[pairs$d],

      p_padj_PASTA = gene_processed$p_val_adj[pairs$p],
      d_padj_PASTA = gene_processed$p_val_adj[pairs$d],

      p_percent1_PASTA = gene_processed$percent.1[pairs$p],
      d_percent1_PASTA = gene_processed$percent.1[pairs$d],

      p_percent2_PASTA = gene_processed$percent.2[pairs$p],
      d_percent2_PASTA = gene_processed$percent.2[pairs$d],

    ) %>%
      # Per-replicate counts added via bind_cols
      dplyr::bind_cols(
        purrr::map_dfc(seq_len(n_trt), ~ tibble::tibble(
          !!paste0("p_nread_treatment_", .x) := gene_processed[[paste0("polyA_treatment_", .x)]][pairs$p]
        )),
        purrr::map_dfc(seq_len(n_ctrl), ~ tibble::tibble(
          !!paste0("p_nread_control_", .x) := gene_processed[[paste0("polyA_control_", .x)]][pairs$p]
        )),
        purrr::map_dfc(seq_len(n_trt), ~ tibble::tibble(
          !!paste0("d_nread_treatment_", .x) := gene_processed[[paste0("polyA_treatment_", .x)]][pairs$d]
        )),
        purrr::map_dfc(seq_len(n_ctrl), ~ tibble::tibble(
          !!paste0("d_nread_control_", .x) := gene_processed[[paste0("polyA_control_", .x)]][pairs$d]
        ))
      )

    pairwise_comparisons <- pairwise_comparisons %>%
      dplyr::mutate( RED_type = dplyr::if_else( p_region == "Intron", "REDi", "REDu") )

    pairwise_comparisons <- pairwise_comparisons %>%
      dplyr::group_by(Ensembl_ID, RED_type) %>%
      dplyr::slice_max(p_percent_change + d_percent_change, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup()

    pairwise_comparisons <- pairwise_comparisons %>%
      dplyr::mutate( RED = log2(d_treatment_prop/p_treatment_prop) - log2(d_control_prop/p_control_prop),
                     RED = ifelse(is.finite(RED), RED, NA_real_))

    pairwise_comparisons <- pairwise_comparisons %>%
      dplyr::mutate(
        RED_fisher_pval = sapply(seq_len(nrow(pairwise_comparisons)), function(i) {
          compute_statistics(p_site_id[i], d_site_id[i], gene_processed )$fisher_pval }),
        RED_ttest_pval = sapply(seq_len(nrow(pairwise_comparisons)), function(i) {
          compute_statistics(p_site_id[i], d_site_id[i], gene_processed )$ttest_pval }),
        RED_LRT_pval = sapply(seq_len(nrow(pairwise_comparisons)), function(i) {
          compute_statistics_LRT(p_site_id[i], d_site_id[i], gene_processed )$lrt_pval }),
        RED_LRT_cov_pval = if (!is.null(covariate_formula)) {
          sapply(seq_len(nrow(pairwise_comparisons)), function(i) {
            compute_statistics_LRT_cov(
              p_site_id[i], d_site_id[i], gene_processed,
              trt_covariates, ctrl_covariates,
              covariate_formula = covariate_formula
            )$lrt_pval })
        } else { NA_real_ },
      )

    pairwise_comparisons <- pairwise_comparisons %>%
      dplyr::mutate( PASTA_pval = pchisq(-2 * (log(p_pval_PASTA) + log(d_pval_PASTA)), df = 4, lower.tail = FALSE) )

    pairwise_comparisons <- pairwise_comparisons %>%
      dplyr::mutate( RED_Direction = ifelse(RED > 0, "Lengthened", "Shortened") )

    pairwise_comparisons <- pairwise_comparisons %>%
      dplyr::mutate(
        RLD  = dplyr::if_else(RED_type == "REDu", RLDu, RLDi),
        RLD_direction = dplyr::if_else(RLD > 0, "Lengthened", "Shortened")
      )

    # Report only the treatment - control difference for each index: gDPAU is
    # gene-level (same across this gene's REDi/REDu rows); DPAU is per-pair
    # (distal/(distal+proximal) of the row's two sites). NA if either group is
    # below gdpau_min_reads.
    pairwise_comparisons <- pairwise_comparisons %>%
      dplyr::mutate(
        gDPAU = gDPAU_diff_g,
        DPAU  = vapply(seq_len(dplyr::n()),
                       function(i) .dpau(trt_count_cols_g,  p_site_id[i], d_site_id[i]) -
                                   .dpau(ctrl_count_cols_g, p_site_id[i], d_site_id[i]),
                       numeric(1))
      )

    return(pairwise_comparisons)
  })

  # --- Guard: si no hay resultados, avisar claramente en vez de fallar en group_by ---
  Matrix <- Matrix[!sapply(Matrix, is.null)]

  if (length(Matrix) == 0) {
    warning("No valid gene-level results for CellType = '", CellType,
            "'. This CellType may have too few peaks/samples after filtering. Returning empty tibble.")
    return(tibble::tibble())
  }

  Matrix <- Matrix %>%
    dplyr::bind_rows() %>%
    dplyr::group_by(RED_type) %>%
    dplyr::bind_rows() %>%
    # BH correction across all genes, separately per test and RED type
    dplyr::group_by(RED_type) %>%
    dplyr::mutate(
      RED_fisher_padj    = p.adjust(RED_fisher_pval,     method = "BH"),
      RED_LRT_padj       = p.adjust(RED_LRT_pval,         method = "BH"),
      RED_LRT_cov_padj   = p.adjust(RED_LRT_cov_pval,     method = "BH"),  # <-- add
      PASTA_padj         = p.adjust(PASTA_pval,           method = "BH"),
      # Bug fix: `Condition` was accepted as an argument and used to FILTER
      # the input `Peaks`, but was never carried into the output -- so
      # `deg_list`s built from this function (e.g. Red_scores_All2 in G3)
      # never actually had a `Condition` column, even though PolyAPlot()
      # unconditionally does `deg_list %>% dplyr::select(p_peak, d_peak,
      # Condition)`, which errors ("Condition doesn't exist") the moment
      # highlight/shared-peak labeling runs. Stamping the argument's value
      # onto every row here closes that gap; NA if Condition wasn't supplied.
      Condition          = if (is.null(Condition)) NA_character_ else Condition
    ) %>%
    dplyr::ungroup() %>%
    # --- Final column ordering ---
    dplyr::select(
      Gene_symbol, Ensembl_ID, chr, strand, npas,
      RED_type, RLD, RLD_direction,
      RED, RED_Direction,
      DPAU, gDPAU,
      RED_LRT_pval,     RED_LRT_padj,
      RED_LRT_cov_pval, RED_LRT_cov_padj,
      PASTA_pval,       PASTA_padj,
      Condition,
      p_site_id, d_site_id,
      p_peak, d_peak,
      p_pos, d_pos,
      p_region, d_region,
      p_PAS_hexamer, d_PAS_hexamer,
      p_PolyaStrength, d_PolyaStrength,
      p_estimate_PASTA, d_estimate_PASTA,
      p_pval_PASTA, d_pval_PASTA,
      p_padj_PASTA, d_padj_PASTA,
      p_percent1_PASTA, d_percent1_PASTA,
      p_percent2_PASTA, d_percent2_PASTA,
      p_treatment_prop, d_treatment_prop,
      p_control_prop,   d_control_prop,
      p_percent_change, d_percent_change
      #dplyr::all_of(p_nread_trt_cols),
      #dplyr::all_of(p_nread_ctrl_cols),
      #dplyr::all_of(d_nread_trt_cols),
      #dplyr::all_of(d_nread_ctrl_cols)
    )

  return(Matrix)
}



#' Compute statistical tests for a proximal–distal polyA site pair
#'
#' For a given proximal (p) and distal (d) polyA site pair within a single gene,
#' computes two statistical tests quantifying differential site usage between
#' treatment and control conditions: a Fisher's exact test on pooled replicate
#' counts, and an unpaired two-sample t-test on per-replicate log2 proportion
#' ratios.
#'
#' @param p_idx Integer. Row index of the proximal polyA site in
#'   \code{gene_processed}
#' @param d_idx Integer. Row index of the distal polyA site in
#'   \code{gene_processed}
#' @param gene_processed \code{data.frame} containing polyA site data for a
#'   single gene. Must contain columns matching \code{^polyA_treatment_} and
#'   \code{^polyA_control_} for per-replicate counts
#'
#' @return A named \code{list} with two elements:
#'   \describe{
#'     \item{\code{fisher_pval}}{P-value from Fisher's exact test applied to a
#'       2×2 contingency matrix of pooled proximal/distal counts across
#'       treatment and control replicates. Returns \code{NA_real_} if the test
#'       errors or produces a warning}
#'     \item{\code{ttest_pval}}{P-value from an unpaired two-sample t-test
#'       comparing per-replicate \eqn{\log_2(d\_prop / p\_prop)} ratios between
#'       treatment and control. Returns \code{NA_real_} if either group has
#'       fewer than 2 non-\code{NA} observations, or if the test errors or
#'       produces a warning}
#'   }
#'
#' @keywords internal
#'
compute_statistics <- function(p_idx, d_idx, gene_processed) {
  trt_cols  <- grep("^polyA_treatment_", colnames(gene_processed), value = TRUE)
  ctrl_cols <- grep("^polyA_control_",   colnames(gene_processed), value = TRUE)
  p_counts_trt  <- unlist(gene_processed[p_idx, trt_cols])
  d_counts_trt  <- unlist(gene_processed[d_idx, trt_cols])
  p_counts_ctrl <- unlist(gene_processed[p_idx, ctrl_cols])
  d_counts_ctrl <- unlist(gene_processed[d_idx, ctrl_cols])
  # --- Fisher test on pooled counts across replicates ---
  nr <- matrix(
    c(sum(p_counts_trt, na.rm = TRUE), sum(p_counts_ctrl, na.rm = TRUE),
      sum(d_counts_trt, na.rm = TRUE), sum(d_counts_ctrl, na.rm = TRUE)),
    nrow = 2,
    dimnames = list(c("proximal", "distal"), c("treatment", "control"))
  )
  fisher_pval <- tryCatch(
    stats::fisher.test(nr)$p.value,
    error   = function(e) NA_real_,
    warning = function(w) NA_real_
  )
  # --- Unpaired two-sample t-test on per-replicate log2 ratios ---
  total_trt  <- colSums(gene_processed[, trt_cols], na.rm = TRUE)
  total_ctrl <- colSums(gene_processed[, ctrl_cols], na.rm = TRUE)
  p_prop_trt  <- p_counts_trt  / total_trt
  d_prop_trt  <- d_counts_trt  / total_trt
  p_prop_ctrl <- p_counts_ctrl / total_ctrl
  d_prop_ctrl <- d_counts_ctrl / total_ctrl
  log2ratio_trt  <- log2(d_prop_trt  / p_prop_trt)
  log2ratio_ctrl <- log2(d_prop_ctrl / p_prop_ctrl)
  ttest_pval <- tryCatch({
    if (sum(!is.na(log2ratio_trt)) >= 2 & sum(!is.na(log2ratio_ctrl)) >= 2) {
      stats::t.test(log2ratio_trt, log2ratio_ctrl, paired = FALSE)$p.value
    } else { NA_real_ }
  },
  error   = function(e) NA_real_,
  warning = function(w) NA_real_
  )
  return(list(
    fisher_pval = fisher_pval,
    ttest_pval  = ttest_pval
  ))
}

#' Compute a covariate-adjusted LRT p-value for a proximal–distal polyA site pair
#'
#' Same as compute_statistics_LRT(), but allows arbitrary per-replicate
#' covariates to be included in both the full and null models, so that the
#' likelihood ratio test isolates the condition effect while adjusting for
#' covariates (e.g., Age, Pmi, batch).
#'
#' @param p_idx,d_idx,gene_processed As in compute_statistics_LRT()
#' @param trt_covariates data.frame, one row per treatment replicate, in the
#'   same order as the treatment columns in gene_processed
#' @param ctrl_covariates data.frame, one row per control replicate, in the
#'   same order as the control columns in gene_processed
#' @param covariate_formula Character. RHS terms to add, e.g. "Age + Pmi"
#'
#' @keywords internal
compute_statistics_LRT_cov <- function(p_idx, d_idx, gene_processed,
                                       trt_covariates, ctrl_covariates,
                                       covariate_formula) {
  trt_cols  <- grep("^polyA_treatment_", colnames(gene_processed), value = TRUE)
  ctrl_cols <- grep("^polyA_control_",   colnames(gene_processed), value = TRUE)

  p_counts_trt  <- unlist(gene_processed[p_idx, trt_cols])
  d_counts_trt  <- unlist(gene_processed[d_idx, trt_cols])
  p_counts_ctrl <- unlist(gene_processed[p_idx, ctrl_cols])
  d_counts_ctrl <- unlist(gene_processed[d_idx, ctrl_cols])

  lrt_pval <- tryCatch({

    df <- data.frame(
      d_count   = c(d_counts_trt, d_counts_ctrl),
      p_count   = c(p_counts_trt, p_counts_ctrl),
      condition = c(rep("treatment", length(d_counts_trt)),
                    rep("control",   length(d_counts_ctrl))),
      rbind(trt_covariates, ctrl_covariates)
    )
    df$condition <- factor(df$condition, levels = c("control", "treatment"))
    df <- df[stats::complete.cases(df) & (df$d_count + df$p_count) > 0, ]

    if (nrow(df) < 3 || length(unique(df$condition)) < 2) {
      NA_real_
    } else {
      full_f <- stats::as.formula(paste("cbind(d_count, p_count) ~ condition +", covariate_formula))
      null_f <- stats::as.formula(paste("cbind(d_count, p_count) ~",            covariate_formula))

      fit_full <- stats::glm(full_f, data = df, family = stats::binomial())
      fit_null <- stats::glm(null_f, data = df, family = stats::binomial())

      if (!fit_full$converged || !fit_null$converged) {
        NA_real_
      } else {
        lrt <- stats::anova(fit_null, fit_full, test = "LRT")
        lrt[["Pr(>Chi)"]][2]
      }
    }
  },
  error   = function(e) NA_real_,
  warning = function(w) NA_real_
  )

  return(list(lrt_pval = lrt_pval))
}

#' Compute a likelihood ratio test p-value for a proximal–distal polyA site pair
#'
#' For a given proximal (p) and distal (d) polyA site pair within a single gene,
#' fits a binomial GLM of site usage (\code{cbind(d_counts, p_counts) ~ condition})
#' using per-replicate counts, and compares it against an intercept-only null model
#' via a likelihood ratio test. Unlike \code{compute_statistics()}'s Fisher's exact
#' test (which pools replicates into a single 2x2 table) or its t-test (which
#' operates on derived log2 ratios), this test uses the raw per-replicate counts
#' directly in a likelihood framework.
#'
#' @param p_idx Integer. Row index of the proximal polyA site in
#'   \code{gene_processed}
#' @param d_idx Integer. Row index of the distal polyA site in
#'   \code{gene_processed}
#' @param gene_processed \code{data.frame} containing polyA site data for a
#'   single gene. Must contain columns matching \code{^polyA_treatment_} and
#'   \code{^polyA_control_} for per-replicate counts
#'
#' @return A named \code{list} with one element:
#'   \describe{
#'     \item{\code{lrt_pval}}{P-value from the likelihood ratio test comparing
#'       a binomial GLM with condition as a predictor against the null model.
#'       Returns \code{NA_real_} if there are fewer than 2 usable replicates
#'       total, if all counts are zero, or if model fitting errors/warns}
#'   }
#'
#' @keywords internal
#'
compute_statistics_LRT <- function(p_idx, d_idx, gene_processed) {
  trt_cols  <- grep("^polyA_treatment_", colnames(gene_processed), value = TRUE)
  ctrl_cols <- grep("^polyA_control_",   colnames(gene_processed), value = TRUE)

  p_counts_trt  <- unlist(gene_processed[p_idx, trt_cols])
  d_counts_trt  <- unlist(gene_processed[d_idx, trt_cols])
  p_counts_ctrl <- unlist(gene_processed[p_idx, ctrl_cols])
  d_counts_ctrl <- unlist(gene_processed[d_idx, ctrl_cols])

  lrt_pval <- tryCatch({

    df <- data.frame(
      d_count   = c(d_counts_trt, d_counts_ctrl),
      p_count   = c(p_counts_trt, p_counts_ctrl),
      condition = c(rep("treatment", length(d_counts_trt)),
                    rep("control",   length(d_counts_ctrl)))
    )
    df$condition <- factor(df$condition, levels = c("control", "treatment"))

    # Drop replicates with no reads at either site (uninformative)
    df <- df[stats::complete.cases(df) & (df$d_count + df$p_count) > 0, ]

    if (nrow(df) < 2 || length(unique(df$condition)) < 2) {
      NA_real_
    } else {
      fit_full <- stats::glm(
        cbind(d_count, p_count) ~ condition,
        data   = df,
        family = stats::binomial()
      )
      fit_null <- stats::glm(
        cbind(d_count, p_count) ~ 1,
        data   = df,
        family = stats::binomial()
      )

      lrt <- stats::anova(fit_null, fit_full, test = "LRT")
      lrt[["Pr(>Chi)"]][2]
    }
  },
  error   = function(e) NA_real_,
  warning = function(w) NA_real_
  )

  return(list(
    lrt_pval = lrt_pval
  ))
}

#' Export cell barcodes by condition and cell type
#'
#' Iterates over all combinations of condition and cell type in a Seurat object
#' and writes the corresponding barcodes to plain-text files, one file per
#' combination. Output files are organized into subdirectories by condition
#' under \code{output_dir}. Intended as a preprocessing step before
#' barcode-level BAM subsetting (e.g., with \code{subset-bam}).
#'
#' @param seurat_obj A \code{Seurat} object containing metadata columns for
#'   condition, cell type, and barcode
#' @param output_dir Character. Path to the root output directory. Created
#'   recursively if it does not exist. Default \code{"barcodes"}
#' @param celltype_col Character. Name of the metadata column containing cell
#'   type labels. Default \code{"CellType"}
#' @param barcode_col Character. Name of the metadata column containing cell
#'   barcodes. Default \code{"Barcode"}
#' @param group Character. Name of the metadata column used to split cells into
#'   conditions (e.g., \code{"Treatment_Strain"}). One subdirectory is created
#'   per unique value. Default \code{"Treatment_Strain"}
#'
#' @return Invisibly returns a nested \code{list} of the form
#'   \code{exported[[condition]][[cell_type]]}, where each element contains:
#'   \describe{
#'     \item{\code{file}}{Full path to the written barcode file}
#'     \item{\code{n}}{Number of barcodes written}
#'   }
#'   Condition–cell type combinations with zero barcodes are skipped with a
#'   \code{warning()} and are absent from the returned list.
#'
#' @details
#' Output files are named \code{<condition>_<cell_type>_barcodes.txt} and
#' written with \code{write.table()} (no row names, no column names, no
#' quotes), one barcode per line. Progress is reported via \code{message()}
#' for each successfully written file.
#'
#' @seealso \code{\link{ExportBarcodes}} is typically followed by a shell call
#'   to \code{subset-bam} to extract per-cell-type BAM files for tools such as
#'   MAAPER or PASTA.
#'
#' @export
#'
ExportBarcodes <- function(seurat_obj,
                           output_dir   = "barcodes",
                           celltype_col = "CellType",
                           barcode_col  = "Barcode",
                           group        = "Treatment_Strain") {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  exported  <- list()
  conditions <- unique(seurat_obj@meta.data[[group]])
  for (cond in conditions) {
    # Subset the Seurat object by condition
    obj <- subset(seurat_obj, subset = !!rlang::sym(group) == cond)
    study_dir <- file.path(output_dir, cond)
    if (!dir.exists(study_dir)) dir.create(study_dir)
    cell_types <- levels(as.factor(obj@meta.data[[celltype_col]]))
    for (ct in cell_types) {
      obj_sub  <- subset(obj, subset = !!rlang::sym(celltype_col) == ct)
      barcodes <- obj_sub@meta.data[[barcode_col]]
      if (length(barcodes) == 0) {
        warning(sprintf("[%s | %s] No barcodes found, skipping.", cond, ct))
        next
      }
      fname <- file.path(study_dir, paste0(cond, "_", ct, "_barcodes.txt"))
      utils::write.table(
        data.frame(barcodes),
        file      = fname,
        row.names = FALSE,
        col.names = FALSE,
        quote     = FALSE
      )
      exported[[cond]][[ct]] <- list(file = fname, n = length(barcodes))
      message(sprintf("[%s | %s] %d barcodes -> %s", cond, ct, length(barcodes), fname))
    }
  }
  invisible(exported)
}
