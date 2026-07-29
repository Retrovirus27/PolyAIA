library(dplyr)

# Suppress specific, known-benign warnings (e.g. GenomicRanges Seqinfo merge
# notices, BiocParallel multicore serialization notices) while still letting
# any other, unexpected warning through. Used instead of blanket
# suppressWarnings() so real problems aren't silently hidden.
# @noRd
.suppress_known_warnings <- function(expr, patterns) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (any(vapply(patterns, function(p) grepl(p, conditionMessage(w)), logical(1)))) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

#' Load 10X data and annotate with metadata
#'
#' @param data_dir Directory containing .h5 files
#' @param file_pattern Pattern to identify .h5 files
#' @param metadata_file Path to CSV metadata file
#' @param species Choose between species \code{"human"} or \code{"mouse"}
#' @param reference According to the species choose the reference.
#'   Options: \code{"humancortexref"} or \code{"mousecortexref"}
#' @param column_id Column ID for samples to match metadata. Ensure no spaces in values
#' @param genome_build Genome build to use for annotation.
#'   Options: \code{"hg38"}, \code{"hg19"}, \code{"mm10"}, or \code{"mm39"}
#' @param IDtype Gene identifier type for annotation.
#'   Options: \code{"SYMBOL"}, \code{"ENTREZID"}, or \code{"ENSEMBL"}
#' @param columns Columns to include in gene annotation.
#'   Options: \code{"SYMBOL"}, \code{"GENETYPE"}, \code{"ENTREZID"}, \code{"ENZYME"}, \code{"ENSEMBL"}
#' @param include_ranges Logical. Whether to include genomic ranges in the annotation. Default \code{TRUE}
#' @param remove_doublets Logical. Whether to remove doublets using scDblFinder. Default \code{TRUE}
#' @param apply_tags Logical. Whether to apply sample tags to barcodes. Default \code{TRUE}
#' @param run_azimuth Logical. Whether to run Azimuth cell type annotation. If it fails,
#'   the pipeline continues without \code{Population}/\code{Subpopulation}/\code{CellType}
#'   columns rather than aborting. Default \code{TRUE}
#' @param ncores Integer. Number of cores for parallel processing during doublet detection. Default \code{4}
#' @param verbose Logical. Whether to print progress messages. Default \code{TRUE}
#'
#' @return A \code{Seurat} object with annotated metadata
#'
#' @export
UploadSce <- function(data_dir = "Data/SCE_1/",
                      file_pattern = "_feature_bc_matrix.h5",
                      metadata_file = "Data/SCE_1.csv",
                      column_id = "Sample_accession",
                      species = c("human", "mouse"),
                      reference = c("humancortexref", "mousecortexref"),
                      genome_build = c("hg38", "hg19", "mm10", "mm39"),
                      IDtype = c("SYMBOL", "ENTREZID", "ENSEMBL"),
                      columns = c("SYMBOL", "GENETYPE", "ENTREZID", "ENZYME", "ENSEMBL"),
                      include_ranges = TRUE,
                      remove_doublets = TRUE,
                      apply_tags = TRUE,
                      run_azimuth = TRUE,
                      ncores = 4,
                      verbose = TRUE) {

  # Resolve multiple-choice arguments -- do this first so mistakes are caught
  # immediately instead of after the expensive doublet-removal step.
  species      <- match.arg(species)
  reference    <- match.arg(reference)
  genome_build <- match.arg(genome_build)
  IDtype       <- match.arg(IDtype)

  if (species == "human" && reference != "humancortexref") {
    stop("species = 'human' requires reference = 'humancortexref' (got '", reference, "').")
  }
  if (species == "mouse" && reference != "mousecortexref") {
    stop("species = 'mouse' requires reference = 'mousecortexref' (got '", reference, "').")
  }

  # Check Azimuth
  if (!requireNamespace("Azimuth", quietly = TRUE)) {
    stop("Package 'Azimuth' is required. Install with: remotes::install_github('satijalab/azimuth')")
  }

  # Check Signac
  if (!requireNamespace("Signac", quietly = TRUE)) {
    stop("Package 'Signac' is required. Install with: BiocManager::install('Signac')")
  }

  # Validate inputs early with clear, actionable messages
  if (!dir.exists(data_dir)) {
    stop("data_dir does not exist: ", data_dir)
  }
  if (!file.exists(metadata_file)) {
    stop("metadata_file does not exist: ", metadata_file)
  }

  # Load 10X data -- always runs
  file_names <- dir(path = data_dir, pattern = file_pattern)
  if (length(file_names) == 0) {
    stop("No files matching pattern '", file_pattern, "' found in ", data_dir)
  }
  all_files <- paste0(data_dir, "/", file_names)
  names(all_files) <- gsub(file_pattern, "", file_names)
  sce <- DropletUtils::read10xCounts(all_files, sample.names = names(all_files), col.names = TRUE)

  message("Loaded ", ncol(sce), " cells across ", length(all_files), " samples")

  # Load and annotate metadata
  meta_data <- utils::read.csv(metadata_file)
  if (!column_id %in% colnames(meta_data)) {
    stop("column_id '", column_id, "' not found in metadata_file columns: ",
         paste(colnames(meta_data), collapse = ", "))
  }
  rownames(meta_data) <- meta_data[[column_id]]

  # Add annotations to sce
  for (x in colnames(meta_data)) {SummarizedExperiment::colData(sce)[[x]] <- meta_data[[x]][match(SummarizedExperiment::colData(sce)$Sample, meta_data[[column_id]])] }

  message("Annotated with ", ncol(meta_data), " metadata columns")

  # Apply tags if requested
  if (apply_tags == TRUE) {
    message("Applying sample tags to barcodes...")

    # DropletUtils::read10xCounts() (with col.names = TRUE, multiple samples)
    # names columns "<sample index>_<barcode>", where <sample index> is a plain
    # 1-based integer -- NOT the sample.names we supplied. So the map here goes
    # index (as a string, e.g. "1") -> real sample name, to translate that
    # numeric prefix into something readable.
    sample_levels <- levels(as.factor(SummarizedExperiment::colData(sce)$Sample))
    gsm.map <- stats::setNames(sample_levels, seq_along(sample_levels))

    # Apply tags to column names
    new_names <- sapply(
      colnames(sce),
      function(x) {
        parts <- strsplit(x, "_")[[1]]
        prefix <- gsm.map[parts[1]]
        paste(prefix, parts[2], sep = "_")
      }
    )

    if (any(is.na(new_names))) {
      stop("Failed to apply sample tags for ", sum(is.na(new_names)), " barcode(s) -- ",
           "the column-name format from read10xCounts() didn't match the expected ",
           "'<sample index>_<barcode>' pattern. Check that colnames(sce) look like ",
           "'1_AAACCTGAGAAACCAT-1' before tagging.")
    }
    colnames(sce) <- new_names

    message("Tags applied successfully")
  }

  if (remove_doublets == TRUE) {
    message("Cells before doublet removal: ", ncol(sce))
    message("Removing doublets...")
    if (!requireNamespace("scDblFinder", quietly = TRUE)) {
      stop("Package 'scDblFinder' is required. Install with: BiocManager::install('scDblFinder')")
    }

    if (!requireNamespace("BiocParallel", quietly = TRUE)) {
      stop("Package 'BiocParallel' is required. Install with: BiocManager::install('BiocParallel')")
    }

    # Run doublet detection
    # (MulticoreParam forking on some platforms emits a benign
    # "'package:stats' may not be available when loading" notice during
    # result serialization -- suppressed here, everything else still shows.)
    sce <- .suppress_known_warnings(
      scDblFinder::scDblFinder(sce, samples="Sample", BPPARAM = BiocParallel::MulticoreParam(workers = ncores)),
      patterns = "may not be available when loading"
    )
    message("Doublets removed: ", sum(sce$scDblFinder.class == "doublet"), " doublets identified.")

    # Filter out doublets
    sce <- sce[, sce$scDblFinder.class == "singlet"]
    message("Remaining cells after doublet removal: ", ncol(sce))
  }

  # Convert to Seurat object
  message("Converting to Seurat object...")

  # Ensure counts are in dgCMatrix format
  if (!inherits(BiocGenerics::counts(sce), "dgCMatrix")) {
    SummarizedExperiment::assay(sce, "counts") <- methods::as(SummarizedExperiment::assay(sce, "counts"), "dgCMatrix")
  }

  # Convert SCE to Seurat
  seurat_obj <- Seurat::as.Seurat(sce, counts = "counts", data = NULL)

  # Rename the original assay to "RNA"
  seurat_obj <- SeuratObject::RenameAssays(object = seurat_obj, originalexp = "RNA")

  # Convert the RNA assay to the Assay5 class
  seurat_obj[["RNA"]] <- methods::as(object = seurat_obj[["RNA"]], Class = "Assay5")

  # Preserve the row information
  seurat_obj@meta.data$row_names <- rownames(seurat_obj@meta.data)

  message("Seurat object created successfully")

  # Calculate basic QC metrics (nCount_RNA and nFeature_RNA)
  message("Calculating QC metrics...")
  seurat_obj[["nCount_RNA"]] <- Matrix::colSums(seurat_obj[["RNA"]]$counts)
  seurat_obj[["nFeature_RNA"]] <- Matrix::colSums(seurat_obj[["RNA"]]$counts > 0)

  # Calculate Mito and Ribo metrics
  if (!requireNamespace("scCustomize", quietly = TRUE)) {
    stop("Package 'scCustomize' is required. Install with: install.packages('scCustomize')")
  }

  if (species == "human") {
    tryCatch({
      seurat_obj <- scCustomize::Add_Cell_QC_Metrics(object = seurat_obj, species = "human", ensembl_ids = TRUE, add_cell_cycle = FALSE)
      message("Mito and ribo metrics calculated for human and added to metadata")
    }, error = function(e) {
      warning("Add_Cell_QC_Metrics failed and was skipped -- QC metrics will be absent from metadata.\n  Reason: ", conditionMessage(e))
    })
  }

  if (species == "mouse") {
    tryCatch({
      seurat_obj <- scCustomize::Add_Cell_QC_Metrics(object = seurat_obj, species = "mouse", ensembl_ids = TRUE, add_cell_cycle = FALSE)
      message("Mito and ribo metrics calculated for mouse and added to metadata")
    }, error = function(e) {
      warning("Add_Cell_QC_Metrics failed and was skipped -- QC metrics will be absent from metadata.\n  Reason: ", conditionMessage(e))
    })
  }

  anno <- AnnotateGenes(rownames(seurat_obj), IDtype = IDtype, species = species, genome_build = genome_build, columns = columns, include_ranges = include_ranges, verbose = verbose)
  colnames(anno)[1]<- "ID"

  # Drop the identifier columns we don't need on the Seurat side, but only if
  # they're actually present -- keeps this safe when `columns` is customized.
  anno_clean <- anno %>%
    dplyr::select(-tidyselect::any_of(c("ENTREZID", "ENZYME", "GENETYPE", "SYMBOL")))

  seurat_obj@assays[["RNA"]]@meta.data <- seurat_obj@assays[["RNA"]]@meta.data %>%
    dplyr::left_join(anno_clean, by = "ID", suffix = c("", ".anno")) %>%
    dplyr::select(-tidyselect::ends_with(".anno"))  # drop any remaining duplicates from anno side


  # Run Azimuth cell type annotation
  if (run_azimuth) {
    message("Running Azimuth annotation with ", reference, "...")

    # Azimuth internally calls some SeuratObject/Seurat generics (e.g. `Key<-`)
    # unqualified, assuming they've been attached via library() rather than
    # accessed via `::`. Since this package only Imports them, attach them to
    # the search path here (if not already) so RunAzimuth can find them.
    for (pkg in c("SeuratObject", "Seurat")) {
      if (!paste0("package:", pkg) %in% search() && requireNamespace(pkg, quietly = TRUE)) {
        suppressPackageStartupMessages(attachNamespace(pkg))
      }
    }

    seurat_obj <- tryCatch({
      seurat_obj_dummy <- Azimuth::RunAzimuth(seurat_obj, reference = reference)

      # Rename Azimuth predictions to standardized names
      seurat_obj$Population <- seurat_obj_dummy$predicted.class[match(seurat_obj_dummy$row_names, seurat_obj$row_names)]
      seurat_obj$Population_Score <- seurat_obj_dummy$predicted.class.score[match(seurat_obj_dummy$row_names, seurat_obj$row_names)]
      seurat_obj$Subpopulation <- seurat_obj_dummy$predicted.subclass[match(seurat_obj_dummy$row_names, seurat_obj$row_names)]
      seurat_obj$Subpopulation_Score <- seurat_obj_dummy$predicted.subclass.score[match(seurat_obj_dummy$row_names, seurat_obj$row_names)]

      message("Azimuth annotation completed")

      # Create broad cell type categories
      message("Creating CellType categories...")
      seurat_obj$CellType <- dplyr::case_when(
        seurat_obj$Subpopulation %in% c("L2/3 IT", "L5 IT", "L5 ET", "L5/6 NP", "L6 IT", "L6 IT Car3", "L6 CT", "L6b") ~ "Exc",
        seurat_obj$Subpopulation %in% c("Meis2","Lamp5", "Vip", "Sst", "Sst Chodl", "Pvalb", "Sncg") ~ "Inh",
        seurat_obj$Subpopulation %in% c("Astro") ~ "Astro",
        seurat_obj$Subpopulation %in% c("Oligo") ~ "Oligo",
        seurat_obj$Subpopulation %in% c("OPC") ~ "OPC",
        seurat_obj$Subpopulation %in% c("Micro-PVM") ~ "Mic",
        seurat_obj$Subpopulation %in% c("VLMC","Endo","Peri") ~ "End"
      )

      seurat_obj
    }, error = function(e) {
      warning("Azimuth annotation failed and was skipped -- Population/Subpopulation/CellType ",
              "will be absent from metadata.\n  Reason: ", conditionMessage(e))
      seurat_obj
    })
  } else {
    message("Skipping Azimuth annotation (run_azimuth = FALSE).")
  }

  message("Pipeline completed successfully!")
  gc()
  return(seurat_obj)
}

#' Annotate gene features with genomic information
#'
#' @param data Gene names/IDs to annotate (character vector or \code{rownames}
#'   from a Seurat object)
#' @param IDtype Type of gene identifier.
#'   Options: \code{"SYMBOL"}, \code{"ENTREZID"}, or \code{"ENSEMBL"}
#' @param species Species for annotation. Options: \code{"human"} or \code{"mouse"}
#' @param genome_build Genome build for range annotation.
#'   Options: \code{"hg38"}, \code{"hg19"}, \code{"mm10"}, or \code{"mm39"}
#' @param columns Columns to retrieve from the annotation database.
#'   Options: \code{"SYMBOL"}, \code{"GENETYPE"}, \code{"ENTREZID"},
#'   \code{"ENZYME"}, \code{"ENSEMBL"}.
#'   Default: \code{c("SYMBOL", "GENETYPE", "ENTREZID", "ENZYME", "ENSEMBL")}
#' @param include_ranges Logical. Whether to include genomic range information.
#'   Default \code{TRUE}
#' @param verbose Logical. Whether to print progress messages. Default \code{TRUE}
#'
#' @return A \code{data.frame} with gene annotations. If \code{include_ranges = TRUE},
#'   genomic coordinates are appended from the corresponding \code{TxDb} package
#'   (\code{genome_build} must be specified)
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Annotate human genes by SYMBOL
#' annot <- AnnotateGenes(rownames(seurat_obj),
#'                        IDtype = "SYMBOL",
#'                        species = "human")
#'
#' # Annotate mouse genes by ENSEMBL without ranges
#' annot <- AnnotateGenes(gene_ids,
#'                        IDtype = "ENSEMBL",
#'                        species = "mouse",
#'                        include_ranges = FALSE)
#' }
AnnotateGenes <- function(
    data,
    IDtype     = c("SYMBOL", "ENTREZID", "ENSEMBL"),
    species    = c("human", "mouse"),
    genome_build = c("hg38", "hg19", "mm10", "mm39"),
    columns    = c("SYMBOL", "GENETYPE", "ENTREZID", "ENZYME", "ENSEMBL"),
    include_ranges = TRUE,
    verbose    = TRUE
) {

  # Match arguments
  IDtype <- match.arg(IDtype)
  species <- match.arg(species)
  genome_build <- match.arg(genome_build)

  # Validate genome_build vs species
  human_builds <- c("hg38", "hg19")
  mouse_builds <- c("mm10", "mm39")

  if (species == "human" && !genome_build %in% human_builds) {
    stop("For species = 'human', genome_build must be one of: ", paste(human_builds, collapse = ", "))
  }
  if (species == "mouse" && !genome_build %in% mouse_builds) {
    stop("For species = 'mouse', genome_build must be one of: ", paste(mouse_builds, collapse = ", "))
  }

  # Load appropriate species library
  if (verbose) message("Loading ", species, " annotation library (", genome_build, ")...")

  if (species == "human") {
    if (!requireNamespace("Homo.sapiens", quietly = TRUE)) {
      stop("Package 'Homo.sapiens' is required. Install with: BiocManager::install('Homo.sapiens')")
    }
    org_db <- Homo.sapiens::Homo.sapiens
  }

  if (species == "mouse") {
    if (!requireNamespace("Mus.musculus", quietly = TRUE)) {
      stop("Package 'Mus.musculus' is required. Install with: BiocManager::install('Mus.musculus')")
    }
    org_db <- Mus.musculus::Mus.musculus
  }

  # Only load a TxDb/EnsDb if genomic ranges were actually requested -- this
  # also means genome builds without ranges support (e.g. mm39) still work
  # fine as long as include_ranges = FALSE.
  txdb <- NULL
  if (include_ranges) {
    if (species == "human") {
      txdb <- switch(genome_build,
                     "hg38" = {
                       if (!requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE)) {
                         stop("Package 'EnsDb.Hsapiens.v86' is required for hg38. Install with: BiocManager::install('EnsDb.Hsapiens.v86')")
                       }
                       EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
                     },
                     "hg19" = {
                       if (!requireNamespace("EnsDb.Hsapiens.v75", quietly = TRUE)) {
                         stop("Package 'EnsDb.Hsapiens.v75' is required for hg19. Install with: BiocManager::install('EnsDb.Hsapiens.v75')")
                       }
                       EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75
                     }
      )
    }

    if (species == "mouse") {
      txdb <- switch(genome_build,
                     "mm10" = {
                       if (!requireNamespace("EnsDb.Mmusculus.v79", quietly = TRUE)) {
                         stop("Package 'EnsDb.Mmusculus.v79' is required for mm10. Install with: BiocManager::install('EnsDb.Mmusculus.v79')")
                       }
                       EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79
                     },
                     "mm39" = {
                       stop("genome_build = 'mm39' does not yet have genomic ranges support. ",
                            "Use genome_build = 'mm10' for ranges, or set include_ranges = FALSE ",
                            "to annotate mm39 genes without ranges.")
                     }
      )
    }
  }

  # Validate columns
  available_cols <- AnnotationDbi::columns(org_db)
  invalid_cols <- dplyr::setdiff(columns, available_cols)

  if (length(invalid_cols) > 0) {
    warning("Invalid columns: ", paste(invalid_cols, collapse = ", "), "\nAvailable columns: ", paste(available_cols, collapse = ", "))
    columns <- intersect(columns, available_cols)
  }

  # Ensure IDtype is in columns
  if (!IDtype %in% columns) {
    columns <- c(IDtype, columns)
  }

  if (verbose) message("Retrieving annotation data for ", length(data), " genes...")

  # Get annotation data
  annotation <- AnnotationDbi::select(org_db, keys = data, columns = columns, keytype = IDtype)

  # Remove duplicates based on keytype
  if (verbose) message("Removing duplicates based on ", IDtype, "...")

  id_column <- annotation[[IDtype]]
  annotation <- annotation[!duplicated(id_column), ]

  if (verbose) message("Retained ", nrow(annotation), " unique genes")

  # Add gene length if requested
  if (include_ranges) {
    if (verbose) message("Calculating ranges using ", genome_build, " build...")

    # Extract genes and transcripts from EnsDb
    # (Merging per-chromosome Seqinfo internally triggers a benign "no sequence
    # levels in common" notice -- suppressed here, everything else still shows.)
    genes_txdb <- .suppress_known_warnings(
      Signac::GetGRangesFromEnsDb(txdb),
      patterns = "no sequence levels in common"
    )
    genes_df <- data.frame(genes_txdb)

    # Aggregate transcripts by gene
    if (verbose) message("Aggregating transcripts by gene...")

    genes_aggregated <- genes_df %>%
      dplyr::group_by(gene_id) %>%
      dplyr::summarise(
        chromosome = dplyr::first(seqnames),
        gene_start = min(start),
        gene_end = max(end),
        gene_width = gene_end - gene_start + 1,
        strand = dplyr::first(strand),
        gene_name = dplyr::first(gene_name),
        gene_biotype = dplyr::first(gene_biotype),
        n_transcripts = dplyr::n_distinct(tx_id),
        transcript_info = list(data.frame(
          tx_id = tx_id,
          tx_chromosome = seqnames,
          tx_start = start,
          tx_end = end,
          tx_strand = strand,
          tx_width = width,
          tx_type = type
        )),
        .groups = "drop"
      ) %>%
      dplyr::rename(ENSEMBL = gene_id)

    # Join with annotation
    annotation_final <- dplyr::left_join(annotation, genes_aggregated, by = "ENSEMBL")

    if (verbose) message("Ranges information added with ", sum(!is.na(annotation_final$n_transcripts)), " genes having transcript information")
  } else {
    annotation_final <- annotation
  }

  if (verbose) {
    message("\n=== Annotation Summary ===")
    message("Total genes annotated: ", nrow(annotation_final))
    message("Genome build: ", genome_build)
    message("Columns: ", paste(names(annotation_final), collapse = ", "))
    if (include_ranges && "gene_width" %in% names(annotation_final)) {
      message("Gene length range: ", min(annotation_final$gene_width, na.rm = TRUE), " - ", max(annotation_final$gene_width, na.rm = TRUE), " bp")
    }
  }

  return(annotation_final)
}


#' Annotation of RNA data
#'
#' @param data A Seurat object.
#' @param idtype ID type used for gene annotation. Default is \code{"ENSEMBL"}.
#' @param Txdb Optional \code{TxDb}/\code{EnsDb} object used to add genomic
#'   range information. If \code{NULL} (default), ranges are skipped.
#'
#' @return A Seurat object with annotated RNA assay.
#' @noRd
#'
Annot <- function(data, idtype, Txdb = NULL) {

  annotation <- AnnotationDbi::select(
    Homo.sapiens::Homo.sapiens,
    keys    = row.names(data),
    columns = c("SYMBOL", "GENETYPE", "ENTREZID", "ENZYME", "ENSEMBL"),
    keytype = idtype
  )

  # Deduplicate on the key column
  id_column  <- annotation[[idtype]]
  annotation <- annotation[!duplicated(id_column), ]

  # Safe row-matched assignment — handles missing/unmatched ENSEMBL IDs
  features <- rownames(data)
  rownames(annotation) <- annotation[[idtype]]

  data@assays$RNA@meta.data$SYMBOL   <- annotation[features, "SYMBOL"]
  data@assays$RNA@meta.data$GENETYPE <- annotation[features, "GENETYPE"]
  data@assays$RNA@meta.data$ENTREZID <- annotation[features, "ENTREZID"]
  data@assays$RNA@meta.data$ENSEMBL  <- annotation[features, "ENSEMBL"]

  if (!is.null(Txdb)) {
    genes_df <- data.frame(Txdb)

    genes_aggregated <- genes_df %>%
      dplyr::group_by(gene_id) %>%
      dplyr::summarise(
        chromosome   = dplyr::first(seqnames),
        gene_start   = min(start),
        gene_end     = max(end),
        gene_width   = gene_end - gene_start + 1,
        strand       = dplyr::first(strand),
        gene_name    = dplyr::first(gene_name),
        gene_biotype = dplyr::first(gene_biotype),
        n_transcripts = dplyr::n_distinct(tx_id),
        transcript_info = list(data.frame(
          tx_id        = tx_id,
          tx_chromosome = seqnames,
          tx_start     = start,
          tx_end       = end,
          tx_strand    = strand,
          tx_width     = width,
          tx_type      = type
        )),
        .groups = "drop"
      ) %>%
      dplyr::rename(ENSEMBL = gene_id)

    meta <- data@assays$RNA@meta.data
    meta <- dplyr::left_join(meta, genes_aggregated, by = "ENSEMBL")
    rownames(meta) <- rownames(data@assays$RNA@meta.data)
    data@assays$RNA@meta.data <- meta
  }

  return(data)
}

#' Create an annotated Seurat object
#'
#' Loads raw count data, removes doublets, and performs cell type annotation
#' using one or two reference datasets.
#'
#' @param data A raw count matrix or path to an \code{.h5} file.
#' @param DataName Character string used to name the resulting Seurat object.
#' @param Condition Optional condition/group label attached to the sample.
#'   Default \code{NULL}.
#' @param RefAnnotation1 First reference dataset for SingleR annotation. Default is \code{NULL}.
#' @param RefAnnotation2 Second reference dataset for SingleR annotation. Default is \code{NULL}.
#' @param Type Input data type. One of \code{"h5"}, \code{"Ambient.h5"}, or \code{"counts"}. Default is \code{"h5"}.
#' @param RDoublets Logical. Whether to remove doublets using \code{scDblFinder}. Default is \code{TRUE}.
#' @param ApplyTags Logical. Whether to prefix cell barcodes with the sample
#'   tag (\code{DataName}). Default is \code{TRUE}.
#' @param IDType Gene identifier type used for annotation (e.g. \code{"ENSEMBL"}).
#'   Default \code{"ENSEMBL"}.
#' @param ncores Integer. Number of cores for parallel processing during doublet detection. Default \code{4}
#'
#' @return An annotated Seurat object.
#'
#' @examples
#' # UploadData(data, DataName = "Sample1", Type = "h5")
#' # UploadData(data, DataName = "Sample1", RefAnnotation1 = ImmGenData(), RDoublets = TRUE)
UploadData <- function(data, DataName, Condition = NULL, RefAnnotation1 = NULL, RefAnnotation2 = NULL,
                       Type = c("h5", "Ambient.h5","counts"), RDoublets = TRUE, ApplyTags = TRUE,
                       IDType = "ENSEMBL", ncores = 4){

  Type <- match.arg(Type)

  if (Type == "h5") {
    # Read the count data
    count_data <- Seurat::Read10X_h5(data, use.names = FALSE, unique.features = TRUE)
    # single cell experiment
    sce <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=count_data))
  }

  if (Type == "Ambient.h5") {
    # Extract the specific files for this sample from the data vector
    filtered_path <- data[grepl(paste0(DataName, "_filtered_feature_bc_matrix\\.h5"), data)]
    raw_path      <- data[grepl(paste0(DataName, "_raw_feature_bc_matrix\\.h5"),      data)]
    cluster_path  <- data[grepl(paste0(DataName, "_clusters\\.csv"),                  data)]

    if (!requireNamespace("SoupX", quietly = TRUE)) {
      stop("Package 'SoupX' is required for Type = 'Ambient.h5'. ",
           "Install with: install.packages('SoupX')")
    }
    if (length(raw_path) == 0)     stop("Raw h5 file not found for: ",     DataName)
    if (length(cluster_path) == 0) stop("Cluster CSV file not found for: ", DataName)

    # Read raw (unfiltered) droplets
    raw_data <- Seurat::Read10X_h5(raw_path, use.names = TRUE, unique.features = TRUE)
    count_data <- Seurat::Read10X_h5(filtered_path, use.names = TRUE, unique.features = TRUE)

    # Build SoupChannel (tod = table of droplets, toc = table of cells)
    sc <- SoupX::SoupChannel(tod = raw_data, toc = count_data)

    # Load pre-computed cluster labels
    cluster_df <- utils::read.csv(cluster_path, row.names = 1)
    clusters   <- stats::setNames(as.character(cluster_df[[1]]), rownames(cluster_df))
    sc         <- SoupX::setClusters(sc, clusters)

    # Estimate and correct contamination
    sc        <- SoupX::autoEstCont(sc)
    corrected <- SoupX::adjustCounts(sc, roundToInt = TRUE)

    message(sprintf("SoupX: estimated contamination fraction = %.1f%%", sc$fit$rhoEst * 100))

    # Rebuild SCE with corrected counts
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = corrected))
  }

  if (Type == "counts") {
    # Count Read
    sce <- DropletUtils::read10xCounts(data)
  }

  if (!is.null(Condition)) {
    SummarizedExperiment::colData(sce)$Condition <- Condition
  }

  # Log normalization
  sce <- scuttle::logNormCounts(sce)

  if (!is.null(RefAnnotation1)) {
    # Annotation
    Annot1 <- SingleR::SingleR(test = sce, ref = RefAnnotation1, labels = RefAnnotation1$label.main, assay.type.test = "logcounts")
    sce$RefAnnotation1 <- Annot1$labels # Aggregate labels
  }

  if (!is.null(RefAnnotation2)) {
    # Annotation
    Annot2 <- SingleR::SingleR(test = sce, ref = RefAnnotation2, labels = RefAnnotation2$label.main, assay.type.test = "logcounts")
    sce$RefAnnotation2  <- Annot2$labels # Aggregate labels
  }

  if (RDoublets) {
    suppressWarnings( sce <- scDblFinder::scDblFinder(sce, BPPARAM=BiocParallel::MulticoreParam(workers = ncores)) )
    n_doublets <- sum(sce$scDblFinder.class == "doublet")
    n_singlets <- sum(sce$scDblFinder.class == "singlet")
    sce <- sce[, sce$scDblFinder.class == "singlet"]
    message(sprintf("scDblFinder: removed %d doublets, %d singlets retained", n_doublets, n_singlets))
  }

  if (ApplyTags == TRUE) {
    message("Applying sample tag to barcodes...")
    colnames(sce) <- paste(DataName, colnames(sce), sep = "_")
    message("Tags applied successfully")
  }

  counts_matrix <- sce@assays@data@listData[["counts"]]
  colnames(counts_matrix) <- colnames(sce)

  # Create a Seurat object
  seurat_obj <- SeuratObject::CreateSeuratObject(counts = counts_matrix, project = DataName)
  seurat_obj <- SeuratObject::AddMetaData(seurat_obj, metadata = as.data.frame(SummarizedExperiment::colData(sce)))


  # Return object
  return(seurat_obj)
}

#' Normalize, integrate, and compute a UMAP embedding for a Seurat object
#'
#' Runs the standard Seurat v5 multi-sample workflow used in this pipeline's
#' example scripts (see the "Normalization and integration" section of
#' \code{G1_Universal_SC_analysis.Rmd}): split the requested assay's layers by
#' sample, log-normalize, find variable features, scale, run PCA, integrate
#' the per-sample layers (Harmony by default), merge the layers back, and
#' compute a 2D UMAP from the integrated space -- all in one call.
#'
#' @param seurat_obj A Seurat object, normally already QC-filtered.
#' @param sample_col Metadata column identifying each cell's sample, used to
#'   split assay layers before integration -- required by
#'   \code{IntegrateLayers()} for per-sample methods like Harmony (without it,
#'   \code{IntegrateLayers(method = HarmonyIntegration, ...)} fails with
#'   \code{Error in names(groups) <- "group" : attempt to set an attribute on
#'   NULL}, since Harmony builds its groups from the per-sample layer names).
#'   Default \code{"Sample"}.
#' @param assay Name of the assay to split/normalize/integrate. Default \code{"RNA"}.
#' @param nfeatures Number of variable features to select. Default \code{2000}.
#' @param method Integration method passed to \code{Seurat::IntegrateLayers()}.
#'   Default \code{Seurat::HarmonyIntegration}.
#' @param new.reduction Name of the integrated reduction to create. Default
#'   \code{"harmony"}.
#' @param reduction.name Name of the 2D UMAP reduction to create. Default
#'   \code{"umap_harmony"} (this pipeline's convention, expected by
#'   \code{UmapPlot()}).
#' @param reduction.key Key prefix for the UMAP reduction's columns. Default
#'   \code{"Harmony_"}.
#' @param custom_pca Number of PCA/integration dimensions to use for
#'   \code{IntegrateLayers()} and \code{RunUMAP()}. If \code{NULL} (default),
#'   this is chosen automatically as the elbow ("codo") of the PCA
#'   standard-deviation curve -- the PC lying farthest from the straight line
#'   joining the first and last PC's stdev (the point of maximum curvature on
#'   the ElbowPlot). Supply a number (1-30) to override the automatic choice.
#'   PCA always computes 30 PCs, so \code{custom_pca} cannot exceed 30.
#' @param verbose Logical. Whether to print progress messages, including which
#'   \code{n_dims} was selected. Default \code{TRUE}.
#'
#' @return The input \code{seurat_obj}, with its assay layers split and
#'   re-joined, and \code{"pca"}, \code{new.reduction}, and
#'   \code{reduction.name} reductions computed. The \code{n_dims} used is
#'   stored at \code{seurat_obj@misc$n_dims} for later reference.
#'
#' @examples
#' \dontrun{
#' seurat_obj <- SeuratPipeline(seurat_obj)
#' seurat_obj <- SeuratPipeline(seurat_obj, custom_pca = 10)
#' }
#'
#' @export
SeuratPipeline <- function(
    seurat_obj,
    sample_col = "Sample",
    assay = "RNA",
    nfeatures = 2000,
    method = Seurat::HarmonyIntegration,
    new.reduction = "harmony",
    reduction.name = "umap_harmony",
    reduction.key = "Harmony_",
    custom_pca = NULL,
    verbose = TRUE
) {
  if (!sample_col %in% colnames(seurat_obj@meta.data)) {
    stop("sample_col '", sample_col, "' not found in seurat_obj metadata.")
  }
  if (!assay %in% SeuratObject::Assays(seurat_obj)) {
    stop("assay '", assay, "' not found in seurat_obj. Available assays: ",
         paste(SeuratObject::Assays(seurat_obj), collapse = ", "))
  }

  # PCA always computes 30 PCs; the elbow selection considers this same range.
  max_dims <- 30L

  # Validate the numeric arguments up front. A string like custom_pca = "10"
  # would otherwise reach RunUMAP()/IntegrateLayers()'s `dims` (or, historically,
  # RunPCA()'s `npcs`) and fail deep inside irlba with the opaque "non-numeric
  # argument to binary operator" (from `nv + 7`). Reject non-numeric input here
  # -- strings are NOT silently coerced.
  if (!is.null(custom_pca)) {
    if (!is.numeric(custom_pca) || length(custom_pca) != 1 || is.na(custom_pca) || custom_pca < 1) {
      stop("`custom_pca` must be a single positive number (numeric, unquoted -- ",
           "e.g. 10, not \"10\"), or NULL to auto-select by the elbow method.")
    }
    if (custom_pca > max_dims) {
      stop("`custom_pca` (", custom_pca, ") cannot exceed ", max_dims,
           " (PCA computes ", max_dims, " PCs).")
    }
  }

  # Seurat v5 needs the assay split into one layer per sample *before*
  # IntegrateLayers() -- see @param sample_col above for why.
  seurat_obj[[assay]] <- split(seurat_obj[[assay]], f = seurat_obj[[sample_col, drop = TRUE]])

  # Internal Seurat steps are always run with verbose = FALSE so their
  # per-layer "Performing log-normalization" progress bars don't flood the
  # console; this function's own high-level messages (gated on `verbose`) are
  # kept as the progress signal instead.
  if (verbose) message("Normalizing, finding variable features, scaling, and running PCA...")
  seurat_obj <- Seurat::NormalizeData(seurat_obj, assay = assay, verbose = FALSE)
  seurat_obj <- Seurat::FindVariableFeatures(
    seurat_obj, selection.method = "vst", nfeatures = nfeatures, assay = assay, verbose = FALSE
  )
  seurat_obj <- Seurat::ScaleData(
    seurat_obj, features = Seurat::VariableFeatures(seurat_obj), assay = assay, verbose = FALSE
  )
  seurat_obj <- Seurat::RunPCA(
    seurat_obj, features = Seurat::VariableFeatures(seurat_obj), npcs = max_dims, assay = assay, verbose = FALSE
  )

  # Pick the number of dimensions to carry into integration/UMAP. Default is the
  # geometric elbow of the PCA stdev curve (.suggest_n_dims: the PC farthest
  # from the chord joining the first and last PC -- the actual "codo" of the
  # ElbowPlot) over the 30 PCs; `custom_pca` overrides it with a fixed value.
  if (is.null(custom_pca)) {
    n_dims <- .suggest_n_dims(seurat_obj, reduction.method = "pca", n.pcs = max_dims)
    if (verbose) {
      message("Auto-selected n_dims = ", n_dims, " (elbow / maximum curvature of the PCA stdev curve).")
    }
  } else {
    n_dims <- custom_pca
    if (verbose) message("Using custom_pca = ", n_dims, " dimensions.")
  }

  seurat_obj <- Seurat::IntegrateLayers(
    seurat_obj,
    method         = method,
    orig.reduction = "pca",
    assay          = assay,
    new.reduction  = new.reduction,
    dims           = 1:n_dims,
    verbose        = FALSE
  )

  # Merge the per-sample layers back into one now that integration is done
  seurat_obj[[assay]] <- SeuratObject::JoinLayers(seurat_obj[[assay]])

  seurat_obj <- Seurat::RunUMAP(
    seurat_obj,
    reduction      = new.reduction,
    dims           = 1:n_dims,
    reduction.name = reduction.name,
    reduction.key  = reduction.key,
    verbose        = FALSE
  )

  seurat_obj@misc$n_dims <- n_dims
  seurat_obj
}
