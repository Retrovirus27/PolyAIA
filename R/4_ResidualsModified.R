#' Calculate PolyA Residuals
#'
#'
#' @param object Seurat object containing a polyAsiteAssay
#' @param assay Name of polyAsiteAssay to be used in calculating polyAresiduals
#' @param features Features to include in calculation of polyA residuals.
#' Default is to use all features.
#' @param background Identity of cells to use as background.
#' Default is to use all cells as a background.
#' @param gene.names Column containing the gene where each polyA site is annotated.
#' Default is symbol.
#' @param min.counts.background Features with at least this many counts in the background cells are included in calculation
#' @param min.variance Sets minimum variance. Default is 0.1.
#' @param sample.n Max number of observations to sample in each bin when performing regularization. Default is 1000.
#' @param do.center Return the centered residuals. Default is FALSE.
#' @param do.scale Return the scaled residuals. Default is FALSE.
#' @param residuals.max Clip residuals above this value. Default is NULL (no clipping).
#' @param residuals.min Clip residuals below this value. Default is NULL (no clipping).
#' @param number.bins Number of bins used to perform regularization. Default is 30.
#' @param verbose Print messages.
#'
#'
#' @return Returns a Seurat object with polyAresiduals assay
#'
#' @examples
#' \dontrun{
#' seurat_obj <- CalcPolyAResidualsPolyA(
#'   seurat_obj,
#'   assay = "polyA",
#'   gene.names = "Gene_Symbol",
#'   background = "Control",
#'   min.counts.background = 5
#' )
#' }
#'
#' @export
#' @concept residuals
#'
CalcPolyAResidualsPolyA <- function(object,
                                assay = "polyA",
                                features = NULL,
                                background = NULL,
                                gene.names = "Gene_Symbol",
                                min.counts.background = 5,
                                min.variance = 0.1,
                                sample.n = 1000,
                                do.scale = FALSE,
                                do.center = FALSE,
                                residuals.max = NULL,
                                residuals.min = NULL,
                                number.bins = 30,
                                verbose = TRUE)
{
  if(verbose) {
    message("Calculating background distribution")
  }

  #if features in NULL, then specify all features in polyA assay
  if (is.null(features)) {
    features <- rownames(SeuratObject::LayerData(object, assay=assay, layer="counts"))
  }

  #if background is NULl, then make a dummy variable
  if (is.null(background)) {
    message("Using all cells in order to estimate background distribution")
    object$dummy <- "all"
    Seurat::Idents(object) <- object$dummy
    background.use = "all"
  }  else {
    background.use = background
    if (!(background %in% unique(Seurat::Idents(object)))) {
      stop("background must be one of the Idents of seurat object")
    }
    message(paste0("Using ", background, " as background distribution"))
  }


  #check if symbols are contained in meta features
  if (!(gene.names %in% colnames(object[[assay]]@meta.features))) {
    stop("Gene.names column not found in meta.features, please make sure
         you are specific gene.names correctly")
  }

  if (sum(is.na(object[[assay]]@meta.features[features,gene.names])) > 0) {
    features.no.anno <- features[is.na(object[[assay]]@meta.features[features,gene.names])]
    message(paste0("Removing ", length(features.no.anno), " sites without a gene annotation"))
    features <- setdiff(features, features.no.anno)
  }

  ##############################################################################
  #get pseudobulked fraction of reads from background
  background.dist <- GetBackgroundDist(object = object, features = features,
                                       background = background.use,
                                       gene.names = gene.names,
                                       assay = assay,
                                       min.counts.background = min.counts.background)
  # Remove features without a gene annotation: `gene` above is built as
  # paste0(<gene name>, "_", <strand>), so a missing/blank gene name (NA or
  # "") collapses to the sentinel values "_-" or "_+" depending on strand --
  # this drops exactly those unannotated sites.
  background.dist <- subset(background.dist, gene!="_-")
  background.dist <- subset(background.dist, gene!="_+")
  features.use <- background.dist$peak

  ##############################################################################
  #calculate sum of counts within each gene
  m <- SeuratObject::LayerData(object = object, assay=assay, layer="counts")
  m <- m[background.dist$peak,]
  #m <- m[order(match(rownames(m), background.dist$peak)), ]
  gene.sum <- rowsum(m, group=background.dist$gene)
  genes <- rownames(gene.sum)

  ##############################################################################
  #fit dirichlet multinomial distribution
  if(verbose) {
    message("Running Dirichlet Multinomial Regression")
  }

  ncells = dim(object)[2]
  background.cells <- Seurat::WhichCells(object, idents=background)
  m.background <- as.matrix(m[,background.cells], nrow = nrow(m))

  #fit dirichlet multinomial for each gene
  res <- lapply(genes, DirichletMultinomial, background.dist=background.dist,
                m.background = m.background, gene.sum=gene.sum,  ncells = ncells)

  # ============================================================
  # RETURN OBJECT
  # ============================================================
  res <- list(
    residuals = res,
    counts = m
  )

  AR <- AssembleResiduals(res = res)

  rm(res)

  # Note: min.variance/number.bins/sample.n/verbose are forwarded here so that
  # customizing them on CalcPolyAResidualsPolyA() actually takes effect --
  # previously this call relied entirely on RegDMVarPolyA()'s own defaults,
  # silently ignoring whatever the caller passed in.
  RDM <- RegDMVarPolyA(ec = AR$ec,
                       var = AR$var,
                       m = AR$m,
                       m.background = m.background,
                       background.dist = background.dist,
                       gene.sum = gene.sum,
                       background.cells = background.cells,
                       min.variance = min.variance,
                       number.bins = number.bins,
                       sample.n = sample.n,
                       verbose = verbose)

  # Same note as above: assay/do.center/do.scale/residuals.max/residuals.min/
  # verbose are forwarded from this function's own parameters instead of
  # being hardcoded, so they're actually respected.
  Residuals <- CalcResiduals(object = object,
                             assay = assay,
                             m = AR$m,
                             ec = AR$ec,
                             var.reg = RDM,
                             do.center = do.center,
                             do.scale = do.scale,
                             residuals.max = residuals.max,
                             residuals.min = residuals.min,
                             verbose = verbose)

  # Logged here rather than inside CalcResiduals() -- see the note there --
  # so the command log captures this function's own small scalar/vector
  # arguments instead of CalcResiduals()'s large m/ec/var.reg matrices.
  Residuals <- Seurat::LogSeuratCommand(object = Residuals)

  return(Residuals)
}

#' Get Background Distribution
#'
#' Calculated Pseudobulk Ratios of Each Isoform within a gene for background distribution.
#' @param object Seurat object containing a polyAsiteAssay
#' @param assay Name of polyAsiteAssay to be used in calculating polyAresiduals
#' @param features Features to include in calculation of polyA residuals.
#' If NULL, use all features.
#' @param background Identity of cells to use as background
#' If NULL, uses all cells combined as a background
#' @param gene.names Name of column containing gene annotations
#' @param min.counts.background Features with at least this many counts in the background cells are included in calculation
#'
#' @return Returns a data frame containing all peaks within genes that have multiple polyA sites that meet min.counts.background criteria
#'
#' @concept residuals
#' @noRd
#'
GetBackgroundDist <- function(object, features, background, gene.names, assay,  min.counts.background) {
  # returns the pseudobulked background distribution for peaks specified
  # must contain gene information in meta data

  suppressMessages(nt.pseudo <- Seurat::AverageExpression(object, features = features, assays = assay, slot="counts"))
  nt.pseudo <- data.frame(background = nt.pseudo[[1]][,background]) #subset just the background
  nt.pseudo$background <- nt.pseudo$background * sum(Seurat::Idents(object)==background)
  nt.pseudo$gene <- paste0(object[[assay]]@meta.features[features, gene.names], "_", object[[assay]]@meta.features[features, "strand"])
  nt.pseudo$peak <- rownames(nt.pseudo)

  nt.pseudo <- nt.pseudo[nt.pseudo$background>min.counts.background,] #subset to peaks with min number of counts
  genes.use <- nt.pseudo$gene[duplicated(nt.pseudo$gene)] #only use genes with at least 2 peaks per gene
  nt.pseudo <- nt.pseudo[nt.pseudo$gene %in% genes.use,]

  if ( length(genes.use)  ==  0) {
    stop("Found no genes with more than 2 features within a gene. Please make sure you are including
         all peaks within a gene you would like to include.")
  }

  tmp <- stats::aggregate(nt.pseudo$background, list(nt.pseudo$gene), FUN=sum)
  colnames(tmp) <- c("gene", "sum")
  nt.pseudo <- merge(nt.pseudo, tmp, by="gene")
  nt.pseudo$frac <- nt.pseudo$background/ nt.pseudo$sum
  return(nt.pseudo)
}

#' Run Dirichlet Multinomial Distribution
#'
#' Fit dirichlet multinomial distribution on each peak within a gene using background cells.
#' Then calculate expected value and variance for each cell based on estimates from dirichlet multinomial regression.
#'
#' @param gene.test which gene to use
#' @param background.dist dataframe containing the isoform ratios for each
#' @param m.background matrix of background distribution
#' @param gene.sum sum of count within each gene for each cell
#' @param ncells number of cells
#'
#' @return Returns a list where first element is matrix of expected values for each peak within the genes,
#' second value is matrix of variance for each peak within the gene
#'
#' @importFrom MGLM MGLMfit
#' @concept residuals
#' @noRd
#'

DirichletMultinomial <- function(
    gene.test,
    background.dist,
    m.background,
    gene.sum,
    ncells
) {
  peaks <- background.dist$peak[background.dist$gene == gene.test]
  t <- m.background[rownames(m.background) %in% peaks,]
  t <- t(t)
  fit <- try(compareFit <- suppressWarnings(MGLM::MGLMfit(t, dist="DM")), silent=TRUE)

  if (!inherits(fit, "try-error")) {
    param <- compareFit@estimate
    sum.p <- sum(param)
    n <-  as.numeric(gene.sum[gene.test,])
    expect.tmp <- data.frame(matrix(nrow=ncells, ncol=length(peaks)))
    var.tmp <- data.frame(matrix(nrow=ncells, ncol=length(peaks)))
    #calculate expected and variance for each peak
    for (i in 1:length(peaks)) {
      expect.x <- n*param[i]/sum.p
      var.x <- expect.x*(1- param[i]/sum.p)*(n + sum.p)/(1+sum.p)
      expect.tmp[,i] <- expect.x
      var.tmp[,i] <- var.x
    }
    colnames(expect.tmp) <- peaks
    colnames(var.tmp) <- peaks
    return(list(ec = expect.tmp, var = var.tmp))
  }
  # No else / no explicit return here (bug fix): this used to return an
  # all-zero-filled ec/var list on a failed fit, which kept that gene's
  # peaks in the final residual matrix with meaningless zero expected
  # counts/variance (floored to `min.variance` downstream) instead of
  # excluding them. Returning nothing (NULL, implicitly) matches PASTA's
  # original behavior: AssembleResiduals()'s `is.null(res$residuals[[i]])`
  # check then correctly drops this gene's peaks entirely.
}

#' Assemble Residuals Matrices
#'
#' Combines Dirichlet-Multinomial results across all genes into
#' expected counts (ec), variance (var), and observed counts (m) matrices.
#'
#' @param res List output from CalcPolyAResiduals2 containing residuals,
#'   counts, and other components
#' @param verbose Print progress messages. Default is TRUE
#'
#' @return List with three components:
#'   \itemize{
#'     \item ec: Expected counts matrix (features x cells)
#'     \item var: Variance matrix (features x cells)
#'     \item m: Observed counts matrix (features x cells)
#'   }
#'
#' @export
#' @concept residuals
#'
AssembleResiduals <- function(res, verbose = TRUE) {

  if (verbose) {
    message("Combine Dirichlet-Multinomial...")
  }

  # Combine Dirichlet-Multinomial results across all genes.
  # `2:length(res$residuals)` would silently become the descending sequence
  # c(2, 1) if only one gene were processed (R's `:` doesn't check that its
  # left side is <= its right side), causing an out-of-bounds
  # res$residuals[[2]] lookup. seq_along()[-1] returns integer(0) instead in
  # that case, which correctly skips the loop since there's nothing left to
  # combine.
  ec <- res$residuals[[1]]$ec   # expected counts
  var <- res$residuals[[1]]$var # variance
  for (i in seq_along(res$residuals)[-1]) {
    if (!is.null(res$residuals[[i]])) {
      ec <- cbind(ec, res$residuals[[i]]$ec)
      var <- cbind(var, res$residuals[[i]]$var)
    }
  }

  if (verbose) {
    message("Transpose...")
  }

  ec <- t(ec)
  var <- t(var)

  # Bug fix: the per-gene data.frames built in DirichletMultinomial() never
  # had real colnames -- they default to "1","2",...,"ncells" (a data.frame's
  # row index), which become ec/var's COLUMN names after this transpose.
  # PASTA's original explicitly restores real cell barcodes here
  # (`colnames(ec) <- colnames(object)`); this port dropped that step. It
  # didn't error before because CalcResiduals()'s `m - ec` subtraction
  # silently inherits `m`'s (correct) dimnames as the first operand -- but
  # any code that indexes ec/var BY cell name (e.g. subsetting to
  # background cells) needs this to actually be right, so set it from
  # res$counts (== m), which already carries the object's real barcodes in
  # the same cell order.
  colnames(ec) <- colnames(res$counts)
  colnames(var) <- colnames(res$counts)

  gc()

  res$counts <- res$counts[rownames(ec), ]

  return(list(
    ec = ec,
    var = var,
    m =  res$counts
  ))
}


#' Regularized Variance Estimation for Dirichlet Multinomial Model
#'
#' Performs kernel regression-based variance regularization on expected counts
#' and variance matrices from a Dirichlet-Multinomial model. Uses stratified
#' sampling and binning to efficiently estimate smoothed variance across the
#' expected count and gene sum space.
#'
#' @param ec Expected counts matrix (features x cells) from Dirichlet-Multinomial model
#' @param var Variance matrix (features x cells) from Dirichlet-Multinomial model
#' @param m Observed counts matrix (features x cells)
#' @param m.background Background counts matrix for cells used in background distribution
#' @param background.dist Data frame with peak-to-gene mapping and background distribution info
#' @param gene.sum Matrix of total counts per gene per cell
#' @param background.cells Vector of column indices or names for background cells
#' @param min.variance Minimum variance threshold. Values below this are set to this minimum. Default is 0.1
#' @param number.bins Number of bins for regularization grid in each dimension. Default is 30
#' @param sample.n Maximum number of observations to sample per bin for kernel regression. Default is 1000
#' @param verbose Print progress messages. Default is TRUE
#'
#' @return Regularized variance matrix (features x cells) with same dimensions as input
#'
#'
#' @export
#'
RegDMVarPolyA <- function(ec,var, m,
                      m.background,
                      background.dist,
                      gene.sum,
                      background.cells,
                      min.variance = 0.1,
                      number.bins = 30,
                      sample.n = 1000,
                      verbose = TRUE) {

  if(verbose) cat("\n=== RegDMVarPolyA: Regularized Variance Estimation ===\n")

  # ============================================================
  # STEP 1: Prepare gene-to-peak mapping
  # ============================================================

  peaks_in_ec <- rownames(ec)
  background.dist.tmp <- background.dist[background.dist$peak %in% peaks_in_ec, ]

  # Map peaks to genes
  peak_to_gene <- stats::setNames(background.dist.tmp$gene, background.dist.tmp$peak)
  gene_for_peak <- peak_to_gene[peaks_in_ec]

  # Extract gene sums for corresponding peaks
  n_matrix <- gene.sum[gene_for_peak, , drop = FALSE]
  rownames(n_matrix) <- peaks_in_ec

  if(verbose) {
    cat(sprintf("Dataset: %s peaks \u00d7 %s cells = %s total elements\n",
                format(nrow(ec), big.mark=","),
                format(ncol(ec), big.mark=","),
                format(length(ec), big.mark=",")))
  }

  # ============================================================
  # STEP 2: Extract background cells subset
  # Memory optimization: work only with background cells for modeling
  #
  # NOTE (bug fix): this used to do `ec_bg <- ec` (the FULL ec/var/n
  # matrices, all cells) instead of actually subsetting to
  # `background.cells` -- `background.cells` was accepted as a parameter,
  # printed in the message below, and then discarded without ever indexing
  # anything with it. That meant the kernel-regression grid (trained below
  # on `ec_bg`/`var_bg`/`n_bg`) was fit on ALL cells across every treatment
  # group instead of the intended clean background/null distribution,
  # contaminating the variance regularization with real biological signal.
  # PASTA's original (`ec.background <- ec[,background.cells]`) does this
  # subsetting correctly; restored that here.
  # ============================================================

  ec_bg <- ec[, background.cells, drop = FALSE]
  var_bg <- var[, background.cells, drop = FALSE]
  n_bg <- n_matrix[, background.cells, drop = FALSE]

  rm(n_matrix, background.dist.tmp, peak_to_gene)
  gc()

  if(verbose) {
    cat(sprintf("Background subset: %s cells (%s elements)\n",
                format(length(background.cells), big.mark=","),
                format(length(ec_bg), big.mark=",")))
  }


  rm(background.cells)

  # ============================================================
  # STEP 3: Filter valid observations (n > 0)
  # Critical: vectorize only AFTER filtering to reduce memory
  # ============================================================

  valid_mask <- n_bg > 0
  n_valid <- sum(valid_mask)

  if(verbose) {
    cat(sprintf("Valid observations (n>0): %s of %s (%.1f%%)\n",
                format(n_valid, big.mark=","),
                format(length(n_bg), big.mark=","),
                100 * n_valid / length(n_bg)))
  }

  # Vectorize only valid elements
  ec_vec <- as.vector(ec_bg[valid_mask])
  n_vec <- as.vector(n_bg[valid_mask])
  var_vec <- as.vector(var_bg[valid_mask])

  # Free memory
  #rm(ec_bg, m_bg, var_bg, n_bg, valid_mask)
  rm(ec_bg, var_bg, n_bg, valid_mask)
  gc(verbose = FALSE)

  # ============================================================
  # STEP 4: Remove outliers using 99th percentile cutoff
  # ============================================================

  cutoff.ec <- stats::quantile(ec_vec, 0.99, na.rm = TRUE)
  cutoff.n <- stats::quantile(n_vec, 0.99, na.rm = TRUE)

  valid_range <- ec_vec < cutoff.ec & n_vec < cutoff.n

  if(verbose) {
    cat(sprintf("After outlier removal: %s observations (%.1f%% retained)\n",
                format(sum(valid_range), big.mark=","),
                100 * sum(valid_range) / length(valid_range)))
  }

  ec_vec <- ec_vec[valid_range]
  n_vec <- n_vec[valid_range]
  var_vec <- var_vec[valid_range]

  # Calculate data range
  max.n <- max(n_vec, na.rm = TRUE)
  max.ec <- max(ec_vec, na.rm = TRUE)
  min.n <- min(n_vec, na.rm = TRUE)
  min.ec <- min(ec_vec, na.rm = TRUE)

  # ============================================================
  # STEP 5: Bin data into regular grid
  # ============================================================

  lx <- number.bins
  ly <- number.bins

  n_step <- (max.n - min.n) / lx
  ec_step <- (max.ec - min.ec) / ly

  n_grid <- min.n + n_step * 0:lx
  ec_grid <- min.ec + ec_step * 0:ly

  # Assign bins to each observation
  ec_bin <- findInterval(ec_vec, ec_grid)
  n_bin <- findInterval(n_vec, n_grid)
  ec_n <- paste0(ec_bin, "_", n_bin)

  # ============================================================
  # STEP 6: Stratified sampling within bins
  # Reduces data size while preserving distribution
  # ============================================================

  if(verbose) cat("Performing stratified sampling within bins...\n")

  sample_idx <- unlist(lapply(split(seq_along(ec_n), ec_n), function(idx) {
    if (length(idx) <= sample.n) return(idx)
    sampled <- sample(idx, sample.n)
    idx[idx %in% sampled]
  }), use.names = FALSE)

  if(verbose) {
    cat(sprintf("Sampled dataset for kernel regression: %s observations\n",
                format(length(sample_idx), big.mark=",")))
    cat(sprintf("  \u2192 Average ~%.0f observations per bin (from %d\u00d7%d grid)\n",
                length(sample_idx) / (lx * ly), lx, ly))
  }

  # ============================================================
  # STEP 7: Kernel regression on sampled data
  # Estimates smooth variance function on regular grid
  # ============================================================

  x_matrix <- cbind(n_vec[sample_idx], ec_vec[sample_idx])
  y_vector <- var_vec[sample_idx]

  # Create regular grid for kernel estimates
  n_grid_midpoints <- calculate_midpoints(min.n, max.n, lx)
  ec_grid_midpoints <- calculate_midpoints(min.ec, max.ec, ly)

  grid <- matrix(NA, nrow = lx * ly, ncol = 2)
  grid[, 1] <- rep(n_grid_midpoints, length(ec_grid_midpoints))
  grid[, 2] <- rep(ec_grid_midpoints, each = length(n_grid_midpoints))

  if(verbose) cat("Running kernel regression...\n")
  if (!requireNamespace("gplm", quietly = TRUE)) {
    stop("Package 'gplm' is required but not installed. Please install it using: install.packages('gplm')")
  }
  mh <- gplm::kreg(x = x_matrix, y = y_vector, grid = grid)

  # ============================================================
  # STEP 8: Create lookup matrix for regularized variance
  # ============================================================

  # Store regularized variance in 2D matrix: [n_bin, ec_bin]
  reg_var_matrix <- matrix(mh$y, nrow = lx, ncol = ly)

  if(verbose) {
    cat(sprintf("Regularized variance grid: %d\u00d7%d = %d values\n",
                lx, ly, lx * ly))
  }

  # ============================================================
  # STEP 9: Apply regularized variance to full dataset
  # Process in chunks to avoid long vector limitations (>2.1B elements)
  # ============================================================

  if(verbose) cat("\nApplying regularized variance to full dataset...\n")

  n_total_elements <- length(ec)
  chunk_size <- 10e6  # 10 million elements per chunk
  n_chunks <- ceiling(n_total_elements / chunk_size)

  if(verbose) {
    cat(sprintf("Processing %s elements in %d chunks of %.0fM each\n",
                format(n_total_elements, big.mark=","),
                n_chunks,
                chunk_size/1e6))
  }

  # Pre-allocate result matrix
  var_fit <- matrix(NA, nrow = nrow(ec), ncol = ncol(ec))

  # Extract gene sums for corresponding peaks
  n_matrix <- gene.sum[gene_for_peak, , drop = FALSE]
  rownames(n_matrix) <- peaks_in_ec

  rm(
    ec_vec, n_vec, var_vec, n_valid,
    sample_idx, n_step, ec_step,
    x_matrix, y_vector, mh, grid,
    ec_bin, n_bin, ec_n,
    cutoff.ec, cutoff.n, valid_range,
    n_grid_midpoints, ec_grid_midpoints,
    gene_for_peak
  )



  gc()


  # Process in chunks
  for(chunk_i in 1:n_chunks) {
    start_idx <- (chunk_i - 1) * chunk_size + 1
    end_idx <- min(chunk_i * chunk_size, n_total_elements)

    if(verbose) {
      cat(sprintf("  Chunk %d/%d: elements %s to %s",
                  chunk_i, n_chunks,
                  format(start_idx, big.mark=","),
                  format(end_idx, big.mark=",")))
    }

    # Extract chunk
    ec_chunk <- as.vector(ec)[start_idx:end_idx]
    n_chunk <- as.vector(n_matrix)[start_idx:end_idx]
    var_chunk <- as.vector(var)[start_idx:end_idx]

    # Assign bins (sin clipear - permitir out-of-range)
    ec_bin_chunk <- findInterval(ec_chunk, ec_grid)
    n_bin_chunk <- findInterval(n_chunk, n_grid)

    # Identificar bins v\u00e1lidos (dentro del rango 1:lx y 1:ly)
    valid_bins <- n_bin_chunk >= 1 & n_bin_chunk <= lx &
      ec_bin_chunk >= 1 & ec_bin_chunk <= ly

    # Inicializar resultado con NAs
    reg_var_chunk <- rep(NA_real_, length(ec_chunk))

    # Calcular \u00edndice lineal solo para bins v\u00e1lidos
    if(sum(valid_bins) > 0) {
      linear_idx <- (n_bin_chunk[valid_bins] - 1) * ly + ec_bin_chunk[valid_bins]
      reg_var_chunk[valid_bins] <- reg_var_matrix[linear_idx]
    }

    # Fallback a varianza original si regularizada es NA
    missing_idx <- is.na(reg_var_chunk)
    n_missing <- sum(missing_idx)
    if(n_missing > 0) {
      reg_var_chunk[missing_idx] <- var_chunk[missing_idx]
    }

    # Apply minimum variance threshold
    below_threshold <- reg_var_chunk < min.variance
    n_below <- sum(below_threshold)
    if(n_below > 0) {
      reg_var_chunk[below_threshold] <- min.variance
    }

    # Store results
    var_fit[start_idx:end_idx] <- reg_var_chunk

    if(verbose) {
      cat(sprintf(" \u2713 (missing: %s, below threshold: %s)\n",
                  format(n_missing, big.mark=","),
                  format(n_below, big.mark=",")))
    }

    # Clean up memory periodically
    rm(ec_chunk, n_chunk, var_chunk, ec_bin_chunk, n_bin_chunk,
       linear_idx, reg_var_chunk, missing_idx, below_threshold)

    if(chunk_i %% 5 == 0) gc(verbose = FALSE)
  }

  # Restore row and column names
  rownames(var_fit) <- rownames(ec)
  colnames(var_fit) <- colnames(ec)

  if(verbose) {
    cat("\n=== Regularization completed successfully ===\n")
    cat(sprintf("Output: %s \u00d7 %s variance matrix\n", nrow(var_fit), ncol(var_fit)))
    cat(sprintf("Range: [%.4f, %.4f]\n", min(var_fit, na.rm=TRUE), max(var_fit, na.rm=TRUE)))
  }

  return(var_fit)
}


#' Calculate midpoints for regular grid
#' @param a Lower bound
#' @param b Upper bound
#' @param n_intervals Number of intervals
#' @return Vector of midpoint values
#' @noRd
calculate_midpoints <- function(a, b, n_intervals) {
  width <- (b - a) / n_intervals
  first_point <- a + width / 2
  midpoints <- seq(first_point, by = width, length.out = n_intervals)
  return(midpoints)
}

#' Calculate Residuals from Dirichlet Multinomial Model
#'
#' @param object Seurat object
#' @param assay Name of assay to use
#' @param m Observed counts matrix
#' @param ec Expected counts matrix
#' @param var.reg Regularized variances matrix
#' @param do.center Center the residuals? Default FALSE
#' @param do.scale Scale the residuals? Default FALSE
#' @param residuals.max Clip residuals above this value. Default NULL (no clipping)
#' @param residuals.min Clip residuals below this value. Default NULL (no clipping)
#' @param verbose Print messages? Default TRUE
#'
#' @return Seurat object with residuals in scale.data layer
#'
#' @concept residuals
#' @noRd
#'
CalcResiduals <- function(
    object,
    assay = "polyA",
    m,          # observed counts
    ec,         # expected counts
    var.reg,    # regularized variances
    do.center = FALSE,
    do.scale = FALSE,
    residuals.max = NULL,
    residuals.min = NULL,
    verbose = TRUE
) {

  if (verbose) message("Calculating residuals matrix...")

  # Calculate residual matrix
  residual.matrix <- (as.matrix(m) - as.matrix(ec)) / sqrt(var.reg)
  residual.matrix <- as.matrix(residual.matrix, nrow = nrow(residual.matrix))

  if (verbose) message("Residuals calculated for ", nrow(residual.matrix),
                      " features and ", ncol(residual.matrix), " cells")

  # Reorder to match counts slot
  features.order <- rownames(residual.matrix)[order(
    match(
      rownames(residual.matrix),
      rownames(SeuratObject::LayerData(object = object, assay = assay, layer = "counts"))
    )
  )]
  residual.matrix <- residual.matrix[features.order, ]

  # Scale and center if requested
  if (do.center || do.scale) {
    if (verbose) message("Scaling residuals (center=", do.center,
                        ", scale=", do.scale, ")...")
    residual.matrix <- scale(residual.matrix, center = do.center, scale = do.scale)
  }

  # Clip residuals if thresholds provided
  if (!is.null(residuals.max)) {
    if (verbose) message("Clipping residuals at max = ", residuals.max)
    residual.matrix[residual.matrix > residuals.max] <- residuals.max
  }

  if (!is.null(residuals.min)) {
    if (verbose) message("Clipping residuals at min = ", residuals.min)
    residual.matrix[residual.matrix < residuals.min] <- residuals.min
  }

  # Set default assay
  Seurat::DefaultAssay(object = object) <- assay

  # Store residuals in scale.data layer
  if (verbose) message("Storing residuals in scale.data layer...")
  SeuratObject::LayerData(object, assay = assay, layer = "scale.data") <- residual.matrix

  # NOTE (bug fix): Seurat::LogSeuratCommand() is intentionally NOT called
  # here anymore. It captures every formal argument of its *caller* --
  # verbatim, with no size-based trimming for plain matrices (confirmed in
  # SeuratObject's own source, R/command.R: everything except an argument
  # that `inherits(x, "Seurat")` gets stored as-is in
  # object@commands[[...]]@params, with a literal
  # "#TODO Institute some check of object size?" left in the code). This
  # helper's own formals include the full `m`/`ec`/`var.reg` matrices
  # (each features x cells, tens to hundreds of MB) -- calling
  # LogSeuratCommand() from inside CalcResiduals() silently stored all
  # three of those matrices a second time inside the command log, which is
  # the most likely explanation for the object being several hundred MB
  # larger than PASTA's. It's now logged once, at the end of
  # CalcPolyAResidualsPolyA() itself, whose formals are just small
  # scalar/vector arguments -- mirroring PASTA's original structure, where
  # CalcPolyAResiduals() is a single flat function and the log call
  # naturally happens in that frame.
  if (verbose) message("Residuals calculation completed successfully")

  return(object)
}
