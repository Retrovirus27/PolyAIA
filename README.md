# PolyAIA — PolyA Isoform Analysis

## Overview

PolyAIA is an R package for **single-cell alternative polyadenylation (APA)**
analysis. It wraps an end-to-end workflow around
[Seurat](https://github.com/satijalab/seurat),
[Signac](https://github.com/stuart-lab/signac) and
[PASTA](https://github.com/satijalab/PASTA): loading and annotating single-cell
data, quality control with per-cell-type thresholds, normalization and
integration, building polyA site assays from peak counts, and differential
polyA usage analysis with RED scores.

The pipeline runs in three stages:

| Stage | What it does | Key functions |
|-------|--------------|---------------|
| **1. Single cell** | Load 10X data, annotate genes + cell types, QC, normalize/integrate | `UploadSce()`, `QCPlot()`, `QCFilter()`, `SeuratPipeline()`, `UmapPlot()`, `DotPlot()` |
| **2. PolyA assay** | Build per-sample polyA assays, unify peaks, requantify, annotate | `UploadPolyAAssays()`, `FilterPeaks()`, `RequantifyPolyA()`, `Peaksdb()` |
| **3. Differential APA** | Residuals, differential polyA usage, RED scores, plots | `DEPsMatrix()`, `DEPolyAPeaks()`, `PolyAPlot()` |

## Installation

```r
# install remotes and BiocManager if necessary
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

remotes::install_github("satijalab/PASTA")
remotes::install_github("Retrovirus27/PolyAIA")
```

Several dependencies come from Bioconductor and are not installed
automatically — install any that the command above reports as missing with
`BiocManager::install()`.

## Quick start

`DEPsMatrix()` runs the whole differential APA analysis in one call. It
requires:

* A **Seurat object** with an annotated `polyA` assay and an `RNA` assay
  (stages 1–2 above build it from raw 10X and polyA-peak files);
* A **comparisons table**, one row per contrast;
* Optionally, a **RepeatMasker table** to drop transposable-element sites.

```r
library(PolyAIA)

comparisons <- data.frame(
  ident1    = c("Alcohol_B6", "Control_3xTg"),
  ident2    = c("Control_B6", "Control_B6"),
  treatment = c("B6_Alcohol", "3xTg_Control"),
  control   = c("B6_Control", "B6_Control"),
  Condition = c("B6_Alcohol_vs_B6_Control", "3xTg_Control_vs_B6_Control"),
  stringsAsFactors = FALSE
)

deps <- DEPsMatrix(
  seu             = seurat_polyA,   # Seurat object with RNA + polyA assays
  comparisons     = comparisons,
  celltype_col    = "Subpopulation",
  condition_cols  = c("Treatment", "Strain"),
  group.by        = c("Strain", "Treatment", "Subpopulation", "Number"),
  background      = "ControlB6",    # identity used as the null distribution
  background_cols = c("Treatment", "Strain"),
  background_sep  = "",
  Pasta           = TRUE,           # FALSE = skip residuals/DE, RED from counts only
  min_cells       = 10,             # drop tiny subpopulations
  repeat_masker   = "Remove",       # drop intronic TE-overlapping sites
  rmsk            = rmsk,
  return_seurat   = TRUE
)

deps$red_scores   # RED score table
deps$polyAdb      # analyzed peak-level annotation
deps$seu          # residual-bearing Seurat object, for PolyAPlot()
```

Please note the following options in the `DEPsMatrix` function:

* By default, `Pasta = TRUE` computes polyA residuals and runs the
  differential test. Set `Pasta = FALSE` to skip both and derive RED scores
  from pseudobulk counts only.
* `repeat_masker` controls how polyA sites overlapping transposable elements
  are treated: `"none"`, `"Remove"` or `"Keep"` (flagged but retained).
  `rmsk` is required for the last two.
* `features` is derived automatically from `genomic_positions` (3'-most exon
  and intron by default); pass it explicitly to restrict the analysis.

Please refer to `?DEPsMatrix` for the full list of arguments and detailed
usage.

## Example

A complete walkthrough on a mouse cortex dataset (3xTg / B6, control vs.
alcohol) — QC, integration, cell-type visualization, differential APA and
coverage plots — is available in
[`examples/PolyAIA_example.Rmd`](examples/PolyAIA_example.Rmd).

## Key features

* **`QCFilter()`** — QC thresholds that vary **per cell type** rather than one
  global cutoff, since RNA content differs a lot between populations.
* **`SeuratPipeline()`** — normalize → variable features → scale → PCA →
  Harmony → UMAP in one call, picking the number of dimensions at the **elbow**
  (point of maximum curvature) of the PCA standard-deviation curve.
* **`DEPsMatrix()`** — integrated differential APA: optional PASTA residuals
  and differential testing, pseudobulk aggregation, RepeatMasker handling and
  per-gene RED scores.
* **`PolyAPlot()`** — coverage tracks with polyA sites, peaks and gene
  annotation, filterable by cell population and comparison.

## Links

[Seurat](https://github.com/satijalab/seurat): R package for the analysis,
integration, and exploration of scRNA-seq data

[Signac](https://github.com/stuart-lab/signac): R package for the analysis of
single-cell chromatin and other genomic-range assays

[PASTA](https://github.com/satijalab/PASTA): PolyA Site analysis using relative
Transcript Abundance, the framework PolyAIA's residual and differential models
build on

[MAAPER](https://github.com/Vivianstats/MAAPER): Model-based analysis of
alternative polyadenylation using 3' end-linked reads

## License

MIT — see [LICENSE](LICENSE).
