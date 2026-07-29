utils::globalVariables(c(
  "%>%", "..count..", ".data", ".x",
  "gene", "gene_id", "gene_name", "gene_biotype",
  "tx_id", "seqnames", "start", "end", "strand", "width",
  "Sample", "Condition", "study",
  "Metric", "Stage", "Value", "Value_no_outlier",
  "Q1", "Q3", "IQR",
  "Percent_Removed", "Num_Cells", "Total_Cells", "Cell_Count",
  "detected", "qc_fail", "subsets_Mito_percent",
  "x_mean", "y_mean", "label_text",
  "gene_start", "gene_end",
  "value", "cells", "Variable", "total",
  "ID...1", "ID...4",
  "ENTREZID", "ENZYME", "GENETYPE", "SYMBOL",
  "TRUE", "FALSE",
  "expr", "avg_exp", "pct_exp", "avg_xp_scaled", "GeneGroup", "fisher_p",
  "PC1", "PC2", "colData", "Gene_symbol", "Cells",
  "Sex", "peak", "p_val_adj", "symbol",
  "CellType", ".", "Intron.exon.location", "polyA_treatment", "polyA_control",
  "prop_treatment", "prop_change", "stranded", "pos", "desc",
  "prop_control", "hg38_Position", "fisher_pval", "n",
  # DEPsMatrix
  "d", "p", "p_region", "RED_type", "slice_max",
  "p_percent_change", "d_percent_change",
  "p_treatment_count", "d_treatment_count",
  "p_control_count", "d_control_count",
  "RED", "p_peak", "d_peak",
  # FilterPeaks
  "n_peaks",
  # DEPsMatrix
  ".env", "RLD", "d_treatment_prop", "p_treatment_prop", "d_control_prop", "p_control_prop",
  "pchisq", "p_padj_PASTA", "d_padj_PASTA",
  # Default-argument package objects
  "Homo.sapiens", "TxDb.Hsapiens.UCSC.hg19.knownGene",
  # AnnotateGenes / plot_annotation_from_granges
  "type", "x", "xend", "y",
  # DEPsMatrix — polyA site columns
  "p_site_id", "d_site_id",
  "p_pos", "d_pos",
  "p_region", "d_region",
  "p_PAS_hexamer", "d_PAS_hexamer",
  "p_PolyaStrength", "d_PolyaStrength",
  "p_estimate_PASTA", "d_estimate_PASTA",
  "p_pval_PASTA", "d_pval_PASTA",
  "p_percent1_PASTA", "d_percent1_PASTA",
  "p_percent2_PASTA", "d_percent2_PASTA",
  "RED_fisher_pval", "RED_fisher_padj",
  "PASTA_pval", "PASTA_padj",
  "RED_Direction", "RLD_direction",
  "chr", "npas",
  # DEPsMatrix — stats
  "p.adjust",
  # UploadSce
  "ID"
))
