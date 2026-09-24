# Register column names used in non-standard evaluation (dplyr/ggplot2 pipelines
# and data.frame column references) so R CMD check doesn't flag them as
# "no visible binding for global variable". These are not real global variables;
# they are column names resolved at runtime inside data masks.
utils::globalVariables(c(
  # RED-score / DEPsMatrix pipeline (R/3_PolyAAnalysis.R)
  "Ensembl_ID", "ensembl_gene_id",
  "RED_LRT_pval", "RED_LRT_padj",
  "RED_LRT_cov_pval", "RED_LRT_cov_padj",
  "gDPAU", "DPAU",
  # RED pair-type assignment (R/3_PolyAAnalysis.R)
  "RED_type", "p_region", "d_region",
  # ClusterPlots elbow/variance plots (R/2_ScePlots.R)
  "PC", "Stdev", "CumulativeStdev", "color",
  # plot_annotation_from_granges gene track (R/3_PolyAAnalysis.R)
  "has_cds", "block_h",
  # RED quantifiability columns (R/3_PolyAAnalysis.R)
  "p_reads_treatment", "p_reads_control", "d_reads_treatment", "d_reads_control",
  "p_usage", "d_usage", "rep_pairs_frac",
  "ok_reads", "ok_usage", "ok_replicates", "quantifiable"
))
