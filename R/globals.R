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
  # ClusterPlots elbow/variance plots (R/2_ScePlots.R)
  "PC", "Stdev", "CumulativeStdev", "color",
  # plot_annotation_from_granges gene track (R/3_PolyAAnalysis.R)
  "has_cds", "block_h"
))
