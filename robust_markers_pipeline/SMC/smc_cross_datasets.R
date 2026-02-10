suppressPackageStartupMessages({
  .libPaths(c("~/seurat_v5/", .libPaths())) # using seurat v5
  library(Seurat)
  library(CellChat)
  library(patchwork)
  library(tidyverse)
  library(glue)
})
setwd("/home/data/fhz/project/reference_map_PH")
source("~/code-hub/code-hub/single-cell/io.R")
source("utils/smc_cross_dataset_utils.R")

# seu <- LoadH5ad("data/lung/13_smc_all_celltypes.h5ad")
# qs::qsave(seu, file="data/lung/35_smc.qs")
seu <- qs::qread("data/lung/35_smc.qs")
seu
# An object of class Seurat 
# 19007 features across 18860 samples within 1 assay 
# Active assay: RNA (19007 features, 0 variable features)
# 1 layer present: counts

# =========================
# Run full pipeline and save outputs
# =========================
output_dir <- "results/smc_cross_datasets"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

params <- list(
  assay = "RNA",
  sample_col = "sample",
  cluster_col = "final_celltype",
  study_col = "study",
  min_cells = 30,
  min_samples_cluster = 3,
  min_cpm = 1,
  min_samples_cpm = 3,
  fdr_max = 0.05,
  min_median_logFC = 0.25,
  min_prop_same_sign = 0.6,
  min_expr_prop_in_A = 0.5
)

full_out <- run_full_marker_pipeline(
  obj = seu,
  assay = params$assay,
  sample_col = params$sample_col,
  cluster_col = params$cluster_col,
  study_col = params$study_col,
  min_cells = params$min_cells,
  min_samples_cluster = params$min_samples_cluster,
  min_cpm = params$min_cpm,
  min_samples_cpm = params$min_samples_cpm,
  fdr_max = params$fdr_max,
  min_median_logFC = params$min_median_logFC,
  min_prop_same_sign = params$min_prop_same_sign,
  min_expr_prop_in_A = params$min_expr_prop_in_A
)

summary_path <- file.path(output_dir, "phase4_summary.csv")
write.csv(full_out$phase4$summary, summary_path, row.names = FALSE)

safe_name <- function(x) {
  gsub("[^A-Za-z0-9._-]", "_", x)
}

status_ok <- full_out$phase4$summary$cluster[full_out$phase4$summary$status == "ok"]
for (cluster_id in status_ok) {
  marker_df <- full_out$phase4$results[[cluster_id]]$markers
  out_path <- file.path(output_dir, paste0("markers__", safe_name(cluster_id), ".csv"))
  write.csv(marker_df, out_path, row.names = FALSE)
}

qs::qsave(full_out, file = file.path(output_dir, "full_out.qs"))
saveRDS(params, file = file.path(output_dir, "run_params.rds"))

cat("Pipeline finished.\n")
cat("Output directory:", normalizePath(output_dir), "\n")
cat("Summary file:", normalizePath(summary_path), "\n")
