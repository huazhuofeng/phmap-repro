suppressPackageStartupMessages({
  library(jsonlite)
})

source(file.path(dirname(dirname(normalizePath(sys.frame(1)$ofile %||% "scripts/run_pipeline.R"))), "R", "robust_marker_utils.R"))

`%||%` <- function(x, y) if (is.null(x)) y else x

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[i]
    if (!startsWith(key, "--")) stop("Invalid arg: ", key)
    key <- sub("^--", "", key)
    if (i == length(args)) stop("Missing value for --", key)
    out[[key]] <- args[i + 1]
    i <- i + 2
  }
  out
}

read_seurat_object <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "qs") {
    if (!requireNamespace("qs", quietly = TRUE)) stop("Input is .qs but package qs is not installed")
    return(qs::qread(path))
  }
  if (ext == "rds") return(readRDS(path))
  stop("Unsupported input format: ", ext, ". Use .qs or .rds")
}

save_outputs <- function(full_out, outdir, params) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  summary_path <- file.path(outdir, "phase4_summary.csv")
  write.csv(full_out$phase4$summary, summary_path, row.names = FALSE)

  safe_name <- function(x) gsub("[^A-Za-z0-9._-]", "_", x)

  ok_clusters <- full_out$phase4$summary$cluster[full_out$phase4$summary$status == "ok"]
  for (cluster_id in ok_clusters) {
    marker_df <- full_out$phase4$results[[cluster_id]]$markers
    consistency_df <- full_out$phase4$results[[cluster_id]]$consistency

    write.csv(marker_df, file.path(outdir, paste0("markers__", safe_name(cluster_id), ".csv")), row.names = FALSE)
    write.csv(consistency_df, file.path(outdir, paste0("consistency__", safe_name(cluster_id), ".csv")), row.names = FALSE)
  }

  saveRDS(full_out, file = file.path(outdir, "full_out.rds"))
  if (requireNamespace("qs", quietly = TRUE)) {
    qs::qsave(full_out, file = file.path(outdir, "full_out.qs"))
  }

  saveRDS(params, file = file.path(outdir, "run_params.rds"))
  write_json(params, path = file.path(outdir, "run_params.json"), auto_unbox = TRUE, pretty = TRUE)
}

main <- function() {
  args <- parse_args()

  required <- c("input", "outdir")
  miss <- setdiff(required, names(args))
  if (length(miss) > 0) stop("Missing required args: ", paste(miss, collapse = ", "))

  params <- list(
    assay = args$assay %||% "RNA",
    sample_col = args$sample_col %||% "sample",
    cluster_col = args$cluster_col %||% "final_celltype",
    study_col = args$study_col %||% "study",
    min_cells = as.numeric(args$min_cells %||% "30"),
    min_samples_cluster = as.numeric(args$min_samples_cluster %||% "3"),
    min_cpm = as.numeric(args$min_cpm %||% "1"),
    min_samples_cpm = as.numeric(args$min_samples_cpm %||% "3"),
    fdr_max = as.numeric(args$fdr_max %||% "0.05"),
    min_median_logFC = as.numeric(args$min_median_logFC %||% "0.25"),
    min_prop_same_sign = as.numeric(args$min_prop_same_sign %||% "0.6"),
    min_expr_prop_in_A = as.numeric(args$min_expr_prop_in_A %||% "0.5")
  )

  message("[INFO] reading object: ", args$input)
  obj <- read_seurat_object(args$input)

  message("[INFO] running pipeline...")
  full_out <- run_full_marker_pipeline(
    obj = obj,
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

  message("[INFO] writing outputs to: ", args$outdir)
  save_outputs(full_out, args$outdir, params)
  message("[INFO] done")
}

main()
