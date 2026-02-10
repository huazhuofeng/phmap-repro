suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

phase0_validate_inputs <- function(
  obj,
  sample_col = "sample",
  cluster_col = "final_celltype",
  study_col = "study",
  min_cells = 30
) {
  stopifnot(inherits(obj, "Seurat"))
  md <- obj@meta.data

  required_cols <- c(sample_col, cluster_col, study_col)
  miss_cols <- setdiff(required_cols, colnames(md))
  if (length(miss_cols) > 0) {
    stop("Missing required metadata columns: ", paste(miss_cols, collapse = ", "))
  }

  meta_clean <- data.frame(
    cell = rownames(md),
    sample = as.character(md[[sample_col]]),
    cluster = as.character(md[[cluster_col]]),
    study = as.character(md[[study_col]]),
    stringsAsFactors = FALSE
  )

  bad_na <- is.na(meta_clean$sample) | is.na(meta_clean$cluster) | is.na(meta_clean$study)
  if (any(bad_na)) stop("Found NA in sample/cluster/study.")

  bad_empty <- meta_clean$sample == "" | meta_clean$cluster == "" | meta_clean$study == ""
  if (any(bad_empty)) stop("Found empty string in sample/cluster/study.")

  sample_study_n <- tapply(meta_clean$study, meta_clean$sample, function(x) length(unique(x)))
  bad_samples <- names(sample_study_n)[sample_study_n != 1]
  if (length(bad_samples) > 0) {
    stop("Inconsistent sample->study mapping: ", paste(bad_samples, collapse = ", "))
  }

  pair_tab <- table(meta_clean$sample, meta_clean$cluster)
  pair_cell_counts <- as.data.frame(pair_tab, stringsAsFactors = FALSE)
  colnames(pair_cell_counts) <- c("sample", "cluster", "n_cells")

  list(
    meta_clean = meta_clean,
    pair_cell_counts = pair_cell_counts,
    thresholds = list(min_cells = min_cells),
    ok = TRUE
  )
}

phase1_build_pseudobulk <- function(
  obj,
  phase0_out,
  assay = "RNA",
  min_cells = 30,
  pb_sep = "___PB___",
  sample_col = "sample",
  cluster_col = "final_celltype",
  study_col = "study"
) {
  stopifnot(inherits(obj, "Seurat"))
  if (!assay %in% names(obj@assays)) stop("Assay not found: ", assay)

  meta_clean <- phase0_out$meta_clean
  pair_tab <- table(meta_clean$sample, meta_clean$cluster)
  keep_pairs_idx <- which(pair_tab >= min_cells, arr.ind = TRUE)
  if (nrow(keep_pairs_idx) == 0) stop("No (sample, cluster) pairs pass min_cells")

  keep_samples <- rownames(pair_tab)[keep_pairs_idx[, "row"]]
  keep_clusters <- colnames(pair_tab)[keep_pairs_idx[, "col"]]
  keep_pair_key <- paste0(keep_samples, pb_sep, keep_clusters)

  cell_pair_key <- paste0(meta_clean$sample, pb_sep, meta_clean$cluster)
  keep_cells <- meta_clean$cell[cell_pair_key %in% keep_pair_key]
  obj_filtered <- subset(obj, cells = keep_cells)

  md_f <- obj_filtered@meta.data
  std_meta <- data.frame(
    cell = rownames(md_f),
    sample = as.character(md_f[[sample_col]]),
    cluster = as.character(md_f[[cluster_col]]),
    study = as.character(md_f[[study_col]]),
    stringsAsFactors = FALSE
  )

  std_meta$pair_key <- paste0(std_meta$sample, pb_sep, std_meta$cluster)
  pair_levels <- unique(std_meta[, c("pair_key", "sample", "cluster")])
  pair_levels$pb_id <- sprintf("PB%06d", seq_len(nrow(pair_levels)))
  rownames(pair_levels) <- pair_levels$pair_key

  std_meta$pb_id <- unname(pair_levels[std_meta$pair_key, "pb_id"])
  obj_filtered$.__pb_id__ <- std_meta$pb_id

  DefaultAssay(obj_filtered) <- assay
  pb_list <- AggregateExpression(
    obj_filtered,
    assays = assay,
    slot = "counts",
    group.by = ".__pb_id__",
    return.seurat = FALSE
  )
  counts_pb <- pb_list[[assay]]

  pb_levels <- pair_levels[, c("pb_id", "sample", "cluster")]
  rownames(pb_levels) <- pb_levels$pb_id
  sample_to_study <- tapply(std_meta$study, std_meta$sample, function(x) unique(x)[1])

  coldata_pb <- data.frame(
    pb_id = colnames(counts_pb),
    sample = pb_levels[colnames(counts_pb), "sample"],
    cluster = pb_levels[colnames(counts_pb), "cluster"],
    study = unname(sample_to_study[pb_levels[colnames(counts_pb), "sample"]]),
    stringsAsFactors = FALSE
  )
  rownames(coldata_pb) <- coldata_pb$pb_id

  list(
    obj_filtered = obj_filtered,
    counts_pb = counts_pb,
    coldata_pb = coldata_pb,
    qc = list(
      n_cells_input = ncol(obj),
      n_cells_filtered = ncol(obj_filtered),
      n_pb_columns = ncol(counts_pb),
      min_cells = min_cells,
      assay = assay
    )
  )
}

run_phase01_minimal <- function(
  obj,
  assay = "RNA",
  sample_col = "sample",
  cluster_col = "final_celltype",
  study_col = "study",
  min_cells = 30
) {
  p0 <- phase0_validate_inputs(
    obj = obj,
    sample_col = sample_col,
    cluster_col = cluster_col,
    study_col = study_col,
    min_cells = min_cells
  )
  p1 <- phase1_build_pseudobulk(
    obj = obj,
    phase0_out = p0,
    assay = assay,
    min_cells = min_cells,
    sample_col = sample_col,
    cluster_col = cluster_col,
    study_col = study_col
  )
  list(phase0 = p0, phase1 = p1)
}

phase2_make_A_vs_rest <- function(counts_pb, coldata_pb, A) {
  if (!all(colnames(counts_pb) == rownames(coldata_pb))) stop("count/meta mismatch")
  if (!A %in% coldata_pb$cluster) stop("A not in cluster list: ", A)

  samples_with_A <- unique(coldata_pb$sample[coldata_pb$cluster == A])
  out_counts <- list()
  out_meta <- list()

  for (sample_id in samples_with_A) {
    idx_s <- which(coldata_pb$sample == sample_id)
    idx_A <- idx_s[coldata_pb$cluster[idx_s] == A]
    idx_rest <- idx_s[coldata_pb$cluster[idx_s] != A]
    if (length(idx_A) != 1 || length(idx_rest) < 1) next

    count_A <- counts_pb[, idx_A, drop = FALSE]
    count_rest <- Matrix::Matrix(Matrix::rowSums(counts_pb[, idx_rest, drop = FALSE]), sparse = TRUE)
    colnames(count_A) <- paste0(sample_id, "__A")
    colnames(count_rest) <- paste0(sample_id, "__rest")
    out_counts[[length(out_counts) + 1]] <- cbind(count_A, count_rest)

    study_s <- unique(coldata_pb$study[idx_s])
    if (length(study_s) != 1) next
    out_meta[[length(out_meta) + 1]] <- data.frame(
      pb_id = c(paste0(sample_id, "__A"), paste0(sample_id, "__rest")),
      sample = sample_id,
      study = study_s,
      group = c("A", "rest"),
      stringsAsFactors = FALSE
    )
  }

  if (length(out_counts) == 0) stop("No usable samples for A vs rest")

  counts_ar <- do.call(cbind, out_counts)
  meta_ar <- do.call(rbind, out_meta)
  rownames(meta_ar) <- meta_ar$pb_id
  counts_ar <- counts_ar[, meta_ar$pb_id, drop = FALSE]

  list(
    counts_ar = counts_ar,
    meta_ar = meta_ar,
    qc = list(
      target_cluster = A,
      n_samples_with_A = length(samples_with_A),
      n_samples_usable = length(unique(meta_ar$sample))
    )
  )
}

phase2_run_edger_A_vs_rest <- function(counts_ar, meta_ar, min_cpm = 1, min_samples_cpm = 3) {
  if (!requireNamespace("edgeR", quietly = TRUE)) stop("Please install edgeR")
  if (!all(colnames(counts_ar) == rownames(meta_ar))) stop("count/meta mismatch")

  y <- edgeR::DGEList(counts = counts_ar)
  y <- edgeR::calcNormFactors(y)
  keep <- rowSums(edgeR::cpm(y) > min_cpm) >= min_samples_cpm
  y <- y[keep, , keep.lib.sizes = FALSE]
  if (nrow(y) == 0) stop("No genes passed CPM filter")

  meta_ar$group <- factor(meta_ar$group, levels = c("rest", "A"))
  meta_ar$study <- factor(meta_ar$study)
  design <- model.matrix(~ 0 + group + study, data = meta_ar)
  colnames(design) <- make.names(colnames(design))

  y <- edgeR::estimateDisp(y, design)
  fit <- edgeR::glmQLFit(y, design, robust = TRUE)

  coef_A <- grep("^groupA$", colnames(design), value = TRUE)
  coef_rest <- grep("^grouprest$", colnames(design), value = TRUE)
  if (length(coef_A) != 1 || length(coef_rest) != 1) stop("Invalid design columns")

  contrast <- rep(0, ncol(design)); names(contrast) <- colnames(design)
  contrast[coef_A] <- 1; contrast[coef_rest] <- -1

  qlf <- edgeR::glmQLFTest(fit, contrast = contrast)
  de_table <- edgeR::topTags(qlf, n = Inf)$table
  de_table$gene <- rownames(de_table)

  list(
    de_table = de_table,
    fit_info = list(
      n_genes_input = nrow(counts_ar),
      n_genes_tested = nrow(de_table),
      n_samples = length(unique(meta_ar$sample))
    )
  )
}

phase3_compute_consistency <- function(counts_ar, meta_ar) {
  if (!requireNamespace("edgeR", quietly = TRUE)) stop("Please install edgeR")

  y <- edgeR::DGEList(counts = counts_ar)
  y <- edgeR::calcNormFactors(y)
  logcpm <- edgeR::cpm(y, log = TRUE, prior.count = 1)

  samples <- unique(meta_ar$sample)
  logFC_by_sample <- matrix(NA_real_, nrow = nrow(logcpm), ncol = length(samples), dimnames = list(rownames(logcpm), samples))
  exprA_by_sample <- matrix(NA_real_, nrow = nrow(counts_ar), ncol = length(samples), dimnames = list(rownames(counts_ar), samples))

  for (sample_id in samples) {
    col_A <- rownames(meta_ar)[meta_ar$sample == sample_id & meta_ar$group == "A"]
    col_R <- rownames(meta_ar)[meta_ar$sample == sample_id & meta_ar$group == "rest"]
    if (length(col_A) != 1 || length(col_R) != 1) next
    logFC_by_sample[, sample_id] <- logcpm[, col_A] - logcpm[, col_R]
    exprA_by_sample[, sample_id] <- as.numeric(counts_ar[, col_A] > 0)
  }

  median_logFC <- apply(logFC_by_sample, 1, median, na.rm = TRUE)
  direction <- sign(median_logFC)

  prop_same_sign <- vapply(seq_len(nrow(logFC_by_sample)), function(i) {
    v <- logFC_by_sample[i, ]
    v <- v[is.finite(v)]
    if (length(v) == 0) return(NA_real_)
    if (direction[i] == 0) return(0)
    mean(sign(v) == direction[i])
  }, numeric(1))

  expr_prop_in_A <- apply(exprA_by_sample, 1, mean, na.rm = TRUE)

  data.frame(
    gene = rownames(logFC_by_sample),
    median_logFC_by_sample = median_logFC,
    prop_same_sign = prop_same_sign,
    expr_prop_in_A = expr_prop_in_A,
    stringsAsFactors = FALSE
  )
}

phase3_rank_markers <- function(
  de_table,
  consistency_table,
  fdr_max = 0.05,
  min_median_logFC = 0.25,
  min_prop_same_sign = 0.6,
  min_expr_prop_in_A = 0.5
) {
  merged <- merge(de_table, consistency_table, by = "gene", all.x = TRUE)
  merged <- merged[!is.na(merged$FDR), , drop = FALSE]
  merged <- merged[merged$FDR <= fdr_max, , drop = FALSE]
  merged <- merged[merged$median_logFC_by_sample >= min_median_logFC, , drop = FALSE]
  merged <- merged[merged$prop_same_sign >= min_prop_same_sign, , drop = FALSE]
  merged <- merged[merged$expr_prop_in_A >= min_expr_prop_in_A, , drop = FALSE]

  if (nrow(merged) == 0) {
    merged$stability_score <- numeric(0)
    return(merged)
  }

  merged$stability_score <- merged$median_logFC_by_sample * merged$prop_same_sign
  merged <- merged[order(-merged$stability_score, -merged$median_logFC_by_sample, -merged$logFC, merged$FDR), , drop = FALSE]
  merged
}

run_phase4_all_clusters <- function(
  phase1_out,
  min_samples_cluster = 3,
  min_cpm = 1,
  min_samples_cpm = 3,
  fdr_max = 0.05,
  min_median_logFC = 0.25,
  min_prop_same_sign = 0.6,
  min_expr_prop_in_A = 0.5
) {
  clusters <- sort(unique(phase1_out$coldata_pb$cluster))
  results <- vector("list", length(clusters)); names(results) <- clusters
  summary_rows <- vector("list", length(clusters))

  for (i in seq_along(clusters)) {
    cluster_id <- clusters[i]
    ar <- tryCatch(
      phase2_make_A_vs_rest(phase1_out$counts_pb, phase1_out$coldata_pb, cluster_id),
      error = function(e) e
    )

    if (inherits(ar, "error")) {
      results[[cluster_id]] <- list(status = "failed", reason = conditionMessage(ar), markers = NULL)
      summary_rows[[i]] <- data.frame(cluster = cluster_id, status = "failed", n_samples_usable = 0, n_markers = 0, note = conditionMessage(ar), stringsAsFactors = FALSE)
      next
    }

    n_samples_usable <- ar$qc$n_samples_usable
    if (n_samples_usable < min_samples_cluster) {
      note_txt <- paste0("Skipped: n_samples_usable=", n_samples_usable, " < ", min_samples_cluster)
      results[[cluster_id]] <- list(status = "skipped", reason = note_txt, n_samples_usable = n_samples_usable, markers = NULL)
      summary_rows[[i]] <- data.frame(cluster = cluster_id, status = "skipped", n_samples_usable = n_samples_usable, n_markers = 0, note = note_txt, stringsAsFactors = FALSE)
      next
    }

    de <- phase2_run_edger_A_vs_rest(ar$counts_ar, ar$meta_ar, min_cpm, min_samples_cpm)
    cons <- phase3_compute_consistency(ar$counts_ar, ar$meta_ar)
    ranked <- phase3_rank_markers(de$de_table, cons, fdr_max, min_median_logFC, min_prop_same_sign, min_expr_prop_in_A)

    results[[cluster_id]] <- list(status = "ok", n_samples_usable = n_samples_usable, n_genes_tested = de$fit_info$n_genes_tested, markers = ranked, consistency = cons)
    summary_rows[[i]] <- data.frame(cluster = cluster_id, status = "ok", n_samples_usable = n_samples_usable, n_markers = nrow(ranked), note = "", stringsAsFactors = FALSE)
  }

  summary_table <- do.call(rbind, summary_rows)
  rownames(summary_table) <- NULL

  list(
    summary = summary_table,
    results = results,
    params = list(
      min_samples_cluster = min_samples_cluster,
      min_cpm = min_cpm,
      min_samples_cpm = min_samples_cpm,
      fdr_max = fdr_max,
      min_median_logFC = min_median_logFC,
      min_prop_same_sign = min_prop_same_sign,
      min_expr_prop_in_A = min_expr_prop_in_A
    )
  )
}

run_full_marker_pipeline <- function(
  obj,
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
) {
  out01 <- run_phase01_minimal(
    obj = obj,
    assay = assay,
    sample_col = sample_col,
    cluster_col = cluster_col,
    study_col = study_col,
    min_cells = min_cells
  )

  out04 <- run_phase4_all_clusters(
    phase1_out = out01$phase1,
    min_samples_cluster = min_samples_cluster,
    min_cpm = min_cpm,
    min_samples_cpm = min_samples_cpm,
    fdr_max = fdr_max,
    min_median_logFC = min_median_logFC,
    min_prop_same_sign = min_prop_same_sign,
    min_expr_prop_in_A = min_expr_prop_in_A
  )

  list(phase0 = out01$phase0, phase1 = out01$phase1, phase4 = out04)
}
