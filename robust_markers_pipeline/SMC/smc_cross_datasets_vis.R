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

seu <- qs::qread("data/lung/35_smc.qs")
seu

seu <- NormalizeData(seu)

files_to_use <- list.files("results/smc_cross_datasets")
files_to_use <- files_to_use[grepl("^markers_", files_to_use)]

markers.list <- list()
for (f in files_to_use) {
  sample_name = str_extract(f, "SMC_c\\d+")
  print(glue("🚀 {sample_name}"))
  marker <- read.csv(glue("results/smc_cross_datasets/{f}"))
  markers.list[[sample_name]] <- marker
}

gene_vec <- unlist(lapply(markers.list, function(x) x$gene[1:10]))

ht.mtx.list <- list()
studies <- levels(seu$study)
for (s in studies) {
   mtx.avg <- AverageExpression(
     seu[, seu$study == s], 
     features = gene_vec, 
     group.by = "final_celltype"
   )
   mtx.avg <- as.matrix(mtx.avg[[1]])
   mtx.avg <- mtx.avg[gene_vec,]
   mtx.scaled <- t(scale(t(mtx.avg)))
   ht.mtx.list[[s]] <- mtx.scaled
}



# functions ---------------------------------------------------------------

plot_simple_heatmap <- function(
    mat, 
    title = NULL, 
    show_row_names = TRUE,
    palette = "viridis"
) {
  library(ComplexHeatmap)
  library(circlize)
  library(viridisLite)
  library(grid)
  
  ## 颜色映射
  cols <- switch(
    palette,
    viridis = viridisLite::viridis(256),
    magma   = viridisLite::magma(256),
    rwb     = c("#B2182B", "#FFFFFF", "#2166AC"),
    stop('"palette" should be one of "viridis", "rwb", "magma" ')
  )
  
  col_fun <- circlize::colorRamp2(
    seq(min(mat, na.rm = TRUE), max(mat, na.rm = TRUE), length.out = length(cols)),
    cols
  )
  
  ## column annotation：直接使用列名作为分组
  group <- factor(colnames(mat))
  
  ha <- HeatmapAnnotation(
    group = group,
    col = list(
      group = structure(
        c('#788d66', '#885578', '#fad09f', '#ff8a9a', '#d157a0'),
        names = levels(group)
      )
    ),
    show_annotation_name = FALSE
  )
  
  ## 绘图
  Heatmap(
    mat,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = show_row_names, 
    col = col_fun,
    row_names_gp = grid::gpar(fontface = "italic"),
    top_annotation = ha,
    column_title = title,
    column_title_gp = grid::gpar(fontsize = 8)
  )
}


# plots -------------------------------------------------------------------

p.list <- lapply(names(ht.mtx.list), function(ds){
  plot_simple_heatmap(ht.mtx.list[[ds]], title=ds, show_row_names = F)
})

p.list2 <- lapply(names(ht.mtx.list), function(ds){
  plot_simple_heatmap(ht.mtx.list[[ds]], title=ds, show_row_names = T)
})

p <- p.list[[1]] + p.list[[2]] + p.list[[3]] + p.list[[4]] + p.list2[[5]]
p

export::graph2pdf(p, file="figures/lung/35_smc_corss_datasets.pdf", width = 15, height = 4)


