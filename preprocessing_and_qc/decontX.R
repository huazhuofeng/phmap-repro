library(Seurat)
library(sceasy)
library(tidyverse)
library(decontX)

# data input
scobj <- sceasy::convertFormat(obj = 'data/00_merged_all.h5ad',from = 'anndata', to = 'seurat')
scobj

# get the counts matrix
counts = GetAssayData(scobj, layer = 'counts')

# decontX
decontX_results <- decontX(counts) 

df <- data.frame(decontX_contamination = decontX_results$contamination)

write.csv(df, file = 'data/02_human_decontX_contamination.csv', quote = F, row.names = F)

