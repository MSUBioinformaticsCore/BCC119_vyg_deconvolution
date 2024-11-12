library(tidyverse)
library(Seurat)

# Capture command line arguments
args <- commandArgs(TRUE)

# Define paths and parameters from arguments
ref.seurat.path = args[1]                  # Path to Seurat reference data
annotation.col = args[2]                   # Column with cell type annotations
remove.types = strsplit(args[3], ",")[[1]] # cell types to discard 

remove.types

# Load Seurat reference object (single-cell data)
seuratOb = readRDS(ref.seurat.path)
seuratOb

# filter seurat object for cell types of interest
cells_oi =
  seuratOb@meta.data %>%
  filter(!(!!sym(annotation.col) %in% remove.types)) 

seuratOb_oi <- subset(seuratOb, cells = rownames(cells_oi))
seuratOb_oi

new_file = gsub("^(.*)[.].*", "\\1", ref.seurat.path)
new_file = paste0(new_file, "_cells_oi.rds")
saveRDS(seuratOb_oi, file = new_file)
