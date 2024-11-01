# process single cell for scaden
# https://scaden.readthedocs.io/en/latest/usage.html#file-formats

library(Seurat)
library(tidyverse)

args <- commandArgs(TRUE)

# Define paths and parameters from arguments
ref.seurat.path = args[1]            # Path to Seurat reference data
annotation.col = args[2]             # Column with cell type annotations
results.dir = args[3]                # Directory for results output

# Subsample reference data ------------------------------------------------
set.seed(1)  # Set random seed for reproducibility
seuratOb = readRDS(ref.seurat.path)

max_cells_per_celltype = 200  # Maximum cells to sample per cell type

# Sample cells to a maximum per cell type using specified annotation column
sampled.metadata <- seuratOb@meta.data %>%
  rownames_to_column('barcode') %>%
  group_by_at(annotation.col) %>% 
  nest() %>%            
  mutate(n = map_dbl(data, nrow)) %>%
  mutate(n = min(n, max_cells_per_celltype)) %>%
  ungroup() %>% 
  mutate(samp = map2(data, n, sample_n)) %>% 
  select(-data) %>%
  unnest(samp)

# Subset Seurat object to include only sampled cells
single.cell.data.sampled <- subset(seuratOb, cells = sampled.metadata$barcode)

# normalize your count data to library size
single.cell.data.sampled = NormalizeData(single.cell.data.sampled, 
                                         normalization.method = "RC",
                                         scale.factor = 1e6)

cpm = t(single.cell.data.sampled[["RNA"]]$data)
cpm.mat = as.matrix(cpm)

write.table(cpm.mat, 
            file = paste0(dirname(ref.seurat.path), "/subsampled_cpm.txt"),
            row.names = FALSE,
            sep = "\t",
            quote= F)

cell_types = 
  single.cell.data.sampled@meta.data %>%
  select(all_of(annotation.col))

colnames(cell_types) = "Celltype"


write.table(cell_types, 
            file = paste0(dirname(ref.seurat.path), "/subsampled_celltypes.txt"),
            row.names = FALSE,
            sep = "\t",
            quote= F)

