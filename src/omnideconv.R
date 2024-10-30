#https://github.com/omnideconv/omnideconv/blob/main/vignettes/omnideconv_example.Rmd

library(tidyverse)
library(Seurat)
library(omnideconv)

args <- commandArgs(TRUE)

ref.seurat.path = args[1] 
bulk.mtx.path = args[2]
tool = args[3]
annotation.col = args[4]
batch.col = args[5]
results.dir = args[6]

setwd(results.dir)
print(tool)

# Set up ------------------------------------------------------------------

# load bulk data
bulk.mtx = readRDS(bulk.mtx.path)

# load reference data
seuratOb = readRDS(ref.seurat.path)


# Subsample reference -----------------------------------------------------

max_cells_per_celltype = 200

sampled.metadata <- seuratOb@meta.data %>%
  rownames_to_column('barcode') %>%
  group_by_at(annotation.col) %>% 
  nest() %>%            
  mutate(n =  map_dbl(data, nrow)) %>%
  mutate(n = min(n, max_cells_per_celltype)) %>%
  ungroup() %>% 
  mutate(samp = map2(data, n, sample_n)) %>% 
  select(-data) %>%
  unnest(samp)

single.cell.data.sampled <- subset(seuratOb, cells = sampled.metadata$barcode)

# DWLS --------------------------------------------------------------------

if(tool == "dwls"){
  
  single_cell_object <- as.matrix(single.cell.data.sampled[["RNA"]]$counts)
  cell.type.annotations <- sampled.metadata %>% pull(all_of(annotation.col))
  batch.ids <- sampled.metadata %>% pull(all_of(batch.col))
  
  signature <- omnideconv::build_model(single_cell_object = single_cell_object,
                                       cell_type_annotations = cell.type.annotations,
                                       batch_ids = batch.ids, 
                                       method = "dwls", 
                                       dwls_method = 'mast_optimized')
  
  saveRDS(signature, file = paste0(tool, "_signature.rds"))
  
  deconvolution <- omnideconv::deconvolute(bulk_gene_expression = bulk.mtx, 
                                           signature,
                                           method='dwls',
                                           dwls_submethod = 'DampenedWLS')
  
  saveRDS(deconvolution, file = paste0(tool, "_deconvolution.rds"))
  
}


# scaden ------------------------------------------------------------------

if(tool == "scaden"){
  
  single_cell_object <- as.matrix(single.cell.data.sampled[["RNA"]]$counts)
  cell.type.annotations <- sampled.metadata %>% pull(all_of(annotation.col))
  batch.ids <- sampled.metadata %>% pull(all_of(batch.col))
  
  signature <- omnideconv::build_model(single_cell_object = single_cell_object,
                                       cell_type_annotations = cell.type.annotations,
                                       batch_ids = batch.ids, 
                                       method = "scaden", 
                                       bulk_gene_expression = bulk.mtx,
                                       verbose =TRUE)
  
  saveRDS(signature, file = paste0(tool, "_signature.rds"))
  
  deconvolution <- omnideconv::deconvolute(bulk_gene_expression = bulk.mtx, 
                                           signature,
                                           method = "scaden")
  
  saveRDS(deconvolution, file = paste0(tool, "_deconvolution.rds"))
  
}


# SCDC --------------------------------------------------------------------

if(tool == "scdc"){
  
  single_cell_object <- as.matrix(seuratOb[["RNA"]]$counts)
  cell.type.annotations <- seuratOb@meta.data %>% pull(all_of(annotation.col))
  batch.ids <- seuratOb@meta.data %>% pull(all_of(batch.col))
  

  deconvolution <- omnideconv::deconvolute(bulk_gene_expression = bulk.mtx,
                                           method = "scdc",
                                           single_cell_object = single_cell_object, 
                                           cell_type_annotations = cell.type.annotations, 
                                           batch_ids = batch.ids,
                                           quality_control = TRUE,
                                           ct_sub = colnames(seuratOb))
  
  saveRDS(deconvolution, file = paste0(tool, "_deconvolution_withQC.rds"))
  
  deconvolution <- omnideconv::deconvolute(bulk_gene_expression = bulk.mtx,
                                           method = "scdc",
                                           single_cell_object = single_cell_object, 
                                           cell_type_annotations = cell.type.annotations, 
                                           batch_ids = batch.ids,
                                           quality_control = FALSE)
  
  saveRDS(deconvolution, file = paste0(tool, "_deconvolution_noQC.rds"))
  
}

# music --------------------------------------------------------------------

if(tool == "music"){
  
  single_cell_object <- as.matrix(seuratOb[["RNA"]]$counts)
  cell.type.annotations <- seuratOb@meta.data %>% pull(all_of(annotation.col))
  batch.ids <- seuratOb@meta.data %>% pull(all_of(batch.col))
  
  deconvolution <- omnideconv::deconvolute(bulk_gene_expression = bulk.mtx,
                                           method = "music",
                                           single_cell_object = single_cell_object, 
                                           cell_type_annotations = cell.type.annotations, 
                                           batch_ids = batch.ids)
  
  saveRDS(deconvolution, file = paste0(tool, "_deconvolution.rds"))
  
}

