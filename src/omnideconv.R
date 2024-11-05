#https://github.com/omnideconv/omnideconv/blob/main/vignettes/omnideconv_example.Rmd
# Load necessary libraries
library(tidyverse)
library(Seurat)
library(omnideconv)

# Capture command line arguments
args <- commandArgs(TRUE)

# Define paths and parameters from arguments
ref.seurat.path = args[1]            # Path to Seurat reference data
bulk.mtx.path = args[2]              # Path to bulk matrix data
tool = args[3]                       # Deconvolution tool (e.g., 'dwls', 'music')
annotation.col = args[4]             # Column with cell type annotations
batch.col = args[5]                  # Column with batch information
results.dir = args[6]                # Directory for results output

# Set working directory to results directory
setwd(results.dir)
print(tool)  # Print the tool being used

# Set up data -------------------------------------------------------------

# Load bulk gene expression data
bulk.mtx = readRDS(bulk.mtx.path)

# Load Seurat reference object (single-cell data)
seuratOb = readRDS(ref.seurat.path)

# # Subsample reference data ------------------------------------------------
# set.seed(1)  # Set random seed for reproducibility
# 
# max_cells_per_celltype = 200  # Maximum cells to sample per cell type
# 
# # Sample cells to a maximum per cell type using specified annotation column
# sampled.metadata <- seuratOb@meta.data %>%
#   rownames_to_column('barcode') %>%
#   group_by_at(annotation.col) %>% 
#   nest() %>%            
#   mutate(n = map_dbl(data, nrow)) %>%
#   mutate(n = min(n, max_cells_per_celltype)) %>%
#   ungroup() %>% 
#   mutate(samp = map2(data, n, sample_n)) %>% 
#   select(-data) %>%
#   unnest(samp)
# 
# # Subset Seurat object to include only sampled cells
# single.cell.data.sampled <- subset(seuratOb, cells = sampled.metadata$barcode)

# DWLS -------------------------------------------------------------

if (tool == "dwls") {
  
  # Extract RNA count matrix and metadata for cell type and batch
  single_cell_object <- as.matrix(seuratOb[["RNA"]]$counts)
  cell.type.annotations <- seuratOb@meta.data %>% pull(all_of(annotation.col))
  batch.ids <- seuratOb@meta.data %>% pull(all_of(batch.col))
  
  # Build DWLS model using omnideconv
  signature <- omnideconv::build_model(single_cell_object = single_cell_object,
                                       cell_type_annotations = cell.type.annotations,
                                       batch_ids = batch.ids, 
                                       method = "dwls", 
                                       dwls_method = 'mast_optimized')
  
  # Save model signature
  saveRDS(signature, file = paste0(tool, "_signature.rds"))
  
  # Perform deconvolution using DWLS and save results
  deconvolution <- omnideconv::deconvolute(bulk_gene_expression = bulk.mtx, 
                                           signature,
                                           method = 'dwls',
                                           dwls_submethod = 'DampenedWLS')
  
  saveRDS(deconvolution, file = paste0(tool, "_deconvolution.rds"))
  
}

# scaden ------------------------------------------------------------------
# can't get numpy to load properly
# if(tool == "scaden"){
#   
#   single_cell_object <- as.matrix(single.cell.data.sampled[["RNA"]]$counts)
#   cell.type.annotations <- sampled.metadata %>% pull(all_of(annotation.col))
#   batch.ids <- sampled.metadata %>% pull(all_of(batch.col))
#   
#   signature <- omnideconv::build_model(single_cell_object = single_cell_object,
#                                        cell_type_annotations = cell.type.annotations,
#                                        batch_ids = batch.ids, 
#                                        method = "scaden", 
#                                        bulk_gene_expression = bulk.mtx,
#                                        verbose =TRUE)
#   
#   saveRDS(signature, file = paste0(tool, "_signature.rds"))
#   
#   deconvolution <- omnideconv::deconvolute(bulk_gene_expression = bulk.mtx, 
#                                            signature,
#                                            method = "scaden")
#   
#   saveRDS(deconvolution, file = paste0(tool, "_deconvolution.rds"))
#   
# }


# SCDC --------------------------------------------------------------------
# issue with ct_sub vs ct.sub argument
#
# if(tool == "scdc"){
#   
#   single_cell_object <- as.matrix(seuratOb[["RNA"]]$counts)
#   cell.type.annotations <- seuratOb@meta.data %>% pull(all_of(annotation.col))
#   batch.ids <- seuratOb@meta.data %>% pull(all_of(batch.col))
#   
# 
#   deconvolution <- omnideconv::deconvolute(bulk_gene_expression = bulk.mtx,
#                                            method = "scdc",
#                                            single_cell_object = single_cell_object, 
#                                            cell_type_annotations = cell.type.annotations, 
#                                            batch_ids = batch.ids,
#                                            quality_control = TRUE,
#                                            ct_sub = NULL)
#   
#   saveRDS(deconvolution, file = paste0(tool, "_deconvolution_withQC.rds"))
#   
#   deconvolution <- omnideconv::deconvolute(bulk_gene_expression = bulk.mtx,
#                                            method = "scdc",
#                                            single_cell_object = single_cell_object, 
#                                            cell_type_annotations = cell.type.annotations, 
#                                            batch_ids = batch.ids,
#                                            quality_control = FALSE,
#                                            ct_sub = NULL)
#   
#   saveRDS(deconvolution, file = paste0(tool, "_deconvolution_noQC.rds"))
# 
#}

# music --------------------------------------------------------------------

if (tool == "music") {
  
  # Extract RNA count matrix and metadata for cell type and batch
  single_cell_object <- as.matrix(seuratOb[["RNA"]]$counts)
  cell.type.annotations <- seuratOb@meta.data %>% pull(all_of(annotation.col))
  batch.ids <- seuratOb@meta.data %>% pull(all_of(batch.col))
  
  # Perform deconvolution using MUSIC and save results
  deconvolution <- omnideconv::deconvolute(bulk_gene_expression = bulk.mtx,
                                           method = "music",
                                           single_cell_object = single_cell_object, 
                                           cell_type_annotations = cell.type.annotations, 
                                           batch_ids = batch.ids)
  
  saveRDS(deconvolution, file = paste0(tool, "_deconvolution.rds"))
  
}
