# Load required libraries
library(Biobase)      # For creating ExpressionSet objects
library(SCDC)         # For single-cell deconvolution
library(Seurat)       # For handling Seurat single-cell objects
library(tidyverse)    # For data manipulation

# Capture command line arguments
args <- commandArgs(TRUE)

# Define paths and parameters from command-line arguments
ref.seurat.path = args[1]        # Path to Seurat reference data
bulk.mtx.path = args[2]          # Path to bulk RNA-seq data
tool = args[3]                   # Tool to be used (in this case, 'scdc')
annotation.col = args[4]         # Column with cell type annotations
batch.col = args[5]              # Column with batch information
results.dir = args[6]            # Directory for output results

# Load the bulk RNA-seq data
bulk.mtx = readRDS(bulk.mtx.path)

# Load the single-cell Seurat reference object
seuratOb = readRDS(ref.seurat.path)

# Convert bulk RNA-seq data to ExpressionSet ------------------------------

# `SCDC` requires an ExpressionSet for bulk data, which includes sample and feature metadata.
bulk.eset <- getESET(
  bulk.mtx,
  pdata = data.frame(Sample = colnames(bulk.mtx)),  # Sample metadata
  fdata = rownames(bulk.mtx)                        # Feature metadata (gene names)
)

# Convert reference Seurat object to ExpressionSet ------------------------

# Convert the single-cell RNA count data from Seurat to ExpressionSet
ref.eset <- getESET(
  seuratOb[["RNA"]]$counts,       # RNA count data from Seurat object
  pdata = seuratOb@meta.data,     # Metadata for cells
  fdata = rownames(seuratOb)      # Feature data (gene names)
)

# Reference Quality Control (QC) ------------------------------------------

# Extract cell type annotations from the metadata and identify unique cell types
cell.type.annotations <- seuratOb@meta.data %>% pull(all_of(annotation.col))
all_types = unique(cell.type.annotations)

# Perform QC on the reference data using SCDC's quality control function
ref.qc <- SCDC_qc(
  ref.eset,
  ct.varname = annotation.col,  # Column with cell type information
  sample = batch.col,           # Column with batch information
  qcthreshold = 0.5,            # QC threshold for filtering
  ct.sub = as.character(all_types) # List of cell types to include
)

# Save QC-processed reference data
save(ref.qc, file = paste0(results.dir, "/SCDC_QCdata.Rdata"))

# Estimate Cell Type Proportions ------------------------------------------

# Run SCDC to estimate cell type proportions in the bulk data
scdc.res = SCDC_prop(
  bulk.eset = bulk.eset,                # Bulk ExpressionSet
  sc.eset = ref.qc$sc.eset.qc,          # Quality-controlled reference ExpressionSet
  ct.varname = annotation.col,          # Cell type column in metadata
  sample = batch.col,                   # Batch information column
  ct.sub = as.character(all_types)      # Cell types to include in estimation
)

# Save the estimated cell type proportions
save(scdc.res, file = paste0(results.dir, "/SCDC_res.Rdata"))
