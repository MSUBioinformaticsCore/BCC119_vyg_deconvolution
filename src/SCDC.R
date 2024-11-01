## Deconvolution with SCDC

library(Biobase)
library(SCDC)
library(Seurat)
library(tidyverse)

args <- commandArgs(TRUE)

ref.seurat.path = args[1] 
bulk.mtx.path = args[2]
tool = args[3]
annotation.col = args[4]
batch.col = args[5]
results.dir = args[6]

# Set paths
data.dir = "/mnt/research/bioinformaticsCore/projects/Yuzbasiyan-Gurkan/BCC119_vyg_deconvolution/data"
results.dir = "/mnt/research/bioinformaticsCore/projects/Yuzbasiyan-Gurkan/BCC119_vyg_deconvolution/results"


# Load the bulk RNA-seq data
bulk.mtx = readRDS(bulk.mtx.path)

# Load the single-cell reference
seuratOb = readRDS(ref.seurat.path)

# `SCDC` takes `ExpressionSet` objects with raw read counts as input. 
# It needs an `ExpressionSet` for each reference individual.
# Convert the bulk data to an `ExpressionSet`

bulk.eset <- getESET(bulk.mtx,
                     pdata = data.frame(Sample = colnames(bulk.mtx)),
                     fdata = rownames(bulk.mtx))

# Convert reference  `SeuratObject` into an `ExpressionSet` object.

ref.eset <- getESET(seuratOb[["RNA"]]$counts,
                    pdata = seuratOb@meta.data,
                    fdata = rownames(seuratOb))

# Reference QC
cell.type.annotations <- seuratOb@meta.data %>% pull(all_of(annotation.col))
all_types = unique(cell.type.annotations)

ref.qc <- SCDC_qc(ref.eset,
                  ct.varname = annotation.col,
                  sample = batch.col,
                  qcthreshold = 0.5,
                  ct.sub = as.character(all_types))

save(ref.qc, file = paste0(results.dir, "/SCDC_QCdata.Rdata"))

### Estimate cell type proportions
scdc.res = SCDC_prop(
  bulk.eset = bulk.eset,
  sc.eset = ref.qc$sc.eset.qc,
  ct.varname = annotation.col,
  sample = batch.col,
  ct.sub = as.character(all_types))

save(scdc.res, file = paste0(results.dir, "/SCDC_res.Rdata"))
