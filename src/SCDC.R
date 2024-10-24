## Deconvolution with SCDC

library(Biobase)
library(SCDC)
library(Seurat)
library(tidyverse)

# Set paths
data.dir = "/mnt/research/bioinformaticsCore/projects/Yuzbasiyan-Gurkan/BCC119_vyg_deconvolution/data"
results.dir = "/mnt/research/bioinformaticsCore/projects/Yuzbasiyan-Gurkan/BCC119_vyg_deconvolution/results"


# Load the bulk RNA-seq data
bulk = read.csv(paste0(data.dir, "/HS_MSU_RNAseq Analysis_rawCountsWithAnnotations.csv"))


# Load the single-cell reference
seuratOb = readRDS(paste0(data.dir,"/Ammons_scrna/canine_naive_n6_annotated.rds"))

# identify rows with duplicated external_gene_name
dups = bulk$external_gene_name[duplicated(bulk$external_gene_name)]

# how many duplicated row names?  
length(dups)

# remove rows with duplicated external_gene_name
bulk.nodups = 
  bulk %>%
  filter(!external_gene_name %in% dups)

# select external_gene_name as rownames
bulk.expr = bulk.nodups[,c(12:34)]
rownames(bulk.expr) = bulk.nodups$external_gene_name

bulk.mtx = as.matrix(bulk.expr)


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
all_types = unique(seuratOb@meta.data$celltype.l1)

ref.qc <- SCDC_qc(ref.eset,
                  ct.varname = "celltype.l1",
                  sample = "orig.ident",
                  qcthreshold = 0.7,
                  ct.sub = as.character(all_types))

save(ref.qc, file = paste0(results.dir, "/SCDC_QCdata.Rdata"))

### Estimate cell type proportions
scdc.res = SCDC_prop(
  bulk.eset = bulk.eset,
  sc.eset = ref.qc$sc.eset.qc,
  ct.varname = "celltype.l1",
  sample = "orig.ident",
  ct.sub = as.character(all_types))

save(scdc.res, file = paste0(results.dir, "/SCDC_res.Rdata"))
