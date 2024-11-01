library(tidyverse)
#Set paths
data.dir = "/mnt/research/bioinformaticsCore/projects/Yuzbasiyan-Gurkan/BCC119_vyg_deconvolution/data"
results.dir = "/mnt/research/bioinformaticsCore/projects/Yuzbasiyan-Gurkan/BCC119_vyg_deconvolution/results"

#Load the bulk RNA-seq data
bulk = read.csv(paste0(data.dir, "/HS_MSU_RNAseq Analysis_rawCountsWithAnnotations.csv"))
dim(bulk)

bulk.mat = as.matrix(bulk[,c(12:34)])
rownames(bulk.mat) = bulk$ensembl_gene_id

# Load the gene lengths
gene_lengths = read.delim(paste0(data.dir, "/gene_lengths_Canis_lupus_familiaris.CanFam3.1.100.txt"))

# Normalize counts
keep_ids = intersect(bulk$ensembl_gene_id, gene_lengths$gene)
bulk.filt = bulk.mat[keep_ids,]
bulk.filt = bulk.filt[order(rownames(bulk.filt)),]

keep_lengths = 
  gene_lengths %>%
  filter(gene %in% keep_ids) %>%
  arrange(gene) 

calculate_tpm <- function(count_matrix, gene_lengths) {
  # Step 1: Calculate Reads Per Kilobase (RPK)
  # Convert gene lengths from base pairs to kilobases for TPM calculation
  gene_lengths_kb <- gene_lengths / 1000
  
  # Divide each count by gene length in kilobases
  rpk <- sweep(count_matrix, 1, gene_lengths_kb, "/")
  
  # Step 2: Calculate the scaling factor (sum of RPKs per sample)
  scaling_factor <- colSums(rpk)
  
  # Step 3: Calculate TPM
  tpm <- sweep(rpk, 2, scaling_factor, "/") * 1e6
  
  return(tpm)
}

tpm.mtx = calculate_tpm(bulk.filt, keep_lengths$merge)

# change to external_gene_name 
# identify rows with duplicated external_gene_name
dups = bulk$external_gene_name[duplicated(bulk$external_gene_name)]

# remove rows with duplicated external_gene_name
bulk.nodups = 
  bulk %>%
  filter(!external_gene_name %in% dups)

keep_ids = intersect(bulk.nodups$ensembl_gene_id, rownames(tpm.mtx))
tpm.nodups = tpm.mtx[keep_ids,]

keep_gene_names = 
  bulk.nodups %>%
  filter(ensembl_gene_id %in% keep_ids) %>%
  select(ensembl_gene_id, external_gene_name) %>%
  arrange(ensembl_gene_id) 

identical(keep_gene_names$ensembl_gene_id, rownames(tpm.nodups))

tpm.nodups.symbol = tpm.nodups
rownames(tpm.nodups.symbol) = keep_gene_names$external_gene_name

saveRDS(tpm.nodups.symbol,  file = paste0(data.dir,"bulk_tpm_mat.Rds"))




