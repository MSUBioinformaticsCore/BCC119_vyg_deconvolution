# Load necessary libraries
library(tidyverse)

# Set paths for data and results directories
data.dir = "/mnt/research/bioinformaticsCore/projects/Yuzbasiyan-Gurkan/BCC119_vyg_deconvolution/data"
results.dir = "/mnt/research/bioinformaticsCore/projects/Yuzbasiyan-Gurkan/BCC119_vyg_deconvolution/results"

# Load the bulk RNA-seq data
bulk = read.csv(paste0(data.dir, "/HS_MSU_RNAseq Analysis_rawCountsWithAnnotations.csv"))
bulk_cellline = read.delim(paste0(data.dir, "/HS_BD_OD_PJ_RNAseq_rawCountsWithAnnotations.txt"))

# Check the dimensions of the bulk RNA-seq data
dim(bulk)
dim(bulk_cellline)

# join the two datasets
bulk_all = left_join(bulk, bulk_cellline)

# Extract relevant columns for RNA-seq data and convert to matrix
bulk.mat = as.matrix(bulk_all[,c(12:38)])
rownames(bulk.mat) = bulk_all$ensembl_gene_id

# Load the gene lengths for normalization
gene_lengths = read.delim(paste0(data.dir, "/gene_lengths_Canis_lupus_familiaris.CanFam3.1.100.txt"))

# Filter the data to include only genes present in both the bulk data and gene length file
keep_ids = intersect(bulk$ensembl_gene_id, gene_lengths$gene)
bulk.filt = bulk.mat[keep_ids,]
bulk.filt = bulk.filt[order(rownames(bulk.filt)),]

# Filter and sort gene lengths to match filtered bulk data
keep_lengths = 
  gene_lengths %>%
  filter(gene %in% keep_ids) %>%
  arrange(gene) 

# Define a function to calculate Transcripts Per Million (TPM)
calculate_tpm <- function(count_matrix, gene_lengths) {
  # Step 1: Calculate Reads Per Kilobase (RPK)
  # Convert gene lengths from base pairs to kilobases
  gene_lengths_kb <- gene_lengths / 1000
  
  # Divide each count by gene length in kilobases
  rpk <- sweep(count_matrix, 1, gene_lengths_kb, "/")
  
  # Step 2: Calculate scaling factor (sum of RPKs per sample)
  scaling_factor <- colSums(rpk)
  
  # Step 3: Calculate TPM by normalizing RPKs
  tpm <- sweep(rpk, 2, scaling_factor, "/") * 1e6
  
  return(tpm)
}

# Calculate TPM matrix from filtered bulk data and gene lengths
tpm.mtx = calculate_tpm(bulk.filt, keep_lengths$merge)

# Handle duplicate gene names
# Identify rows with duplicated gene names
dups = bulk$external_gene_name[duplicated(bulk$external_gene_name)]

# Remove rows with duplicated gene names from the bulk data
bulk.nodups = 
  bulk %>%
  filter(!external_gene_name %in% dups)

# Retain only common genes between filtered data and TPM matrix
keep_ids = intersect(bulk.nodups$ensembl_gene_id, rownames(tpm.mtx))
tpm.nodups = tpm.mtx[keep_ids,]

# Map Ensembl IDs to external gene names and arrange by Ensembl ID
keep_gene_names = 
  bulk.nodups %>%
  filter(ensembl_gene_id %in% keep_ids) %>%
  select(ensembl_gene_id, external_gene_name) %>%
  arrange(ensembl_gene_id) 

# Confirm if order of IDs in TPM matrix matches external gene names
identical(keep_gene_names$ensembl_gene_id, rownames(tpm.nodups))

# Replace Ensembl IDs with external gene names in the TPM matrix
tpm.nodups.symbol = tpm.nodups
rownames(tpm.nodups.symbol) = keep_gene_names$external_gene_name

# Save the processed TPM matrix as an RDS file
saveRDS(tpm.nodups.symbol,  file = paste0(data.dir,"/bulk_tpm_mat.Rds"))

# save as text for scaden
write.table(tpm.nodups.symbol,
            file = paste0(data.dir,"/bulk_tpm_mat.txt"),
            sep = "\t",
            quote = FALSE,
            row.names = TRUE)

# filter pseudobulk genes for genes present in bulk
# don't tmp norm because 10x data is 3 prime biased
pseudobulk_counts_matrix = readRDS(paste0(data.dir,"/Ammons_scrna/pseudobulk/pseudobulk_counts_matrix.rds"))
keep_pseudobulk_genes = intersect(keep_gene_names$external_gene_name, rownames(pseudobulk_counts_matrix))

pseudobulk_filt_matrix = pseudobulk_counts_matrix[keep_pseudobulk_genes,]

# Save the matrix as an RDS file
saveRDS(pseudobulk_filt_matrix,  file = paste0(data.dir,"/Ammons_scrna/pseudobulk/pseudobulk_filtered_genes_counts_matrix.rds"))

# save as text for scaden
write.table(pseudobulk_filt_matrix,
            file = paste0(data.dir,"/Ammons_scrna/pseudobulk/pseudobulk_filtered_genes_counts_matrix.txt"),
            sep = "\t",
            quote = FALSE,
            row.names = TRUE)