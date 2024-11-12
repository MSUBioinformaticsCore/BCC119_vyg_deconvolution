# Load required libraries
library(Seurat)
library(tidyverse)

# Capture command line arguments
args <- commandArgs(TRUE) 

ref.seurat.path = args[1]              # Path to Seurat reference data
annotation.col = args[2]               # Column with cell type annotations
pseudo.dir = args[3]                   # Directory for results output
num_datasets = as.numeric(args[4])     # Number of pseudobulk datasets to create
cells_per_sample = as.numeric(args[5]) # Number of cells per pseudobulk dataset

# ref.seurat.path = "/mnt/research/bioinformaticsCore/projects/Yuzbasiyan-Gurkan/BCC119_vyg_deconvolution/data/Ammons_scrna/can"
# annotation.col = "celltype.l1"
# pseudo.dir = "/mnt/research/bioinformaticsCore/projects/Yuzbasiyan-Gurkan/BCC119_vyg_deconvolution/data/Ammons_scrna/pseudobulk"

# Assume `seurat_obj` is your Seurat object with:
# - RNA counts in the "RNA" assay
# - Cell types labeled in the `CellType` metadata column


seurat_obj = readRDS(ref.seurat.path)

# Function to generate a pseudobulk sample with cell type proportions
generate_pseudobulk_with_proportions <- function(seurat_obj, num_cells) {
  
  # Randomly select `num_cells` from the Seurat object
  # sampled_cells <- sample(Cells(seurat_obj), num_cells, replace = FALSE)
  # Get all cell types and their cell counts in the Seurat object
  cell_types <- as.character(unique(seurat_obj@meta.data %>% pull(all_of(annotation.col))))
  cell_types = cell_types[!is.na(cell_types)]
  
  # Generate random proportions for each cell type, then scale to sum to 1
  random_proportions <- runif(length(cell_types))
  random_proportions <- random_proportions / sum(random_proportions)
  
  # Determine number of cells to sample from each cell type based on proportions
  cells_to_sample <- round(random_proportions * num_cells)
  
  # Sample cells from each cell type according to the generated proportions
  seurat_obj@meta.data$barcodes = Cells(seurat_obj)
  
  sampled_cells <- unlist(lapply(seq_along(cell_types), function(i) {
    cell_type <- cell_types[i]
    cell_type_cells <- 
      seurat_obj@meta.data %>%
      filter(.[[annotation.col]] == cell_type) %>%
      pull(barcodes)
    
    sample(cell_type_cells, size = cells_to_sample[i], replace = TRUE)
  }))
  
  # Subset the Seurat object to include only the sampled cells
  subset_obj <- subset(seurat_obj, cells = sampled_cells)
  
  # Sum counts across the sampled cells to create a pseudobulk sample
  pseudobulk_counts <- Matrix::rowSums(subset_obj[["RNA"]]@counts)
  
  # Calculate the proportion of each cell type in the sampled cells
  cell.type.annotations = subset_obj@meta.data %>% pull(all_of(annotation.col))
  cell_type_counts <- table(as.character(cell.type.annotations))
  cell_type_proportions <- cell_type_counts / num_cells
  
  # Return a list containing pseudobulk counts and cell type proportions
  return(list(pseudobulk_counts = pseudobulk_counts, cell_type_proportions = cell_type_proportions))
}

# Create a list to store pseudobulk datasets and cell type proportions
pseudobulk_data <- vector("list", num_datasets)

# Generate pseudobulk datasets and store cell type proportions

for (i in 1:num_datasets) {
  print(i)
  set.seed(i)
  pseudobulk_data[[i]] <- generate_pseudobulk_with_proportions(seurat_obj, cells_per_sample)
}

# Extract pseudobulk count matrices and cell type proportions
pseudobulk_counts_matrix <- do.call(cbind, lapply(pseudobulk_data, function(x) x$pseudobulk_counts))
colnames(pseudobulk_counts_matrix) <- paste0("Pseudobulk_", 1:num_datasets)

# Extract cell type proportions into a separate data frame
# fix missing values
cell_type_proportions_df <- lapply(pseudobulk_data, function(x) 
  as.data.frame((x$cell_type_proportions)) %>%
    colnames())

cell_type_proportions_list = list()

for(i in 1:num_datasets){
  
  cell_type_proportions_list[[i]] = as.data.frame((pseudobulk_data[[i]]$cell_type_proportions))
  colnames(cell_type_proportions_list[[i]]) = c("CellType", "Proportion")
  cell_type_proportions_list[[i]]$Sample = paste0("Pseudobulk_", i)
}

cell_type_proportions_df =
  do.call(rbind, cell_type_proportions_list) %>%
  pivot_wider(id_cols = Sample,
              names_from = CellType, 
              values_from = Proportion,
              values_fill = 0)

# Save results as RDS files or text objects if needed
saveRDS(pseudobulk_counts_matrix, file = paste0(pseudo.dir, "/pseudobulk_counts_matrix.rds"))

write.table(pseudobulk_counts_matrix, 
            file = paste0(pseudo.dir, "/pseudobulk_counts_matrix.txt"),
            sep = "/t",
            quote = F,
            row.names = TRUE)

saveRDS(cell_type_proportions_df, file = paste0(pseudo.dir, "/pseudobulk_real_cell_type_proportions.rds"))


