install.packages(c("tidyverse", "Seurat", "pak"))

if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("meichendong/SCDC")

library(omnideconv)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biobase")

# minimal installation
pak::pkg_install("omnideconv/omnideconv")