# Estimate immune cell proportions in tumor samples
## Stephanie Hickey, Ph.D.

# Project Summary

Here we use single-cell RNA-seq data from primary osteosarcoma in 6 treatment-naïve dogs ([Ammons et al, Commun Biol. 2024 Apr 24](https://www.nature.com/articles/s42003-024-06182-w) as a reference data set to estimate immune cell proportions in bulk tumor samples using parallel methods [MuSiC](https://xuranw.github.io/MuSiC/articles/MuSiC.html), [SCDC](https://meichendong.github.io/SCDC/articles/SCDC.html#scdc-pre-process-of-scrna-seq-data-1), and [DWLS](https://cran.r-project.org/web/packages/DWLS/index.html).

# Set up

## Make project directories
```{bash, eval=FALSE}
PROJECT=/mnt/research/bioinformaticsCore/projects/Yuzbasiyan-Gurkan/BCC119_vyg_deconvolution
mkdir $PROJECT/src
mkdir $PROJECT/data
mkdir $PROJECT/results
```

## Download scRNA-seq reference 
```{bash, eval=FALSE}
mkdir $PROJECT/data/Ammons_scrna
cd $PROJECT/data/Ammons_scrna
wget https://zenodo.org/records/10891255/files/canine_naive_n6_annotated.rds
```

## Install R packages

Install the R packages required for this analysis:

* [pak](https://pak.r-lib.org)
* [BiocManager](https://www.bioconductor.org/install/)
* [tidyverse](https://www.tidyverse.org)
* [Biobase](https://www.bioconductor.org/packages/release/bioc/html/Biobase.html)
* [Seurat](https://satijalab.org/seurat/articles/install_v5)
* [SCDC](https://meichendong.github.io/SCDC/articles/SCDC.html)
* [omnideconv](https://github.com/omnideconv/omnideconv)

```{bash}
module purge
module load R-bundle-CRAN/2023.12-foss-2023a

Rscript $PROJECT/src/packages.R
```

## Get gene lengths

Using [gfttools](http://www.genemine.org/gtftools.php). 
```{bash}
module purge
module load Conda/3

pip install gtftools

cd $PROJECT/data

gtftools -l  gene_lengths_Canis_lupus_familiaris.CanFam3.1.100.txt Canis_lupus_familiaris.CanFam3.1.100.gtf
```

## Make random pseuodbulk samples

Creates n pseudobulk samples by randomly selecting m cells from each of the listed cell types in varying random proportions. 

**script:** `random_pseudobulk.R`

**arguments:**

1. Path to Seurat reference data    
2. Seurat meta data column with cell type annotations   
3. Directory for results output   
4. Number of pseudobulk datasets to create (n)    
5. Number of cells to sample from the single cell data set to make each pseuobulk sample (m)   

**output**

* `pseudobulk_counts_matrix.rds`: an `rds` file with a genes by n data sets matrix with summed counts from m cells   
* `pseudobulk_counts_matrix.txt`: a `txt` file with a genes by n data sets matrix with summed counts from m cells   
* `pseudobulk_real_cell_type_proportions.rds`:  the ground truth cell type proportions for each pseudobulk sample     

```{bash}
mkdir $PROJECT/data/Ammons_scrna/pseudobulk 

module purge
module load R-bundle-CRAN/2023.12-foss-2023a

Rscript $PROJECT/src/random_pseudobulk.R \
  $PROJECT/data/Ammons_scrna/canine_naive_n6_annotated.rds \
  celltype.l1 \
  $PROJECT/data/Ammons_scrna/pseudobulk \
  100 \
  10000

```

## Normalize bulk 

Uses the gene length file from gtftools to calculate tpm from the raw count data. Includes only genes with ensembl IDs present in both files. Gene symbols are used in the final matrix to match the single-cell reference data. Duplicated gene symbols have been removed. 

Filters the pseudobulk matrix for genes present in the bulk matrix. Length normalization was not performed on the pseudobulk matrix because it was generated using cells processed with a Chromium Next GEM Single Cell 3ʹ v3.1 kit.

**script:** `format_bulk.R`

**arguments:** None. Hard coded to be specific for this project.

**output:** 

* `$PROJECT/data/bulk_tpm_mat.Rds`    
* `$PROJECT/data/bulk_tpm_mat.txt`    
* `$PROJECT/data/Ammons_scrna/pseudobulk/pseudobulk_filtered_genes_counts_matrix.rds`   
* `$PROJECT/data/Ammons_scrna/pseudobulk/pseudobulk_filtered_genes_counts_matrix.txt`   

```{bash}
module purge
module load R-bundle-CRAN/2023.12-foss-2023a

Rscript $PROJECT/src/format_bulk.R
```

# Deconvolution analysis

## Deconvolution with DWLS and MuSiC using onmideconv

*onmideconv* is an R package for unified access to computational methods for estimating immune cell fractions from bulk RNA sequencing data. I was only able to run *DWLS* and *MuSiC* on hpcc with this tool.

**script:** `$PROJECT/src/submit_omnideconv.sh`

**arguments:**

1. Path to reference data seurat object   
2. Path to bulk data to be deconvoluted    
3. Tool to use    
4. Seurat meta data column with cell type annotations   
5. Seurat meta data column with sample/batch annotations    
6. Results directory    

**output:**

* `<tool>_deconvolution.rds` containing the predicted cell type proportions.

```{bash}
cd $PROJECT/run

for tool in dwls #music  
do

  echo $tool
  
  sbatch $PROJECT/src/submit_omnideconv.sh \
    $PROJECT/src/omnideconv.R \
    $PROJECT/data/Ammons_scrna/canine_naive_n6_annotated.rds \
    $PROJECT/data/bulk_tpm_mat.Rds \
    $tool \
    celltype.l1 \
    orig.ident \
    $PROJECT/results
  
done

# pseudobulk
mkdir $PROJECT/results/pseudobulk

for tool in dwls #music  
do

  echo $tool
  
  sbatch $PROJECT/src/submit_omnideconv.sh \
    $PROJECT/src/omnideconv.R \
    $PROJECT/data/Ammons_scrna/canine_naive_n6_annotated.rds \
    $PROJECT/data/Ammons_scrna/pseudobulk/pseudobulk_filtered_genes_counts_matrix.rds\
    $tool \
    celltype.l1 \
    orig.ident \
    $PROJECT/results/pseudobulk
  
done
  
```

## Deconvolution with SCDC

**script:** `$PROJECT/src/SCDC.R`

**arguments:**

1. Path to reference data seurat object   
2. Path to bulk data to be deconvoluted    
3. Tool to use    
4. Seurat meta data column with cell type annotations   
5. Seurat meta data column with sample/batch annotations    
6. Results directory    

**output:**

* `SCDC_QCdata.Rdata`: a list including: 1) a probability matrix for each single cell input; 2) a clustering QCed ExpressionSet object; 3) a heatmap of QC result.    
* `SCDC_res.Rdata`: a list including: Estimated proportions, basis matrix, predicted gene expression levels for bulk samples    

```{bash}

sbatch $PROJECT/src/submitSCDC.sh \
  $PROJECT/src/SCDC.R \
  $PROJECT/data/Ammons_scrna/canine_naive_n6_annotated.rds \
  $PROJECT/data/bulk_tpm_mat.Rds \
  SCDC \
  celltype.l1 \
  orig.ident \
  $PROJECT/results

# pseudobulk
sbatch $PROJECT/src/submitSCDC.sh \
  $PROJECT/src/SCDC.R \
  $PROJECT/data/Ammons_scrna/canine_naive_n6_annotated.rds \
  $PROJECT/data/Ammons_scrna/pseudobulk/pseudobulk_filtered_genes_counts_matrix.rds \
  SCDC \
  celltype.l1 \
  orig.ident \
  $PROJECT/results/pseudobulk
    
```