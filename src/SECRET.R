## Deconvolution with SCDC

library(presto)
library(metap)
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

# SECRET requires a list of cell-type marker genes

Idents(object = seuratOb) <- "celltype.l1"

markers_list = list()  

for (ct in as.character(all_types)){
  
  print(ct)
  markers_list[[ct]] = 
    FindConservedMarkers(seuratOb, 
                         ident.1 = ct,
                         min.pct = 0.25,
                         only.pos = TRUE,
                         grouping.var = "orig.ident")
  
  markers_list[[ct]]$CellType = ct
}

save(markers_list, file = paste0(results.dir, "/ref_markers.Rdata"))

markers_list = lapply(markers_list, function(df){
  df$Gene = rownames(df);
  df =  
    df %>%
    select(Gene, CellType, minimump_p_val) %>%
    filter(minimump_p_val < .01)
})

markers_df = do.call(rbind, markers_list)

# count for marker gene frequency
w=as.data.frame(table(markers_df$Gene))
colnames(w)[1]='genes'
rownames(w)=w$genes

# make a matrix from the seurat object
refdat = as.matrix(seuratOb[["RNA"]]$counts)

### find the shared genes and subset the data
glist=Reduce(intersect,list(rownames(bulk.mtx),
                            rownames(refdat),
                            (markers_df$Gene)))

w=w[glist,]

bulkdat2<-bulk.mtx[as.character(w$genes),]
bulkdat2 = bulkdat2[order(row.names(bulkdat2)),]

refdat2<-refdat[as.character(w$genes),]
refdat2 = refdat2[order(row.names(refdat2)),]

print(identical(rownames(bulkdat2),rownames(refdat2)))
#> [1] TRUE
w = w[order(rownames(w)),]
print(identical(rownames(bulkdat2),rownames(w)))
#> [1] TRUE

w$weight=1/w$Freq
### Estimate cell type propotion using SECRET

library(SECRET)

wt=w$weight
secret.res=SECRET(bulkdat2,
                  refdat2,
                  withUnknown = T,
                  w=wt,
                  yNorm = 'cpm',
                  bNorm='cpm')[[1]]

save(secret.res, file = paste0(results.dir, "/SECRET_res.Rdata"))
