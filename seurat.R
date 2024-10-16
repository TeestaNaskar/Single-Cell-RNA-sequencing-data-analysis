### this script is for reading and running the seurat object from placental single seq data
###from WashU
setwd("/Users/teestanaskar/Dropbox/Teesta/Placenta/WashU.SCdata.humanplacenta")
seuratObj = readRDS("sc.NormByLocationRep.Harmony.final.rds")

library(scCustomize)
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(patchwork)
library(readxl)
library(openxlsx)
library(tibble)
library(magrittr)
library(dplyr)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("glmGamPoi")
library(glmGamPoi)
ens <- mapIds(EnsDb.Hsapiens.v79, keys =reference@assays$RNA@data@Dimnames[[1]], column ='SYMBOL', keytype ='GENEID')
###the rownames are in ensemble so, we need to transform them
require(EnsDb.Hsapiens.v79)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EnsDb.Hsapiens.v79")
BiocManager::install("EnsDb.Hsapiens.v87")
library(EnsDb.Hsapiens.v86)
keys = seuratObj@assays$RNA@data@Dimnames[[1]]
keys = rownames(seuratObj)
keys = sub("\\..*$", "", keys)
ens <- mapIds(EnsDb.Hsapiens.v86, keys = keys, column ='SYMBOL', keytype ='GENEID')
length(ens)==length(keys) #it shoud Return 'TRUE'
keep <- !is.na(ens) 
ens <- ens[keep]

seuratObj <- seuratObj[keep,] 
seuratObj@assays$RNA@data@Dimnames[[1]]=ens 
names(seuratObj@assays$RNA@data@Dimnames[[1]])=c()
rownames(seuratObj) <- ens
all(rownames(seuratObj) == names(ens))

ElbowPlot(object = seuratObj, 
          ndims = 40)

#QC step
# Plot the elbow plot
ElbowPlot(object = seuratObj, 
          ndims = 50)

DimPlot(seuratObj, label = TRUE)
seuratObj = FindClusters(seuratObj, resolution = 0.8)
DimPlot(seuratObj, label = TRUE)
#finding number of cluster that has majority cell types
seuratObj = FindClusters(seuratObj, resolution = 0.2)
#its giving 16 cell types
metadata = seuratObj@meta.data
write.csv(metadata, "metadata.csv")
##save count data 
count = seuratObj@assays$RNA@counts 
count.frame = count %>% as.matrix %>% t %>% as.data.frame
write.csv(count, "count.csv")
write.table(seuratObj@assays[["RNA"]]@counts, file='count.tsv', quote=FALSE, sep='\t', col.names = TRUE)


#subsetting
subset.seuratobject <- subset(seuratObj, subset = "Old_Sample_ID" == c("C_HUZ_718148_F", "C_HUZ_718148_M", "HUZ-718134-2407", "HUZ-718134-2407"))
#getting gene name corresponding to ensemble ids
value = rownames(data)
value = sub("\\..*$", "", value)
ensembl <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol", "chromosome_name"),
                  #filters = "ensembl_gene_id",   # with filters it does not work
                  values = value,mart = ensembl)
idx <- match(value, genemap$ensembl_gene_id)
data$entrez <- genemap$entrezgene_id[idx]
data$chr <- genemap$chromosome_name[idx]
data$gene_name <- genemap$hgnc_symbol[idx]
