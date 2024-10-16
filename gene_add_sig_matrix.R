#adding genes corresponding to the ensemble ids in the signature matrix
## this script is for adding gene symbol in the signature matrix
setwd("/Users/teestanaskar/Dropbox/Teesta/Placenta/WashU.SCdata.humanplacenta")
library(DESeq2)
library(openxlsx)
library(dplyr)
library(apeglm)
library(biomaRt)
library(annotables)
install.packages("devtools")
devtools::install_github("stephenturner/annotables")
library(tidyverse)
library(org.Rn.eg.db)
library(EnsDb.Hsapiens.v86)
install.packages("EnsDb.Hsapiens.v86")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EnsDb.Hsapiens.v86")
#load signature matrix data file
signature_matrix = read.csv("signaturematrix.washu.csv")
rownames(signature_matrix) = signature_matrix[,1]
#method 1 by using biomaRt
listEnsembl()
useEnsembl(biomart = "genes")
ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)
ensembl.con <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl") #for rat
ensembl.con <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") # for human
attr <- listAttributes(ensembl.con)
filters <- listFilters(ensembl.con)

genes <- getBM(attributes = c('external_gene_name','ensembl_gene_id'),
               filters = "ensembl_gene_id",
               values = rownames(signature_matrix),
               mart = ensembl.con)
Added.genelist <- write.csv(genes, "genes.ensembleid.biomart.csv")
#rename colname of vst counts
names(signature_matrix)[1]<-paste('ensembl_gene_id')
addedgenes.counts <- left_join(signature_matrix,genes, by = "ensembl_gene_id")
#save
write.csv(addedgenes.counts, "signatureMatrix.with.gene.names.csv")
