#snRNA-seq data analysis
library(tidyverse)
library(Seurat)
library(sctransform)

#Read the count matrix
count=Read10X("./snRNA.TestData/",gene.column=1)

#Creat Seurat object for each sample
snRNA=CreateSeuratObject(counts=count,project = "snRNA",min.cells=3,min.features=2000)

#Normalization
snRNA=SCTransform(snRNA)

#Dimension reduction
snRNA=RunPCA(snRNA,features = VariableFeatures(object = snRNA))
snRNA=RunUMAP(snRNA,dims=1:30)

#Clustering
snRNA=FindNeighbors(snRNA,dims=1:30)
snRNA=FindClusters(snRNA,resolution = 0.5)

#Identification of cell type-specific genes
marker=FindAllMarkers(snRNA,assay="SCT",only.pos=T,min.pct = 0.1,logfc.threshold=0.25)
write_tsv(marker,"test.snRNA.celltype.specific.genes.txt")

#Save results
saveRDS(snRNA,file="test.snRNA.Seurat.rds")
