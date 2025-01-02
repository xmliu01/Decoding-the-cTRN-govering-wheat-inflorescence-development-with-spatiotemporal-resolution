library(tidyverse)
library(Seurat)
library(RColorBrewer)
library(sctransform)

#Read the count matrix
W2.5_1.count=Read10X("./data/W2.5_1/",gene.column=1)
W2.5_2.count=Read10X("./data/W2.5_2/",gene.column=1)
W3_1.count=Read10X("./data/W3_1/",gene.column=1)
W3_2.count=Read10X("./data/W3_2/",gene.column=1)
W3.25_1.count=Read10X("./data/W3.25_1/",gene.column = 1)
W3.25_2.count=Read10X("./data/W3.25_2/",gene.column = 1)
W3.5_1.count=Read10X("./data/W3.5_1/",gene.column = 1)
W3.5_2.count=Read10X("./data/W3.5_2/",gene.column = 1)
W4_1.count=Read10X("./data/W4_1/",gene.column = 1)
W4_2.count=Read10X("./data/W4_2/",gene.column = 1)
W5_1.count=Read10X("./data/W5_1/",gene.column = 1)
W5_2.count=Read10X("./data/W5_2/",gene.column = 1)

#Creat Seurat object for each sample
W2.5_1=CreateSeuratObject(counts=W2.5_1.count,project = "W2.5_1",min.cells=3,min.features=2000)
W2.5_2=CreateSeuratObject(counts=W2.5_2.count,project = "W2.5_2",min.cells=3,min.features=2000)
W3_1=CreateSeuratObject(counts=W3_1.count,project = "W3_1",min.cells=3,min.features=2000)
W3_2=CreateSeuratObject(counts=W3_2.count,project = "W3_2",min.cells=3,min.features=2000)
W3.25_1=CreateSeuratObject(counts=W3.25_1.count,project = "W3.25_1",min.cells=3,min.features=2000)
W3.25_2=CreateSeuratObject(counts=W3.25_2.count,project = "W3.25_2",min.cells=3,min.features=2000)
W3.5_1=CreateSeuratObject(counts=W3.5_1.count,project = "W3.5_1",min.cells=3,min.features=2000)
W3.5_2=CreateSeuratObject(counts=W3.5_2.count,project = "W3.5_2",min.cells=3,min.features=2000)
W4_1=CreateSeuratObject(counts=W4_1.count,project = "W4_1",min.cells=3,min.features=2000)
W4_2=CreateSeuratObject(counts=W4_2.count,project = "W4_2",min.cells=3,min.features=2000)
W5_1=CreateSeuratObject(counts=W5_1.count,project = "W5_1",min.cells=3,min.features=2000)
W5_2=CreateSeuratObject(counts=W5_2.count,project = "W5_2",min.cells=3,min.features=2000)

#Merge all object into one Seurat object
all=merge(W2.5_1,list(W2.5_2,W3_1,W3_2,W3.25_1,W3.25_2,W3.5_1,W3.5_2,W4_1,W4_2,W5_1,W5_2),add.cell.ids=c("W2.5_1","W2.5_2","W3_1","W3_2","W3.25_1","W3.25_2","W3.5_1","W3.5_2","W4_1","W4_2","W5_1","W5_2"),project="all")

#Normalization
all=SCTransform(all)

#Dimension reduction
all=RunPCA(all,features = VariableFeatures(object = all))
all=RunUMAP(all,dims=1:30)

#Clustering
all=FindNeighbors(all,dims=1:30)
all=FindClusters(all,resolution = 0.5)

#Identification of cell type-specific genes
marker=FindAllMarkers(all,assay="SCT",only.pos=T,min.pct = 0.1,logfc.threshold=0.25)
write_tsv(marker,"snRNA.celltype.specific.genes.txt")

#Save results
saveRDS(all,file="snRNA.Seurat.rds")