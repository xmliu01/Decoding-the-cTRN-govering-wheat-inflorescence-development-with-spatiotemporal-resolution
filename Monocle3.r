library(Seurat)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(RColorBrewer)

#Partition for all cells
snrna=readRDS("snRNA.Seurat.rds")
data=GetAssayData(snrna,assay="RNA",slot="counts")
cell_metadata=snrna@meta.data
gene_annotation=data.frame(gene_short_name=rownames(data))
rownames(gene_annotation)=rownames(data)
cds=new_cell_data_set(data,cell_metadata=cell_metadata,gene_metadata=gene_annotation)
cds=preprocess_cds(cds,num_dim=100)
cds=reduce_dimension(cds)
cds=cluster_cells(cds)
cluster=as.data.frame(cds@clusters$UMAP$clusters)
partitions=as.data.frame(cds@clusters$UMAP$partitions)
write.table(partitions,"snrna.partition.txt",quote=F,sep="\t")

#Trajectory inference and pseudotime analysis for spikelet formation
par4=subset(snrna,subset=R7=="R7.P4" | seurat_clusters==5 | seurat_clusters==15 )
data=GetAssayData(par4,assay="RNA",slot="counts")
cell_metadata=par4@meta.data
gene_annotation=data.frame(gene_short_name=rownames(data))
rownames(gene_annotation)=rownames(data)
par4.cds=new_cell_data_set(data,cell_metadata=cell_metadata,gene_metadata=gene_annotation)
par4.cds=preprocess_cds(par4.cds,num_dim=50,method="PCA")
par4.cds=reduce_dimension(par4.cds)
par4.cds=cluster_cells(par4.cds)
par4.cds=learn_graph(par4.cds)
par4.cds=order_cells(par4.cds)
pseudotime.gene=graph_test(par4.cds,neighbor_graph="principal_graph",cores=24)
write.table(pr_deg_ids,"snrna.partition4.pseudotime.gene.txt",quote=F,sep="\t")

#Trajectory inference and pseudotime analysis for floret formation
par1=subset(snrna,subset=R7=="R7.P1" | seurat_clusters==10 | seurat_clusters==14 )
data=GetAssayData(par1,assay="RNA",slot="counts")
cell_metadata=par1@meta.data
gene_annotation=data.frame(gene_short_name=rownames(data))
rownames(gene_annotation)=rownames(data)
par1.cds=new_cell_data_set(data,cell_metadata=cell_metadata,gene_metadata=gene_annotation)
par1.cds=preprocess_cds(par1.cds,num_dim=50,method="PCA")
par1.cds=reduce_dimension(par1.cds)
par1.cds=cluster_cells(par1.cds)
par1.cds=learn_graph(par1.cds)
par1.cds=order_cells(par1.cds)
pseudotime.gene=graph_test(par1.cds,neighbor_graph="principal_graph",cores=24)
write.table(pr_deg_ids,"snrna.partition1.pseudotime.gene.txt",quote=F,sep="\t")
