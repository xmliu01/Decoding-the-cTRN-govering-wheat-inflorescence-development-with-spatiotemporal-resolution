library(ArchR)
library(ggplot2)
library(tidyverse)
library(Seurat)
addArchRThreads(threads=24)
library("BSgenome.Triticum.aestivum.IWGSC")
genomeAnnotation=createGenomeAnnotation(genome=BSgenome.Triticum.aestivum.IWGSC,filter=F)
library("Triticum.aestivum.IWGSC.TxDb")
TSS=promoters(Triticum.aestivum.IWGSC.TxDb,upstream=0,downstream=0)
exons=exons(Triticum.aestivum.IWGSC.TxDb)
genes=genes(Triticum.aestivum.IWGSC.TxDb)
genes$symbol=genes$gene_id
geneAnnotation=createGeneAnnotation(TSS=TSS,exons=exons,genes=genes)

#Read input data
ArrowFiles=createArrowFiles(inputFiles="test.snATAC.fragment.tsv.gz",sampleNames="snATAC",geneAnnotation=geneAnnotation,genomeAnnotation=genomeAnnotation,minTSS=4,minFrags=500,maxFrags=1e+5,promoterRegion=c(3000,500),addGeneScoreMat=T,addTileMat=T)

#Create ArchR object
projSpike1=ArchRProject(ArrowFiles=ArrowFiles,outputDirectory="./res",copyArrows=F,geneAnnotation=geneAnnotation,genomeAnnotation=genomeAnnotation)

#Dimension reduction
projSpike1=addIterativeLSI(ArchRProj=projSpike1,useMatrix="TileMatrix",name="IterativeLSI",iterations=2,clusterParams=list(resolution=c(0.2),sampleCells=10000,n.start=10),varFeatures=25000,dimsToUse=1:30)

#Clustering
projSpike1=addClusters(input=projSpike1,reducedDims="IterativeLSI",method="Seurat",name="Clusters",resolution=0.6)
projSpike1=addUMAP(ArchRProj=projSpike1,reducedDims="IterativeLSI",name="UMAP",nNeighbors=30,minDist=0.5,metric="cosine",force=T)

#Save results
saveArchRProject(ArchRProj=projSpike1,outputDirectory="test.snATAC.output",load=F)


