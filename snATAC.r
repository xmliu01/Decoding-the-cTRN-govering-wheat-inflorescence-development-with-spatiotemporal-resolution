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
ArrowFiles=createArrowFiles(inputFiles=c("./data/W2.5_1.fragment.tsv.gz","./data/W2.5_2.fragment.tsv.gz","./data/W3_1.fragment.tsv.gz","./data/W3_2.fragment.tsv.gz","./data/W3_3.fragment.tsv.gz","./data/W3_4.fragment.tsv.gz","./data/W3.25_1.fragment.tsv.gz","./data/W3.25_2.fragment.tsv.gz","./data/W3.25_3.fragment.tsv.gz","./data/W3.25_4.fragment.tsv.gz","./data/W3.5_1.fragment.tsv.gz","./data/W3.5_2.fragment.tsv.gz","./data/W3.5_3.fragment.tsv.gz","./data/W3.5_4.fragment.tsv.gz","./data/W4_1.fragment.tsv.gz","./data/W4_2.fragment.tsv.gz","./data/W5_1.fragment.tsv.gz","./data/W5_2.fragment.tsv.gz"),sampleNames=c("W2.5_1","W2.5_2","W3_1","W3_2","W3_3","W3_4","W3.25_1","W3.25_2","W3.25_3","W3.25_4","W3.5_1","W3.5_2","W3.5_3","W3.5_4","W4_1","W4_2","W5_1","W5_2"),geneAnnotation=geneAnnotation,genomeAnnotation=genomeAnnotation,minTSS=5,minFrags=3000,maxFrags=1e+5,promoterRegion=c(3000,500),addGeneScoreMat=T,addTileMat=T)

#Create ArchR object
doubScores=addDoubletScores(input=ArrowFiles,k=10,knnMethod="UMAP",LSIMethod = 1)
projSpike1=ArchRProject(ArrowFiles=ArrowFiles,outputDirectory="./res",copyArrows=F,geneAnnotation=geneAnnotation,genomeAnnotation=genomeAnnotation)
projSpike1=filterDoublets(projSpike1,filterRatio=1)

#Dimension reduction
projSpike1=addIterativeLSI(ArchRProj=projSpike1,useMatrix="TileMatrix",name="IterativeLSI",iterations=2,clusterParams=list(resolution=c(0.2),sampleCells=10000,n.start=10),varFeatures=25000,dimsToUse=1:30)

#Clustering
projSpike1=addClusters(input=projSpike1,reducedDims="IterativeLSI",method="Seurat",name="Clusters",resolution=0.6)
projSpike1=addUMAP(ArchRProj=projSpike1,reducedDims="IterativeLSI",name="UMAP",nNeighbors=30,minDist=0.5,metric="cosine",force=T)

#Peak calling
projSpike1=addGroupCoverages(ArchRProj=projSpike1,groupBy="Clusters")
pathToMacs2=findMacs2()
projSpike1=addReproduciblePeakSet(ArchRProj=projSpike1,groupBy="Clusters",pathToMacs2=pathToMacs2,genomeSize=14600000000)
projSpike1=addPeakMatrix(projSpike1)

#Identification of cell type-specific genes
markersGenes=getMarkerFeatures(ArchRProj=projSpike1,useMatrix="GeneScoreMatrix",groupBy="Clusters",bias=c("TSSEnrichment","log10(nFrags)"),testMethod="wilcoxon",maxCells=2000)
markerList=getMarkers(markersGenes,cutOff="FDR <= 0.05 & Log2FC >= 1")

#Identification of cell type-specific peaks
markersPeaks=getMarkerFeatures(ArchRProj=projSpike1,useMatrix="PeakMatrix",groupBy="Clusters",bias=c("TSSEnrichment","log10(nFrags)"),testMethod="wilcoxon",maxCells=2000)
markerList=getMarkers(markersPeaks,cutOff="FDR <= 0.05 & Log2FC >= 1")

#Save results
saveArchRProject(ArchRProj=projSpike1,outputDirectory="Save-all",load=F)



