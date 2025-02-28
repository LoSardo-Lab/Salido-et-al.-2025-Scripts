#this code was adapted from the Seurat vignette on single-cell data analysis, publicly available here:https://satijalab.org/seurat/articles/pbmc3k_tutorial
#and from 

#load packages
library(Seurat)
library(readr)
library(dplyr)
library(plyr)
library(patchwork)
library(data.table)
library(ggplot2)
library(cowplot)
library(pheatmap)
library(tidyverse)
library(RCurl)

#load individual samples, make them seurat objects, edit 'orig.ident' variable to match the sample name, and add a 'genotype' variable
sample1 <- Read10X(data.dir = "filePath/sample1 directory")
sample1 <- CreateSeuratObject(counts = sample1, project = "myProject", min.cells = 3, min.features = 200)
sample1$orig.ident = rep("sample1",length(sample1$orig.ident))
sample1$genotype = rep("risk_int",length(sample1$orig.ident))

sample2 <- Read10X(data.dir = "filePath/sample2 directory")
sample2 <- CreateSeuratObject(counts = sample2, project = "myProject", min.cells = 3, min.features = 200)
sample2$orig.ident = rep("sample2",length(sample2$orig.ident))
sample2$genotype = rep("risk_ko",length(sample2$orig.ident))

sample3 <- Read10X(data.dir = "filePath/sample3 directory")
sample3 <- CreateSeuratObject(counts = sample3, project = "myProject", min.cells = 3, min.features = 200)
sample3$orig.ident = rep("sample3",length(sample3$orig.ident))
sample3$genotype = rep("nonrisk_int",length(sample3$orig.ident))

sample4 <- Read10X(data.dir = "filePath/sample4 directory")
sample4 <- CreateSeuratObject(counts = sample4, project = "myProject", min.cells = 3, min.features = 200)
sample4$orig.ident = rep("sample4",length(sample4$orig.ident))
sample4$genotype = rep("nonrisk_ko",length(sample4$orig.ident))

#combine samples into a single seurat object
allSamples = merge(sample1,y=c(sample2, sample3, sample4),project="myProject", add.cell.ids=c("sample1","sample2","sample3","sample4"))

#calculate percent mitochondrial RNA
allSamples[["percent.mt"]] <- PercentageFeatureSet(allSamples, pattern = "^MT-")

#Generate violin plots of number of genes per cell (nFeature), number of UMIs per cell (nCount), and percent mt RNA per cell (percent.mt)
#to get a sense of the data and where the filters should be set
VlnPlot(allSamples, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0, ncol = 3)

#filter out low UMI and high mt
cleanDat <- subset(allSamples, subset = nFeature_RNA > 2500 & nFeature_RNA < 8500 & percent.mt < 20 & nCount_RNA > 5000 & nCount_RNA < 120000)

#SCTransform the data
myData <- SCTransform(cleanDat, vars.to.regress = "percent.mt", verbose = FALSE)

#Principle component analysis
myData <- RunPCA(object = myData)

#plot variance by principle components to decide how may dimensions to include moving forward
ElbowPlot(myData, ndims = 50, reduction = "pca")

#cluster the data
myData <- FindNeighbors(object = myData, dims = 1:5)
myData <- FindClusters(object = myData, resolution = 0.4)

#Create a UMAP of the clustered data
myData<- RunUMAP(myData, dims = 1:5, reduction = "pca")

#save cleaned and clustered data
saveRDS(myData, "filePath/myData.rds")

#Plot the UMAP colored by cluster
DimPlot(myData , pt.size = 1.6, reduction = "umap", group.by = 'seurat_clusters',
        cols = c("darkorange2","gold2","darkolivegreen2","olivedrab","aquamarine1","deepskyblue1",
                 "blue3","mediumpurple2","mediumorchid2","plum1",
                 "firebrick2","honeydew3","burlywood3","azure4","cornsilk2","lightpink","lightcyan1"))

#Plot the UMAP colored by genotype
DimPlot(mySCT, pt.size = 1.6, reduction = "umap", group.by = 'genotype',
        cols = c("firebrick","gold2","turquoise4","darkolivegreen3"))

#Differential expression analysis between all clusters
allMarkers <- FindAllMarkers(myData, min.pct = 0.1)

#Differential expression analysis with all risk-associated clusters as one group compared to all other clusters as the second group
#numbers listed in the 'ident.1' and 'ident.2' portions of the command correspond to the risk-associated clusters and all other clusters
RRvAllElseMarkers <- FindMarkers(myData, ident.1 = c(3,4,6), ident.2 = c(0,1,7,8,9,2,5), min.pct = 0.1)

#Heatmap 
DoHeatmap(myData, features = geneList, group.by = "seurat_clusters", angle = 0) + scale_fill_distiller(palette = "PRGn") + theme(axis.text.y = element_text(size = 10))

