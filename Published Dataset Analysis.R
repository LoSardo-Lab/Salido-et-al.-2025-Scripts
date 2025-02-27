#load packages
library(Seurat)
library(readr)
library(dplyr)
library(patchwork)
library(data.table)
library(ggplot2)
library(cowplot)
library(pheatmap)
library(sctransform)

#read in data
publishedData = read.table("filepath/publishedData.txt.gz", header = T, row.names = 1)

#create a seurat object
publishedData = CreateSeuratObject(counts = publishedData, project = "publishedDataAnalysis", min.cells = 3, min.features = 200)

#calculate percent mitochondrial RNA
publishedData[["percent.mt"]] <- PercentageFeatureSet(publishedData, pattern = "^MT-")

#Generate violin plots of number of genes per cell (nFeature), number of UMIs per cell (nCount), and percent mt RNA per cell (percent.mt)
VlnPlot(publishedData, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0, ncol = 3)

#filter out low UMI and high mt
cleanDat <- subset(cleanDat, subset = nFeature_RNA > 1000 & nFeature_RNA < 4000 & percent.mt < 6 & nCount_RNA > 100 & nCount_RNA < 20000)

#SCTransform the data
cleanDat <- SCTransform(cleanDat, vars.to.regress = "percent.mt", verbose = FALSE)

#Principle component analysis
pubDat <- RunPCA(object = cleanDat)

#plot variance explained by each principle component
ElbowPlot(pubDat, ndims = 50, reduction = "pca")

#cluster the data
pubDat <- FindNeighbors(object = pubDat, dims = 1:5)
pubDat <- FindClusters(object = pubDat, resolution = 0.5)
pubDat<- RunUMAP(pubDat, dims = 1:5, reduction = "pca")

#save the cleaned and clustered data
saveRDS(pubDat, "filePath/pubDat.rds")

#generate featureplots of the data colored by gene markers of various cell types to identify which clusters are comprised of what type of cell
FeaturePlot(pubDat, pt.size = 0.9, features = c("gene1","gene2","gene3"))

#label the clusters according to the cell types they represent
new.cluster.ids <- c("Fibroblasts & SMC", "SMC", "SMC", "EC", "SMC", "Macrophage","T & NK","Macrophage","Macrophage","Macrophage","Mast, B, & plasma","Other")
pubDat <- RenameIdents(pubDat, new.cluster.ids)

