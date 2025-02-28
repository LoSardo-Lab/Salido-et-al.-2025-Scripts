#this code was adapted from the UC Davis Bioinformatics core tutorial, publicly available here: https://ucdavis-bioinformatics-training.github.io/2021-August-Advanced-Topics-in-Single-Cell-RNA-Seq-Trajectory-and-Velocity/data_analysis/monocle_fixed

#load packages
library(monocle3)
library(monocle)
library(Seurat)
library(SeuratWrappers)
library(patchwork)
library(dplyr)
library(Signac)
library(Matrix)
library(ggplot2)
library(pheatmap)

#Read in clustered seurat object
myData <- readRDS("filePath/myData.rds")

#subset by genotype
risk_int <- subset(myData, genotype == 'risk_int')

#make a monocle object out of the seurat object
cds <- as.cell_data_set(risk_int)

#cluster cells
cds <- cluster_cells(cds, resolution=1e-3)

#plot clusters and compare with whole
p1 <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)

#trajectory calculation
cds <- learn_graph(cds, use_partition = TRUE, verbose = FALSE)

#plot the trajectory
plot_cells(cds,
           color_cells_by = "cluster",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

#choose a starting point for the trajectory
cds <- order_cells(cds)

#plot pseudotime colored by trajectory
plot_cells(cds,
           cell_size = 1,
           color_cells_by = "pseudotime",
           group_cells_by = "cluster",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           trajectory_graph_color = "grey60")

#find genes that vary over the trajectory
myTrajectoryGenes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 8)

#set the gene names as rownames in the graph test results
myTrajectoryGenes$gene_name <- rownames(myTrajectoryGenes)

#filter by corrected p-value and arrange in ascending order
significant_genes_df = cds_graph_test_results %>% filter(q_value <= 0.05) %>% arrange(q_value)

#create gene index
significant_genes_df$gene_index = 1:nrow(significant_genes_df)

#add pseudotime and pseudotime_value_round variables to the seurat object used for the pseudotime analysis
risk_int@meta.data$pseudotime = pseudotime(cds)
risk_int$pseudotime_value_round = round(risk_int$pseudotime, digits = 1)
risk_int$pseudotime_value_round = factor(risk_int$pseudotime_value_round,
                                   levels = unique(risk_int$pseudotime_value_round)[order(unique(risk_int$pseudotime_value_round))])

#set Idents
Idents(risk_int) = "pseudotime_value_round"

#create average expression table
avgexp = AverageExpression(risk_int, 
                           assays = c("SCT"), 
                           layer = "data",
                           return.seurat = TRUE)

#create top_de_genes object
top_de_genes = significant_genes_df$gene_name

#extract expression values
expression.values = FetchData(avgexp, vars = cds_graph_test_results$gene_name, 
                              layer = "scale.data")

#create gene clustering function
cluster_pseudotime_genes = function(expression.values, k.treecut=5, keep.hclust=F){
  #expression.values: expression values from FetchData
  #k.treecut: cutoff value for the clustertree, based on this there will be more or less groups of genes
  expression.values = as.data.frame(t(expression.values))
  d = dist(expression.values, method = "euclidean")
  clust = hclust(d, method = "complete")
  plot(clust, labels = FALSE)
  clust = cutree(clust, k = k.treecut) %>%  data.frame()
  names(clust)="hcluster"
  expression.values = cbind(expression.values,clust)
  expression.values = arrange(expression.values,hcluster)
  
  order.df=NULL
  for (k in levels(factor(expression.values$hcluster))) {
    k_sub = expression.values[expression.values$hcluster==k,]
    k_sub$hcluster=NULL
    k_sub.means = colMeans(k_sub)
    df = data.frame("ps"=names(k_sub.means[k_sub.means == max(k_sub.means)]),"k"=k)
    order.df=rbind(order.df,df)
  }
  order.df = arrange(order.df,ps)
  expression.values$hcluster = factor(expression.values$hcluster,levels=order.df$k )
  expression.values = arrange(expression.values,hcluster)
  if (keep.hclust) {
    return(expression.values) 
  }else{
    expression.values$hcluster = sapply(expression.values$hcluster,function(k) strrep("_",k) )
    rownames(expression.values) = paste0(expression.values$hcluster,rownames(expression.values))
    expression.values$hcluster=NULL
    return(expression.values)}
}

#Cluster genes
riskMatrix <- cluster_pseudotime_genes(expression.values)

#plot trajectory genes in a heatmap
expression.heatmap = pheatmap(riskMatrix,
                              labels_col = "", fontsize_row = 3,
                              cluster_cols = FALSE, cluster_rows = FALSE,
                              angle_col = 0,border_color = 0,
                              treeheight_row = 0,
                              fontsize = 9, 
                              scale = "column",
                              breaks = seq(-4, 4, length.out = 100),
                              color = colorRampPalette(c("#053061", "#2171B5", "white", "#FFD92F", "#A50F15"))(100)
) 

#save the heatmap
save_plot(filename = "filePath/trajectoryHeatmap.png", plot = expression.heatmap, limitsize = FALSE,
          base_height = 30, base_width = 15)