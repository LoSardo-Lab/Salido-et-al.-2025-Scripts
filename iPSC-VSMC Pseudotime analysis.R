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
cds <- as.cell_data_set(risk_int
                        )
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
cds_graph_test_results <- graph_test(cds, neighbor_graph = "principal_graph", cores = 8)