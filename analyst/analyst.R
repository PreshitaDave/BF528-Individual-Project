##########################################################################################
# Author: Preshita Dave
# Description: This script analyses and identifies marker genes for clusters that were
# created. 

##########################################################################################

library(Seurat)
library(patchwork)
library(dplyr)
library(data.table)

setwd('/projectnb/bf528/users/swiss_cheese_2022/project_4/analyst')

cells<- readRDS('/projectnb/bf528/users/swiss_cheese_2022/project_4/github_project-4-swisscheese/pbmc_prelim.rds')


# Identifies DEGs between two groups of cells using a Wilcoxon Rank Sum test
all.markers <- FindAllMarkers(object = cells, min.pct = 0.25, only.pos = TRUE)

#filtering for p-adj<0.05 for retaining significant genes
cells.markers.padj <- all.markers[all.markers$p_val_adj<0.05,]

write.csv(cells.markers.padj, 'Cell_markers.csv')
#genes of interest: goi based on the paper supplemental information
goi <- c("GCG", "INS",  "SST", "PPY", "GHRL", "KRT19",
                         "CPA1", "PDGFRB", "VWF", "PECAM1", "CD34",
                         "CD163", "CD68", "IgG","CD3","CD8", "TPSAB1", "KIT", "CPA3")

violin <- VlnPlot(cells, features = goi, pt.size=0)
violin
# "IgG","CD3","CD8" are not found when the violin plot is plotted, hence we can remove them from goi 
# Hence Cytotoxic T cells can be eliminated from the analysis. 
goi <- c("GCG", "INS",  "SST", "PPY", "GHRL", "KRT19",
         "CPA1", "PDGFRB", "VWF", "PECAM1", "CD34",
         "CD163", "CD68", "TPSAB1", "KIT", "CPA3")

violin2 <- VlnPlot(cells, features = goi, pt.size=0)
png(filename = "violin.png", width = 1200, height = 800)
violin2
dev.off()
# based on the violin plot, we can see that 
# epsilon - GHRL not there
# Vascular - VWF,PECAM1, CD34 not at all 
# Macrophage - CD68 only, CD163 and IgG not there
# Cytotoxic T - CD3 and CD8 not at all 
# Mast - TPSAB1, KIT, CPA3 not at all 


# to create a feature plot to understand the distribution of cells 
feature_plot <- FeaturePlot(cells, features = c("GCG","INS","SST","PPY","GHRL","KRT19",
                                               "CPA1","PDGFRB","VWF","PECAM1","CD34",
                                               "CD163","CD68","TPSAB1",
                                               "KIT","CPA3"))
png(filename = "feature.png", width = 1200, height = 800)
feature_plot
dev.off()
# finding top 10 markers for every cluster with the highest positive log2fc score
top_10 <- cells.markers.padj %>% group_by(cluster) %>% top_n(10, avg_log2FC)


# to find the presence of specific markers in each cluster
alpha_clust <- subset(cells.markers.padj, gene %like% "GCG")
beta_clust <- subset(cells.markers.padj, gene %like% "INS")
delta_clust <- subset(cells.markers.padj, gene %like% "SST")
gamma_clust <- subset(cells.markers.padj, gene %like% "PPY")
epsilon_clust <- subset(cells.markers.padj, gene %like% "GHRL")
ductal_clust <- subset(cells.markers.padj, gene %like% "KRT19")
acinar_clust <- subset(cells.markers.padj, gene %like% "CPA1")
stellate_clust <- subset(cells.markers.padj, gene %like% "PDGFRB")
vascular_clust <- subset(cells.markers.padj, gene %like% "VWF") # or PECAM1 
macrophage_clust <- subset(cells.markers.padj, gene %like% "CD68")

# to find all the markers in each cluster
cluster0.markers <- FindMarkers(cells, ident.1 = 0, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.5)
cluster1.markers <- FindMarkers(cells, ident.1 = 1, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.5)
cluster2.markers <- FindMarkers(cells, ident.1 = 2, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.5)
cluster3.markers <- FindMarkers(cells, ident.1 = 3, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.5)
cluster4.markers <- FindMarkers(cells, ident.1 = 4, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.5)
cluster5.markers <- FindMarkers(cells, ident.1 = 5, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.5)
cluster6.markers <- FindMarkers(cells, ident.1 = 6, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.5)
cluster7.markers <- FindMarkers(cells, ident.1 = 7, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.5)
cluster8.markers <- FindMarkers(cells, ident.1 = 8, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.5)
cluster9.markers <- FindMarkers(cells, ident.1 = 9, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.5)
cluster10.markers <- FindMarkers(cells, ident.1 = 10, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.5)
cluster11.markers <- FindMarkers(cells, ident.1 = 11, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.5)
cluster12.markers <- FindMarkers(cells, ident.1 = 12, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.5)

### Assigning cell types to clusters based on expression of marker genes

## Alpha 
# Clusters 1, 2
alpha.markers <- FindMarkers(cells, ident.1 = 1, ident.2 = 2)
head(alpha.markers, n = 5)
# GCG most significantly expressed in cluster 1 


## Beta
# Clusters 4, 6
beta.markers <- FindMarkers(cells, ident.1 = 4, ident.2 = 6)
head(beta.markers, n = 5)
# INS most significantly expressed in cluster 4

## Delta 
# Cluster 0 significantly expresses SST 

## Gamma
# Cluster 5 significantly expresses PPY 

## Ductal 
# Clusters 8, 9
ductal.markers <- FindMarkers(cells, ident.1 = 9, ident.2 = 8)
head(ductal.markers, n = 5)
# KRT19 significantly expressed in cluster 9 

## Acinar 
# Clusters 3, 10 
acinar.markers <- FindMarkers(cells, ident.1 = 3, ident.2 = 10)
head(acinar.markers, n = 5)
# CPA1 significantly expressed in cluster 3
# SERPINA3, REG1B are markers for acinar cells significantly expressed in cluster 10 
# (https://www.panglaodb.se/markers.html?cell_type=%27Acinar%20cells%27)

## Stellate 
# Cluster 7 expresses PDGFRB significantly 

## Vascular 
# Cluster 12 expresses VWF and PECAM1 significantly 

## Macrophage 
# Cluster 10, 11
macrophage.markers <- FindMarkers(cells, ident.1 = 11, ident.2 = 10)
head(macrophage.markers, n = 5)
# Cluster 11 expresses CD68 significantly 

# Based on this, we can assign identities to the clusters based on the markers
# cluster 0 - SST -> Delta
# cluster 1 - GCG -> Alpha1
# cluster 2 - GCG -> Alpha2
# cluster 3 - CPA1 -> Acinar1
# cluster 4 - INS -> Beta1
# cluster 5 - PPY -> Gamma
# cluster 6 - INS -> Beta2
# cluster 7 - PDGFRB -> Stellate
# cluster 8 - KRT19 -> Ductal2
# cluster 9 - KRT19 -> Ductal1
# cluster 10 - CPA1 -> Acinar2
# cluster 11 - CD68 -> Macrophage
# cluster 12 - PECAM1 -> Vascular

# rename cluster with cell type 
cluster_cell_types <- c("Delta","Alpha1","Alpha2","Acinar1","Beta1","Gamma","Beta2", "Stellate","Ductal2","Ductal1","Acinar2","Macrophage","Vascular")
names(cluster_cell_types) <- levels(cells)
cells <- RenameIdents(cells, cluster_cell_types)


# visualize the top marker genes per cluster
top_5 <- cells.markers.padj %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

# heat map
heatmap_5 <- DoHeatmap(cells, features = top_5$gene, group.bar = TRUE, angle=70, size = 4)
png(filename = "heatmap.png", width = 1200, height = 800)
heatmap_5
dev.off()

# visualize the clustered cells using a projection method - umap
ElbowPlot(cells, ndims=50)
umap <- RunUMAP(cells, dims = 1:30)
umap <- DimPlot(cells, reduction = "umap", label=TRUE)
png(filename = "umap.png", width = 1200, height = 800)
umap
dev.off()

### Finding novel marker genes with stricter min percent
novel_genes <- FindAllMarkers(cells, only.pos = TRUE, min.pct = 0.8, logfc.threshold = 1)
# stricter adjusted p value
signficant_novel_genes <- novel_genes[novel_genes$p_val_adj<0.005,]
write.csv(signficant_novel_genes, "novel_genes.csv")

