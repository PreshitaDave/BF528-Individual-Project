# Author: Preshita Dave
# File Description - Preprocessing the UMI Counts Matrix and saving the rds object with cluster information

# loading required libraries
library(dplyr)
library(fishpond)
library(patchwork)
library(Seurat)
library(tximport)
library(Matrix)
library(tidyverse)
#BiocManager::install("EnsDb.Hsapiens.v79")
library(EnsDb.Hsapiens.v79)

setwd('/projectnb/bf528/students/preshita/project-4/programmer')

data_dir <- '/projectnb/bf528/students/preshita/project-4/data/Alevin_Output/alevin/quants_mat.gz'

txi <- tximport(data_dir, type="alevin")

# getting the ensembl ids
ensembl <-  rownames(txi$counts) 
ensembl <- sub("[.][0-9]*","",ensembl) #remove ends
rownames(txi$counts) <- ensembl

# convert ensembl ids to gene names
symbols <- select(EnsDb.Hsapiens.v79, keys= ensembl, keytype = "GENEID", columns = c("SYMBOL","GENEID")) 

# rename the ensembl ids to gene names
txi$counts <- txi$counts[rownames(txi$counts) %in% symbols$GENEID,]
rownames(txi$counts) <- symbols$SYMBOL


############################################################################################

# creating a seurat object
pc <- CreateSeuratObject(counts = txi$counts, project="panc")

# subsetting mitchondrial genes
pc[["percent.mt"]] <- PercentageFeatureSet(pc, pattern = "^MT-")

VlnPlot(pc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1+plot2

#subsetting the matrix 
pc <- subset(pc, subset = nFeature_RNA > 200 & percent.mt < 10)

pc

#normalize the data
pc <- NormalizeData(pc)

#Identification of highly variable features 
pc <- FindVariableFeatures(pc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(pc), 10)

#variable plot
plot1 <- VariableFeaturePlot(pc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

#Scaling the data
all.genes <- rownames(pc)
pc <- ScaleData(pc, features = all.genes)

#Perform linear dimensional reduction
pc <- RunPCA(pc, features = VariableFeatures(object = pc))
# Examine and visualize PCA results a few different ways
print(pc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pc, dims = 1:2, reduction = "pca")

DimPlot(pc, reduction = "pca")
DimHeatmap(pc, dims = 1, cells = 500, balanced = TRUE)

#Determine the ‘dimensionality’ of the dataset
ElbowPlot(pc, ndims = 50)

#Cluster the cells
pc <- FindNeighbors(pc, dims = 1:30)
pc <- FindClusters(pc, resolution = 0.5)
head(Idents(pc), 5)

pc <- RunUMAP(pc, dims = 1:30)
DimPlot(pc, reduction = "umap")

# plot counts, proportions
counts <- as.vector(table(Idents(pc)))
names <- names(table(Idents(pc)))

png('barplot.png')
bar <- barplot(height = table(Idents(pc)), names.arg = (names), xlab = 'Cluster', ylab = 'Total cells')
text(x = bar, y = table(Idents(pc))-60, label = table(Idents(pc)), pos = 3, cex = 0.8)
dev.off()

png('pie.png')
slices <- as.vector(prop.table(table(Idents(pc))))
lbls <- names(prop.table(table(Idents(pc))))
pie(slices, labels = lbls, main="Relative proportions")
dev.off()

# saving the RDS file
saveRDS(pc, file = "panc.rds")








