# ==============================================================================
#          Setup Seurat Object
# ==============================================================================

library(dplyr)
library(Seurat)
library(patchwork)

# load snail dataset
snail.data <- readRDS("/Users/andreahassig/Desktop/lyons lab/scRNA-seq/WTAATAC_BD_files/snail-full-v2ref_Seurat.rds")


# =============================================================================
#          NORMALIZING THE DATA
# =============================================================================

library(ggplot2)
library(sctransform)
library(BiocManager)
library(glmGamPoi)

# run sctransform normalization
snail.data <- SCTransform(snail.data, verbose = FALSE)


# =============================================================================
#          RUN PCA & UMAP
# =============================================================================

snail.data <- RunPCA(snail.data, verbose = FALSE)
snail.data <- RunUMAP(snail.data, dims = 1:50, verbose = FALSE)

snail.data <- FindNeighbors(snail.data, dims = 1:50, verbose = FALSE)
snail.data <- FindClusters(snail.data, verbose = FALSE)
DimPlot(snail.data, label = TRUE)


# =============================================================================
#          IDENTIFY HIGHLY VARIABLE FEATURES
# =============================================================================

snail.data <- FindVariableFeatures(snail.data, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(snail.data), 10)
# plot variable features with labels
plot1 <- VariableFeaturePlot(snail.data)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


# =============================================================================
#          SCALING THE DATA
# =============================================================================

# shift the expression of each gene, so that the mean expression across cells is 0
# and the variance across cells is 1
all.genes <- rownames(snail.data)
snail.data <- ScaleData(snail.data, features = all.genes)


# ==============================================================================
#          LINEAR DIMENSIONAL REDUCTION & CLUSTERING
# ==============================================================================

snail.data <- RunPCA(snail.data, features = VariableFeatures(object = snail.data))
DimHeatmap(snail.data, dims = 50, cells = 500, balanced = TRUE)

# Cluster the cells
snail.data <- FindNeighbors(snail.data, dims = 1:50)
snail.data <- FindClusters(snail.data, resolution = 0.5)

snail.data <- RunUMAP(snail.data, dims = 1:50)
DimPlot(snail.data, reduction = "umap", label = TRUE)


# ==============================================================================
#          FINDING DIFFERENTIALLY EXPRESSED FEATURES
# ==============================================================================

# find all markers of cluster 2
cluster2.markers <- FindMarkers(snail.data, ident.1 = 2)
head(cluster2.markers, n = 5)

VlnPlot(snail.data, features = c("STRG.73049", "STRG.60049", "STRG.38977", "STRG.79864", "STRG.3200"))
