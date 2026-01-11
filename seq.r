# ==============================================================================
#          PART 1
# ==============================================================================
# 1. Libraries
library(Seurat)
library(dplyr)
library(ggplot2)

# 2. Load Data
pbmc <- readRDS("/Users/andreahassig/Desktop/lyons lab/scRNA-seq/WTAATAC_BD_files/snail-full-v2ref_Seurat.rds")
anno_df <- read.csv("/Users/andreahassig/Desktop/lyons lab/scRNA-seq/scRNAseq counts and annotations/annotated_all_markers_50dim.csv")

# 3. Initial Prep
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc)


# ==============================================================================
#          PART 2
# ==============================================================================
# 1. Calculate Mitochondrial % 
mito_genes <- unique(as.character(anno_df$gene[grep("mitochondrial|MT-", anno_df$Protein_Name, ignore.case = TRUE)]))
genes_present <- intersect(mito_genes, rownames(pbmc))
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, features = genes_present)

# 2. Run the 50-Dimension 
pbmc <- RunPCA(pbmc, npcs = 50, verbose = FALSE)
pbmc <- FindNeighbors(pbmc, dims = 1:50)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# 3. Create the UMAP 
pbmc <- RunUMAP(pbmc,
    dims = 1:50,
    seed.use = 42,
    min.dist = 0.2,
    spread = 2.0
)

# Plot A: The Violins (Genes, Counts, and Soup)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
# with points
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Plot B: The Diagonal (Technical Consistency)
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Plot C: The Elbow Plot (Justifying 50 Dimensions)
ElbowPlot(pbmc, ndims = 50)

# Plot D: The Map
DimPlot(pbmc, reduction = "umap", label = TRUE) + NoLegend()


# ==============================================================================
#          PART 3
# ==============================================================================

# Variable Gene Comparison: 2000 vs 3000
print("Running Variable Gene Comparison: 2000 vs 3000...")

# 1. Standard 2000 Genes
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
top10_2000 <- head(VariableFeatures(pbmc), 10)
plot_2000 <- VariableFeaturePlot(pbmc) + ggtitle("Top 2000 Genes")
plot_2000_labeled <- LabelPoints(plot = plot_2000, points = top10_2000, repel = TRUE)

# 2. Expanded 3000 Genes
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
top10_3000 <- head(VariableFeatures(pbmc), 10)
plot_3000 <- VariableFeaturePlot(pbmc) + ggtitle("Top 3000 Genes")
plot_3000_labeled <- LabelPoints(plot = plot_3000, points = top10_3000, repel = TRUE)

# 3. Compare Side-by-Side
library(patchwork)
print(plot_2000_labeled + plot_3000_labeled)


# PCA Heatmaps for Dimensions 1, 45, and 50
print("Generating PCA Heatmaps for Dimensions 1, 45, and 50...")

# This shows the top genes driving each of these PCs
# We look at 500 cells to get a good representative sample
DimHeatmap(pbmc, dims = c(1, 45, 50), cells = 500, balanced = TRUE)


# Dimensionality Showdown: 40 vs 50
print("Running the final Dimensionality Showdown: 40 vs 50...")

# --- Test 1: 40 Dimensions ---
pbmc <- FindNeighbors(pbmc, dims = 1:40, verbose = FALSE)
pbmc <- RunUMAP(pbmc, dims = 1:40, seed.use = 42, verbose = FALSE)
plot_40 <- DimPlot(pbmc, reduction = "umap", label = TRUE) +
    ggtitle("40 Dimensions (Cleaner Math)") + NoLegend()

# --- Test 2: 50 Dimensions ---
pbmc <- FindNeighbors(pbmc, dims = 1:50, verbose = FALSE)
pbmc <- RunUMAP(pbmc, dims = 1:50, seed.use = 42, verbose = FALSE)
plot_50 <- DimPlot(pbmc, reduction = "umap", label = TRUE) +
    ggtitle("50 Dimensions (Higher Sensitivity)") + NoLegend()

# View them side-by-side in external window
library(patchwork)
print(plot_40 + plot_50)


# Refined Variable Gene Selection with Custom Cutoffs
print("Refining Variable Gene Selection with Custom Cutoffs...")

# 1. Re-run the selection but look at the actual values
pbmc <- FindVariableFeatures(pbmc,
    selection.method = "vst",
    nfeatures = 3000
)

# 2. Define thresholds to 'cut' the bottom and left
# can adjust further if needed
min_mean <- 0.1 # Cuts the left side (low expression)
min_disp <- 0.5 # Cuts the bottom (low variance)

# 3. Visualize the 'Cut'
v_plot <- VariableFeaturePlot(pbmc)
v_plot_refined <- v_plot +
    geom_vline(xintercept = min_mean, linetype = "dashed", color = "blue") +
    geom_hline(yintercept = min_disp, linetype = "dashed", color = "blue") +
    ggtitle("Variable Gene Selection: Identifying the 'High Signal' Zone")

print(v_plot_refined)

h

# Identify the genes in this new window
vf_info <- pbmc@assays$RNA@meta.features
new_high_signal <- rownames(vf_info[vf_info$vst.mean > min_mean & vf_info$vst.variance.standardized > min_var, ])

# Update and Re-run
VariableFeatures(pbmc) <- new_high_signal
pbmc <- ScaleData(pbmc, features = new_high_signal, verbose = FALSE)
pbmc <- RunPCA(pbmc, npcs = 50, verbose = FALSE)
pbmc <- FindNeighbors(pbmc, dims = 1:50, verbose = FALSE)
pbmc <- FindClusters(pbmc, resolution = 0.5, verbose = FALSE)
pbmc <- RunUMAP(pbmc, dims = 1:50, verbose = FALSE)

DimPlot(pbmc, reduction = "umap", label = TRUE)


# =============================================================================
#          PART 4
# =============================================================================

# This will take a minute...
print("Finding markers for all clusters... sit tight.")
all_markers <- FindAllMarkers(pbmc, 
                               only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.25)

# Group to see the top 5 genes for every cluster
top5_markers <- all_markers %>%
    group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC)

View(top5_markers)

# no thresholding
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1)

pbmc.top5 <- pbmc.markers %>%
    group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC)

View(pbmc.top5)


# 1. Replace 'Inf' with a high but readable number (like 10)
all_markers_clean <- pbmc.markers
all_markers_clean$avg_log2FC[is.infinite(all_markers_clean$avg_log2FC)] <- 10

# 2. Filter for quality: High pct.1 (>0.4) and a big gap between pct.1 and pct.2
# This removes genes that are only in 2 cells (which often cause the 'Inf' error)
reliable_markers <- all_markers_clean %>%
    filter(p_val_adj < 0.05) %>%
    filter(pct.1 > 0.40 & pct.2 < 0.20) %>% 
    group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC)

View(reliable_markers)




library(patchwork)

# Replace 'GENE' with a top marker from our list 
target_gene <- "STRG.90096" 

plot_feature <- FeaturePlot(pbmc, features = target_gene) + 
    theme(legend.position = "none") +
    ggtitle(paste("Location:", target_gene))

plot_violin <- VlnPlot(pbmc, features = target_gene, pt.size = 0.1) + 
    theme(legend.position = "none") +
    ggtitle(paste("Expression Level:", target_gene))

# This puts them side-by-side
print(plot_feature | plot_violin)

library(ggplot2)


# just violin
VlnPlot(pbmc, features = c("STRG.110237"))

VlnPlot(pbmc, 
        features = "STRG.110237", 
        slot = "data",      # Forces use of log-normalized expression levels
        pt.size = 0.2,      # Adds the dots back
        cols = NULL) +      # Keeps standard cluster colors
  ylim(0, 7) +              # CAPS the y-axis so we can see the distribution




