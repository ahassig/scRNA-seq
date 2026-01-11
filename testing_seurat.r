# ==============================================================================
#          Setup Seurat Object
# ==============================================================================

library(dplyr)
library(Seurat)
library(patchwork)

# load snail dataset
snail.data <- readRDS("/Users/andreahassig/Desktop/lyons lab/scRNA-seq/WTAATAC_BD_files/snail-full-v2ref_Seurat.rds")

# Visualize QC metrics as a violin plot
VlnPlot(snail.data, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

# FeatureScatter is typically used to visualize feature-feature relationships
feature.plot <- FeatureScatter(snail.data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


# =============================================================================
#          NORMALIZING THE DATA
# =============================================================================

library(ggplot2)
library(sctransform)
library(BiocManager)
library(glmGamPoi)

# run sctransform normalization
snail.data <- SCTransform(
    snail.data,
    method = "glmGamPoi", # Forces the fast math engine
    variable.features.n = 3000, # Increase to 3,000 variable genes
    vst.flavor = "v2", # The most modern/stable SCT version
    verbose = TRUE # Keep this TRUE so we can see it moving!
)


# =============================================================================
#          RE-REFINED VARIABLE FEATURE SELECTION
# =============================================================================

# 1. Visualize thresholds (SCTransform version)
# We use the SCT assay here
v_plot <- VariableFeaturePlot(snail.data, assay = "SCT")
v_plot_refined <- v_plot +
    geom_vline(xintercept = 0.001, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 1.2, linetype = "dashed", color = "red") +
    ggtitle("Refined Thresholds (SCT Data)")

print(v_plot_refined)

# 2. Get the data table specifically from the SCT assay
var_info <- HVFInfo(snail.data, assay = "SCT")

# 3. Filter the genes based on your thresholds
# (SCT uses 'gmean' for the x-axis and 'residual_variance' for the y-axis)
refined_genes <- rownames(var_info)[
    var_info$gmean >= 0.001 &
        var_info$residual_variance >= 1.2
]
length(refined_genes)

# 4. Update the object with your hand-picked genes
VariableFeatures(snail.data) <- refined_genes
message("Genes remaining after cut: ", length(VariableFeatures(snail.data)))


# =============================================================================
#          RUN PCA & UMAP
# =============================================================================

# 1. Run PCA
snail.data <- RunPCA(snail.data, assay = "SCT", verbose = FALSE)

# 2. Visualize the top genes driving the PCs
DimHeatmap(snail.data, dims = 1:15, cells = 500, balanced = TRUE)

# 3. Continue with UMAP and Clustering
snail.data <- RunUMAP(snail.data, dims = 1:30, verbose = FALSE)
snail.data <- FindNeighbors(snail.data, dims = 1:30, verbose = FALSE)
# snail.data <- FindClusters(snail.data, verbose = FALSE)
snail.data <- FindClusters(snail.data, resolution = 0.5, verbose = FALSE)

DimPlot(snail.data, label = TRUE)

# umap with counts
FeaturePlot(snail.data, features = "nCount_RNA")
# lower the max cutoff to better visualize
FeaturePlot(snail.data, features = "nCount_RNA", max.cutoff = 1500)


# ==============================================================================
#          FINDING DIFFERENTIALLY EXPRESSED FEATURES
# ==============================================================================

snail.data <- PrepSCTFindMarkers(snail.data)

# find top 5 markers of each cluster
all.markers <- FindAllMarkers(
    snail.data,
    only.pos = TRUE,
    min.pct = 0.25, # gene must be expressed in ≥25% of cells
    logfc.threshold = 0.5 # gene must be enriched in cluster
)

# pull top 5 markers for each cluster
top5 <- all.markers %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 5)

top5

# get the gene list (in right order, removing duplicates)
top5.genes <- top5 %>%
    arrange(cluster, desc(avg_log2FC)) %>%
    distinct(gene, .keep_all = TRUE) %>%
    pull(gene) %>%
    unique()

# make dotplot
DotPlot(
    snail.data,
    features = top5.genes,
    assay = "SCT",
    dot.scale = 3 # smaller dots
) +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5) # smaller text
    )

head(top5, n = 100)
print(top5, n = 100)

# cluster2.markers <- FindMarkers(snail.data, ident.1 = 2)
# head(cluster2.markers, n = 5)
# VlnPlot(snail.data, features = c("STRG.73049", "STRG.60049", "STRG.38977", "STRG.79864", "STRG.3200"))


fasta_of_interest_blast <- read.csv("/Users/andreahassig/Desktop/lyons lab/scRNA-seq/scRNAseq counts and annotations/full_annotations_50dim.csv")
unique_fasta_blast <- fasta_of_interest_blast %>%
    distinct(Query_ID_truncated, .keep_all = TRUE)
write.csv(unique_fasta_blast, "unique_annotations_30dim.csv")

# combine with the marker table, with only unique hits
unique_fasta_blast$gene <- unique_fasta_blast$Query_ID_truncated
markers_annotated_50dim <- merge(unique_fasta_blast, markers.txt, by = "gene", all = TRUE)

write.csv(markers_annotated_50dim, "annotated_all_markers_50dim.csv")
