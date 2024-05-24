# Load necessary libraries
if (!requireNamespace("Seurat", quietly = TRUE)) {
  install.packages("Seurat")
}
library(Seurat)
library(dplyr)
library(ggplot2)

# Read the data and create Seurat object
data.dir <- "/Users/zamilakarasart/Desktop/Master thesis/GSM6879661_C7_filtered_feature_bc_"
GSM6879661_C7_data <- Read10X(data.dir = data.dir)
seuratObject <- CreateSeuratObject(counts = GSM6879661_C7_data, 
                                   project = "Kaeda_green_myeloid_cells", 
                                   min.cells = 3, 
                                   min.features = 200)

# Calculate percentage of mitochondrial genes
seuratObject[["percent.mt"]] <- PercentageFeatureSet(seuratObject, pattern = "^mt-")
head(seuratObject@meta.data, 5)

# Visualize QC metrics
VlnPlot(seuratObject, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(seuratObject, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seuratObject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter cells based on QC metrics
seuratObject <- subset(seuratObject, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalize data
seuratObject <- NormalizeData(seuratObject, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
seuratObject <- FindVariableFeatures(seuratObject, selection.method = "vst", nfeatures = 2000)

# Visualize top 10 highly variable genes
top10 <- head(VariableFeatures(seuratObject), 10)
plot1 <- VariableFeaturePlot(seuratObject)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Scale the data
seuratObject <- ScaleData(seuratObject, features = VariableFeatures(object = seuratObject))

# Perform PCA
seuratObject <- RunPCA(seuratObject, features = VariableFeatures(object = seuratObject))
DimPlot(seuratObject, reduction = "pca")
DimHeatmap(seuratObject, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(seuratObject)

# Find neighbors and clusters
seuratObject <- FindNeighbors(seuratObject, dims = 1:20)
seuratObject <- FindClusters(seuratObject, resolution = 0.5)

# Run UMAP for visualization
seuratObject <- RunUMAP(seuratObject, dims = 1:20)
DimPlot(seuratObject, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


seuratObject_markers_clusters <- FindAllMarkers(seuratObject, 
                                                 only.pos = TRUE, 
                                                 min.pct = 0.25, 
                                                 logfc.threshold = 0.25)


# Extract top 10 markers for each cluster
top_markers <- seuratObject_markers_clusters %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)

# Check the top markers
head(top_markers)

# Visualize top markers using DoHeatmap
top10_markers <- top_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC) %>% 
  ungroup() %>% 
  pull(gene) %>% 
  unique()

DoHeatmap(seuratObject, features = top10_markers) + NoLegend()

# Visualize top markers for specific clusters using VlnPlot
# Example: Top markers for cluster 0
cluster_10_markers <- top_markers %>% 
  filter(cluster == 10) %>% 
  pull(gene)

VlnPlot(seuratObject, features = cluster_10_markers, pt.size = 0.1) + NoLegend()
print(cluster_10_markers)



# Markers for myeloid cells

macrophage_markers <- c("Cd68", "Cd14", "Itgam", "Mrc1", "Fcgr3", "Adgre1", "Csf1r", "Cd163", "Trem2", "Maf")
dc_markers <- c("Itgax", "Hla-Dra", "Hla-Drb1", "Cd80", "Cd86", "Cd40", "Xcr1", 
                "Clec9a", "Batf3", "Cd1c", "Fcer1a", "Sirpa", "Ly6c2", 'Fcgr1',
                "Ccr2", "Cebpb", "Itgax", "H2-Ab1", "Clec9a", "Xcr1")
langerhans_markers <- c("Cd1a", "Langerin", "Hla-Dra", "Hla-Drb1", "Itgax", "Cd1c", "Epcam", "Clec4k")
mast_cell_markers <- c("Kit", "Fcεri", "Mcpt1", "Cma1", "Tpsb2", "Hdc", "Car4", "Il1rl1", "Cd34")
pDC_markers <- c("Clec4c", "Il3ra", "Tcf4", "Lyl1", "Gzmb", "Irf7", "Lilra4", "Pdha2", "SpiB")
mono_DCs_markers <- c("Ly6c2", "Fcgr1", "Ccr2", "Cebpb")
DCs_marker <- c("H2-Ab1", "Itgax")
cDC1_markers <- c("Clec9a","Xcr1","Itgae", "Batf3", "Irf8")
cDC2_markers <- c("Clec10a", "Cd14", "Cd209a", "Itgam", "Irf4")
activation_marker <- c("Ccr7", "Cd274", "Pdcd1lg2", "Fscn1", "Cd200")
proliferation_marker <- c("Mki67", "Birc5")

all_markers <- c(mono_DCs_markers, 
                 DCs_marker, 
                 cDC1_markers, 
                 cDC2_markers, 
                 activation_marker, 
                 proliferation_marker)

# Create DotPlot for all markers
DotPlot(seuratObject, features = all_markers) + 
  RotatedAxis() +
  ggtitle("Dot Plot of All Markers") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Visualize the expression of macrophage markers
FeaturePlot(seuratObject, features = macrophage_markers)
VlnPlot(seuratObject, features = macrophage_markers, pt.size = 0.1)

# Visualize the expression of dendritic cell markers
FeaturePlot(seuratObject, features = dc_markers)
VlnPlot(seuratObject, features = dc_markers, pt.size = 0.1)

# Visualize the expression of Langerhans cell markers
FeaturePlot(seuratObject, features = langerhans_markers)
VlnPlot(seuratObject, features = langerhans_markers, pt.size = 0.1)

# Visualize the expression of mast cell markers
FeaturePlot(seuratObject, features = mast_cell_markers)
VlnPlot(seuratObject, features = mast_cell_markers, pt.size = 0.1)

# Visualize the expression of plasmacytoid dendritic cell markers
FeaturePlot(seuratObject, features = pDC_markers)
VlnPlot(seuratObject, features = pDC_markers, pt.size = 0.1)

# Visualize the expression of conventional dendritic cell a markers
FeaturePlot(seuratObject, features = cDC1_markers)
VlnPlot(seuratObject, features = cDC1_markers, pt.size = 0.1)

FeaturePlot(seuratObject, features = proliferation_marker)
VlnPlot(seuratObject, features = proliferation_marker, pt.size = 0.1)

dot_plot_all <- DotPlot(seuratObject, features = c("Jaml", "Irf8", "Irf4")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "JAML Expression in Different Clusters",
       x = "Clusters",
       y = "Expression Level") +
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  scale_size(range = c(2, 8))

# Display the plot
print(dot_plot_all)


cluster_2_markers <- top_markers %>% 
  filter(cluster == 2) %>% 
  pull(gene)

VlnPlot(seuratObject, features = cluster_2_markers, pt.size = 0.1) + NoLegend()
print(cluster_2_markers)

cluster_3_markers <- top_markers %>% 
  filter(cluster == 3) %>% 
  pull(gene)

VlnPlot(seuratObject, features = cluster_3_markers, pt.size = 0.1) + NoLegend()
print(cluster_3_markers)


cluster_4_markers <- top_markers %>% 
  filter(cluster == 4) %>% 
  pull(gene)

VlnPlot(seuratObject, features = cluster_4_markers, pt.size = 0.1) + NoLegend()
print(cluster_4_markers)


cluster_5_markers <- top_markers %>% 
  filter(cluster == 5) %>% 
  pull(gene)

VlnPlot(seuratObject, features = cluster_5_markers, pt.size = 0.1) + NoLegend()
print(cluster_5_markers)


# Load necessary libraries
library(Seurat)
library(dplyr)
library(ggplot2)




# Subset the Seurat object to include only the specified clusters
cDC_clusters <- subset(seuratObject, idents = c(2, 5, 6))

# Recalculate the UMAP for the subsetted data
cDC_clusters <- ScaleData(cDC_clusters)
cDC_clusters <- RunPCA(cDC_clusters, npcs = 50, verbose = FALSE)
ElbowPlot(cDC_clusters)
cDC_clusters <- FindNeighbors(cDC_clusters, dims = 1:10)
cDC_clusters <- FindClusters(cDC_clusters, resolution = 0.5)
cDC_clusters <- RunUMAP(cDC_clusters, dims = 1:10)

# Visualize the UMAP
DimPlot(cDC_clusters, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()




# Identify markers for each cluster in the subsetted data
cluster_markers <- FindAllMarkers(cDC_clusters, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

mono_DCs_markers <- c("Ly6c2", "Fcgr1", "Ccr2", "Cebpb")
DCs_marker <- c("H2-Ab1", "Itgax")
cDC1_markers <- c("Clec9a","Xcr1","Itgae", "Batf3", "Irf8")
cDC2_markers <- c("Clec10a", "Cd14", "Cd209a", "Itgam", "Irf4")
pDC_markers <- c("Clec4c", "Il3ra", "Tcf4", "Lyl1", "Gzmb", "Irf7", "Lilra4", "Pdha2", "SpiB", 
                 "Bst2", "Ccr9", "Siglech", "Tlr7", "Tlr9", "Blr1", "Ptpn13", "Zbtb46")
activation_marker <- c("Ccr7", "Cd274", "Pdcd1lg2", "Fscn1", "Cd200")
proliferation_marker <- c("Mki67", "Birc5")

all_markers <- c(mono_DCs_markers, 
                 DCs_marker, 
                 cDC1_markers, 
                 cDC2_markers,
                 activation_marker, 
                 proliferation_marker)

# Create DotPlot for all markers
DotPlot(cDC_clusters, features = all_markers) + 
  RotatedAxis() +
  ggtitle("Dot Plot of All Markers") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Visualize expression of these markers
FeaturePlot(cDC_clusters, features = mono_DCs_markers)
VlnPlot(cDC_clusters, features = mono_DCs_markers, pt.size = 0.1)

FeaturePlot(cDC_clusters, features = DCs_marker)
VlnPlot(cDC_clusters, features = DCs_marker, pt.size = 0.1)

FeaturePlot(cDC_clusters, features = cDC1_markers)
VlnPlot(cDC_clusters, features = cDC1_markers, pt.size = 0.1)

FeaturePlot(cDC_clusters, features = cDC2_markers)
VlnPlot(cDC_clusters, features = cDC2_markers, pt.size = 0.1)

FeaturePlot(cDC_clusters, features = pDC_markers)
VlnPlot(cDC_clusters, features = pDC_markers, pt.size = 0.1)

FeaturePlot(cDC_clusters, features = activation_marker)
VlnPlot(cDC_clusters, features = activation_marker, pt.size = 0.1)

FeaturePlot(cDC_clusters, features = prolifiration_marker)
VlnPlot(cDC_clusters, features = prolifiration_marker, pt.size = 0.1)



cluster_name <- c("0" = "CCR7_DC.3",
                  "1" = "cDC2.2",
                  "2" = "cDC2.1",
                  "3" = "CCR7_DC.2",
                  "4" = "CCR7_DC.1")

# Rename clusters
cDC_clusters <- RenameIdents(cDC_clusters, cluster_name)

current_clusters <- Idents(cDC_clusters)

# Order clusters alphabetically
ordered_clusters <- factor(current_clusters, levels = sort(levels(current_clusters)))

# Set the ordered clusters back to the Seurat object
cDC_clusters <- SetIdent(cDC_clusters, value = ordered_clusters)

# Create the Dot Plot with ordered clusters
dot_plot <- DotPlot(cDC_clusters, features = c("Jaml", "Irf8", "Irf4", "Tbet")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "JAML Expression",
       x = "Clusters",
       y = "Expression Level") +
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  scale_size(range = c(2, 8))

# Display the plot
print(dot_plot)
