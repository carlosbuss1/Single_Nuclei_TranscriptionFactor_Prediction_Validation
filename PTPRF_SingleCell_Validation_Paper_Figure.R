# Load necessary libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(harmony)

# Define paths
data_dir <- "/home/ulb/Desktop/carlos/Human Liver Atlas/rawData_human/countTable_human"
annotations_file <- "/home/ulb/Desktop/carlos/Human Liver Atlas/rawData_human/countTable_human/annot_humanAll.csv"
output_dir <- "/home/ulb/Desktop/carlos/Human Liver Atlas/output"

# Load the dataset and create Seurat object
liveratlas.data <- Read10X(data.dir = data_dir, gene.column = 1)
liveratlas <- CreateSeuratObject(counts = liveratlas.data, project = "liveratlas", min.cells = 3, min.features = 100)

# Add mitochondrial percentage to metadata
liveratlas[["percent.mt"]] <- PercentageFeatureSet(liveratlas, pattern = "^MT-")

# Quality control and filtering
qc_plot <- VlnPlot(liveratlas, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) + ggtitle("Quality Control Metrics")
liveratlas <- subset(liveratlas, subset = nFeature_RNA < 5000 & percent.mt < 7.5)

# Data normalization, scaling, and PCA
liveratlas <- NormalizeData(liveratlas, normalization.method = "LogNormalize", scale.factor = 10000)
liveratlas <- FindVariableFeatures(liveratlas, selection.method = "vst", nfeatures = 5000)
liveratlas <- ScaleData(liveratlas)
liveratlas <- RunPCA(liveratlas)

# Load and add annotation data
HumanLiverAtlas.anno <- read.csv(annotations_file)
rownames(HumanLiverAtlas.anno) <- HumanLiverAtlas.anno$cell

# Filter for complete metadata and subset liveratlas to include only matching cells
HumanLiverAtlas.anno <- HumanLiverAtlas.anno %>% filter(!is.na(diet) & !is.na(annot))
matching_cells <- intersect(colnames(liveratlas), HumanLiverAtlas.anno$cell)
liveratlas <- subset(liveratlas, cells = matching_cells)

# Add cell type annotations and Condition to metadata
liveratlas <- AddMetaData(liveratlas, metadata = HumanLiverAtlas.anno[, c("annot", "diet")], col.name = c("cellType", "Condition"))

# Verify that the metadata is complete and proceed with Harmony
if (!all(c("Condition", "cellType") %in% colnames(liveratlas@meta.data))) {
  stop("Error: Required metadata columns 'Condition' or 'cellType' are missing.")
}

# Run Harmony to integrate and correct batch effects
liveratlas <- RunHarmony(liveratlas, group.by.vars = "Condition", dims.use = 1:20)

# Use Harmony embeddings for UMAP and clustering
liveratlas <- RunUMAP(liveratlas, reduction = "harmony", dims = 1:20)
liveratlas <- FindNeighbors(liveratlas, reduction = "harmony", dims = 1:20)
liveratlas <- FindClusters(liveratlas, resolution = 0.5)

# UMAP Plot for Cell Type (without labels on the plot, only in the legend)
umap_plot_clusters <- DimPlot(liveratlas, reduction = "umap", group.by = "cellType") +
  ggtitle("UMAP Plot with Cell Type Annotations") +
  theme_minimal(base_size = 14) +  # Adjusted for a more visible publication figure
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 12),
    panel.border = element_rect(color = "black", fill = NA, size = 0.7),
    axis.line = element_line(size = 0.5)
  )

# UMAP Plot for Diet Condition (using bright pink and bright blue for vibrant contrast)
umap_plot_diet <- DimPlot(liveratlas, reduction = "umap", group.by = "Condition") +
  scale_color_manual(values = c("Lean" = "#4169E1", "Obese" = "#FF1493")) +  # Bright blue and bright pink
  ggtitle("UMAP Plot with Diet Annotations") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 12),
    panel.border = element_rect(color = "black", fill = NA, size = 0.7),
    axis.line = element_line(size = 0.5)
  )

# Save plots with 600 DPI resolution
ggsave(filename = file.path(output_dir, "UMAP_plot_clusters_600dpi.png"), plot = umap_plot_clusters, width = 8, height = 6, dpi = 600)
ggsave(filename = file.path(output_dir, "UMAP_plot_diet_600dpi.png"), plot = umap_plot_diet, width = 8, height = 6, dpi = 600)

# Define genes of interest and retrieve expression data by condition
selected_genes <- c("PTPRF", "NFYA", "ARNT", "NFATC2", "NR4A2", "SP1")
expression_data <- GetAssayData(liveratlas, slot = "counts")[selected_genes, ]
expression_df <- as.data.frame(t(expression_data))

# Add condition metadata to expression data
expression_df <- cbind(expression_df, Condition = liveratlas@meta.data$Condition)

# Save expression data table
write.csv(expression_df, file = file.path(output_dir, "expression_data_by_condition.csv"), row.names = TRUE)

# Dot Plot for selected genes by cell type
dotplot <- DotPlot(liveratlas, features = selected_genes, group.by = "cellType") +
  scale_color_gradient(low = "lightblue", high = "darkblue") +  # Color gradient for dot intensity
  ggtitle("Dot Plot of Selected Genes by Cell Type") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 12)
  )

# Save Dot Plot with 600 DPI resolution
ggsave(filename = file.path(output_dir, "DotPlot_selected_genes_600dpi.png"), plot = dotplot, width = 10, height = 6, dpi = 600)

