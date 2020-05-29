library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(gridExtra)


load_dataset <- function(file, name) {
    gene_expression_matrix <- read.csv(file, row.names=1)  # Load csv file
    seuratObj <- CreateSeuratObject(counts=gene_expression_matrix, project=name)  # Create Seurat object
    print(cat(name, dim(seuratObj)))
    return(seuratObj)
}

qc_plots <- function(data, name) {
    # Add percentage of mitochondrial genes to the Seurat object
    data[['percent.mito']] <- PercentageFeatureSet(data, pattern = '^MT-')
    # Visualize QC metrics as a violin plot
    violin_plot <- VlnPlot(data, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mito'), ncol = 3)
    # Visualize relationship between nCount_RNA and nFeature_RNA
    scatter_plot <- FeatureScatter(data, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')
    return(list(violin_plot, scatter_plot))
    # while (!is.null(dev.list()))  dev.off()
}

feature_selection <- function(data) {
    data <- FindVariableFeatures(data, selection.method = 'vst', nfeatures = 2000)

    # Identify the 10 most highly variable genes
    top10 <- head(VariableFeatures(data), 10)

    # plot variable features with and without labels
    plot1 <- VariableFeaturePlot(data)
    plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
    plot1 + plot2
    return(data)
}

scaling <- function(data) {
    all.genes <- rownames(data)
    data <- ScaleData(data, features = all.genes)
    return(data)
}

umap_clustering <- function(data, rds_name, dimensions) {
    data <- FindNeighbors(data, dims = 1:dimensions)
    data <- FindClusters(data, resolution = 0.5)
    # Look at cluster IDs of the first 5 cells
    head(Idents(data), 5)
    data <- RunUMAP(data, dims = 1:dimensions)
    umap_plot <- DimPlot(data, reduction = 'umap', label = TRUE)
    saveRDS(data, file = paste0('../plots/umap_', rds_name, '.rds'))
    return(umap_plot)
}


# Load datasets
data_dir <- file.path(dirname(getwd()), 'data')
print(data_dir)
E13_hom_file = file.path(data_dir, 'E13_hom.csv')
E14_hom_file = file.path(data_dir, 'E14_hom.csv')

######################
# Quality Control
######################

# E13_HOM
E13_hom <- load_dataset(E13_hom_file, 'E13_hom')
qc_plots(E13_hom, 'E13_hom')
E13_hom <- subset(E13_hom, subset=nFeature_RNA < 6000 & nCount_RNA < 35000)  # Filter out cells found from quality control
print(cat('E13_hom subset: ', dim(E13_hom)))

# E14_HOM
E14_hom <- load_dataset(E14_hom_file, 'E14_hom')
qc_plots(E14_hom, 'E14_hom')
E14_hom <- subset(E14_hom, subset=nFeature_RNA < 7500 & nCount_RNA < 45000)  # Filter out cells found from quality control
print(cat('E14_hom subset: ', dim(E14_hom)))

# Merge E14_hom and E13_hom only for visualisation purposes - avoid 2 plots
E13_E14 <- merge(E13_hom, E14_hom, add.cell.ids=c('E13_hom', 'E14_hom'))
pdf(file='feature_scatter_plot.pdf')
FeatureScatter(E13_E14, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')


# LogNormalize normalizes the feature expression measurements for each cell by the total expression,
# multiplies this by a scale factor (10,000 by default), and log-transforms the result.
E13_hom <- NormalizeData(E13_hom, normalization.method = 'LogNormalize', scale.factor = 10000)  # Normalise data
E14_hom <- NormalizeData(E14_hom, normalization.method = 'LogNormalize', scale.factor = 10000)  # Normalise data


# Feature Selection
E13_hom <- feature_selection(E13_hom)
E14_hom <- feature_selection(E14_hom)

# Scaling - 0 mean, 1 variance
E13_hom <- scaling(E13_hom)
E14_hom <- scaling(E14_hom)

# Dimensionality reduction
E13_hom <- RunPCA(E13_hom, features = VariableFeatures(object = E13_hom))
E14_hom <- RunPCA(E14_hom, features = VariableFeatures(object = E14_hom))

# Plot 2 principal components of PCA for both datasets side by side
par(mfrow=c(1,2))
pca_13_hom <- plot(DimPlot(E13_hom, reduction = 'pca'))
pca_14_hom <- plot(DimPlot(E14_hom, reduction = 'pca'))
grid.arrange(pca_13_hom, pca_14_hom, nrow=1)

# Plot heatmap
par(mfrow=c(1,2))
heatmap_13_hom <- DimHeatmap(E13_hom, dims = 1, cells = 500, balanced = TRUE)
heatmap_14_hom <- DimHeatmap(E14_hom, dims = 1, cells = 500, balanced = TRUE)
grid.arrange(heatmap_13_hom, heatmap_14_hom, nrow=1)

# Define number of principal components
JackStraw(E13_hom, num.replicate = 100)
ScoreJackStraw(E13_hom, dims = 1:20)

JackStraw(E14_hom, num.replicate = 100)
ScoreJackStraw(E14_hom, dims = 1:20)

par(mfrow=c(1,2))
elbow_13_hom <- ElbowPlot(E13_hom)
elbow_14_hom <- ElbowPlot(E14_hom)
grid.arrange(elbow_13_hom, elbow_14_hom, nrow=1)


# Clustering
umap_clustering(E13_hom, 'E13_hom', 15)
umap_clustering(E14_hom, 'E14_hom', 15)
