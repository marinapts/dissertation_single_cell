library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(gridExtra)
library(optparse)
library(cowplot)


load_dataset = function(file, name) {
    print(paste0('Loading dataset ', file))
    gene_expression_matrix = read.csv(file, row.names=1)  # Load csv file
    seuratObj = CreateSeuratObject(counts=gene_expression_matrix, project=name)  # Create Seurat object
    print(cat(name, dim(seuratObj)))
    return(seuratObj)
}

quality_control_plots = function(data, name) {
    # Add percentage of mitochondrial genes to the Seurat object
    data[['percent.mt']] = PercentageFeatureSet(data, pattern = '^mt-')

    # Visualize QC metrics as a violin plot
    violin_plot = VlnPlot(data, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
    # Visualize relationship between nCount_RNA and nFeature_RNA
    scatter_plot_1 = FeatureScatter(data, feature1 = 'nCount_RNA', feature2 = 'percent.mt')
    scatter_plot_2 = FeatureScatter(data, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')
    scatter_plot = scatter_plot_1 + scatter_plot_2

    ggsave(file=paste0(plots_dir, name, '_qc_violin.pdf'), plot=violin_plot, height=30, width=20, units='cm')
    ggsave(file=paste0(plots_dir, name, '_qc_scatter.pdf'), plot=scatter_plot, height=20, width=30, units='cm')

    return(data)
}

feature_selection = function(data, name, n_features) {
    data = FindVariableFeatures(data, nfeatures = n_features)

    # Identify the 10 most highly variable genes
    top10 = head(VariableFeatures(data), 10)

    # plot variable features with and without labels
    plot1 = VariableFeaturePlot(data)
    plot2 = LabelPoints(plot = plot1, points = top10, repel = TRUE)
    feature_selection_plot = plot1 + plot2
    ggsave(file=paste0(plots_dir, name, '_feature_selection.pdf'), plot=feature_selection_plot, height=30, width=35, units='cm')
    return(data)
}

scaling = function(data) {
    # Regressing out nCount_RNA, percent.mt removes unwanted sources of variation.
    # Centering and scaling the data matrix
    all.genes = rownames(data)
    # data = ScaleData(data, features = all.genes)
    data = ScaleData(data, vars.to.regress = c('nCount_RNA', 'percent.mt'))
    return(data)
}

plot_PCA = function(E13, E14, name) {
    pca_13 = plot(FeaturePlot(object = E13, features = 'nCount_RNA'))
    pca_14 = plot(FeaturePlot(object = E14, features = 'nCount_RNA'))
    pca_plots = pca_13 + pca_14
    ggsave(file=paste0(plots_dir, name, '_pca.pdf'), plot=pca_plots, width=30, units='cm')
}

plot_heatmap = function(E13, E14, name) {
    heatmap_13 = plot(DimHeatmap(E13, dims = 1, cells = 500, balanced = TRUE))
    heatmap_14 = plot(DimHeatmap(E14, dims = 1, cells = 500, balanced = TRUE))
    heatmap = heatmap_13 + heatmap_14
    ggsave(file=paste0(plots_dir, name, '_heatmap.pdf'), plot=heatmap, width=30, units='cm')
}

elbow_plot = function(E13, E14, name, dimensions) {
    elbow_13 = ElbowPlot(E13, ndims=dimensions)
    elbow_14 = ElbowPlot(E13, ndims=dimensions)
    elbow = elbow_13 + elbow_14
    ggsave(file=paste0(plots_dir, name, '_elbow.pdf'), plot=elbow, width=30, units='cm')
}

run_umap = function(data, from_dimension=1, to_dimension) {
    print(from_dimension)
    data = FindNeighbors(data, dims = from_dimension:to_dimension)
    data = FindClusters(data, resolution = 0.5)
    # Look at cluster IDs of the first 5 cells
    # head(Idents(data), 5)
    data = RunUMAP(data, dims = from_dimension:to_dimension)
    return(data)
}

plot_umap_clustering = function(E13, E14, name, from_dimension, to_dimension) {
    E13 = run_umap(E13, from_dimension, to_dimension)
    E14 = run_umap(E14, from_dimension, to_dimension)
    umap_plot_13 = DimPlot(E13, reduction = 'umap', label = TRUE)
    umap_plot_14 = DimPlot(E14, reduction = 'umap', label = TRUE)
    umap_plot = umap_plot_13 + umap_plot_14
    ggsave(file=paste0(plots_dir, name, '_umap.pdf'), plot=umap_plot, width=30, units='cm')
    saveRDS(E13, file = paste0(rds_dir, 'umap_E13', '.rds'))
    saveRDS(E14, file = paste0(rds_dir, 'umap_E14', '.rds'))
}


# ==================================================================

# Set dataset to either het or hom
mouse_model = 'hom'
# mouse_model = 'het'
print(paste0('Mouse model: ', mouse_model))

E13_dataset_name = paste0('E13_', mouse_model)  # E13_hom or E13_het
E14_dataset_name = paste0('E14_', mouse_model)  # E14_hom or E14_het
# plots_dir = paste0('./plots_', mouse_model, '/')
# rds_dir = paste0('./rds_', mouse_model, '/')
plots_dir = paste0('./plots_', mouse_model, '_12_june_3000/')
rds_dir = paste0('./rds_', mouse_model, '_12_june_3000/')
dir.create(plots_dir)
dir.create(rds_dir)

# Load datasets
data_dir = file.path(dirname(getwd()), 'data')
print(data_dir)
E13_csv_file = file.path(data_dir, paste0(E13_dataset_name, '.csv'))
E14_csv_file = file.path(data_dir, paste0(E14_dataset_name, '.csv'))

E13 = load_dataset(E13_csv_file, E13_dataset_name)
E14 = load_dataset(E14_csv_file, E14_dataset_name)


# Quality Control
E13 = quality_control_plots(E13, E13_dataset_name)
E14 = quality_control_plots(E14, E14_dataset_name)

# Filter out cells found from quality control
E13_subset = subset(E13, subset = nFeature_RNA > 0 & nFeature_RNA < 6000 &
                                  nCount_RNA > 1000 & nCount_RNA < 30000 &
                                  percent.mt > 1 & percent.mt < 5)

E14_subset = subset(E14, subset = nFeature_RNA > 1000 & nFeature_RNA < 7000 &
                                  nCount_RNA > 100 & nCount_RNA < 30000 &
                                  percent.mt > 1 & percent.mt < 10)

print(cat('E13 subset: ', dim(E13_subset)))
print(cat('E14 subset: ', dim(E14_subset)))

# Merge E14_subset and E13_subset only for visualisation purposes - avoid 2 plots
E13_E14 = merge(E13_subset, E14_subset, add.cell.ids=c(E13_dataset_name, E14_dataset_name))
feature_scatter_plot = FeatureScatter(E13_E14, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')
ggsave(file=paste0(plots_dir, 'E13_E14_', mouse_model, '_feature_plot.pdf'), plot=feature_scatter_plot)


# LogNormalize normalizes the feature expression measurements for each cell by the total expression,
# multiplies this by a scale factor (10,000 by default), and log-transforms the result.
E13_normalised = NormalizeData(E13_subset, normalization.method = 'LogNormalize', scale.factor = 10000)  # Normalise data
E14_normalised = NormalizeData(E14_subset, normalization.method = 'LogNormalize', scale.factor = 10000)  # Normalise data

saveRDS(E13_normalised, file = paste0(rds_dir, 'E13_normalised.rds'))
saveRDS(E14_normalised, file = paste0(rds_dir, 'E14_normalised.rds'))

# Feature Selection - select highly variable genes for dimensionality reduction
E13_subset = feature_selection(E13_normalised, E13_dataset_name, 2000)
E14_subset = feature_selection(E14_normalised, E14_dataset_name, 2000)

# Scaling - 0 mean, 1 variance
E13_processed = scaling(E13_subset)
E14_processed = scaling(E14_subset)

# Dimensionality reduction
E13_processed = RunPCA(E13_processed, features = VariableFeatures(object = E13_processed))
E14_processed = RunPCA(E14_processed, features = VariableFeatures(object = E14_processed))

# Plot 2 principal components of PCA for both datasets side by side
plot_PCA(E13_processed, E14_processed, paste0('E13_E14_', mouse_model))

# plot_heatmap(E13_processed, E14_processed, paste0('E13_E14_', mouse_model))

# Define number of PCs from the elbow plot --> 20
elbow_plot(E13_processed, E14_processed, paste0('E13_E14_', mouse_model), 50)

# Clustering
plot_umap_clustering(E13_processed, E14_processed, paste0('E13_E14_', mouse_model, '_1_20'), 1, 20)
plot_umap_clustering(E13_processed, E14_processed, paste0('E13_E14_', mouse_model, '_10_20'), 10, 20)
