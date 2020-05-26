library(dplyr)
library(Seurat)
library(patchwork)
# require(ggplot2)

load_dataset <- function(file, name) {
    gene_expression_matrix <- read.csv(file, row.names=1)  # Load csv file
    seuratObj <- CreateSeuratObject(counts=gene_expression_matrix, project=name)  # Create Seurat object
    print(cat(name, dim(seuratObj)))
    return(seuratObj)
}

qc_plots <- function(data, name) {
    dev.off()
    # Add percentage of mitochondrial genes to the Seurat object
    data[["percent.mito"]] <- PercentageFeatureSet(data, pattern = "^MT-")
    # Visualize QC metrics as a violin plot
    violin_plot <- VlnPlot(data, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mito'), ncol = 3)

    # Visualize relationship between nCount_RNA and nFeature_RNA
    scatter <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

    # Remove values from the datasets
    # data_subset <- subset(data, subset=nFeature_RNA < 7500 & nCount_RNA < 50000)
}

data_dir <- file.path(dirname(getwd()), 'data')
print(data_dir)
E14_hom_file = file.path(data_dir, 'E14_hom.csv')
E13_hom_file = file.path(data_dir, 'E13_hom.csv')

# E14_HOM
E14_hom <- load_dataset(E14_hom_file, 'E14_hom')
qc_plots(E14_hom, 'E14_hom')
E14_hom <- subset(E14_hom, subset=nFeature_RNA < 7500 & nCount_RNA < 50000)  # Filter out cells found from quality control
print(cat('E14_hom subset: ', dim(E14_hom)))
# LogNormalize normalizes the feature expression measurements for each cell by the total expression,
# multiplies this by a scale factor (10,000 by default), and log-transforms the result.
E14_hom <- NormalizeData(E14_hom, normalization.method = "LogNormalize", scale.factor = 10000)  # Normalise data


# E13_HOM
E13_hom <- load_dataset(E13_hom_file, 'E13_hom')
qc_plots(E13_hom, 'E13_hom')
E13_hom <- subset(E13_hom, subset=nFeature_RNA < 6000 & nCount_RNA < 35000)  # Filter out cells found from quality control
print(cat('E13_hom subset: ', dim(E13_hom)))


# Merge E14_hom and E13_hom only for visualisation purposes - avoid 2 plots
E13_E14 <- merge(E13_hom, E14_hom, add.cell.ids=c('E13_hom', 'E14_hom'))
scatter <- FeatureScatter(E13_E14, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
