library(dplyr)
library(Seurat)
library(patchwork)

data_dir <- file.path(dirname(getwd()), 'data')
E14_het_file = file.path(data_dir, 'E14_het.csv')
E14_hom_file = file.path(data_dir, 'E14_hom.csv')

# Load the dataset
E14_het_matrix <- read.csv(E14_het_file, row.names=1)
E14_hom_matrix <- read.csv(E14_hom_file, row.names=1)
print(cat('E14 het:', dim(E14_het_matrix)))
print(cat('E14 hom:', dim(E14_hom_matrix)))

# Create Seurat object for each dataset
E14_het <- CreateSeuratObject(counts=E14_het_matrix, project="E14_het", min.cells=3, min.features=200)
E14_hom <- CreateSeuratObject(counts=E14_hom_matrix, project="E14_hom", min.cells=3, min.features=200)

# Merge into one single seurat object. Add cell ids in case you have overlapping barcodes between the datasets
alldata <- merge(E14_het, E14_hom, add.cell.ids=c('E14_het', 'E14_hom'))
table(Idents(alldata))  # Checks number of cells per each sample

# Add percentage of mitochondrial genes
E14_het[["percent.mito"]] <- PercentageFeatureSet(E14_het, pattern = "^MT-")
alldata[["percent.mito"]] <- PercentageFeatureSet(alldata, pattern = "^MT-")
# Show QC metrics for the first 5 cells
head(alldata@meta.data, 10)

VlnPlot(E14_het, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mito'), ncol = 3)
VlnPlot(alldata, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mito'), ncol = 3)
# Visualize relationship between nCount_RNA and nFeature_RNA
FeatureScatter(E14_het, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(alldata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Remove values from the datasets
E14_het <- subset(E14_het, subset=nFeature_RNA > 6500 & nCount_RNA > 40000)
E14_hom <- subset(E14_hom, subset=nFeature_RNA > 7500 & nCount_RNA > 45000)
print(cat('E14 het:', dim(E14_het)))
print(cat('E14 hom:', dim(E14_hom)))
