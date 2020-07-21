library(Seurat)
library(tensorflow)
library(cellassign)
library(pheatmap)
library(here)
library(scater)
library(rhdf5)
library(hdf5r)
library(SeuratDisk)
library(scran)

set.seed(1234)                # set R rng seed
reticulate::py_set_seed(1234) # set python rng seed

data_dir <- here('ann_data/exp_04')

dataset <- 'E13_hom_norm_all_genes'
h5ad_file <- file.path(data_dir, paste0(dataset, '.h5ad'))
Convert(h5ad_file, dest = 'h5seurat', overwrite = TRUE)
E14_hom <- LoadH5Seurat(file.path(data_dir, paste0(dataset, '.h5seurat')))
E14_hom <- as.SingleCellExperiment(E14_hom, assay = 'RNA')  # Convert to SingleCellExperiment class
E14_hom <- computeSumFactors(E14_hom)
sfactors <- sizeFactors(E14_hom)

marker_list <- list(
    `Neural progenitors` = c('Pax6', 'Vim', 'Sox2'),
    `Intermediate progenitors` = c('Eomes', 'Btg2'),
    `Post-mitotic neurons` = c('Tbr1', 'Sox5'),
    # `Ectopic cells` = c('Gsx2', 'Prdm13', 'Dlx1', 'Dlx2', 'Dlx5', 'Gad1', 'Gad2', 'Ptf1a', 'Msx3', 'Helt', 'Olig3')
    `Ectopic cells` = c('Gsx2', 'Prdm13', 'Dlx1', 'Dlx2', 'Dlx5', 'Gad1', 'Gad2', 'Ptf1a', 'Msx3', 'Helt')
)

# Construct binary gene-by-celltype marker gene matrix
mgi <- marker_list_to_mat(marker_list, include_other = FALSE)
pheatmap(mgi)


# Run cellassign
fit_E14_hom <- cellassign(exprs_obj = E14_hom[rownames(mgi),],
                          marker_gene_info = mgi,
                          s = sfactors,
                          # s = colSums(SummarizedExperiment::assay(E14_hom, "counts")),
                          learning_rate = 1e-2,
                          shrinkage = TRUE,
                          verbose = TRUE,
                          num_runs = 10)

E14_hom$cell_type <- fit_E14_hom$cell_type
plotUMAP(E14_hom, colour_by = 'cell_type')

# Convert SCE object to a dataframe and save UMIs and cell types to a csv file
df = data.frame(umi = colnames(E14_hom), cell_type = E14_hom[['cell_type']])
csv_anno = file.path(data_dir, 'annotations/E13_hom_cellassign.csv')
write.csv(df, file=csv_anno, row.names = F)
