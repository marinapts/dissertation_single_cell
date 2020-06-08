# ========================================================================================
#
# Remove small clusters containing cajal-retzius cells that are visible in the umap plots
#
# ========================================================================================

save_umap_without_cajal_retzius = function(E13, E14, name) {
    umap_E13 = DimPlot(E13, reduction='umap', label=TRUE)
    umap_E14 = DimPlot(E14, reduction='umap', label=TRUE)
    umap_plot = umap_E13 + umap_E14
    saveRDS(E13, file = paste0(rds_dir, 'umap_E13_no_cajal_retzius', '.rds'))
    saveRDS(E14, file = paste0(rds_dir, 'umap_E14_no_cajal_retzius', '.rds'))
    ggsave(file=paste0(plots_dir, name, '_umap_no_cajal_retzius.pdf'), plot=umap_plot, width=25, units='cm')
}

set.seed(42)

# Load either het or hom datasets
# mouse_model = 'hom'
mouse_model = 'het'
print(paste0('Mouse model: ', mouse_model))
plots_dir = paste0('./plots_', mouse_model, '/')
rds_dir = paste0('./rds_', mouse_model, '/')


# To load the rds seurat files directly, run only the following:
E13_final = readRDS(file = paste0(rds_dir, 'umap_E13.rds'))
E14_final = readRDS(file = paste0(rds_dir, 'umap_E14.rds'))

# How many cells are in each cluster
table(Idents(E13_final))
table(Idents(E14_final))
# What proportion of cells are in each cluster?
prop.table(table(Idents(E13_final)))
prop.table(table(Idents(E14_final)))

if (mouse_model == 'hom') {
    cajal_retzius_cluster_E13 = '9'
    cajal_retzius_cluster_E14 = '10'
} else {
    cajal_retzius_cluster_E13 = '8'
    cajal_retzius_cluster_E14 = '8'
}

# What are the cell names of all cajal-retzius cells?
WhichCells(E13_final, idents = cajal_retzius_cluster_E13)
WhichCells(E14_final, idents = cajal_retzius_cluster_E14)

# Create a Seurat object of all cells except these cajal-retzius cells
E13_final_subset = subset(E13_final, idents = c(cajal_retzius_cluster_E13), invert = TRUE)
E14_final_subset = subset(E14_final, idents = c(cajal_retzius_cluster_E14), invert = TRUE)

print(dim(E13_final_subset))
print(dim(E14_final_subset))

save_umap_without_cajal_retzius(E13_final_subset, E14_final_subset, paste0('E13_E14_', mouse_model))

# Embeddings(object = object[["umap"]])
# umap_embeddings = Embeddings(object = E13_final_subset[["umap"]])

# umap_embeddings = as.data.frame(E13_final_subset@reductions$umap@cell.embeddings, stringsAsFactors = F)
# # E13_after_deletion = umap_embeddings[umap_embeddings$UMAP_2 < 9,]
# cells_to_remove = umap_embeddings[umap_embeddings$UMAP_2 > 9,]
# cell_names = rownames(cells_to_remove)
# # final = subset(E13_final_subset, cells = cell_names)
# final = subset(E13_final_subset, cells = -c(cell_names))


# drop <- c(cell_names)
# for(i in columns.to.remove) {
#     E13_final_subset[[i]] <- NULL
# }
