import scanpy as sc

E13_all = sc.read('ann_data/E13_hom_all_genes.h5ad')
E13_var = sc.read('ann_data/E13_hom_variable_genes.h5ad')

# When removing Cajal-Retzius cells, there was 1 extra cell in E13_all which we need to remove,
# because the 2 objects should have the same number of cells
rownames_intersection = set(E13_all.obs_names).intersection(set(E13_var.obs_names))
E13_all = E13_all[E13_all.obs_names.isin(rownames_intersection), :]

E13_all.write('ann_data/E13_hom_all_genes.h5ad')
