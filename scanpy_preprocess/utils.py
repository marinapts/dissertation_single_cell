import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors


def get_colormap(color='red'):
    # Define a nice colour map for gene expression
    if color == 'red':
        colors2 = plt.cm.Reds(np.linspace(0, 1, 128))
    elif color == 'blue':
        colors2 = plt.cm.Blues(np.linspace(0, 1, 128))
    colors3 = plt.cm.Greys_r(np.linspace(0.7, 0.8, 20))
    colorsComb = np.vstack([colors3, colors2])
    mymap = colors.LinearSegmentedColormap.from_list('my_colormap', colorsComb)
    return mymap


def get_known_marker_genes(adata):
    marker_genes = dict()
    # marker_genes['neural_progen'] = ['Pax6', 'Vim', 'Sox2']
    # marker_genes['intermediate_progen'] = ['Eomes', 'Btg2']
    # marker_genes['post_mitotic'] = ['Tbr1', 'Sox5']
    # marker_genes['ectopic'] = ['Gsx2', 'Prdm13', 'Dlx1', 'Dlx2', 'Dlx5', 'Gad1', 'Gad2', 'Ptf1a', 'Msx3', 'Helt', 'Olig3']

    # main_cell_types = marker_genes['neural_progen'] + marker_genes['intermediate_progen'] + marker_genes['post_mitotic']

    marker_genes['Neural Progenitors'] = ['Pax6', 'Vim', 'Sox2']
    marker_genes['Intermediate Progenitors'] = ['Eomes', 'Btg2']
    marker_genes['Post-mitotic Neurons'] = ['Tbr1', 'Sox5']
    marker_genes['Ectopic'] = ['Gsx2', 'Prdm13', 'Dlx1', 'Dlx2', 'Dlx5', 'Gad1', 'Gad2', 'Ptf1a', 'Msx3', 'Helt', 'Olig3']
    main_cell_types = marker_genes['Neural Progenitors'] + marker_genes['Intermediate Progenitors'] + marker_genes['Post-mitotic Neurons']

    var_names = set(adata.var_names)
    columns = set(adata.obs.columns)
    gene_names = var_names.union(columns)
    available_ectopic = gene_names.intersection(marker_genes['Ectopic'])

    return marker_genes, main_cell_types, available_ectopic


def probability_distr_of_overlap(gene_overlap_norm):
    """Converts the gene overlap matrix of known marker genes with the DE genes
    so that each cluster is a probability distribution over the cell types.
    Each value is the percentage of cells in that cluster belong to a cell type (assumption based on the overlap)
    """
    gene_overlap_norm_distr = gene_overlap_norm.copy()
    for cl in gene_overlap_norm.columns:
        gene_overlap_norm_distr[cl] = np.nan_to_num(gene_overlap_norm[cl] / sum(gene_overlap_norm[cl]))

    print('Prob. distribution of marker gene overlap:')
    print(gene_overlap_norm_distr)

    return gene_overlap_norm_distr
