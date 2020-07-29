import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors


def get_colormap():
    # Define a nice colour map for gene expression
    colors2 = plt.cm.Reds(np.linspace(0, 1, 128))
    colors3 = plt.cm.Greys_r(np.linspace(0.7, 0.8, 20))
    colorsComb = np.vstack([colors3, colors2])
    mymap = colors.LinearSegmentedColormap.from_list('my_colormap', colorsComb)
    return mymap


def get_known_marker_genes(adata):
    marker_genes = dict()
    marker_genes['neural_progen'] = ['Pax6', 'Vim', 'Sox2']
    marker_genes['intermediate_progen'] = ['Eomes', 'Btg2']
    marker_genes['post_mitotic'] = ['Tbr1', 'Sox5']
    marker_genes['ectopic'] = ['Gsx2', 'Prdm13', 'Dlx1', 'Dlx2', 'Dlx5', 'Gad1', 'Gad2', 'Ptf1a', 'Msx3', 'Helt', 'Olig3']

    main_cell_types = marker_genes['neural_progen'] + marker_genes['intermediate_progen'] + marker_genes['post_mitotic']
    var_names = set(adata.var_names)
    columns = set(adata.obs.columns)
    gene_names = var_names.union(columns)
    available_ectopic = gene_names.intersection(marker_genes['ectopic'])

    return marker_genes, main_cell_types, available_ectopic
