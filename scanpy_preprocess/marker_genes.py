import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import umap
import seaborn as sb
import json
from matplotlib import colors
from pathlib import Path


sc.logging.print_versions()
sc.set_figure_params(facecolor="white", figsize=(6, 4))
sc.settings.verbosity = 3

np.random.seed(2211)

# Define a nice colour map for gene expression
colors2 = plt.cm.Reds(np.linspace(0, 1, 128))
colors3 = plt.cm.Greys_r(np.linspace(0.7, 0.8, 20))
colorsComb = np.vstack([colors3, colors2])
mymap = colors.LinearSegmentedColormap.from_list('my_colormap', colorsComb)



if __name__ == '__main__':
    # dataset = 'E13_hom'
    # adata = sc.read('../ann_data/exp_04/' + dataset + '_norm_variable_genes.h5ad')
    dataset = 'integration_HOMs.h5ad'
    adata = sc.read(Path('../ann_data/', dataset))

    # Construct dictionary of known marker genes
    marker_genes = dict()
    marker_genes['neural_progen'] = ['Pax6', 'Vim', 'Sox2']
    marker_genes['intermediate_progen'] = ['Eomes', 'Btg2']
    marker_genes['post_mitotic'] = ['Tbr1', 'Sox5']
    marker_genes['ectopic'] = ['Gsx2', 'Prdm13', 'Dlx1', 'Dlx2', 'Dlx5', 'Gad1',
                               'Gad2', 'Ptf1a', 'Msx3', 'Helt', 'Olig3']

    var_names = set(adata.var_names)
    columns = set(adata.obs.columns)
    gene_names = var_names.union(columns)
    available_ectopic = gene_names.intersection(marker_genes['ectopic'])
