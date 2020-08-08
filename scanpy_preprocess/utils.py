import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import seaborn as sns
from matplotlib import colors
from pathlib import Path


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


def get_updated_marker_genes(adata):
    """Known marker genes + the ones that have high correlation with the known (>0.5 corr.) - found by common_plots.py
    """
    marker_genes = dict()
    marker_genes['Neural Progenitors'] = ['Arx', 'Ccnd2', 'Dbi', 'Fabp7', 'Mfge8', 'Pax6', 'Sox2', 'Tox3', 'Vim', 'Zfp36l1']
    marker_genes['Intermediate Progenitors'] = ['Btg2', 'Cdh13', 'Chl1', 'Cnr1', 'Mpped1', 'Ppp2r2b', 'Sox5', 'Tbr1']
    marker_genes['Post-mitotic Neurons'] = ['Arpp21', 'Cdh13', 'Dab1', 'Fam49a', 'Gsx2', 'Mpped1', 'Neurod2',
                                            'Neurod6', 'Ppp2r2b', 'Rabgap1l', 'Sox5', 'Tbr1', 'Tnik']
    marker_genes['Ectopic'] = ['Arx', 'Ccdc88a', 'Ccnd2', 'Cdca7', 'Cdk14', 'Dlx1', 'Dlx2', 'Dlx5', 'Dlx6', 'Dlx6os1',
                               'Etv1', 'Gad1', 'Gad2', 'Gm13889', 'Gsx2', 'Helt', 'Msx3', 'Nlk', 'Nrxn3', 'Pfn2',
                               'Prdm13', 'Rnd3', 'Slain1', 'Sp8', 'Sp9', 'Tiam2']

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


def marker_gene_overlap(adata, marker_genes, fig_dir, dataset_name, key_added='rank_genes_groups', updated=False):
    """Calculates the overlap of the known marker_genes with the DE top_n_genes (+ some plots)
    Args:
        adata (AnnData)
        marker_genes (dict): Dictionary of known marker genes
        key_added (str, optional): Key in adata.uns to find the ranked DE genes
        updated (bool, optional): If the method is called again - for plotting purposes

    Returns:
        (pd.DataFrame): A table with clusters as columns and cell types as rows
    """
    gene_overlap_norm = sc.tl.marker_gene_overlap(adata, reference_markers=marker_genes,
                                                  key=key_added, normalize='reference')
    print(gene_overlap_norm)

    # Plot the overlap on a heatmap
    if updated == True:
        fig_title = 'gene_overlap_heatmap_updated'
        rotation = 25
        size = (14, 9)
    else:
        fig_title = 'gene_overlap_heatmap_'
        rotation = 0
        size = (14, 6)

    plt.figure(figsize=size)
    # ax0 = plt.subplot(111)
    # ax1 = sns.heatmap(gene_overlap_norm, cbar=True, annot=True, square=True)
    ax1 = sns.heatmap(gene_overlap_norm, cbar=True, annot=True, square=False)

    ax1.set_xticklabels(ax1.get_xticklabels(), rotation=rotation)
    ax1.set_title('Marker gene overlap heatmap')
    fig = ax1.get_figure()
    # fig.set_size_inches(12, 8)
    fig.savefig(Path(fig_dir, fig_title + dataset_name + '.eps'))

    return gene_overlap_norm


def annotate_clusters_based_on_overlap(gene_overlap_norm, adata, clusters_key='leiden'):
    """Creates annotations for each cluster based on the gene overlap
    Args:
        gene_overlap_norm (pd.DataFrame): Overlap of known marker genes with DE genes
        clusters_key (str, optional): Key in adata.obs where clusters are stored
    Returns:
        (list): The new annotation names for each cluster
    """
    cluster_annotations = dict()

    # Specify a cell type for each cluster
    for cluster in gene_overlap_norm.columns:
        overlaps = gene_overlap_norm.loc[:, cluster]
        max_overlap = max(overlaps.values)

        if max_overlap == 0:
            cluster_annotations[cluster] = 'Unknown'
        else:
            cluster_annotations[cluster] = overlaps.idxmax()
    print('cluster_annotations', cluster_annotations)

    # Map each cluster to its cell type
    new_cluster_names = list()
    for annotation in adata.obs[clusters_key]:
        cell_type_name = cluster_annotations[str(annotation)]
        new_cluster_names.append(cell_type_name)

    return new_cluster_names
