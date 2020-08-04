import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import argparse
from pathlib import Path
from utils import get_colormap, get_known_marker_genes, probability_distr_of_overlap

sc.logging.print_versions()
sc.set_figure_params(facecolor="white", figsize=(8, 4))
sc.settings.verbosity = 3
np.random.seed(2211)

DATASET_NAME = ''
FIGDIR = './figures/marker_genes/'

sc.settings.figdir = FIGDIR
sc.settings.file_format_figs = 'eps'
sc.settings._vector_friendly = False
sc.settings.autosave = False
sc.settings.autoshow = False

vmin = -5
vmax = 5
use_raw = False


cell_type_mapping = {
    'neural_progen': 'Neural progenitors',
    'intermediate_progen': 'Intermediate progenitors',
    'post_mitotic': 'Post-mitotic neurons',
    'ectopic': 'Ectopic cells',
    'unknown': 'Unknown'
}


def plot_clusters(adata, clusters_key, is_integrated_dataset, main_cell_types, available_ectopic):
    sc.pl.umap(adata, color=[clusters_key], legend_loc='on data', save=DATASET_NAME)
    if is_integrated_dataset is True:
        sc.pl.umap(adata, color=['batch'], save='_batch_' + DATASET_NAME)

    # Plot expression of marker genes
    sc.pl.umap(adata, color=main_cell_types + list(available_ectopic), cmap=get_colormap(),
               vmin=-5, vmax=5, use_raw=False, save='_markergenes_expr_' + DATASET_NAME)


def differential_expression(adata, clusters_key, n_genes, DEGs_file, updated):
    fig_title = '_updated_' + DATASET_NAME if updated else DATASET_NAME
    # Calculate marker genes
    sc.tl.rank_genes_groups(adata, clusters_key, n_genes=n_genes, method='wilcoxon', use_raw=True)  # 100 genes by default
    sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False, use_raw=True, save=fig_title)

    top_ranked_genes_per_cluster = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
    print(top_ranked_genes_per_cluster.head(20))

    # Write top_n_genes for each cluster in a csv file
    if DEGs_file is not None:
        top_ranked_genes_per_cluster.to_csv(DEGs_file)

    sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, use_raw=use_raw, swap_axes=True, vmin=vmin, vmax=vmax,
                                    cmap='bwr', show_gene_labels=True, var_group_rotation=45, save=fig_title)
    sc.tl.dendrogram(adata, groupby=clusters_key, use_rep='X_umap', var_names=marker_genes, use_raw=use_raw)
    sc.pl.dendrogram(adata, groupby=clusters_key, save=fig_title)
    sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, save=fig_title)

    return adata.uns['rank_genes_groups']


def marker_gene_overlap(adata, ranked_genes, updated):
    gene_overlap_norm = sc.tl.marker_gene_overlap(adata, reference_markers=marker_genes,
                                                  key='rank_genes_groups', normalize='reference')
    print(gene_overlap_norm)

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
    fig.savefig(FIGDIR + fig_title + DATASET_NAME + '.eps')

    return gene_overlap_norm


def annotate_clusters_based_on_overlap(gene_overlap_norm):
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
    for leiden in adata.obs[clusters_key]:
        # cell_type_name = cell_type_mapping[cluster_annotations[leiden]]
        cell_type_name = cluster_annotations[leiden]
#         print(cell_type_name)
        new_cluster_names.append(cell_type_name)

    return new_cluster_names


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Preprocessing of single-cell RNA-seq datasets')
    parser.add_argument('--dataset', nargs='+', help='One or more space separated datasets, e.g. E14_hom E13_hom integration')
    parser.add_argument('--dataset_type', type=str, default='variable', help='One of: ["variable", "all"]')
    parser.add_argument('--update_file', action='store_true', help='Update h5ad file')
    args = parser.parse_args()

    alldata = {}

    for dataset in args.dataset:
        DATASET_NAME = dataset

        is_integrated_dataset = True if dataset == 'integration' else False
        if is_integrated_dataset:
            data_filename = dataset + '_HOMs_' + args.dataset_type
        else:
            data_filename = dataset + '_' + args.dataset_type + '_genes'

        filepath = Path('ann_data/', data_filename + '.h5ad')
        adata = sc.read(filepath)
        print(adata)
        clusters_key = 'integr_clusters' if is_integrated_dataset else 'leiden'

        marker_genes, main_cell_types, available_ectopic = get_known_marker_genes(adata)
        if dataset == 'E13_hom' and args.dataset_type == 'variable':
            marker_genes['Neural Progenitors'].remove('Vim')
            main_cell_types.remove('Vim')

        # Show initial clusters and expression of known marker genes - nothing new
        plot_clusters(adata, clusters_key, is_integrated_dataset, main_cell_types, available_ectopic)

        top_n_genes = 100
        DEGs_filename = dataset + '_' + args.dataset_type + '_' + str(top_n_genes)
        DEGs_file = Path('DEGs', DEGs_filename + '.csv')

        # DE to calculate marker genes for the clusters
        ranked_genes = differential_expression(adata, clusters_key, top_n_genes, DEGs_file, updated=False)

        # Write marker genes in a csv for annotation by SCSA
        groups = ranked_genes['names'].dtype.names
        dat = pd.DataFrame({group + '_' + key[:1]: ranked_genes[key][group] for group in groups for key in ['names', 'logfoldchanges', 'scores', 'pvals']})
        dat.to_csv('cellmarker/DE_' + data_filename + '.csv')

        # Get overlap of known marker genes with the marker genes found from DE
        gene_overlap_norm = marker_gene_overlap(adata, ranked_genes, updated=False)
        gene_overlap_norm_distr = probability_distr_of_overlap(gene_overlap_norm)

        # Annotate clusters
        adata.obs[clusters_key + '_annotations'] = annotate_clusters_based_on_overlap(gene_overlap_norm)
        sc.pl.umap(adata, color=[clusters_key + '_annotations'], save='_annotations_' + DATASET_NAME)

        # DE to find marker genes for the updated clusters (cell types)
        updated_ranked_genes = differential_expression(adata, clusters_key + '_annotations', top_n_genes, None, updated=True)
        updated_gene_overlap_norm = marker_gene_overlap(adata, updated_ranked_genes, updated=True)

        if args.update_file:
            adata.write(filepath)
            print('{} file updated'.format(filepath))
