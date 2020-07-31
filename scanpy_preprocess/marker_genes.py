import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sb
import argparse
from pathlib import Path
from utils import get_colormap, get_known_marker_genes, probability_distr_of_overlap

sc.logging.print_versions()
sc.set_figure_params(facecolor="white", figsize=(8, 4))
sc.settings.verbosity = 3
np.random.seed(2211)

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
    sc.pl.umap(adata, color=[clusters_key], legend_loc='on data')
    if is_integrated_dataset is True:
        sc.pl.umap(adata, color=['batch'])

    # Plot expression of marker genes
    sc.pl.umap(adata, color=main_cell_types + list(available_ectopic), cmap=get_colormap(),
               vmin=-5, vmax=5, use_raw=False)


def differential_expression(adata, clusters_key, n_genes, show_extra_de_plots, DEGs_file):
    # Calculate marker genes
    sc.tl.rank_genes_groups(adata, clusters_key, n_genes=n_genes, method='wilcoxon', use_raw=True)  # 100 genes by default
    # sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False, use_raw=True)

    top_ranked_genes_per_cluster = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
    print(top_ranked_genes_per_cluster.head(20))

    # Write top_n_genes for each cluster in a csv file
    if DEGs_file is not None:
        top_ranked_genes_per_cluster.to_csv(DEGs_file)

    if show_extra_de_plots:
        sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, use_raw=use_raw, swap_axes=True,
                                        vmin=vmin, vmax=vmax, cmap='bwr', show_gene_labels=True)
        sc.tl.dendrogram(adata, groupby=clusters_key, use_rep='X_umap', var_names=marker_genes, use_raw=use_raw)
        sc.pl.dendrogram(adata, groupby=clusters_key)
        sc.pl.rank_genes_groups_dotplot(adata, n_genes=5)

    return adata.uns['rank_genes_groups']


def marker_gene_overlap(adata, ranked_genes):
    fig, ax = plt.subplots(figsize=(10, 3))
    gene_overlap_norm = sc.tl.marker_gene_overlap(adata, reference_markers=marker_genes,
                                                  key='rank_genes_groups', normalize='reference')
    print(gene_overlap_norm)
    sb.heatmap(gene_overlap_norm, cbar=True, annot=True)
    return gene_overlap_norm


def annotate_clusters_based_on_overlap(gene_overlap_norm):
    cluster_annotations = dict()

    # Specify a cell type for each cluster
    for cluster in gene_overlap_norm.columns:
        overlaps = gene_overlap_norm.loc[:, cluster]
        max_overlap = max(overlaps.values)

        if max_overlap == 0:
            cluster_annotations[cluster] = 'unknown'
        else:
            cluster_annotations[cluster] = overlaps.idxmax()
    print('cluster_annotations', cluster_annotations)

    # Map each cluster to its cell type
    new_cluster_names = list()
    for leiden in adata.obs[clusters_key]:
        cell_type_name = cell_type_mapping[cluster_annotations[leiden]]
#         print(cell_type_name)
        new_cluster_names.append(cell_type_name)

    return new_cluster_names


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Preprocessing of single-cell RNA-seq datasets')
    parser.add_argument('--dataset', type=str, default='E14_hom', help='One of: ["E13_hom", "E14_hom", "integration"]')
    parser.add_argument('--dataset_type', type=str, default='variable', help='One of: ["variable", "all"]')
    parser.add_argument('--plot_initial_clusters', action='store_true', help='Whether to plot the initial clusters')
    parser.add_argument('--plot_diff_expression', action='store_true', help='Whether to plot the results from DE')
    parser.add_argument('--update_file', action='store_true', help='Update h5ad file')
    args = parser.parse_args()

    is_integrated_dataset = True if args.dataset == 'integration' else False
    if is_integrated_dataset:
        data_filename = args.dataset + '_HOMs_' + args.dataset_type
    else:
        data_filename = args.dataset + '_' + args.dataset_type + '_genes'

    filepath = Path('ann_data/', data_filename + '.h5ad')
    adata = sc.read(filepath)
    print(adata)
    clusters_key = 'integr_clusters' if is_integrated_dataset else 'leiden'

    marker_genes, main_cell_types, available_ectopic = get_known_marker_genes(adata)
    if args.dataset == 'E13_hom' and args.dataset_type == 'variable':
        marker_genes['neural_progen'].remove('Vim')
        main_cell_types.remove('Vim')

    # Show clustering
    if args.plot_initial_clusters:
        plot_clusters(adata, clusters_key, is_integrated_dataset, main_cell_types, available_ectopic)

    top_n_genes = 100
    DEGs_filename = args.dataset + '_' + args.dataset_type + '_' + str(top_n_genes)
    DEGs_file = Path('DEGs', DEGs_filename + '.csv')

    # DE to calculate marker genes for the clusters
    ranked_genes = differential_expression(adata, clusters_key, top_n_genes, args.plot_diff_expression, DEGs_file)
    # Get overlap of known marker genes with the marker genes found from DE
    gene_overlap_norm = marker_gene_overlap(adata, ranked_genes)
    gene_overlap_norm_distr = probability_distr_of_overlap(gene_overlap_norm)

    # Annotate clusters
    adata.obs[clusters_key + '_annotations'] = annotate_clusters_based_on_overlap(gene_overlap_norm)
    sc.pl.umap(adata, color=[clusters_key + '_annotations'])

    # DE to find marker genes for the updated clusters (cell types)
    updated_ranked_genes = differential_expression(adata, clusters_key + '_annotations', top_n_genes,
                                                   args.plot_diff_expression, None)
    updated_gene_overlap_norm = marker_gene_overlap(adata, updated_ranked_genes)

    if args.update_file:
        adata.write(filepath)
        print('{} file updated'.format(filepath))
