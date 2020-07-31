import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import umap
import scanpy as sc
import argparse
import hdbscan
import seaborn as sb
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score, silhouette_score
from sklearn.cluster import KMeans
from pathlib import Path

from scanpy_preprocess.utils import get_colormap, get_known_marker_genes, probability_distr_of_overlap

sc.logging.print_version_and_date()
sns.set(style='white', context='notebook', rc={'figure.figsize': (9, 5)})
MYSEED = 2211
np.random.seed(MYSEED)

vmin = -5
vmax = 5
use_raw = False

label2num = {
    'Neural progenitors': 1,
    'Intermediate progenitors': 2,
    'Post-mitotic neurons': 3,
    'Ectopic cells': 4,
    'Unknown': 5
}
num2label = {v: k for k, v in label2num.items()}  # reverse mapping


def get_all_data(dataset):
    adata = sc.read(Path('ann_data', dataset + '_genes.h5ad'))
    bottleneck = pd.read_csv(Path('scDeepCluster/bottleneck', dataset + '_bottleneck.csv'), header=None)

    return adata, bottleneck


def fit_umap(data, n_neighbors=15, min_dist=0.5, metric='euclidean'):
    np.random.seed(42)
    umap_fit = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist, metric=metric, random_state=42)
    embedding = umap_fit.fit_transform(data)
    return embedding


def plot_clustering(embedding, labels, cmap='tab10', title='Bottleneck clusters'):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    scatter = ax.scatter(embedding[:, 0], embedding[:, 1], c=labels, cmap=cmap, marker='.')

    legend1 = ax.legend(*scatter.legend_elements(), title='Clusters')
    ax.add_artist(legend1)
    plt.title(title, fontsize=18)
    plt.show()


def plot_marker_genes(embedding, orig_data, labels, n_cols=3, title='', cmap='tab10'):

    len_labels = len(labels)
    fig_size = (15, 70)
    fig = plt.figure(1)
    n_rows = len_labels // n_cols
    # n_rows += len_labels % n_cols

    # Create subplots
    subplot_idx = 1
    for label in labels:
        if label in orig_data.var_names:
            ax = fig.add_subplot(n_rows, n_cols, subplot_idx)
            scatter = ax.scatter(embedding[:, 0], embedding[:, 1], c=orig_data[:, label].X, cmap=cmap, marker='.')
            plt.title(label, fontsize=15)
            subplot_idx = subplot_idx + 1

    plt.show()


def differential_expression(adata, bottleneck_anndata, clusters_key, n_genes, key_added, show_extra_de_plots=False, DEGs_file=None):
    # Calculate marker genes in adata based on the kmeans clusters in bottleneck_anndata
    sc.tl.rank_genes_groups(adata, clusters_key, key_added=key_added, n_genes=n_genes, method='wilcoxon', use_raw=True)

    # Save all top ranked genes
    top_ranked_genes_bottleneck = pd.DataFrame(adata.uns[key_added]['names'])
    print(top_ranked_genes_bottleneck.head(20))

    # Add top ranked genes to the bottleneck_anndata object - cannot n plot them otherwise because I can't provide a key different than ranked_genes_groups
    bottleneck_anndata.uns['rank_genes_groups'] = adata.uns[key_added]
    sc.pl.rank_genes_groups(bottleneck_anndata, n_genes=20, sharey=False, use_raw=True)

    # Write top_n_genes for each cluster in a csv file
    if DEGs_file is not None:
        top_ranked_genes_bottleneck.to_csv(DEGs_file)

    if show_extra_de_plots is True:
        sc.pl.rank_genes_groups_heatmap(bottleneck_anndata, n_genes=5, use_raw=use_raw, swap_axes=True,
                                        vmin=vmin, vmax=vmax, cmap='bwr', show_gene_labels=True)
        sc.tl.dendrogram(bottleneck_anndata, groupby=clusters_key, use_rep='X_umap', var_names=marker_genes, use_raw=use_raw)
        sc.pl.dendrogram(bottleneck_anndata, groupby=clusters_key)
        sc.pl.rank_genes_groups_dotplot(bottleneck_anndata, n_genes=5)

    return adata.uns[key_added]


def marker_gene_overlap(bottleneck_anndata, marker_genes, ranked_genes):
    fig, ax = plt.subplots(figsize=(10, 3))
    gene_overlap_norm = sc.tl.marker_gene_overlap(bottleneck_anndata, reference_markers=marker_genes,
                                                  key='rank_genes_bottleneck', normalize='reference')
    print(gene_overlap_norm)
    sb.heatmap(gene_overlap_norm, cbar=True, annot=True)
    plt.show()
    return gene_overlap_norm


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Preprocessing of single-cell RNA-seq datasets')
    parser.add_argument('--dataset', type=str, default='E14_hom', help='One of: ["E13_hom", "E14_hom", "integration"]')
    parser.add_argument('--dataset_type', type=str, default='variable', help='One of: ["variable", "all"]')
    parser.add_argument('--n_neighbors', type=int, default=30)
    parser.add_argument('--min_dist', type=float, default=0.5)
    parser.add_argument('--hdbscan', action='store_true')
    parser.add_argument('--top_n_genes', type=int, default=100)
    parser.add_argument('--show_known_marker_genes_plots', action='store_true')
    args = parser.parse_args()

    # One of: E1*_hom_variable, E1*_hom_all, integration_variable, integration_all (where * is 3 or 4)
    dataset = args.dataset + '_' + args.dataset_type
    adata, bottleneck = get_all_data(dataset)
    n_neighbors = 15

    umap_embedding = fit_umap(bottleneck, n_neighbors=args.n_neighbors, min_dist=args.min_dist)

    # K-Means
    kmeans = KMeans(n_clusters=10, n_init=20, random_state=MYSEED).fit(umap_embedding)
    plot_clustering(umap_embedding, kmeans.labels_, title='K-Means')
    print('K-Means silhouette_score:', silhouette_score(umap_embedding, kmeans.labels_, metric='euclidean'))

    # HDBSCAN
    if args.hdbscan:
        hdbscan_fit = hdbscan.HDBSCAN(min_samples=10, min_cluster_size=500)
        hdbscan_labels = hdbscan_fit.fit_predict(umap_embedding)
        plot_clustering(umap_embedding, hdbscan_labels, title='HDBSCAN')
        clustered = (hdbscan_labels >= 0)
        pct_clustered = np.sum(clustered) / bottleneck.shape[0]  # Percentage of cells that were clustered
        print(pct_clustered)

    # Plot expression of known marker genes
    marker_genes, main_cell_types, available_ectopic = get_known_marker_genes(adata)
    if args.show_known_marker_genes_plots:
        plot_marker_genes(umap_embedding, adata, main_cell_types, cmap='PuBu')
        plot_marker_genes(umap_embedding, adata, marker_genes['ectopic'], cmap='PuRd')

        # Plot the expression of the top 20 marker genes found by DE on the initial dataset
        if args.dataset == 'E13_hom':  # E13_hom doesn't have marker genes for ectopic cells so we use E14_hom
            E14_adata = sc.read(Path('ann_data', 'E14_hom_' + args.dataset_type + '_genes.h5ad'))
            ectopic = pd.DataFrame(E14_adata.uns['rank_genes_groups']['names'])['Ectopic cells'].head(20)
            plot_marker_genes(umap_embedding, adata, ectopic, cmap='PuBu', n_cols=4)

    # Initialise AnnData object for the bottleneck (might not be needed)
    bn = sc.AnnData(bottleneck)
    # bn.obs_names = adata.obs_names
    adata.obs['kmeans'] = pd.Categorical(kmeans.labels_)
    bn.obs = adata.obs
    bn.uns = adata.uns
    print(adata)

    # DE on the bottleneck using the kmeans clusters
    DEGs_file = Path('DEGs/bottleneck', dataset + '_' + str(args.top_n_genes) + '.csv')
    rank_key = 'rank_genes_bottleneck'
    ranked_de_genes = differential_expression(adata, bn, 'kmeans', args.top_n_genes, key_added=rank_key, DEGs_file=DEGs_file)
    bn.uns['rank_genes_bottleneck'] = adata.uns['rank_genes_bottleneck']

    # Get overlap of known marker genes with the marker genes found from DE
    gene_overlap_norm = marker_gene_overlap(bn, marker_genes, ranked_de_genes)
    gene_overlap_norm_distr = probability_distr_of_overlap(gene_overlap_norm)
