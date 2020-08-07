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

from scanpy_preprocess.utils import get_known_marker_genes, probability_distr_of_overlap, marker_gene_overlap, annotate_clusters_based_on_overlap, get_updated_marker_genes

sc.logging.print_version_and_date()
sns.set(style='white', rc={'figure.figsize': (15, 10)})
MYSEED = 2211
np.random.seed(MYSEED)

DATASET_NAME = ''
FIGDIR = './figures/bottleneck/'

sc.settings.figdir = FIGDIR
sc.settings.file_format_figs = 'eps'
sc.settings._vector_friendly = False
sc.settings.autosave = True
sc.settings.autoshow = False
sc.settings._frameon = False

vmin = -5
vmax = 5
use_raw = False


celltypes = ['Neural Progenitors', 'Intermediate Progenitors', 'Post-mitotic Neurons', 'Ectopic', 'Unknown']


def get_all_data(dataset):
    adata = sc.read(Path('ann_data', dataset + '_genes.h5ad'))
    bottleneck = pd.read_csv(Path('scDeepCluster/bottleneck', dataset + '_bottleneck.csv'), header=None)
    return adata, bottleneck


def fit_umap(data, n_neighbors=15, min_dist=0.5, metric='euclidean'):
    np.random.seed(MYSEED)
    umap_fit = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist, metric=metric, random_state=MYSEED)
    embedding = umap_fit.fit_transform(data)
    return embedding


def plot_marker_genes(embedding, adata, labels, n_cols=3, title='', plot_type='', cmap='tab10'):
    fig = plt.figure()
    len_labels = len(labels)
    n_rows = (len_labels // n_cols) + 1

    # Create subplots
    subplot_idx = 1
    for label in labels:
        if label in adata.var_names:
            ax = fig.add_subplot(n_rows, n_cols, subplot_idx)
            ax.scatter(embedding[:, 0], embedding[:, 1], c=adata[:, label].X, cmap=cmap, marker='.')
            ax.set_title(label)
            plt.xticks(color='w')
            plt.yticks(color='w')
            subplot_idx = subplot_idx + 1

    fig_name = '_'.join(['marker_genes', plot_type, DATASET_NAME])
    fig.set_rasterized(True)
    fig.savefig(FIGDIR + fig_name + '.eps')
    # plt.show()


def differential_expression(adata, bottleneck_anndata, clusters_key, n_genes, key_added, DEGs_file=None):
    # Find marker genes in adata based on the kmeans clusters in bottleneck_anndata
    sc.tl.rank_genes_groups(adata, clusters_key, key_added=key_added, n_genes=n_genes, method='wilcoxon', use_raw=True)

    # Save all top ranked genes
    top_ranked_genes_bottleneck = pd.DataFrame(adata.uns[key_added]['names'])
    print(top_ranked_genes_bottleneck.head(20))

    # Add top ranked genes to the bottleneck_anndata object - cannot n plot them otherwise because I can't provide a key different than ranked_genes_groups
    bottleneck_anndata.uns['rank_genes_groups'] = adata.uns[key_added]
    sc.pl.rank_genes_groups(bottleneck_anndata, n_genes=20, sharey=False, use_raw=True, save=DATASET_NAME)

    # Write top_n_genes for each cluster in a csv file
    if DEGs_file is not None:
        top_ranked_genes_bottleneck.to_csv(DEGs_file)

    # Show DE plots
    # sc.pl.rank_genes_groups_heatmap(adata, groupby='kmeans', n_genes=5, use_raw=use_raw, swap_axes=True, vmin=vmin, vmax=vmax,
    #                                 cmap='bwr', show_gene_labels=True, var_group_rotation=45, save=DATASET_NAME)
    sc.tl.dendrogram(bottleneck_anndata, groupby=clusters_key, use_rep='X_umap', var_names=marker_genes, use_raw=use_raw)
    sc.pl.dendrogram(bottleneck_anndata, groupby=clusters_key, save=DATASET_NAME)
    # sc.pl.rank_genes_groups_dotplot(bottleneck_anndata, n_genes=5, save=DATASET_NAME)

    return adata.uns[key_added]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Preprocessing of single-cell RNA-seq datasets')
    parser.add_argument('--dataset', nargs='+', help='One or more space separated datasets, e.g. E14_hom E13_hom integration')
    parser.add_argument('--dataset_type', type=str, default='variable', help='One of: ["variable", "all"]')
    parser.add_argument('--n_neighbors', type=int, default=30)
    parser.add_argument('--min_dist', type=float, default=0.5)
    parser.add_argument('--top_n_genes', type=int, default=100)
    parser.add_argument('--n_clusters', type=int, default=10, help='n_clusters for K-Means')
    parser.add_argument('--fig_dir', type=str, default='./figures/bottleneck/')
    parser.add_argument('--write_to_file', action='store_true', help='Write bottleneck AnnData to h5ad file')
    args = parser.parse_args()

    alldata = {}
    FIGDIR = args.fig_dir
    Path(FIGDIR).mkdir(parents=True, exist_ok=True)
    Path('DEGs/bottleneck/').mkdir(parents=True, exist_ok=True)
    Path('ann_data/bottleneck/').mkdir(parents=True, exist_ok=True)

    for dataset in args.dataset:
        DATASET_NAME = dataset

        # One of: E1*_hom_variable, E1*_hom_all, integration_variable, integration_all (where * is 3 or 4)
        dataset_with_type = dataset + '_' + args.dataset_type
        adata, bottleneck = get_all_data(dataset_with_type)

        umap_embedding = fit_umap(bottleneck, n_neighbors=args.n_neighbors, min_dist=args.min_dist)

        # K-Means
        kmeans = KMeans(n_clusters=args.n_clusters, n_init=20, random_state=MYSEED).fit(umap_embedding)
        print('K-Means silhouette_score:', silhouette_score(umap_embedding, kmeans.labels_, metric='euclidean'))

        # HDBSCAN - not used anywhere
        hdbscan_fit = hdbscan.HDBSCAN(min_samples=10, min_cluster_size=500)
        hdbscan_labels = hdbscan_fit.fit_predict(umap_embedding)
        clustered = (hdbscan_labels >= 0)
        pct_clustered = np.sum(clustered) / bottleneck.shape[0]  # Percentage of cells that were clustered
        print('HDBSCAN percentage clustered:', pct_clustered)

        # Plot expression of known marker genes
        marker_genes, main_cell_types, available_ectopic = get_known_marker_genes(adata)
        # marker_genes, main_cell_types, available_ectopic = get_updated_marker_genes(adata)

        # plot_marker_genes(umap_embedding, adata, main_cell_types, plot_type='main_cell_types', cmap='PuBu')
        # plot_marker_genes(umap_embedding, adata, marker_genes['Ectopic'], plot_type='ectopic', cmap='PuRd')

        # Plot the expression of the top 20 marker genes found by DE on the initial dataset
        # if DATASET_NAME == 'E13_hom':  # E13_hom doesn't have marker genes for ectopic cells so we use E14_hom
        #     E14_adata = sc.read(Path('ann_data', 'E14_hom_' + args.dataset_type + '_genes.h5ad'))
        #     ectopic = pd.DataFrame(E14_adata.uns['rank_genes_groups']['names'])['Ectopic cells'].head(20)
        #     plot_marker_genes(umap_embedding, adata, ectopic, plot_type='DE_markers', cmap='PuBu', n_cols=4)

        #
        # =====================================================================
        #               ANALYSIS OF BOTTLENECK STARTS HERE
        # =====================================================================
        #

        # Initialise AnnData object for the bottleneck (might not be needed)
        bn = sc.AnnData(bottleneck)
        # bn.obs_names = adata.obs_names
        adata.obs['kmeans'] = pd.Categorical(kmeans.labels_)
        adata.obs['hdbscan'] = pd.Categorical(hdbscan_labels)
        bn.obs = adata.obs
        bn.uns = adata.uns
        bn.obsm['X_umap'] = umap_embedding

        print('adata', adata)
        print('bn', bn)

        # DE on the bottleneck using the kmeans clusters
        DEGs_file = Path('DEGs/bottleneck', dataset_with_type + '_' + str(args.top_n_genes) + '.csv')
        rank_key = 'rank_genes_bottleneck'

        ranked_de_genes = differential_expression(adata, bn, 'kmeans', args.top_n_genes, key_added=rank_key, DEGs_file=DEGs_file)
        bn.uns[rank_key] = adata.uns[rank_key]

        # Plot the expression of the top 20 marker genes found by DE on the bottleneck!
        # Select clusters that we believe there are ectopic
        # @TODO: See which clusters express more ectopic genes and use that number instead of 9
        # if DATASET_NAME == 'E14_hom':
        #     bn_DE_genes = pd.DataFrame(bn.uns[rank_key]['names'])['9'].head(20)
        #     plot_marker_genes(umap_embedding, adata, bn_DE_genes, plot_type='DE_NEW_9', cmap='PuBu', n_cols=4)

        #     bn_DE_genes = pd.DataFrame(bn.uns[rank_key]['names'])['0'].head(20)
        #     plot_marker_genes(umap_embedding, adata, bn_DE_genes, plot_type='DE_NEW_0', cmap='PuBu', n_cols=4)

        # Get overlap of known marker genes with the marker genes found from DE
        gene_overlap_norm = marker_gene_overlap(bn, marker_genes, FIGDIR, DATASET_NAME, key_added=rank_key)
        gene_overlap_norm_distr = probability_distr_of_overlap(gene_overlap_norm)

        #
        # ==========================================================================
        #       ANNOTATION - might bot be needed the annotation in the bottleneck
        # ==========================================================================
        bn.obs['bottleneck_annotations'] = annotate_clusters_based_on_overlap(gene_overlap_norm, bn, clusters_key='kmeans')
        sc.pl.umap(bn, color='kmeans', title='Initial Clustering', save=DATASET_NAME + '_kmeans', legend_loc='on data')
        sc.pl.umap(bn, color=['leiden_annotations', 'bottleneck_annotations'],
                   title=['E14_hom annotations', 'Bottleneck annotations'], save=DATASET_NAME + '_both_annotations')

        # Print frequencies of cell types for each annotation
        print('Leiden annotations')
        for ct in celltypes:
            print(ct, len(bn[bn.obs['leiden_annotations'] == ct]))

        print('Bottleneck annotations')
        for ct in celltypes:
            print(ct, len(bn[bn.obs['bottleneck_annotations'] == ct]))

        if args.write_to_file:
            # Write in h5ad file to use as an input to scDeepCluster
            filename = '_'.join(['bottleneck', dataset, args.dataset_type]) + '.h5ad'
            bottleneck_file_path = Path('ann_data', 'bottleneck', filename)
            bn.write(bottleneck_file_path)
            print('{} file saved'.format(bottleneck_file_path))
