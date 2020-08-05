import scanpy as sc
import pandas as pd
import numpy as np
import argparse
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from pathlib import Path
from utils import get_colormap, get_known_marker_genes

sc.logging.print_versions()
# sc.set_figure_params(facecolor="white", figsize=(6, 4))
sc.settings.verbosity = 3
np.random.seed(2211)
sns.set(style='whitegrid')
# Create an array with the colors you want to use
colors = ["#3498db", "#9b59b6", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"]
sns.set_palette(sns.color_palette(colors))  # Set your custom color palette

FIG_DIR = './figures/preprocess/'
sc.settings.figdir = FIG_DIR
sc.settings.file_format_figs = 'eps'
sc.settings._vector_friendly = False
sc.settings.autosave = True
sc.settings.autoshow = False
sc.settings._frameon = False

DATASET_NAME = ''


def plot_qc_measures(adata):
    # sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.7, multi_panel=True,
    #              xlabel='test', name=['1', '2', '3'])
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', save='_mt_' + DATASET_NAME)
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', save='_ngenes_' + DATASET_NAME)


def plot_qc_distplot(adata):
    mpl.rcParams['axes.titlesize'] = 18
    mpl.rcParams['axes.labelsize'] = 14
    mpl.rcParams['xtick.labelsize'] = 12
    mpl.rcParams['ytick.labelsize'] = 12
    mpl.rcParams['figure.titlesize'] = 18

    keys = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']
    key_names = ['Number of genes with counts per cell', 'Total number of genes counts per cell',
                 'Proportion of total mitochondrial counts per cell']

    fig, axs = plt.subplots(1, len(keys), figsize=(10, 3))
    for idx, k in enumerate(keys):
        sns.distplot(adata.obs[k], ax=axs[idx], axlabel=key_names[idx], hist=True)
        axs[idx].set_ylabel('Proportion of cells')
    fig.tight_layout()
    fig.show()

    # Violin plots
    titles_2 = ['Number of genes with positive counts per cell', 'Total number of genes counts per cell',
                'Total mitochondrial counts per cell']
    ylabels_2 = ['Number of genes', 'Number of genes', 'Number of mitochondrial genes']

    fig2, axs2 = plt.subplots(1, len(keys), figsize=(15, 3))
    for idx, k in enumerate(keys):
        sns.violinplot(data=adata.obs[k], ax=axs2[idx])
        sns.stripplot(data=adata.obs[k], ax=axs2[idx], jitter=0.7, color='black', size=2)
        axs2[idx].set_ylabel(ylabels_2[idx])
        axs2[idx].set_title(titles_2[idx])
    fig2.tight_layout()
    fig2.show()

    # Scatter plot
    scatter_keys = ['n_genes_by_counts', 'pct_counts_mt']
    fig3, axs3 = plt.subplots(1, len(scatter_keys), figsize=(15, 3))
    for idx, k in enumerate(scatter_keys):
        sns.scatterplot(data=adata.obs, x='total_counts', y=k, ax=axs3[idx], markers=True)
        # axs3[idx].set_ylabel(ylabels_2[idx])
        # axs3[idx].set_title(titles_2[idx])
    fig3.tight_layout()
    fig3.show()

    # # Scatter
    # fig3, axs3 = plt.subplots(1, 2, figsize=(6, 6))
    # sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', ax=axs3[0])
    # sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', ax=axs3[1])


def normalise_data(adata, keep_only_highly_variable=True):
    sc.pp.normalize_total(adata, target_sum=1e4, key_added='size_factors')
    sc.pp.log1p(adata)

    if keep_only_highly_variable:
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        sc.pl.highly_variable_genes(adata, save=DATASET_NAME)
        adata = adata[:, adata.var.highly_variable]
        print('highly variable genes:', adata.shape)
    adata.raw = raw_counts[:, list(adata.var_names)]  # Update raw counts object with the updated list of genes
    return adata


def run_pca(adata):
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pl.pca(adata, color=['total_counts', 'pct_counts_mt'], save='_1_' + DATASET_NAME)
    sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50, save=DATASET_NAME)


def cell_cycle_scoring(adata):
    # Load cell cycle genes from file and split into S and G2M genes
    cell_cycle_genes = [x.strip().lower().capitalize() for x in open('scanpy_preprocess/regev_lab_cell_cycle_genes.txt')]
    s_genes = cell_cycle_genes[:43]
    g2m_genes = cell_cycle_genes[43:]
    cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]

    # Perform cell cycle scoring
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

    # Plots
    adata_cc_genes = adata[:, cell_cycle_genes]  # Subset of the data using only the cell cycle phase genes

    # PCA
    sc.tl.pca(adata_cc_genes)
    sc.pl.pca_scatter(adata_cc_genes, color='phase', title='PCA', save='_before_' + DATASET_NAME)

    # UMAP on 40 PCs
    sc.pp.neighbors(adata_cc_genes, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata_cc_genes)
    sc.pl.umap(adata_cc_genes, color='phase', title='Cell cycle phases before regression', save='_ccbefore_' + DATASET_NAME)

    return cell_cycle_genes


def plot_cell_cycle_after_regression(adata, cell_cycle_genes):
    adata_cc_genes_regressed = adata[:, cell_cycle_genes]
    print('After regressing out the cell cycle')

    # PCA
    sc.tl.pca(adata_cc_genes_regressed)
    sc.pl.pca_scatter(adata_cc_genes_regressed, color='phase', save='_after_' + DATASET_NAME)

    # UMAP on 40 PCs
    sc.pp.neighbors(adata_cc_genes_regressed, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata_cc_genes_regressed)
    sc.pl.umap(adata_cc_genes_regressed, color='phase', title='Cell cycle phases after regression', save='_ccafter_' + DATASET_NAME)


def clustering(adata, dataset, keep_only_highly_variable):
    # Use only the first 40 PCs for clustering --> more distinct clusters
    sc.pp.neighbors(adata, n_neighbors=50, n_pcs=40)
    sc.tl.leiden(adata)
    sc.tl.umap(adata, min_dist=0.5)
    sc.pl.umap(adata, color=['leiden'], legend_loc='on data', save='_leiden_' + DATASET_NAME)

    # Remove Cajal Retzius cells
    if dataset == 'E14_hom':
        cluster = '12'
    elif dataset == 'E13_hom':
        cluster = '11' if keep_only_highly_variable else '10'
    elif dataset == 'E13_het':
        cluster = '10'
    elif dataset == 'E14_het':
        cluster = '11'

    adata = adata[~adata.obs['leiden'].isin([cluster]), :]
    # adata.raw = adata.copy()
    return adata


def plot_marker_genes(adata, marker_genes, main_cell_types):
    if DATASET_NAME == 'E13_hom' and ('Vim' in marker_genes['Neural Progenitors']):
        marker_genes['Neural Progenitors'].remove('Vim')
        main_cell_types.remove('Vim')
        print('Vim removed')

    print('marker_genes:', marker_genes)
    var_names = set(adata.var_names)
    columns = set(adata.obs.columns)
    gene_names = var_names.union(columns)
    available_ectopic = gene_names.intersection(marker_genes['Ectopic'])

    use_raw = False
    vmin = -5
    vmax = 5

    print('3 cell types:')
    sc.pl.umap(adata, color=main_cell_types, cmap=get_colormap('blue'), legend_loc='on data', size=50, ncols=3,
               vmin=vmin, vmax=vmax, use_raw=use_raw, save='_3types_' + DATASET_NAME)

    print('Ectopic:')
    sc.pl.umap(adata, color=available_ectopic, cmap=get_colormap(), size=50, ncols=3,
               vmin=vmin, vmax=vmax, use_raw=use_raw, save='_ectopic_' + DATASET_NAME)

    # Print all in one image
    sc.pl.umap(adata, color=main_cell_types + list(available_ectopic), cmap=get_colormap('blue'), size=50,
               vmin=vmin, vmax=vmax, use_raw=use_raw, save='_all_' + DATASET_NAME)


def plot_heatmap_dotplot(adata, marker_genes, main_cell_types, available_ectopic):
    marker_genes['Ectopic'] = available_ectopic

    # A dot is plotted for each gene and each cluster.
    # Each dot represents two values:
    # 1. mean expression within each cluster (visualized by color)
    # 2. fraction of cells expressing the gene in the cluster (visualized by the size of the dot)
    # sc.pl.dotplot(adata, var_names=main_cell_types + list(available_ectopic), groupby='leiden', color_map='viridis', save=DATASET_NAME,
    print('marker_genes', marker_genes)
    sc.pl.dotplot(adata, var_names=marker_genes, groupby='leiden', color_map='viridis', save=DATASET_NAME,
                  use_raw=False, log=True, standard_scale='group')

    # Heatmap of the expression values of genes in each cluster
    sc.pl.heatmap(adata, var_names=marker_genes, groupby='leiden', log=True, cmap='magma', standard_scale='var',
                  swap_axes=False, dendrogram=True, save=DATASET_NAME)

    # Heatmap of the mean expression values per cluster of each gene
    sc.pl.matrixplot(adata, var_names=marker_genes, groupby='leiden', cmap='magma', log=True, standard_scale='var',
                     save=DATASET_NAME)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Preprocessing of single-cell RNA-seq datasets')
    parser.add_argument('--dataset', nargs='+', help='One or more space separated datasets, e.g. E14_hom E13_hom')
    parser.add_argument('--keep_only_highly_variable', action='store_true', help='Whether to only keep ONLY highly variable genes or all genes')
    parser.add_argument('--write_to_file', action='store_true', help='Write preprocess data to h5ad file')
    args = parser.parse_args()

    alldata = {}

    for dataset in args.dataset:
        DATASET_NAME = dataset
        dataset_path = Path('data', dataset + '.csv')
        print('Reading single-cell csv file: ', dataset_path)
        data = pd.read_csv(dataset_path, index_col=0)
        adata = sc.AnnData(X=data.T)
        print('adata', adata)

        # Quality Control
        print('---Quality control')
        adata.var['mt'] = adata.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, inplace=True)  # compute QC metrics
        print('# genes with 0 counts:', adata[:, adata.var['total_counts'] == 0].shape[1])

        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=3)

        sns.jointplot('log1p_total_counts', 'log1p_n_genes_by_counts', data=adata.obs, kind='hex')
        plt.tight_layout()

        # Quality control plots
        # plot_qc_distplot(adata)
        plot_qc_measures(adata)
        plt.savefig(FIG_DIR + 'jointplot_' + DATASET_NAME)

        # Remove cells with certain threshold
        adata = adata[adata.obs.total_counts < 40000, :]
        adata = adata[adata.obs.n_genes_by_counts < (7500 if 'E14' in dataset else 6000), :]  # for E14: 7500, for E13: 6000
        adata = adata[adata.obs.pct_counts_mt < (10 if 'E14' in dataset else 5), :]  # for E14: 10, for E13: 5
        raw_counts = adata.copy()

        # Normalisation
        print('---Normalisation')
        adata = normalise_data(adata, args.keep_only_highly_variable)
        run_pca(adata)

        print('---Cell cycle scoring')
        cell_cycle_genes = cell_cycle_scoring(adata)
        sc.pp.regress_out(adata, ['pct_counts_mt', 'S_score', 'G2M_score'])
        sc.pp.scale(adata, max_value=10)
        plot_cell_cycle_after_regression(adata, cell_cycle_genes)

        print('---Clustering')
        adata = clustering(adata, dataset, args.keep_only_highly_variable)

        print('---Known marker genes')
        marker_genes, main_cell_types, available_ectopic = get_known_marker_genes(adata)
        plot_marker_genes(adata, marker_genes, main_cell_types)

        plot_heatmap_dotplot(adata, marker_genes, main_cell_types, available_ectopic)

        if args.write_to_file:
            # Write in h5ad file to use as an input to scDeepCluster
            filename = dataset + ('_variable_genes.h5ad' if args.keep_only_highly_variable else '_all_genes.h5ad')
            processed_file = Path('ann_data', filename)
            adata.write(processed_file)
            print('{} file saved'.format(processed_file))

        alldata[dataset] = adata
