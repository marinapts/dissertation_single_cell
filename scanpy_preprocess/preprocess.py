import scanpy as sc
import pandas as pd
import numpy as np
import argparse
from pathlib import Path
from utils import get_colormap, get_known_marker_genes

sc.logging.print_versions()
sc.set_figure_params(facecolor="white", figsize=(6, 4))
sc.settings.verbosity = 3
np.random.seed(2211)


def plot_qc_measures(adata):
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.7, multi_panel=True)
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')


def normalise_data(adata, keep_only_highly_variable=True):
    sc.pp.normalize_total(adata, target_sum=1e4, key_added='size_factors')
    sc.pp.log1p(adata)

    if keep_only_highly_variable:
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        # sc.pl.highly_variable_genes(adata)
        adata = adata[:, adata.var.highly_variable]

    adata.raw = raw_counts[:, list(adata.var_names)]  # Update raw counts object with the updated list of genes
    return adata


def run_pca(adata, show_plots=False):
    sc.tl.pca(adata, svd_solver='arpack')
    if show_plots:
        sc.pl.pca(adata, color=['total_counts', 'pct_counts_mt'])
        sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50)


def cell_cycle_scoring(adata, show_plots=False):
    # Load cell cycle genes from file and split into S and G2M genes
    cell_cycle_genes = [x.strip().lower().capitalize() for x in open('scanpy_preprocess/regev_lab_cell_cycle_genes.txt')]
    s_genes = cell_cycle_genes[:43]
    g2m_genes = cell_cycle_genes[43:]
    cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]

    # Perform cell cycle scoring
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

    if show_plots:
        adata_cc_genes = adata[:, cell_cycle_genes]  # Subset of the data using only the cell cycle phase genes

        # PCA
        sc.tl.pca(adata_cc_genes)
        sc.pl.pca_scatter(adata_cc_genes, color='phase', title='PCA')

        # UMAP on 40 PCs
        sc.pp.neighbors(adata_cc_genes, n_neighbors=10, n_pcs=40)
        sc.tl.umap(adata_cc_genes)
        sc.pl.umap(adata_cc_genes, color='phase', title='UMAP on 40 PCs')

    return cell_cycle_genes


def plot_cell_cycle_after_regression(adata, cell_cycle_genes):
    adata_cc_genes_regressed = adata[:, cell_cycle_genes]
    print('After regressing out the cell cycle')

    # PCA
    sc.tl.pca(adata_cc_genes_regressed)
    sc.pl.pca_scatter(adata_cc_genes_regressed, color='phase')

    # UMAP on 40 PCs
    sc.pp.neighbors(adata_cc_genes_regressed, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata_cc_genes_regressed)
    sc.pl.umap(adata_cc_genes_regressed, color='phase', title='UMAP on 40 PCs')


def clustering(adata, dataset, keep_only_highly_variable):
    # Use only the first 40 PCs for clustering --> more distinct clusters
    sc.pp.neighbors(adata, n_neighbors=50, n_pcs=40)
    sc.tl.leiden(adata)
    sc.tl.umap(adata, min_dist=0.5)
    sc.pl.umap(adata, color=['leiden'], legend_loc='on data')

    # Remove Cajal Retzius cells
    if dataset == 'E14_hom':
        cluster = '12'
    elif dataset == 'E13_hom':
        cluster = '11' if keep_only_highly_variable else '10'
    elif dataset == 'E13_het':
        cluster = '12'
    elif dataset == 'E14_het':
        cluster = '19'

    adata = adata[~adata.obs['leiden'].isin([cluster]), :]
    # adata.raw = adata.copy()
    return adata


def plot_marker_genes(adata, marker_genes, main_cell_types):
    var_names = set(adata.var_names)
    columns = set(adata.obs.columns)
    gene_names = var_names.union(columns)
    available_ectopic = gene_names.intersection(marker_genes['ectopic'])

    use_raw = False
    vmin = -5
    vmax = 5

    # print('3 cell types:')
    # sc.pl.umap(adata, color=main_cell_types, cmap=get_colormap(), legend_loc='on data', size=50, ncols=3,
    #            vmin=vmin, vmax=vmax, use_raw=use_raw)

    # print('Ectopic:')
    # sc.pl.umap(adata, color=available_ectopic, cmap=get_colormap(), size=50, ncols=3,
    #            vmin=vmin, vmax=vmax, use_raw=use_raw)

    # Print all in one image
    sc.pl.umap(adata, color=main_cell_types + list(available_ectopic), cmap=get_colormap(), size=50,
               vmin=vmin, vmax=vmax, use_raw=use_raw)


def plot_heatmap_dotplot(adata, marker_genes, main_cell_types, available_ectopic):
    marker_genes['ectopic'] = available_ectopic

    # A dot is plotted for each gene and each cluster.
    # Each dot represents two values:
    # 1. mean expression within each cluster (visualized by color)
    # 2. fraction of cells expressing the gene in the cluster (visualized by the size of the dot)
    sc.pl.dotplot(adata, var_names=main_cell_types + list(available_ectopic), groupby='leiden', color_map='viridis',
                  use_raw=False, log=True, standard_scale='group')

    # Heatmap of the expression values of genes
    sc.pl.heatmap(adata, var_names=marker_genes, groupby='leiden', cmap='magma', standard_scale='var')
    sc.pl.heatmap(adata, var_names=marker_genes, groupby='leiden', log=True, cmap='magma', standard_scale='var', swap_axes=False)

    # Heatmap of the mean expression values per cluster of each gene
    sc.pl.matrixplot(adata, var_names=marker_genes, groupby='leiden', cmap='magma', log=True, standard_scale='var')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Preprocessing of single-cell RNA-seq datasets')
    parser.add_argument('--dataset', type=str, default='E13_hom', help='Path to the csv file')
    parser.add_argument('--keep_only_highly_variable', action='store_true', help='Whether to only keep ONLY highly variable genes or all genes')
    parser.add_argument('--show_qc_plots', action='store_true', help='Whether to show quality control plot or not')
    parser.add_argument('--show_pca_plots', action='store_true', help='Whether to show pca plots or not')
    parser.add_argument('--show_cell_cycle_plots', action='store_true', help='Whether to show cell cycle plots or not')
    parser.add_argument('--show_marker_genes_plots', action='store_true', help='Whether to show expression plots for each marker gene')
    parser.add_argument('--show_extra_plots', action='store_true', help='Whether to show extra heatmaps and dotplots')
    parser.add_argument('--write_to_file', action='store_true', help='Write preprocess data to h5ad file')
    args = parser.parse_args()

    dataset_path = Path('data', args.dataset + '.csv')
    print('Reading single-cell csv file: ', dataset_path)
    data = pd.read_csv(dataset_path, index_col=0)
    adata = sc.AnnData(X=data.T)
    print(adata)

    # Quality Control
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata.var['mt'] = adata.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)  # compute QC metrics

    if args.show_qc_plots:
        plot_qc_measures(adata)

    # Remove cells with certain threshold
    adata = adata[adata.obs.total_counts < 40000, :]
    adata = adata[adata.obs.n_genes_by_counts < (7500 if 'E14' in args.dataset else 6000), :]  # for E14: 7500, for E13: 6000
    adata = adata[adata.obs.pct_counts_mt < (10 if 'E14' in args.dataset else 5), :]  # for E14: 10, for E13: 5
    raw_counts = adata.copy()

    # Normalisation
    adata = normalise_data(adata, args.keep_only_highly_variable)
    run_pca(adata, args.show_pca_plots)
    cell_cycle_genes = cell_cycle_scoring(adata, args.show_cell_cycle_plots)
    sc.pp.regress_out(adata, ['pct_counts_mt', 'S_score', 'G2M_score'])
    sc.pp.scale(adata, max_value=10)

    # plot_cell_cycle_after_regression(adata, cell_cycle_genes)
    adata = clustering(adata, args.dataset, args.keep_only_highly_variable)

    marker_genes, main_cell_types, available_ectopic = get_known_marker_genes(adata)
    if args.show_marker_genes_plots:
        plot_marker_genes(adata, marker_genes, main_cell_types)

    if args.show_extra_plots:
        plot_heatmap_dotplot(adata, marker_genes, main_cell_types, available_ectopic)

    if args.write_to_file:
        # Write in h5ad file to use as an input to scDeepCluster
        filename = args.dataset + ('_variable_genes.h5ad' if args.keep_only_highly_variable else '_all_genes.h5ad')
        processed_file = Path('ann_data', filename)
        adata.write(processed_file)
        print('{} file saved'.format(processed_file))
