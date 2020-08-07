import scanpy as sc
import pandas as pd
import numpy as np
import argparse
from pathlib import Path
from utils import get_colormap, get_known_marker_genes, probability_distr_of_overlap, marker_gene_overlap, annotate_clusters_based_on_overlap

sc.logging.print_versions()
sc.set_figure_params(facecolor="white", figsize=(8, 4))
sc.settings.verbosity = 3
np.random.seed(2211)

DATASET_NAME = ''
FIGDIR = './figures/marker_genes/'

sc.settings.figdir = FIGDIR
sc.settings.file_format_figs = 'eps'
sc.settings._vector_friendly = False
sc.settings.autosave = True
sc.settings.autoshow = False
sc.settings._frameon = False

vmin = -5
vmax = 5
use_raw = False


def plot_clusters(adata, clusters_key, is_integrated_dataset, main_cell_types, available_ectopic):
    sc.pl.umap(adata, color=[clusters_key], legend_loc='on data', save=DATASET_NAME, title='Leiden clusters')
    if is_integrated_dataset is True:
        sc.pl.umap(adata, color=['batch'], save='_batch_' + DATASET_NAME, title='Leiden clusters')

    # Plot expression of marker genes
    sc.pl.umap(adata, color=main_cell_types + list(available_ectopic), cmap=get_colormap(),
               vmin=-5, vmax=5, use_raw=False, save='_markergenes_expr_' + DATASET_NAME)


def differential_expression(adata, clusters_key, top_n_genes, DEGs_file, key_added='rank_genes_groups', updated=False):
    """Finds top_n differentially expressed genes in each cluster and ranks them using Wilcoxon rank-sum.
    Then writes the DE genes into a csv and updates the adata with adata.uns[key_added]. (+ some plots)
    Args:
        adata (AnnData)
        clusters_key (str): The key of adata.obs that contains the cluster numbers (or annotations)
        top_n_genes (int): Num of genes to appear in table
        DEGs_file (TYPE): Name of the csv (or latex) file with the DEGs for each cluster - to be saved under /DEGs
        fig_name (TYPE): Name for the figure to be saved
    """
    fig_title = '_updated_' + DATASET_NAME if updated else DATASET_NAME
    # Calculate marker genes
    sc.tl.rank_genes_groups(adata, clusters_key, n_genes=top_n_genes, key_added=key_added, method='wilcoxon', use_raw=True)  # 100 genes by default
    sc.pl.rank_genes_groups(adata, n_genes=20, key=key_added, sharey=False, use_raw=True, save=fig_title)

    print('Top ranked genes for key ', key_added)
    top_ranked_genes_per_cluster = pd.DataFrame(adata.uns[key_added]['names'])
    print(top_ranked_genes_per_cluster.head(20))

    # Write top_n_genes for each cluster in a csv file
    if DEGs_file is not None:
        top_ranked_genes_per_cluster.to_csv(DEGs_file + '.csv')
        top_ranked_genes_per_cluster.to_latex(DEGs_file + '.latex')
        (top_ranked_genes_per_cluster.transpose()).head(10).to_latex(DEGs_file + '.latex.T')

    # Some plots
    sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, key=key_added, use_raw=use_raw, swap_axes=True, vmin=vmin, vmax=vmax,
                                    cmap='bwr', show_gene_labels=True, var_group_rotation=45, save=fig_title)
    sc.tl.dendrogram(adata, groupby=clusters_key, use_rep='X_umap', var_names=marker_genes, use_raw=use_raw)
    sc.pl.dendrogram(adata, groupby=clusters_key, save=fig_title)
    sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, key=key_added, save=fig_title)

    return adata.uns[key_added]


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

        #
        # =====================================================================
        #            DIFFERENTIAL EXPRESSION ANALYSIS STARTS HERE
        # =====================================================================
        #

        top_n_genes = 100
        DEGs_filename = dataset + '_' + args.dataset_type + '_' + str(top_n_genes)
        DEGs_file = 'DEGs/' + DEGs_filename

        # DE to calculate marker genes for the clusters
        ranked_genes = differential_expression(adata, clusters_key, top_n_genes, DEGs_file)

        # Write marker genes in a csv for annotation by SCSA
        groups = ranked_genes['names'].dtype.names
        dat = pd.DataFrame({group + '_' + key[:1]: ranked_genes[key][group] for group in groups for key in ['names', 'logfoldchanges', 'scores', 'pvals']})
        dat.to_csv('cellmarker/DE_' + data_filename + '.csv')

        # Get overlap of known marker genes with the marker genes found from DE
        gene_overlap_norm = marker_gene_overlap(adata, marker_genes, FIGDIR, DATASET_NAME)
        gene_overlap_norm_distr = probability_distr_of_overlap(gene_overlap_norm)

        # Annotate clusters
        adata.obs[clusters_key + '_annotations'] = annotate_clusters_based_on_overlap(gene_overlap_norm, adata, clusters_key=clusters_key)
        sc.pl.umap(adata, color=[clusters_key + '_annotations'], save='_annotations_' + DATASET_NAME, title='Manual annotations')

        # DE to find marker genes for the updated clusters (cell types)
        annotations_key = 'rank_genes_annotations'
        updated_ranked_genes = differential_expression(adata, clusters_key + '_annotations', top_n_genes, None,
                                                       key_added=annotations_key, updated=True)

        updated_gene_overlap_norm = marker_gene_overlap(adata, marker_genes, FIGDIR, DATASET_NAME, key_added=annotations_key, updated=True)
        updated_gene_overlap_norm_distr = probability_distr_of_overlap(updated_gene_overlap_norm)

        if args.update_file:
            adata.write(filepath)
            print('{} file updated'.format(filepath))
