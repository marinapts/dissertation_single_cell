import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path


def plot_correlation_matrix(corr_matrix, square=False, title='', fig_name=''):
    fig, ax = plt.subplots(1)
    sns.heatmap(corr_matrix, annot=False, square=square, cmap="coolwarm", yticklabels=True, ax=ax)

    plt.setp(ax)
    plt.title(title)
    plt.tight_layout()
    ax.set_xlabel('Known ectopic marker genes')
    ax.set_ylabel('Genes highly correlated with the known ectopic')
    if fig_name != '':
        plt.savefig(figdir + 'corr_' + fig_name + '.eps', dpi=300)
    plt.show()


def find_high_correlation(which_marker_genes, corr_value_cutoff=0.5, fig_name=''):
    gene_union = corr_matrix.index.union(which_marker_genes)
    subset = corr_matrix.loc[gene_union, which_marker_genes]

    if corr_value_cutoff is None:
        high_corr = subset
    else:
        high_corr = subset.loc[subset.values > corr_value_cutoff]
        high_corr = high_corr[~high_corr.index.duplicated(keep='first')]  # Remove duplicate rows

    plot_correlation_matrix(high_corr, fig_name=fig_name)

    return high_corr


if __name__ == '__main__':
    # Settings
    figdir = './figures/more_plots/'
    Path(figdir).mkdir(parents=True, exist_ok=True)

    # Load preprocessed datasets
    E14 = sc.read('ann_data/E14_hom_variable_genes.h5ad')
    E14_df = E14.to_df()
    E13 = sc.read('ann_data/E13_hom_variable_genes.h5ad')
    E13_df = E13.to_df()

    known_marker_genes = ['Pax6', 'Vim', 'Sox2', 'Eomes', 'Btg2', 'Tbr1', 'Sox5', 'Gsx2', 'Prdm13', 'Dlx1', 'Dlx2',
                          'Dlx5', 'Gad1', 'Gad2', 'Ptf1a', 'Msx3', 'Helt', 'Olig3']
    neural_prog_markers = ['Pax6', 'Vim', 'Sox2']
    intermediate_prog_markers = ['Btg2', 'Tbr1']
    post_mitotic_markers = ['Sox5', 'Gsx2']
    ectopic_markers = ['Gsx2', 'Prdm13', 'Dlx1', 'Dlx2', 'Dlx5', 'Gad1', 'Gad2', 'Msx3', 'Helt']

    # Remove ectopic markers that arent in the df columns
    for gene in ectopic_markers:
        if gene not in E13_df.columns:
            ectopic_markers.remove(gene)

    # Pearson correlation between genes
    corr_matrix = E14_df.corr()

    # Correlation matrix of all known marker genes
    # ---------------------------------------------
    gene_intersection = corr_matrix.index.intersection(known_marker_genes)
    corr_matrix_subset = corr_matrix.loc[gene_intersection, gene_intersection]
    plot_correlation_matrix(corr_matrix_subset, title='Correlation matrix of all known marker genes', fig_name='E14_hom_all_marker_genes')

    # Correlation matrix of the known ectopic with all other genes > 0.5 correlation value
    # -------------------------------------------------------------------------------------
    ectopic_high_corr = find_high_correlation(ectopic_markers, fig_name='ectopic_high')
    neural_prog_high_corr = find_high_correlation(neural_prog_markers, fig_name='neural_prog_high')
    intermediate_prog_high_corr = find_high_correlation(intermediate_prog_markers, fig_name='intermediate_prog_high')
