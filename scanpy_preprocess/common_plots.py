import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns


def plot_correlation_matrix(corr_matrix, name, title):
    sns.heatmap(corr_matrix, square=True, cmap="coolwarm")
    plt.savefig(figdir + 'corr_' + name + '.eps', dpi=300)
    plt.title(title)
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    # Settings
    figdir = './figures/more_plots/'

    # Load preprocessed datasets
    E14 = sc.read('ann_data/E14_hom_variable_genes.h5ad')
    E13 = sc.read('ann_data/E13_hom_variable_genes.h5ad')

    known_marker_genes = ['Pax6', 'Vim', 'Sox2', 'Eomes', 'Btg2', 'Tbr1', 'Sox5', 'Gsx2', 'Prdm13', 'Dlx1', 'Dlx2',
                          'Dlx5', 'Gad1', 'Gad2', 'Ptf1a', 'Msx3', 'Helt', 'Olig3']

    # Pearson correlation between genes
    corr_matrix = E14.to_df().corr()

    gene_intersection = corr_matrix.index.intersection(known_marker_genes)
    corr_matrix_subset = corr_matrix.loc[gene_intersection, gene_intersection]
    # plot_correlation_matrix(corr_matrix_subset, 'E14_hom', 'Correlation matrix of marker genes - 3 cell types')
    sns.heatmap(corr_matrix_subset, square=True, annot=True, cmap="RdBu")
    plt.title('Correlation matrix of marker genes - 3 cell types')
    plt.tight_layout()
    plt.show()
    # plt.savefig(figdir + 'corr_E14_hom.eps', dpi=300)
