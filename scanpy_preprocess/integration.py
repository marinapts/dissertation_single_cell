import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanorama
import argparse
from pathlib import Path


sc.logging.print_versions()
sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3
np.random.seed(2211)


def visualise_qc_metrics(E13, E14):
    # Visualise QC metrics
    for name, adata in [('E13', E13), ('E14', E14)]:
        fig, axs = plt.subplots(1, 4, figsize=(12, 3))
        fig.suptitle('Covariates for filtering: ' + name)

        sns.distplot(adata.obs['total_counts'], kde=False, ax=axs[0])
        sns.distplot(
            adata.obs['total_counts'][adata.obs['total_counts'] < 20000],
            kde=False,
            bins=40,
            ax=axs[1],
        )
        sns.distplot(adata.obs['n_genes_by_counts'], kde=False, bins=60, ax=axs[2])
        sns.distplot(
            adata.obs['n_genes_by_counts'][adata.obs['n_genes_by_counts'] < 4000],
            kde=False,
            bins=60,
            ax=axs[3],
        )
    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Integration of datasets')
    parser.add_argument('--integrate_homs', type=bool, default=True, help='integrate HOMs or days')
    parser.add_argument('--datasets_path', type=str, default='../ann_data/exp_04', help='the path to the datasets that will be integrated')
    parser.add_argument('--integration_file', type=str, default='../ann_data/integration_HOMs.h5ad',
                        help='path and name of the h5ad integration file')
    parser.add_argument('--plot_qc_metrics', type=bool, default=False, help='visualize QC metrics')
    parser.add_argument('--write_to_file', type=bool, default=False, help='Create an h5ad file for the integrated data')
    args = parser.parse_args()

    # Load E13 and E14 datasets
    E13 = sc.read(Path(args.datasets_path, 'E13_hom_norm_variable_genes.h5ad'))
    E14 = sc.read(Path(args.datasets_path, 'E14_hom_norm_variable_genes.h5ad'))

    if args.plot_qc_metrics:
        visualise_qc_metrics(E13, E14)

    # Use scanorama to integrate the datasets. Scanorama returns two lists, one for the integrated
    # embeddings and one for the corrected counts, for each dataset
    adatas = [E13, E14]
    integrated, corrected = scanorama.correct_scanpy(adatas, return_dimred=True)  # Integration and batch correction

    # Concatenate the two datasets and save the integrated embeddings in adata_integr.obsm['scanorama_embedding'].
    adata_integr = E13.concatenate(E14, uns_merge='unique', batch_key='batch')
    adata_integr.obs.batch = list(map(lambda x: 'E13' if x == '0' else 'E14', list(adata_integr.obs.batch)))

    embedding = np.concatenate(integrated, axis=0)
    adata_integr.obsm['scanorama_embedding'] = embedding
    # Compute UMAP to visualize the results and qualitatively assess the data integration task
    sc.pp.neighbors(adata_integr, use_rep='scanorama_embedding', n_neighbors=15)
    sc.tl.leiden(adata_integr, key_added='integr_clusters')
    sc.tl.umap(adata_integr, min_dist=0.5)
    sc.pl.umap(adata_integr, color='integr_clusters', legend_loc='on data', palette=sc.pl.palettes.default_20)
    sc.pl.umap(adata_integr, color='batch', palette=sc.pl.palettes.default_20)

    if args.write_to_file:
        # Write data to file
        adata_integr.write(args.integration_file)
        print('{} integration file saved'.format(args.integration_file))
