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
    parser = argparse.ArgumentParser(description='Integration of HOMs datasets')
    parser.add_argument('--dataset_type', type=str, default='variable', help='"variable" or "all"')
    parser.add_argument('--plot_qc_metrics', action='store_true', help='visualize QC metrics')
    parser.add_argument('--write_to_file', action='store_true', help='Create an h5ad file for the integrated data')
    args = parser.parse_args()

    # Load E13 and E14 datasets
    E13_data = 'E13_hom_' + args.dataset_type + '_genes.h5ad'
    E14_data = 'E14_hom_' + args.dataset_type + '_genes.h5ad'
    print('E13_data', E13_data)
    E13 = sc.read(Path('ann_data', E13_data))
    E14 = sc.read(Path('ann_data', E14_data))

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
    sc.pp.neighbors(adata_integr, use_rep='scanorama_embedding', n_neighbors=50)
    sc.tl.leiden(adata_integr, key_added='integr_clusters')
    sc.tl.umap(adata_integr, min_dist=0.5)
    sc.pl.umap(adata_integr, color='integr_clusters', legend_loc='on data', palette=sc.pl.palettes.default_20)
    sc.pl.umap(adata_integr, color='batch', palette=sc.pl.palettes.default_20)

    if args.write_to_file:
        # Write data to file
        integration_file = 'integration_HOMs_' + args.dataset_type + '.h5ad'
        integration_file_path = Path('ann_data', integration_file)
        adata_integr.write(integration_file_path)
        print('{} integration file saved'.format(integration_file_path))
