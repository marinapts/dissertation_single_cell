# Single-cell analysis

Single-cell analysis pipeline to detect the ectopic cells that appear after inactivation of the Pax6 gene
in the cerebral cortex of mouse embryos.

We have gene expression values in single cells for the embryonic days E13 and E14 for controls and knockouts.

> Knockouts (*Homozygous mutant*): `E13_hom`, `E14_hom`

> Controls (*Heterozygous control*): `E13_het`, `E14_het`

TODO list:
- [x] Preprocessing
- [x] Assign cell types
- [x] Differential expression analysis
- [x] scDeepCluster - AE
- [ ] sci... variational AE


## Pipeline

### Preprocessing
Preprocess datasets that exist under the data directory. New h5ad files (all genes or only highly variable genes)
are saved in *ann_data*.

```bash
python scanpy_preprocess/preprocess.py --dataset E13_hom --keep_only_highly_variable --write_to_file
```

Setting the flag `--keep_only_highly_variable` produces the files `ann_data/E1*_hom_variable_genes.h5ad` while
if we don't pass the parameter then it produces `ann_data/E1*_hom_all_genes.h5ad`.

### Data Integration - HOMs

```bash
python scanpy_preprocess/integration.py --dataset_type variable --write_to_file
```

### Differential expression and cell-type annotation on initial data

```bash
 python scanpy_preprocess/marker_genes.py --dataset E14_hom --dataset_type variable  --update_file
```

### Analysis of the AE bottleneck
```bash
python analyse_bottleneck.py --dataset E14_hom --dataset_type variable --kmeans
```







### Preprocessing using Scanpy
`scanpy_preprocess/preprocess.ipynb`

1. Quality control
2. Normalisation
3. Regression of unwanted variation (cell cycle)
4. PCA
5. UMAP
6. Marker genes for 3 major cell types (Neural progenitors, Intermediate progenitors, Post-mitotic neuron)
7. Ectopic marker genes

Saves the preprocessed files in `ann_data/exp_04`

### Differential expression analysis
`scanpy_preprocess/marker_genes.ipynb`

Find differentiably expressed genes between the clusters.

### Cell type assignment
RMarkdown: `preprocessing/assign_cell_types.R`

Assign cell types to cells using the R package `cellassign` and use the annotations to plot the bottleneck of the
scDeepCluster autoencoder.


### scDeepCluster
`scDeepCluster/train.sh`
`scDeepCluster/cluster.sh`

We train on the preprocessed data and generate a csv file containing the latent space representation of the scDeepCluster
model. Bottleneck has size 32. We plot it using UMAP. Ectopic cells are still not visible in E13_HOM.

### Analyse bottleneck of scDeepCluster
`scDeepCluster/latent_space_HOMs.ipynb`

Uses the csv files with the annotations that were generated by `preprocessing/assign_cell_types.R`
