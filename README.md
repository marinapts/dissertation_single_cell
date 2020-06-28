# Single-cell analysis

Single-cell analysis pipeline to detect the ectopic cells that appear after inactivation of the Pax6 gene
in the cerebral cortex of mouse embryos.

We have gene expression values in single cells for the embryonic days E13 and E14 for controls and knockouts.

> Knockouts (*Homozygous mutant*): `E13_HOM`, `E14_HOM`

> Controls (*Heterozygous control*): `E13_HET`, `E14_HET`

TODO list:
- [x] Preprocessing
- [x] Assign cell types
- [ ] Differential expression analysis
- [x] scDeepCluster - AE
- [ ] sci... variational AE


## Preprocessing with Seurat
1. Quality control
2. Normalisation
3. Regression of unwanted variation (cell cycle)
4. PCA
5. UMAP
6. Marker genes for 3 major cell types (Neural progenitors, Intermediate progenitors, Post-mitotic neuron)
7. Ectopic marker genes

RMarkdown: `preprocessing/preprocess_and_cell_cycle.Rmd`


## Cell type assignment based on marker genes
Assign cell types to cells using the R package `cellassign`

RMarkdown: `preprocessing/cell_assign.Rmd`


## Differential expression analysis
@TODO: Differential expression analysis to find more marker genes for each cell type.

## Clustering using Deep Learning

### scDeepCluster
#### Plot latent space
Bottleneck has size 32. We plot it using t-SNE and UMAP. Ectopic cells are still not visible in E13_HOM.

RMarkdown: `scDeepCluster/plot_latent_space.Rmd`


