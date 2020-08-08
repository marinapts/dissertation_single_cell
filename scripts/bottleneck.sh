#!/bin/bash

# n_clusters=(4 6 8 10 12 14 16)
n_clusters=(10)
for cl in "${n_clusters[@]}"
do
    # fig_dir="./figures/bottleneck/E14_hom/clusters_${cl}/"
    fig_dir="./figures/bottleneck/integration/"
    echo "$cl clusters:  fig_dif=$fig_dir"

    python analyse_bottleneck.py \
    --dataset integration \
    --dataset_type variable \
    --n_neighbors 30 \
    --min_dist 0.5 \
    --top_n_genes 100 \
    --n_clusters $cl \
    --fig_dir $fig_dir \
    --write_to_file
done
