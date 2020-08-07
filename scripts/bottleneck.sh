#!/bin/bash

n_clusters=(4 6 8 10 12 14 16)
for cl in "${n_clusters[@]}"
do
    fig_dir="./figures/bottleneck/clusters_${cl}/"
    echo "$cl clusters:  fig_dif=$fig_dir"

    python analyse_bottleneck.py \
    --dataset E14_hom \
    --dataset_type variable \
    --n_neighbors 30 \
    --min_dist 0.5 \
    --top_n_genes 100 \
    --n_clusters $cl \
    --fig_dir $fig_dir
    # --write_to_file
done
