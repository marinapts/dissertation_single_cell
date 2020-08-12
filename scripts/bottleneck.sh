#!/bin/bash

# n_clusters=(4 6 8 10 12 14 16)
# n_clusters=(5)

# n_bottleneck_size=('02' '04' '08' '16' '32' '64')
n_bottleneck_size=('02')

for bn in "${n_bottleneck_size[@]}"
do
    fig_dir="./figures/bottleneck/bottleneck_size/bottleneck_${bn}/"
    echo "$bn clusters:  fig_dif=$fig_dir"

    python analyse_bottleneck.py \
    --dataset E14_hom \
    --dataset_type variable \
    --n_neighbors 15 \
    --min_dist 0.1 \
    --top_n_genes 100 \
    --n_clusters 5 \
    --fig_dir $fig_dir \
    --bottleneck_size $bn \
    --write_to_file
done
