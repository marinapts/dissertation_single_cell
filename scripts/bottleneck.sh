#!/bin/bash

# n_clusters=(4 6 8 10 12 14 16)

# for cl in "${n_clusters[@]}"
# do
fig_dir="./figures/bottleneck/E13_hom_predictions/"
echo "$cl clusters:  fig_dif=$fig_dir"

python analyse_bottleneck.py \
--dataset E13_hom \
--dataset_type variable \
--n_neighbors 15 \
--min_dist 0.1 \
--top_n_genes 100 \
--n_clusters 5 \
--fig_dir $fig_dir
# --write_to_file
# done
