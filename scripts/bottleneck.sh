#!/bin/bash

python analyse_bottleneck.py \
--dataset E14_hom E13_hom \
--dataset_type variable \
--n_neighbors 30 \
--min_dist 0.5 \
--top_n_genes 100 \
--write_to_file
