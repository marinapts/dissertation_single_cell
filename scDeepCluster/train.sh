python code/scDeepCluster.py \
--data_file ../ann_data/E14_hom_variable_genes.h5ad \
--pretrain_epochs 400 \
--save_dir experiments/E14_hom_var_clusters_/results/E14_hom \
--ae_weight_file experiments/E14_hom_var_clusters_/weights/E14_hom_weights.h5
# --n_clusters 4
# --labels_file ../data_qc/anno_E13_hom.csv \
# --ae_weights True
