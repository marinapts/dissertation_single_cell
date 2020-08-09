python code/scDeepCluster.py \
--data_file ../ann_data/E14_hom_variable_genes.h5ad \
--pretrain_epochs 400 \
--save_dir experiments/E14_hom_clusters_16/results/ \
--ae_weight_file experiments/E14_hom_clusters_16/integration_weights.h5 \
--n_clusters 4
# --labels_file ../data_qc/anno_E13_hom.csv \
# --ae_weights True


python code/scDeepCluster_latent.py \
--data_file ../ann_data/E14_hom_variable_genes.h5ad \
--scDeepCluster_weights experiments/E14_hom_clusters_16/results/scDeepCluster_model_final.h5 \
--latent_output bottleneck/16_clusters_E14_hom_variable_bottleneck.csv \
--n_clusters 4
