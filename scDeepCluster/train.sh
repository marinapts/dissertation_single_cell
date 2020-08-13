python code/scDeepCluster.py \
--data_file ../ann_data/E13_hom_variable_genes.h5ad \
--pretrain_epochs 400 \
--save_dir experiments/E13_hom_predictions/results/ \
--ae_weight_file experiments/E13_hom_predictions/integration_weights.h5 \
--n_clusters 5 \
--labels_file ../ann_data/annotations/E13_hom_predictions.csv \
# --ae_weights True


python code/scDeepCluster_latent.py \
--data_file ../ann_data/E13_hom_variable_genes.h5ad \
--scDeepCluster_weights experiments/E13_hom_predictions/results/scDeepCluster_model_final.h5 \
--latent_output bottleneck/E13_hom_predictions.csv \
--n_clusters 5
