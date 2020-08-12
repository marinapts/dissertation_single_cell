python code/scDeepCluster.py \
--data_file ../ann_data/E14_hom_variable_genes.h5ad \
--pretrain_epochs 400 \
--save_dir experiments/n_layers_E14/bottleneck_02/results/ \
--ae_weight_file experiments/n_layers_E14/bottleneck_02/size_64_weights.h5 \
--n_clusters 5 \
--bottleneck_size 02 \
--labels_file ../ann_data/annotations/E14_hom_anno_leiden.csv
# --labels_file ../data_qc/anno_E13_hom.csv \
# --ae_weights True


python code/scDeepCluster_latent.py \
--data_file ../ann_data/E14_hom_variable_genes.h5ad \
--scDeepCluster_weights experiments/n_layers_E14/bottleneck_02/results/scDeepCluster_model_final.h5 \
--latent_output bottleneck/num_layers_bottleneck_size/bottleneck_02.csv \
--n_clusters 5 \
--bottleneck_size 02
