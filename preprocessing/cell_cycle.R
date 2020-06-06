library(Hmisc)
library(ggplot2)

mouse_model = 'hom'
# mouse_model = 'het'
rds_dir = paste0('./rds_', mouse_model, '/')
plots_dir = paste0('./plots_', mouse_model, '/')

# Load datasets
E13 = readRDS(file = paste0(rds_dir, 'umap_E13_no_cajal_retzius.rds'))
E14 = readRDS(file = paste0(rds_dir, 'umap_E14_no_cajal_retzius.rds'))


# Load cell cycle genes from a regev publication. These are human genes so need to be
# formatted to mouse genes. Then we separate them to S phase and G2M phase markers.
human.cc <- scan('regev_lab_cell_cycle_genes.txt', what = 'character')
mouse.cc <- human.cc %>% tolower() %>% capitalize() %>% unlist
m.s <- mouse.cc[1:43]  # S: DNA synthesis
m.g2s <- mouse.cc[44:97]  # G2: Growth and preparation for mitosis

# CellCycleScoring stores S and G2/M scores in object meta data, along with the predicted classification of each cell
# in either G2M, S or G1 phase.
E13 = CellCycleScoring(object = E13, s.features = m.s, g2m.features = m.g2m, set.ident = F)
E14 = CellCycleScoring(object = E14, s.features = m.s, g2m.features = m.g2m, set.ident = F)

# Colour cells by cell cycle phase on UMAP plots
plot_E13 = DimPlot(SetIdent(E13, value = E13@meta.data$Phase), reduction = 'umap')
plot_E14 = DimPlot(SetIdent(E14, value = E14@meta.data$Phase), reduction = 'umap')
cell_cycle_phase_plot = plot_E13 + plot_E14
ggsave(file=paste0(plots_dir, 'E13_E14_umap_cell_cycle_phase.pdf'), plot=cell_cycle_phase_plot, width=25, units='cm')


# There are three major cell types in these datasets:
# 1. Neural progenitors, marker Pax6
# 2. Intermediate progenitors, marker Eomes
# 3. Post-mitotic neurons, marker Tbr1.

# See all three markers and the developmental trajectory on UMAP plots
feature_plot_E13 = FeaturePlot(E13, reduction = 'umap', features = c('Pax6', 'Eomes', 'Tbr1'))
ggsave(file=paste0(plots_dir, 'E13_', mouse_model,  '_dev_trajectory.pdf'), plot=feature_plot_E13, width=30, height=30, units='cm')

feature_plot_E14 = FeaturePlot(E14, reduction = 'umap', features = c('Pax6', 'Eomes', 'Tbr1'))
ggsave(file=paste0(plots_dir, 'E14_', mouse_model,  '_dev_trajectory.pdf'), plot=feature_plot_E14, width=30, height=30, units='cm')
