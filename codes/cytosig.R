pacman::p_load(tidyverse, readxl, ComplexHeatmap, progeny, ggsci, scales, circlize, decoupleR, ggpubr, Hmisc, ggdist)


# data preparation ----------------------------------------------------

cytosig_analysis <- list()

cytosig_analysis$cell_line_annotation <- read_tsv("data/hpa_ccle_pathway_combined_annotation.tsv")

cytosig_analysis$cell_line_cytosig_results <- read.delim("data/cytosig_cell_line.Zscore")
cytosig_analysis$cell_line_cytosig_results_p <- read.delim("data/cytosig_cell_line.Pvalue")


colnames(cytosig_analysis$cell_line_cytosig_results) <- cytosig_analysis$cell_line_annotation$cell_line_name
colnames(cytosig_analysis$cell_line_cytosig_results_p) <- cytosig_analysis$cell_line_annotation$cell_line_name

cytosig_analysis$cell_line_annotation <- 
  cytosig_analysis$cell_line_annotation %>% 
  arrange(primary_disease)


cytosig_analysis$cell_line_cytosig_results <- cytosig_analysis$cell_line_cytosig_results[,cytosig_analysis$cell_line_annotation$cell_line_name]
cytosig_analysis$cell_line_cytosig_results_p <- cytosig_analysis$cell_line_cytosig_results_p[,cytosig_analysis$cell_line_annotation$cell_line_name]


identical(colnames(cytosig_analysis$cell_line_cytosig_results), cytosig_analysis$cell_line_annotation$cell_line_name)
identical(colnames(cytosig_analysis$cell_line_cytosig_results_p), cytosig_analysis$cell_line_annotation$cell_line_name)



# CytoSig UMAP ------------------------------------------------------------

cytosig_analysis$pca_cytosig <- 
  cytosig_analysis$cell_line_cytosig_results %>% 
  do_pca() %>% 
  get_pca_scores(use_R2cum_PCselection = T,
                 R2cum_lim = 0.8) 

cytosig_analysis$umap_cytosig <- 
  cytosig_analysis$pca_cytosig %>% 
  do_umap() %>% mutate(sample = rownames(.))

cytosig_analysis$umap_cytosig %>% 
  left_join(cytosig_analysis$cell_line_annotation, by = c("sample" = "cell_line_name")) %>%
  mutate(data_source = factor(data_source, levels = c("CCLE", "HPA", "Combined"))) %>% 
  ggplot(aes(UMAP1, UMAP2, color = primary_disease)) +
  geom_point(aes(shape = data_source), alpha = 0.8, size = 1.5) +
  theme_bw() + scale_color_manual(values = disease_colors_combine) +
  guides(color="none", shape = guide_legend(override.aes = list(size = 3))) +
  ggtitle("UMAP CytoSig") +
  theme(panel.border = element_rect(colour = "black", size = 1), 
        axis.title = element_blank(), axis.ticks = element_blank(),
        axis.text = element_blank(), 
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
        aspect.ratio=1,
        plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), 
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15),
        legend.position = c(0.15, 0.1))
ggsave(file = "results/Figure_S6/UMAP_CytoSig.pdf", width = 5, height = 5.3, units = "in", dpi = 600)
dev.off()


# heatmap -------------------------------------------------------------

cytosig_analysis$annotation_heatmap <- 
  cytosig_analysis$cell_line_annotation %>% 
  filter(!is.na(primary_disease) & primary_disease != "Unknown")


row_ha <- rowAnnotation(Disease = cytosig_analysis$annotation_heatmap$primary_disease,
                        col = list(Disease = disease_colors_combine1))

pdf("results/Figure_6/cytosig_heatmap.pdf", width = 12, height = 10)
draw(Heatmap(t(cytosig_analysis$cell_line_cytosig_results[,cytosig_analysis$annotation_heatmap$cell_line_name]),
             show_column_names = T, show_row_names = F, left_annotation = row_ha,
             width = 43*unit(5, "mm"), height = 35*unit(5, "mm"),
             col = colorRamp2(c(-10, 0, 10), c("#2166AC", "white", "#B2182B")),
             row_names_gp = grid::gpar(fontsize = 10),
             cluster_columns = F, cluster_rows = F,
             heatmap_legend_param = list(title = "Z-score",
                                         at = c(-10, -5, 0, 5, 10),
                                         direction = "horizontal")),
     annotation_legend_side = "bottom")
dev.off()



# cytosig correlation ----------------------------------

cytosig_analysis$cytosig_cor <- cor(cytosig_analysis$cell_line_cytosig_results %>% t, method = "spearman")

pdf("results/Figure_S6/cytosig_cor.pdf", width = 11, height = 11)
Heatmap(cytosig_analysis$cytosig_cor, show_column_names = T, show_row_names = T,
        width = ncol(cytosig_analysis$cytosig_cor)*unit(5, "mm"), 
        height = nrow(cytosig_analysis$cytosig_cor)*unit(5, "mm"),
        col = colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B")),
        heatmap_legend_param = list(title = "Spearman's ρ", at = c(-1, -0.5, 0, 0.5, 1), direction = "horizontal"))
dev.off()



# comparison between progeny and cytosig ----------------------------------

progeny_analysis$cell_line_pathways_res_mat %>% colnames %>% 
  identical(cytosig_analysis$cell_line_cytosig_results %>% colnames)
progeny_cytosig_cor <- cor(progeny_analysis$cell_line_pathways_res_mat %>% t, 
                           cytosig_analysis$cell_line_cytosig_results %>% t, method = "spearman")

pdf("results/Figure_6/progeny_cytosig_cor.pdf", width = 15, height = 6)
Heatmap(progeny_cytosig_cor, show_column_names = T, show_row_names = T,
        height = nrow(progeny_cytosig_cor)*unit(5, "mm"),
        width = ncol(progeny_cytosig_cor)*unit(5, "mm"), 
        col = colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B")),
        heatmap_legend_param = list(title = "Spearman's ρ", at = c(-1, -0.5, 0, 0.5, 1), direction = "horizontal"))
dev.off()



savehistory()
save.image()

