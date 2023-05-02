pacman::p_load(tidyverse, readxl, ComplexHeatmap, progeny, ggsci, scales, circlize, decoupleR, ggpubr, Hmisc, ggdist)

# data preparation --------------------------------------------------------

progeny_analysis <- list()
progeny_analysis$exp <- read_tsv("data/hpa_ccle_combined_vst.tsv.gz") %>%
  as.data.frame %>% column_to_rownames("ensg_id") %>% as.matrix

progeny_analysis$annotation <- read_tsv("data/hpa_ccle_pathway_combined_annotation.tsv")

progeny_analysis$annotation <-
  progeny_analysis$annotation %>%
  arrange(primary_disease)

progeny_analysis$exp <- progeny_analysis$exp[,progeny_analysis$annotation$cell_line_name]

identical(colnames(progeny_analysis$exp), progeny_analysis$annotation$cell_line_name)


# progeny calculation -----------------------------------------------------

progeny_analysis$progeny_model = get_progeny(organism = 'human', top = 100) %>%
  dplyr::rename(mor = weight)


aa <- Sys.time()
progeny_analysis$pathways_results <-
  decouple(progeny_analysis$exp, progeny_analysis$progeny_model, minsize = 0)
Sys.time() - aa


progeny_analysis$pathways_res_mat <-
  progeny_analysis$pathways_results %>%
  filter(statistic=='consensus') %>%
  pivot_wider_profile(id_cols = source, names_from = condition, values_from = score) %>%
  as.matrix

progeny_analysis$pathways_res_mat_z <- abs(progeny_analysis$pathways_res_mat) < 1

sum(colnames(progeny_analysis$pathways_res_mat) %in% progeny_analysis$annotation$cell_line_name)
progeny_analysis$pathways_res_mat <- progeny_analysis$pathways_res_mat[,progeny_analysis$annotation$cell_line_name]
progeny_analysis$pathways_res_mat_z <- progeny_analysis$pathways_res_mat_z[,progeny_analysis$annotation$cell_line_name]
identical(colnames(progeny_analysis$pathways_res_mat), progeny_analysis$annotation$cell_line_name)
identical(colnames(progeny_analysis$pathways_res_mat_z), progeny_analysis$annotation$cell_line_name)



progeny_analysis$annotation <- read_tsv("data/hpa_ccle_pathway_combined_annotation.tsv")
progeny_analysis$annotation <-
  progeny_analysis$annotation %>%
  arrange(primary_disease)

names(progeny_analysis)[-3] <- 
  paste0("cell_line_", names(progeny_analysis)[-3])

progeny_analysis$cell_line_pathways_res_mat <- progeny_analysis$cell_line_pathways_res_mat[,progeny_analysis$cell_line_annotation$cell_line_name]
progeny_analysis$cell_line_pathways_res_mat_z <- progeny_analysis$cell_line_pathways_res_mat_z[,progeny_analysis$cell_line_annotation$cell_line_name]
identical(colnames(progeny_analysis$cell_line_pathways_res_mat), progeny_analysis$cell_line_annotation$cell_line_name)
identical(colnames(progeny_analysis$cell_line_pathways_res_mat_z), progeny_analysis$cell_line_annotation$cell_line_name)

progeny_analysis$cell_line_pathways_res_mat %>% 
  t %>% 
  as.data.frame %>% 
  rownames_to_column("cell_line") %>% 
  arrange(-dplyr::desc(cell_line)) %>% 
write_tsv(file = "results/progeny/progeny_results_cell_line.tsv")



# PROGENy UMAP ------------------------------------------------------------

progeny_analysis$cell_line_pca_progeny <- 
  progeny_analysis$cell_line_pathways_res_mat %>% do_pca() %>% 
  get_pca_scores(use_R2cum_PCselection = T, R2cum_lim = 0.8) 

progeny_analysis$cell_line_umap_progeny <- 
  progeny_analysis$cell_line_pca_progeny %>% 
  do_umap() %>% mutate(sample = rownames(.))

progeny_analysis$cell_line_umap_progeny %>% 
  left_join(progeny_analysis$cell_line_annotation, by = c("sample" = "cell_line_name")) %>%
  mutate(data_source = factor(data_source, levels = c("CCLE", "HPA", "Combined"))) %>% 
  ggplot(aes(UMAP1, UMAP2, color = primary_disease)) +
  geom_point(aes(shape = data_source), alpha = 0.8, size = 1.5) +
  theme_bw() + scale_color_manual(values = disease_colors_combine) +
  guides(color="none", shape = guide_legend(override.aes = list(size = 3))) +
  ggtitle("UMAP PROGENy") +
  theme(panel.border = element_rect(colour = "black", size = 1), 
        axis.title = element_blank(), axis.ticks = element_blank(),
        axis.text = element_blank(), 
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
        aspect.ratio=1,
        plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), 
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15),
        legend.position = c(0.15, 0.1))
ggsave(file = "results/Figure_S6/UMAP_PROGENy.pdf", width = 5, height = 5.3, units = "in", dpi = 600)
dev.off()



# heatmap -------------------------------------------------------------

progeny_analysis$cell_line_annotation_heatmap <- 
  progeny_analysis$cell_line_annotation %>% 
  filter(!is.na(primary_disease) & primary_disease != "Unknown")


row_ha <- rowAnnotation(Disease = progeny_analysis$cell_line_annotation_heatmap$primary_disease,
                        col = list(Disease = disease_colors_combine1))

pdf("results/Figure_6/progeny_heatmap.pdf", width = 8, height = 10)
draw(Heatmap(t(progeny_analysis$cell_line_pathways_res_mat[,progeny_analysis$cell_line_annotation_heatmap$cell_line_name]),
             show_column_names = T, show_row_names = F, left_annotation = row_ha,
             width = 15*unit(5, "mm"), height = 35*unit(5, "mm"),
             col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
             row_names_gp = grid::gpar(fontsize = 10),
             cluster_columns = F, cluster_rows = F,
             heatmap_legend_param = list(title = "Z-score",
                                         at = c(-2, -1, 0, 1, 2),
                                         direction = "horizontal")),
     annotation_legend_side = "bottom")
dev.off()



# progeny correlation ----------------------------------

progeny_analysis$cell_line_progeny_cor <- cor(t(progeny_analysis$cell_line_pathways_res_mat), method = "spearman")

pdf("results/Figure_S6/progeny_cor.pdf", width = 6, height = 4.5)
Heatmap(progeny_analysis$cell_line_progeny_cor, show_column_names = T, show_row_names = T, 
        width = ncol(progeny_analysis$cell_line_progeny_cor)*unit(5, "mm"), 
        height = nrow(progeny_analysis$cell_line_progeny_cor)*unit(5, "mm"),
        col = colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B")),
        heatmap_legend_param = list(title = "Spearman's ρ", at = c(-1, -0.5, 0, 0.5, 1), direction = "horizontal"))
dev.off()



savehistory()
save.image()

