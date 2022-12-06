pacman::p_load(tidyverse, readxl, ComplexHeatmap, progeny, ggsci, scales, circlize, decoupleR, ggpubr)

load("data/progeny_analysis.RData")
# data preparation --------------------------------------------------------

# progeny_analysis <- list()
# progeny_analysis$exp <- read_tsv("data/HPA_CCLE_combined_vst.tsv.gz") %>%
#   as.data.frame %>% column_to_rownames("ensg_id") %>% as.matrix
# 
# progeny_analysis$annotation <- read_tsv("data/HPA_CCLE_pathway_combined_annotation.tsv")
# 
# progeny_analysis$annotation <-
#   progeny_analysis$annotation %>%
#   arrange(primary_disease)
# 
# progeny_analysis$exp <- progeny_analysis$exp[,progeny_analysis$annotation$cell_line_name]
# 
# identical(colnames(progeny_analysis$exp), progeny_analysis$annotation$cell_line_name)


# progeny calculation -----------------------------------------------------

# progeny_analysis$progeny_model = get_progeny(organism = 'human', top = 100) %>% 
#   dplyr::rename(mor = weight)
# 
# 
# # aa = Sys.time()
# progeny_analysis$pathways_results <-
#   decouple(progeny_analysis$exp, progeny_analysis$progeny_model, minsize = 0)
# # Sys.time() - aa


# progeny_analysis$pathways_res_mat <- 
#   progeny_analysis$pathways_results %>% 
#   filter(statistic=='consensus') %>%
#   pivot_wider_profile(id_cols = source, names_from = condition, values_from = score) %>%
#   as.matrix
# 
# progeny_analysis$pathways_res_mat_z <- abs(progeny_analysis$pathways_res_mat) < 1
# 
# sum(colnames(progeny_analysis$pathways_res_mat) %in% progeny_analysis$annotation$cell_line_name)
# progeny_analysis$pathways_res_mat <- progeny_analysis$pathways_res_mat[,progeny_analysis$annotation$cell_line_name]
# progeny_analysis$pathways_res_mat_z <- progeny_analysis$pathways_res_mat_z[,progeny_analysis$annotation$cell_line_name]
# identical(colnames(progeny_analysis$pathways_res_mat), progeny_analysis$annotation$cell_line_name)
# identical(colnames(progeny_analysis$pathways_res_mat_z), progeny_analysis$annotation$cell_line_name)



progeny_analysis$pathways_res_mat %>% t %>% as.data.frame %>% 
  rownames_to_column("Cell_line") %>% arrange(-desc(Cell_line)) %>% 
write_tsv(file = "../data_publish/progeny_results.tsv")


# PROGENy UMAP ------------------------------------------------------------

progeny_analysis$pca_progeny <- 
  progeny_analysis$pathways_res_mat %>% do_pca() %>% 
  get_pca_scores(use_R2cum_PCselection = T, R2cum_lim = 0.8) 

progeny_analysis$umap_progeny <- 
  progeny_analysis$pca_progeny %>% 
  do_umap() %>% mutate(sample = rownames(.))

progeny_analysis$umap_progeny %>% 
  left_join(progeny_analysis$annotation, by = c("sample" = "cell_line_name")) %>%
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
ggsave(file = "results/Figure_S5/UMAP_PROGENy.pdf", width = 5, height = 5.3, units = "in", dpi = 600)



# heatmap -------------------------------------------------------------

progeny_analysis$annotation_heatmap <- 
  progeny_analysis$annotation %>% 
  filter(!is.na(primary_disease) & primary_disease != "Unknown")

column_ha <- HeatmapAnnotation(Disease = progeny_analysis$annotation_heatmap$primary_disease,
                               col = list(Disease = disease_colors_combine1),
                               annotation_legend_param = list(Disease = list(ncol = 7)))
pdf("results/Figure_5/progeny_heatmap.pdf", width = 14, height = 4)
draw(Heatmap(progeny_analysis$pathways_res_mat[,progeny_analysis$annotation_heatmap$cell_line_name], 
             show_column_names = F, show_row_names = T, top_annotation = column_ha, 
             col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
             row_names_gp = grid::gpar(fontsize = 10), 
             cluster_columns = F, cluster_rows = F,
             heatmap_legend_param = list(title = "Z-score", 
                                         at = c(-2, -1, 0, 1, 2), 
                                         direction = "horizontal")),
     annotation_legend_side = "bottom")
dev.off()



row_ha <- rowAnnotation(Disease = progeny_analysis$annotation_heatmap$primary_disease,
                        col = list(Disease = disease_colors_combine1))

pdf("results/Figure_5/progeny_heatmap_2.pdf", width = 8, height = 10)
draw(Heatmap(t(progeny_analysis$pathways_res_mat[,progeny_analysis$annotation_heatmap$cell_line_name]), 
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

progeny_analysis$progeny_cor <- cor(t(progeny_analysis$pathways_res_mat), method = "spearman")

pdf("results/Figure_S5/progeny_cor.pdf", width = 6, height = 4.5)
Heatmap(progeny_analysis$progeny_cor, show_column_names = T, show_row_names = T, 
        width = ncol(progeny_analysis$progeny_cor)*unit(5, "mm"), 
        height = nrow(progeny_analysis$progeny_cor)*unit(5, "mm"),
        col = colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B")),
        heatmap_legend_param = list(title = "Spearman's ρ", at = c(-1, -0.5, 0, 0.5, 1), direction = "horizontal"))
dev.off()


# single cell line visualization -------------------------------------------

plot_progeny_pathway <- function(plot_cell_line){
  plot_progeny = progeny_analysis$pathways_res_mat %>% 
    as.data.frame %>% 
    rownames_to_column("pathway") %>% 
    as.tibble %>% 
    select(plot_cell_line, pathway) %>% 
    dplyr::rename(score = 1) %>% 
    mutate(direction = ifelse(score > 0, "up", "down")) %>% 
    mutate(sig = ifelse(abs(score) > 1, "sig", "nonsig") %>% factor(levels = c("sig", "nonsig"))) %>% 
    mutate(y = ifelse(score > 0, -max(abs(score))/20, max(abs(score))/20), 
           x = 1:14, hjust = ifelse(score > 0, 1, 0)) %>% 
    mutate(direction = paste(direction, sig, sep = "_") %>% 
             factor(levels = c("up_sig", "up_nonsig", "down_sig", "down_nonsig"))) %>% 
    mutate(size = ifelse(sig == "sig", 5, 3))
  
  yrange <- plot_progeny$score %>% abs %>% max %>% ceiling
  
  ggdotchart(plot_progeny, x = "pathway", y = "score", sorting = "none", dot.size = plot_progeny$size,
             add = "segments", xlab="PROGENy pathways", title = plot_cell_line, shape = 21) +
    geom_point(aes(fill = direction, color = sig), shape = 21, size = plot_progeny$size) + 
    geom_text(aes(label = pathway, x = x, y = y, hjust = hjust), vjust = 0.5, size = 4.5, angle = 90) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
          axis.line.x = element_blank(), legend.position = c(0.85, 1.05),
          plot.title = element_text(face = "bold", size = 18)) +
    guides(fill = "none", color = guide_legend(override.aes = list(size = c(5,3)))) +
    scale_y_continuous(limits = c(-yrange, yrange)) +
    scale_color_manual(values = c(sig = "black", nonsig = "grey"), name = NULL, labels = c("Significant", "Non-significant")) +
    scale_fill_manual(values = c(up_sig = "#B2182B", up_nonsig = "#fe8f87",
                                 down_sig = "#2166AC", down_nonsig = "#66b8ff")) +
    geom_hline(yintercept = 0) + ylab("Z-score")
  ggsave(paste0("results/PROGENy/", str_replace_all(plot_cell_line, "/", "_"), ".pdf"), width = 5, height = 3.5)
}

# plot_progeny_pathway("VCaP")
# plot_progeny_pathway("PC-3")
# plot_progeny_pathway("LNCaP clone FGC")
# plot_progeny_pathway("MDA-PCa-2b")
# plot_progeny_pathway("DU145")
# plot_progeny_pathway("22Rv1")
# plot_progeny_pathway("NCI-H660")

# plot_progeny_pathway("Hep-G2")
# plot_progeny_pathway("Huh-7")


for(i in colnames(progeny_analysis$pathways_res_mat))
  plot_progeny_pathway(i)

