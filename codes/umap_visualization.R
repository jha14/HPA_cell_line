pacman::p_load(tidyverse, ggrepel, readxl, ggdist, circlize, dendextend)
source("codes/hpa_functions.R")

ccle_hpa_umap <- list()

# annotation table --------------------------------------------------------

ccle_hpa_umap$hpa_annotation <- read_excel("data/hpa_annotation.xlsx") %>% 
  select(HPA_ID, RRID, cell_line_name, hpa_name, primary_disease) %>% 
  dplyr::rename(uniq_id = HPA_ID) %>% mutate(data_source = "HPA")

ccle_hpa_umap$ccle_annotation <- read_excel("data/ccle_2019_annotation.xlsx") %>% 
  select(sra_Run, RRID, cell_line_name, primary_disease) %>% 
  dplyr::rename(uniq_id = sra_Run) %>% mutate(data_source = "CCLE")

ccle_hpa_umap$combined_annotation <- bind_rows(ccle_hpa_umap$hpa_annotation %>% select(-hpa_name),
                                               ccle_hpa_umap$ccle_annotation)


# nTPM expression -----------------------------------------

ccle_hpa_umap$hpa <- read_tsv("data/hpa_celline_all_samples_gene_tpm_103.tsv.gz") %>% 
  group_by(ensg_id, sample) %>% summarise(ntpm_avg = mean(ntpm)) %>% spread(sample, ntpm_avg) %>% 
  select(ensg_id, ccle_hpa_umap$hpa_annotation$hpa_name)
identical(names(ccle_hpa_umap$hpa)[-1], ccle_hpa_umap$hpa_annotation$hpa_name)
ccle_hpa_umap$hpa <- ccle_hpa_umap$hpa %>% rename_with(~c("ensg_id", ccle_hpa_umap$hpa_annotation$uniq_id))

ccle_hpa_umap$ccle <- read_tsv("data/ccle_2019_exp_103.tsv.gz") %>% 
  select(ensg_id, sample, ntpm) %>% spread(sample, ntpm) %>% 
  select(ensg_id, ccle_hpa_umap$ccle_annotation$uniq_id)
identical(names(ccle_hpa_umap$ccle)[-1], ccle_hpa_umap$ccle_annotation$uniq_id)

identical(ccle_hpa_umap$hpa$ensg_id, ccle_hpa_umap$ccle$ensg_id)

ccle_hpa_umap$combined_exp = ccle_hpa_umap$hpa %>% inner_join(ccle_hpa_umap$ccle, by = "ensg_id")


# check expression and annotation -----------------------------------------

identical(names(ccle_hpa_umap$combined_exp)[-1], ccle_hpa_umap$combined_annotation$uniq_id)

intersect_samples <- ccle_hpa_umap$combined_annotation %>% select(RRID) %>% table()
intersect_samples <- intersect_samples[intersect_samples==2]
label_show <- ccle_hpa_umap$combined_annotation$cell_line_name
label_show[!(ccle_hpa_umap$combined_annotation$RRID %in% names(intersect_samples))] <- NA
ccle_hpa_umap$combined_annotation <- ccle_hpa_umap$combined_annotation %>% 
  mutate(label = label_show) %>% mutate(size = ifelse(is.na(label_show), 1, 3))


# pca and umap calculation ------------------------------------------------

pca_cell_lines <- 
  ccle_hpa_umap$combined_exp %>% 
  column_to_rownames("ensg_id") %>% 
  remove_genes(rm_not_expressed = T, 
               rm_no_variance = T) %>% 
  scale_data(logp1_scale = F, 
             zscore_scale = T) %>% 
  do_pca() %>% 
  get_pca_scores(use_R2cum_PCselection = T,
                 R2cum_lim = 0.8) 

umap_cell_lines <- 
  pca_cell_lines %>% 
  do_umap() %>% mutate(uniq_id = rownames(.))


# visualization ------------------------------------------------

umap_cell_lines %>% 
  left_join(ccle_hpa_umap$combined_annotation, by = "uniq_id") %>% 
  ggplot(aes(UMAP1, UMAP2, color = data_source)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  xlab("UMAP-1") + ylab("UMAP-2") +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(panel.border = element_rect(colour = "black", size = 2), 
        axis.title = element_text(size = 20), 
        axis.ticks = element_blank(),
        axis.text = element_blank(), 
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
        aspect.ratio=1,
        legend.text = element_text(size = 20),
        legend.title = element_blank(),
        legend.position = c(0.2, 0.88))
ggsave(file = "results/Figure_S2/HPA_CCLE_data_source.pdf", width = 4, height = 4, units = "in", dpi = 600)
dev.off()


umap_cell_lines_intersect <- 
  umap_cell_lines %>% 
  left_join(ccle_hpa_umap$combined_annotation, by = "uniq_id")
umap_cell_lines_intersect %>% 
  ggplot(aes(UMAP1, UMAP2, color = label)) +
  geom_point(aes(shape = data_source), size = umap_cell_lines_intersect$size, alpha = ifelse(umap_cell_lines_intersect$size == 1, 0.5, 1)) +
  geom_text_repel(aes(label = label), max.overlaps = 100) +
  xlab("UMAP-1") + ylab("UMAP-2") +
  theme_bw() + 
  guides(color="none", shape = guide_legend(override.aes = list(size = 3))) +
  theme(panel.border = element_rect(colour = "black", size = 1), 
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 13),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(), 
        aspect.ratio=1,
        legend.text = element_text(size = 15),
        legend.title = element_blank(),
        legend.position = c(0.9, 0.1))
ggsave(file = "results/Figure_S2/CCLE_HPA_intersect_samples.pdf", width = 8, height = 8, units = "in", dpi = 600)
dev.off()


umap_cell_lines %>% 
  left_join(ccle_hpa_umap$combined_annotation, by = "uniq_id") %>% 
  filter(data_source == "CCLE") %>% 
  ggplot(aes(UMAP1, UMAP2, color = primary_disease, group = primary_disease)) +
  geom_point(aes(shape = data_source), alpha = 0.8, size = 3) +
  geom_text(data = . %>% group_by(primary_disease) %>% 
              summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2)),
            aes(UMAP1, UMAP2, label = primary_disease, colour = primary_disease),
            inherit.aes = F, size = 6) +
  xlab("UMAP-1") + ylab("UMAP-2") +
  guides(shape = "none", color = guide_legend(ncol=3, override.aes = list(size = 4))) +
  theme_bw() + scale_color_manual(values = disease_colors_ccle) +
  theme(panel.border = element_rect(colour = "black", size = 1), 
        axis.title = element_text(size = 18), 
        axis.text = element_text(size = 15), 
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
        aspect.ratio=1,
        legend.text = element_text(size = 15),
        legend.title = element_blank(),
        legend.position = "left")
ggsave(file = "results/Figure_1/CCLE_HPA_disease.pdf", width = 16, height = 8, units = "in", dpi = 600)
dev.off()



# circular dendrogram -----------------------------------------------------

cir_dendro <- list()

cir_dendro$annot <- 
  ccle_hpa_umap$combined_annotation %>% 
  mutate(label_show = paste0(label, "(", data_source, ")"))
cir_dendro$annot$label_show[is.na(cir_dendro$annot$label)] = paste0("    ", 1:1022)
cir_dendro$annot$primary_disease[is.na(cir_dendro$annot$primary_disease)] = "Non-Cancerous"

cir_dendro$exp <- ccle_hpa_umap$combined_exp

cir_dendro$cor <- 1 - cor(cir_dendro$exp[,-1], method = "spearman")

cir_dendro$hc <- 
  cir_dendro$cor %>% 
  as.dist %>% 
  hclust %>% 
  as.dendrogram

cir_dendro$annot <- 
  cir_dendro$annot[match(labels(cir_dendro$hc), cir_dendro$annot$uniq_id),] %>% 
  mutate(color = disease_colors_combine[primary_disease]) %>% 
  mutate(source_color = case_when(
    data_source == "HPA" ~ "#00BFC4",
    data_source == "CCLE" ~ "#F8766D"))

labels(cir_dendro$hc) <- cir_dendro$annot$label_show


pdf(file = "results/Figure_S3/dendrogram.pdf", height = 12, width = 12)
cir_dendro$hc %>% 
  circlize_dendrogram(labels_track_height = NA,
                      dend_track_height = 0.5,
                      labels = T)  
dev.off()



# cell line correlation for same disease -----------------------------------------------------

cl_cor <- 
  list()

cl_cor$annot <- 
  ccle_hpa_umap$combined_annotation

cl_cor$exp <- 
  ccle_hpa_umap$combined_exp

cl_cor$all_cor <- 
  1 - cir_dendro$cor

cl_cor$all_cor_long <- 
  cl_cor$all_cor

cl_cor$all_cor_long[lower.tri(cl_cor$all_cor_long, diag = T)] = 0

cl_cor$all_cor_long <- 
  cl_cor$all_cor_long %>% 
  as.data.frame %>% 
  rownames_to_column("sample1") %>% 
  pivot_longer(!sample1, names_to = "sample2", values_to = "spearman_corr") %>% 
  filter(spearman_corr != 0) %>% 
  left_join(cl_cor$annot %>% dplyr::select(sample1 = uniq_id, source1 = primary_disease)) %>% 
  left_join(cl_cor$annot %>% dplyr::select(sample2 = uniq_id, source2 = primary_disease)) %>% 
  filter(source1 != "Unknown", source2 != "Unknown") %>% 
  filter(source1 != "Non-Cancerous", source2 != "Non-Cancerous") %>% 
  mutate(source = ifelse(source1 == source2, source1, "Rest"))

cl_cor$all_cor_long %>% 
  dplyr::select(spearman_corr, source) %>% 
  filter(source == "Rest") %>% 
  dplyr::select(spearman_corr) %>% 
  pull %>% 
  median

cl_cor$all_cor_long %>% 
  ggplot(aes(x = source, y = spearman_corr, fill = source)) +
  # stat_halfeye(adjust = .5, width = .7, .width = 0, justification = -.2, point_colour = NA) +
  geom_boxplot(width = .5, outlier.shape = NA) +
  theme_bw() + 
  geom_hline(yintercept = 0.8452428, linetype="dashed", color = "grey60") +
  scale_x_discrete(limits = c(sort(unique(cl_cor$all_cor_long$source))[-22], "Rest")) +
  scale_fill_manual(values = disease_colors_combine3) +
  scale_y_continuous(limits = c(0.735, 1.00), breaks = seq(0.75, 1, 0.05))+
  # scale_y_continuous(breaks = seq(0.7, 1, 0.01))+
  theme(panel.border = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 17, colour = "black"),
        panel.background = element_blank(), 
        axis.text = element_text(size = 13, colour = "black"), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle=315, hjust=0),
        plot.margin = margin(0.3, 0.65, 0.3, 0.3, "cm"),
        legend.position = "none",
        legend.title = element_blank()) + ylab("Spearman's p")
ggsave("results/Figure_S3/box_plot_cell_line_cor_disease.pdf", height = 4.5, width = 11, dpi = 600, units = "in")
dev.off()


savehistory()
save.image()

