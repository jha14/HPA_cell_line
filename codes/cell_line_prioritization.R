pacman::p_load(tidyverse, data.table, ggpubr, xlsx, fgsea, ggpubr, ggrepel, ggdist, ComplexHeatmap, ggalluvial)


# data preparation --------------------------------------------------------

tcga_exp <- 
  fread("unzip -p data/rna_cancer_sample.tsv.zip")

tcga_exp$Gene %>% unique %>% length  # n = 19760
tcga_exp$Sample %>% unique %>% length  # n = 7932
tcga_exp$Cancer %>% unique %>% length  # n = 21

tcga_cancer_exp_avg <- tcga_exp %>% dplyr::select(Gene, Cancer, FPKM) %>%
  group_by(Cancer, Gene) %>% summarise(fpkm = mean(FPKM, na.rm = T))
saveRDS(tcga_cancer_exp_avg, file = "data/tcga_cancer_exp_avg.Rds")
rm(tcga_exp, tcga_cancer_exp_avg)



hpa_tcga <- list()

hpa_tcga$tcga_disease_match <- c(BLCA = "Bladder Cancer",
                                 BRCA = "Breast Cancer",
                                 CESC = "Cervical Cancer",
                                 COAD = "Colon/Colorectal Cancer",
                                 GBM  = "Brain Cancer",
                                 HNSC = "Head and Neck Cancer",
                                 KICH = "Kidney Cancer",
                                 KIRC = "Kidney Cancer",
                                 KIRP = "Kidney Cancer",
                                 LIHC = "Liver Cancer",
                                 LUAD = "Lung Cancer",
                                 LUSC = "Lung Cancer",
                                 OV   = "Ovarian Cancer",
                                 PAAD = "Pancreatic Cancer",
                                 PRAD = "Prostate Cancer",
                                 READ = "Colon/Colorectal Cancer",
                                 SKCM = "Skin Cancer",
                                 STAD = "Gastric Cancer",
                                 TGCT = "Testis Cancer",
                                 THCA = "Thyroid Cancer",
                                 UCEC = "Endometrial/Uterine Cancer")


hpa_tcga$tcga_cancer_exp_avg <- readRDS("data/tcga_cancer_exp_avg.Rds") %>% 
  spread(Cancer, fpkm) %>% column_to_rownames("Gene")



hpa_tcga$cell_line_annotation <- exp_landscape$combined_annotation_disease
hpa_tcga$cell_line_annotation$primary_or_metastasis[is.na(hpa_tcga$cell_line_annotation$primary_or_metastasis)] <- "NA"
hpa_tcga$cell_line_annotation$primary_or_metastasis <- factor(hpa_tcga$cell_line_annotation$primary_or_metastasis, 
                                                              levels = c("Primary", "Metastasis", "NA"))

hpa_tcga$cell_line_exp <- 
  exp_landscape$combine %>% 
  spread(sample, ntpm) %>% 
  as.data.frame %>% 
  column_to_rownames("ensg_id")

# align TCGA and cell line genes
sum(rownames(hpa_tcga$tcga_cancer_exp_avg) %in% rownames(hpa_tcga$cell_line_exp))
hpa_tcga$cell_line_exp_aligned <- hpa_tcga$cell_line_exp[rownames(hpa_tcga$tcga_cancer_exp_avg),]
identical(rownames(hpa_tcga$tcga_cancer_exp_avg), rownames(hpa_tcga$cell_line_exp_aligned))



# calculate TCGA cell line correlation ------------------------------------

# TGCT Teratoma
cor(hpa_tcga$tcga_cancer_exp_avg[,"TGCT"], hpa_tcga$cell_line_exp_aligned[,"CVCL_L280"], method = "spearman")

hpa_tcga$corr_results <- list()
for(i in names(hpa_tcga$tcga_disease_match)){
  print(i)
  temp_samples <- hpa_tcga$cell_line_annotation %>% 
    filter(primary_disease == hpa_tcga$tcga_disease_match[i]) %>% 
    dplyr::select(RRID) %>% pull
  cor_temp <- cor(hpa_tcga$tcga_cancer_exp_avg[,i], hpa_tcga$cell_line_exp_aligned[, temp_samples, drop = F], method = "spearman")
  hpa_tcga$corr_results[[i]] <- data.frame(correlation = cor_temp[1,], RRID = colnames(cor_temp)) %>% 
    arrange(-correlation) %>% mutate(rank = 1:length(cor_temp)) %>% 
    left_join(hpa_tcga$cell_line_annotation, by = "RRID")
}


temp_low <- c()
temp_high <- c()
hpa_tcga$corr_results$BLCA
for(i in names(hpa_tcga$corr_results)){
  temp_low <- c(temp_low, range(hpa_tcga$corr_results[[i]]$correlation)[1])
  temp_high <- c(temp_high, range(hpa_tcga$corr_results[[i]]$correlation)[2])
}
min(temp_low)
max(temp_high)
sapply(hpa_tcga$corr_results, nrow)



# calculate TCGA signature GSEA ------------------------------------

# load tcga cancer gene classification
hpa_tcga$tcga_cancer_gene_classification <- 
  hpa_gene_classification(readRDS("data/tcga_cancer_exp_avg.Rds"), "fpkm", "Cancer", "Gene", 4, 10) %>% 
  filter(spec_category %in% c("tissue enriched", "group enriched", "tissue enhanced"))


hpa_tcga$tcga_disease_signature <- list()

for(i in 1:nrow(hpa_tcga$tcga_cancer_gene_classification)){
  cohort_temp <- str_split(hpa_tcga$tcga_cancer_gene_classification$enriched_tissues[i], ";")[[1]]
  for(j in cohort_temp)
    hpa_tcga$tcga_disease_signature[[j]] <- c(hpa_tcga$tcga_disease_signature[[j]], hpa_tcga$tcga_cancer_gene_classification$gene[i])
}

sapply(hpa_tcga$tcga_disease_signature, length)


#prepare logfc of each cell line relative to disease mean exp
hpa_tcga$cl_exp_gsea <- 
  exp_landscape$combine %>% 
  spread(sample, ntpm) %>% 
  as.data.frame %>% 
  column_to_rownames("ensg_id")

sd_genes <- apply(hpa_tcga$cl_exp_gsea, MARGIN = 1, sd)
sum(sd_genes == 0)
hpa_tcga$cl_exp_gsea <- hpa_tcga$cl_exp_gsea[sd_genes != 0,]

zeor_exp <- apply(hpa_tcga$cl_exp_gsea, MARGIN = 1, function(x) sum(x==0))
sum(zeor_exp == 985)

cl_mean_exp <- exp_landscape$combine %>% 
  left_join(hpa_tcga$cell_line_annotation, 
            by = c("sample" = "RRID")) %>% 
  dplyr::select(ensg_id, ntpm, primary_disease) %>% 
  group_by(ensg_id, primary_disease) %>% 
  summarise(ntpm = mean(ntpm, na.rm = T)) %>% 
  group_by(ensg_id) %>% 
  summarise(ntpm = mean(ntpm, na.rm = T)) %>% 
  filter(ensg_id %in% rownames(hpa_tcga$cl_exp_gsea))

identical(rownames(hpa_tcga$cl_exp_gsea), cl_mean_exp$ensg_id)
sum(cl_mean_exp$ntpm == 0)

hpa_tcga$cl_exp_gsea <- hpa_tcga$cl_exp_gsea/cl_mean_exp$ntpm
hpa_tcga$cl_exp_gsea <- hpa_tcga$cl_exp_gsea %>% t
hpa_tcga$cl_exp_gsea <- log2(hpa_tcga$cl_exp_gsea + 1)


#calculate GSEA
hpa_tcga$gsea_res <- list()

for(i in names(hpa_tcga$tcga_disease_match)){
  print(i)
  disease_enriched_genes <- hpa_tcga$tcga_disease_signature[[i]]
  
  print(disease_enriched_genes[!(disease_enriched_genes %in% rownames(hpa_tcga$cell_line_exp))])
  temp_samples <- hpa_tcga$cell_line_annotation %>% 
    filter(primary_disease == hpa_tcga$tcga_disease_match[i]) %>% 
    dplyr::select(RRID) %>% pull
  
  temp_cl_exp <- hpa_tcga$cl_exp_gsea[temp_samples, , drop = F]
  
  disease_signature <- list()
  disease_signature[[i]] <- disease_enriched_genes
  hpa_tcga$gsea_res[[i]] <- c()
  for(j in rownames(temp_cl_exp)){
    gene_rank <- sort(temp_cl_exp[j,])
    hpa_tcga$gsea_res[[i]] <- rbind(hpa_tcga$gsea_res[[i]], 
                                    fgsea(pathways = disease_signature, 
                                          stats = gene_rank, 
                                          eps = 1e-200, 
                                          minSize = 1, 
                                          maxSize = 10000))
  }
  hpa_tcga$gsea_res[[i]] <- hpa_tcga$gsea_res[[i]] %>% 
    mutate(RRID = temp_samples, 
           padj = p.adjust(padj, method = "BH")) %>% 
    arrange(-NES) %>% 
    mutate(rank = 1:nrow(.)) %>% 
    relocate(RRID, rank, .before = pathway) %>% 
    left_join(hpa_tcga$cell_line_annotation, by = "RRID")
}



View(hpa_tcga$gsea_res$BLCA)


# combine correlation and GSEA results ------------------------------------

hpa_tcga$combined_rank <- list()

for(i in names(hpa_tcga$tcga_disease_match)){
  hpa_tcga$combined_rank[[i]] <- 
    hpa_tcga$corr_results[[i]] %>% 
    dplyr::select(uniq_id, RRID, cell_line_name, correlation, correlation_rank = rank) %>% 
    left_join(hpa_tcga$gsea_res[[i]] %>% 
                dplyr::select(uniq_id, gsea_rank = rank, NES, pval, padj, primary_disease, 
                              primary_or_metastasis, sample_collection_site, data_source), by = "uniq_id") %>% 
    mutate(overall_rank = correlation_rank + gsea_rank) %>% 
    mutate(overall_rank = rank(overall_rank, ties.method = "min")) %>% 
    arrange(overall_rank) %>% relocate(overall_rank, .after = cell_line_name) %>% 
    relocate(gsea_rank, .after = padj)
}



# visualization dotplot ---------------------------------------------------

correlation_range <- c()
nes_range <- c()
for(i in names(hpa_tcga$tcga_disease_match)){
  correlation_range <- c(correlation_range, hpa_tcga$combined_rank[[i]] %>% filter(overall_rank <= 10) %>% dplyr::select(correlation) %>% pull)
  nes_range <- c(nes_range, hpa_tcga$combined_rank[[i]] %>% filter(overall_rank <= 10) %>% dplyr::select(NES) %>% pull)
}
range(correlation_range)
range(nes_range)


plot_cell_line <- function(cohort, xbreaks, ybreaks){
  alpha_dot <- ifelse(hpa_tcga$combined_rank[[cohort]]$overall_rank < 6, 1, 0.4)
  p <- hpa_tcga$combined_rank[[cohort]] %>% 
    mutate(text_show = ifelse(overall_rank < 6, cell_line_name, NA)) %>%
    ggplot(aes(x = NES, y = correlation)) +
    geom_point(aes(colour = overall_rank < 6), alpha = alpha_dot, size = 3, shape = 16) +
    # geom_text_repel(aes(label = text_show), max.overlaps = 100, size = 5) +
    ggtitle(cohort) + theme_bw() + 
    scale_x_continuous(breaks = xbreaks) +
    scale_y_continuous(breaks = ybreaks) +
    scale_color_manual(values = c("TRUE" = "#F8766D", "FALSE" = "grey60")) +
    theme(panel.border = element_rect(colour = "black", size = 0.75), 
          # panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size = 15, colour = "black"),
          panel.background = element_blank(), 
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          axis.text = element_text(size = 13, colour = "black"), 
          legend.position = "none",
          legend.title = element_blank()) + 
    ylab("Spearman's ρ") + xlab("NES")
  return(p)
}



hpa_tcga$cell_line_plot <- list()
hpa_tcga$cell_line_plot[["BLCA"]] <- plot_cell_line("BLCA", seq(0.9, 2.1, 0.4), seq(0.74, 0.84, 0.02))
hpa_tcga$cell_line_plot[["BRCA"]] <- plot_cell_line("BRCA", seq(0.9, 2.1, 0.4), seq(0.72, 0.8, 0.02))
hpa_tcga$cell_line_plot[["CESC"]] <- plot_cell_line("CESC", seq(1.6, 2.2, 0.2), seq(0.78, 0.81, 0.01))
hpa_tcga$cell_line_plot[["COAD"]] <- plot_cell_line("COAD", seq(0.8, 2.3, 0.5), seq(0.73, 0.85, 0.03))
hpa_tcga$cell_line_plot[["GBM"]] <- plot_cell_line("GBM", seq(0.8, 2.0, 0.4), seq(0.73, 0.83, 0.02))
hpa_tcga$cell_line_plot[["HNSC"]] <- plot_cell_line("HNSC", seq(1.0, 2.5, 0.5), seq(0.76, 0.82, 0.02))
hpa_tcga$cell_line_plot[["KICH"]] <- plot_cell_line("KICH", seq(0.8, 1.4, 0.2), seq(0.72, 0.78, 0.02))
hpa_tcga$cell_line_plot[["KIRC"]] <- plot_cell_line("KIRC", seq(1.2, 2.0, 0.2), seq(0.69, 0.78, 0.03))
hpa_tcga$cell_line_plot[["KIRP"]] <- plot_cell_line("KIRP", seq(1.1, 2.3, 0.4), seq(0.72, 0.82, 0.02))
hpa_tcga$cell_line_plot[["LIHC"]] <- plot_cell_line("LIHC", seq(0.8, 2.3, 0.5), seq(0.68, 0.78, 0.02))
hpa_tcga$cell_line_plot[["LUAD"]] <- plot_cell_line("LUAD", seq(0.6, 2.4, 0.6), seq(0.64, 0.8, 0.04))
hpa_tcga$cell_line_plot[["LUSC"]] <- plot_cell_line("LUSC", seq(0.6, 2.1, 0.5), seq(0.66, 0.81, 0.05))
hpa_tcga$cell_line_plot[["OV"]] <- plot_cell_line("OV", seq(0.9, 2.1, 0.4), seq(0.74, 0.82, 0.02))
hpa_tcga$cell_line_plot[["PAAD"]] <- plot_cell_line("PAAD", seq(0.5, 2, 0.5), seq(0.7, 0.78, 0.02))
hpa_tcga$cell_line_plot[["PRAD"]] <- plot_cell_line("PRAD", seq(0.9, 2.1, 0.4), seq(0.75, 0.80, 0.01))
hpa_tcga$cell_line_plot[["READ"]] <- plot_cell_line("READ", seq(0.8, 2.4, 0.4), seq(0.73, 0.85, 0.03))
hpa_tcga$cell_line_plot[["SKCM"]] <- plot_cell_line("SKCM", seq(0.9, 2.1, 0.4), seq(0.76, 0.84, 0.02))
hpa_tcga$cell_line_plot[["STAD"]] <- plot_cell_line("STAD", seq(0.6, 2.1, 0.3), seq(0.7, 0.78, 0.02))
hpa_tcga$cell_line_plot[["THCA"]] <- plot_cell_line("THCA", seq(0.8, 2.3, 0.3), seq(0.7, 0.8, 0.02))
hpa_tcga$cell_line_plot[["UCEC"]] <- plot_cell_line("UCEC", seq(0.8, 2, 0.3), seq(0.7, 0.8, 0.02))


ggarrange(plotlist = hpa_tcga$cell_line_plot, nrow = 4, ncol = 5)
ggsave("results/Figure_4/cell_line_prioritization.pdf", width = 15, height = 12)



# TCGA cohort gene classification -----------------------------------------

hpa_tcga$tcga_gene_short <- hpa_tcga$tcga_cancer_gene_classification %>% dplyr::select(spec_category, enriched_tissues) %>% 
  filter(spec_category %in% c("group enriched", "tissue enhanced", "tissue enriched"))

hpa_tcga$tcga_gene_barplot <- data.frame(matrix(0, length(hpa_tcga$tcga_disease_match), 3), 
                                         row.names = names(hpa_tcga$tcga_disease_match))
colnames(hpa_tcga$tcga_gene_barplot) <- c("group enriched", "tissue enhanced", "tissue enriched")

for(i in 1:nrow(hpa_tcga$tcga_gene_short)){
  temp <- str_split(hpa_tcga$tcga_gene_short$enriched_tissues[i], pattern = ";", simplify = T)
  hpa_tcga$tcga_gene_barplot[temp,hpa_tcga$tcga_gene_short$spec_category[i]] <- hpa_tcga$tcga_gene_barplot[temp,hpa_tcga$tcga_gene_short$spec_category[i]] + 1
}

hpa_tcga$tcga_gene_barplot_long <- hpa_tcga$tcga_gene_barplot %>% rownames_to_column("cohort") %>% 
  pivot_longer(!cohort, names_to = "Specificity", values_to = "value") %>% 
  mutate(Specificity = factor(Specificity %>% str_to_sentence %>% str_replace("Tissue", "Cohort"), 
                              levels = c("Cohort enhanced", "Group enriched", "Cohort enriched")), cohort = cohort)


hpa_tcga$tcga_gene_barplot_long %>% 
  ggplot(aes(fill=Specificity, x=cohort, y=value)) + 
  geom_bar(stat="identity", width = 0.8) + 
  xlab("TCGA Cohort") + ylab("Number of genes") +
  theme_bw() + scale_fill_manual(values = gene_category_pal2) +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  theme(panel.border = element_blank(), 
        # panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 10, colour = "black"), 
        panel.background = element_blank(), axis.text.y = element_text(size = 10, colour = "black"), 
        axis.text.x = element_text(size = 10, colour = "black", angle = 90, hjust = 1, vjust = 0.5), 
        axis.line = element_line(), axis.ticks.length = unit(.1, "cm"),
        legend.position = "right", 
        legend.title = element_text(size = 11, colour = "black"),
        legend.text = element_text(size = 10, colour = "black")) + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 900), breaks = c(0, 200, 400, 600, 800)) + ylab("Number of genes")
ggsave("results/Figure_S6/TCGA_cohort_bar.pdf", height = 4, width = 10, dpi = 600, units = "in")



# TCGA cohort gene signature in cell line -----------------------------------------

lihc_signature = hpa_tcga$tcga_cancer_gene_classification %>% dplyr::select(gene, spec_category, enriched_tissues) %>% 
  filter(spec_category %in% c("group enriched", "tissue enhanced", "tissue enriched"))
lihc_signature <- lihc_signature[grep(pattern = "LIHC", x = lihc_signature$enriched_tissues),]

lihc_cell_line <- hpa_tcga$combined_rank$LIHC$RRID

lihc_signature_dot <- hpa_tcga$cell_line_exp_aligned[lihc_signature$gene, lihc_cell_line]
cl_mean_exp_lihc_signature <- cl_mean_exp %>% as.data.frame %>% column_to_rownames("ensg_id")
cl_mean_exp_lihc_signature <- cl_mean_exp_lihc_signature[lihc_signature$gene, , drop = F]
identical(rownames(lihc_signature_dot), rownames(cl_mean_exp_lihc_signature))

# lihc_signature_dot <- lihc_signature_dot/cl_mean_exp_lihc_signature$ntpm
lihc_signature_dot <- log2(lihc_signature_dot + 1)

lihc_signature_dot_long <- lihc_signature_dot %>% rownames_to_column("ensg_id") %>% 
  pivot_longer(!ensg_id, names_to = "lihc_cell_line", values_to = "value") %>% 
  left_join(lihc_signature, by = c("ensg_id" = "gene")) %>% 
  mutate(spec_category = factor(spec_category %>% str_to_sentence %>% str_replace("Tissue", "Cohort"), 
                                levels = c("Cohort enhanced", "Group enriched", "Cohort enriched"))) %>% 
  left_join(hpa_tcga$combined_rank$LIHC %>% dplyr::select(RRID, cell_line_name), by = c("lihc_cell_line" = "RRID")) %>% 
  mutate(cell_line_name = factor(cell_line_name, levels = hpa_tcga$combined_rank$LIHC$cell_line_name))


lihc_signature_dot_long %>% 
  ggplot(aes(x=cell_line_name, y=value)) + 
  # stat_halfeye(adjust = .5, width = .7, .width = 0, justification = -.2, point_colour = NA) +
  geom_boxplot(width = .5, outlier.shape = NA, fill = "#D1CBE5") +
  geom_jitter(color = "#D1CBE5", size = 0.4, alpha = 0.5, shape = 16) +
  # geom_jitter(aes(color = spec_category), size = 0.3, alpha = 0.2, shape = 16) +
  scale_color_manual(values = gene_category_pal2) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  scale_y_continuous(expand = c(0,0), limits = c(-0.1, 15.5)) +
  ylab(expression(paste("Log"[2], "(nTPM+1)"))) +
  xlab("Liver cancer cell lines") +
  ggtitle("TCGA-LIHC signature") +
  theme(panel.border = element_blank(), 
        panel.background = element_blank(),
        # panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 15, colour = "black", face = "bold"), 
        axis.title = element_text(size = 15, colour = "black"), 
        axis.text.y = element_text(size = 13, colour = "black"), 
        axis.text.x = element_text(size = 13, colour = "black", angle = 90, hjust = 1, vjust = 0.5), 
        axis.line = element_line(), axis.ticks.length = unit(.1, "cm"),
        legend.position = "right", 
        legend.title = element_text(size = 11, colour = "black"),
        legend.text = element_text(size = 10, colour = "black"))
ggsave("results/Figure_4/LIHC_cell_line.pdf", height = 4, width = 8.3, dpi = 600, units = "in")



# ALB in liver cancer cell line -----------------------------------------

lihc_signature = hpa_tcga$tcga_cancer_gene_classification %>% dplyr::select(gene, spec_category, enriched_tissues) %>% 
  filter(spec_category %in% c("group enriched", "tissue enhanced", "tissue enriched"))
lihc_signature <- lihc_signature[grep(pattern = "LIHC", x = lihc_signature$enriched_tissues),]

lihc_cell_line <- hpa_tcga$combined_rank$LIHC$RRID
lihc_signature_bar <- 
  hpa_tcga$cell_line_exp_aligned["ENSG00000163631", lihc_cell_line] %>% 
  t %>% as.data.frame

lihc_signature_bar <- log2(lihc_signature_bar + 1)

lihc_signature_bar$cell_line <- rownames(lihc_signature_bar)
lihc_signature_bar <- 
  lihc_signature_bar %>% 
  left_join(hpa_tcga$combined_rank$LIHC %>% dplyr::select(RRID, cell_line_name), by = c("cell_line" = "RRID")) %>% 
  mutate(cell_line_name = factor(cell_line_name, levels = hpa_tcga$combined_rank$LIHC$cell_line_name))
colnames(lihc_signature_bar)[1] = "value"

ggplot(lihc_signature_bar, aes(y=value, x=cell_line_name)) + 
  geom_bar(position = "dodge", stat = "identity", fill = "#D1CBE5", width = 0.7) + theme_bw() + 
  ylab(expression(paste("Log"[2], "(nTPM+1) ALB"))) +
  ggtitle("ALB expression") +
  scale_y_continuous(expand = c(0,0), limits = c(0,15.5)) +
  xlab("Liver cancer cell lines") +
  theme(panel.border = element_blank(), 
        # panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 15, colour = "black", face = "bold"), 
        axis.title = element_text(size = 15, colour = "black"), 
        panel.background = element_blank(), axis.text.y = element_text(size = 13, colour = "black"), 
        axis.text.x = element_text(size = 13, colour = "black", angle = 90, hjust = 1, vjust = 0.5), 
        axis.line = element_line(), axis.ticks.length = unit(.1, "cm"),
        legend.position = "right", 
        legend.title = element_text(size = 11, colour = "black"),
        legend.text = element_text(size = 10, colour = "black"))
ggsave("results/Figure_4/LIHC_cell_line_ALB.pdf", height = 4, width = 8.3, dpi = 600, units = "in")



# TCGA cohort liver cancer cell line correlation -----------------------------------------

liver_cell_line_exp <- 
  hpa_tcga$cell_line_exp_aligned[,hpa_tcga$cell_line_annotation %>% filter(primary_disease == "Liver Cancer") %>% dplyr::select(RRID) %>% pull]
liver_cell_line_exp <- log2(liver_cell_line_exp + 1)
liver_cell_line_exp <- liver_cell_line_exp %>% 
  rownames_to_column("Gene")

lihc_corr_plot <- 
  readRDS("data/tcga_cancer_exp_avg.Rds") %>% 
  filter(Cancer == "LIHC") %>% 
  mutate(LIHC = log2(fpkm+1)) %>% 
  # rename(LIHC = fpkm) %>% 
  left_join(liver_cell_line_exp) %>% 
  left_join(lihc_signature %>% dplyr::select(Gene = gene, spec_category))
lihc_corr_plot$spec_category[is.na(lihc_corr_plot$spec_category)] = "not elevated"
lihc_corr_plot <- lihc_corr_plot %>% 
  mutate(spec_category = factor(spec_category %>% str_to_sentence %>% str_replace("Tissue", "Cohort"), 
                                levels = c("Cohort enhanced", "Group enriched", "Cohort enriched", "Not elevated"))) %>% 
  mutate(size = ifelse(spec_category == "Not elevated", 1, 2)) %>% 
  mutate(transparent = ifelse(spec_category == "Not elevated", 0.3, 0.75))
View(hpa_tcga$combined_rank$LIHC)
# range(lihc_corr_plot$SRR8616135)

ggplot(lihc_corr_plot, aes(x = LIHC, y = CVCL_0336)) + 
  geom_point(aes(color = spec_category), size = lihc_corr_plot$size,
             alpha = lihc_corr_plot$transparent, shape = 16) +
  theme_bw(base_size = 14) +
  scale_color_manual(values = gene_category_pal3) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 14.3)) + scale_y_continuous(expand = c(0,0), limits = c(0, 16.8)) +
  theme(panel.background = element_rect(colour = "black",size = 1),aspect.ratio=1,
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 14, colour = "black"), axis.title.x = element_text(size = 15, colour = "black"), 
        axis.title.y = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 14, colour = "black"),
        legend.position = "none") + 
  xlab(expression(paste("Log"[2], "(FPKM+1) LIHC"))) + ylab(expression(paste("Log"[2], "(nTPM+1) Huh-7")))
ggsave(file = "results/Figure_4/LIHC_cell_line_huh7.pdf", dpi = 600, width = 4, height = 4)

ggplot(lihc_corr_plot, aes(x = LIHC, y = CVCL_0027)) + 
  geom_point(aes(color = spec_category), size = lihc_corr_plot$size,
             alpha = lihc_corr_plot$transparent, shape = 16) +
  theme_bw(base_size = 14) +
  scale_color_manual(values = gene_category_pal3) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 14.3)) + scale_y_continuous(expand = c(0,0), limits = c(0, 16.8)) +
  theme(panel.background = element_rect(colour = "black",size = 1),aspect.ratio=1,
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 14, colour = "black"), axis.title.x = element_text(size = 15, colour = "black"), 
        axis.title.y = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 14, colour = "black"),
        legend.position = "none") + 
  xlab(expression(paste("Log"[2], "(FPKM+1) LIHC"))) + ylab(expression(paste("Log"[2], "(nTPM+1) Hep-G2")))
ggsave(file = "results/Figure_4/LIHC_cell_line_hepg2.pdf", dpi = 600, width = 4, height = 4)


ggplot(lihc_corr_plot, aes(x = LIHC, y = CVCL_0077)) + 
  geom_point(aes(color = spec_category), size = lihc_corr_plot$size,
             alpha = lihc_corr_plot$transparent, shape = 16) +
  theme_bw(base_size = 14) +
  scale_color_manual(values = gene_category_pal3) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 14.3)) + scale_y_continuous(expand = c(0,0), limits = c(0, 16.8)) +
  theme(panel.background = element_rect(colour = "black",size = 1),aspect.ratio=1,
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 14, colour = "black"), axis.title.x = element_text(size = 15, colour = "black"), 
        axis.title.y = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 14, colour = "black"),
        legend.position = "none") + 
  xlab(expression(paste("Log"[2], "(FPKM+1) LIHC"))) + ylab(expression(paste("Log"[2], "(nTPM+1) SNU-398")))
ggsave(file = "results/Figure_4/LIHC_cell_line_SNU398.pdf", dpi = 600, width = 4, height = 4)


# Cell line TCGA correlation plot -------------------------------------------------------------

cl_tcga_cor <- list()
cl_tcga_cor$cell_line_disease_exp <- 
  exp_landscape$combine_disease %>% 
  spread(primary_disease, ntpm) %>% 
  as.data.frame %>% 
  column_to_rownames("ensg_id")

cl_tcga_cor$tcga_disease_exp <- 
  readRDS("data/tcga_cancer_exp_avg.Rds") %>% 
  spread(Cancer, fpkm) %>% column_to_rownames("Gene")

sum(rownames(cl_tcga_cor$tcga_disease_exp) %in% rownames(cl_tcga_cor$cell_line_disease_exp))
cl_tcga_cor$cell_line_disease_exp <- cl_tcga_cor$cell_line_disease_exp[rownames(cl_tcga_cor$tcga_disease_exp),]

cl_tcga_cor$cor <- cor(cl_tcga_cor$cell_line_disease_exp, cl_tcga_cor$tcga_disease_exp, method = "spearman")

cl_tcga_cor$tcga_order <- 
  c("COAD", "READ", "STAD", "SKCM", "HNSC", "BLCA", "CESC", "LUSC", "LUAD",
    "BRCA", "OV", "UCEC", "GBM", "PRAD", "KIRP", "KIRC", "KICH", 
    "LIHC", "PAAD",  "THCA", "TGCT")

cl_tcga_cor$cell_line_order <- 
  c(unique(hpa_tcga$tcga_disease_match[cl_tcga_cor$tcga_order][1:19]), 
    "Bile Duct Cancer", "Esophageal Cancer", "Thyroid Cancer", 
    "Sarcoma", "Bone Cancer", "Rhabdoid", "Gallbladder Cancer", "Lymphoma",
    "Leukemia", "Myeloma", "Neuroblastoma", "Testis Cancer")

cl_tcga_cor$cor <- 
  cl_tcga_cor$cor[cl_tcga_cor$cell_line_order, cl_tcga_cor$tcga_order]

rownames(cl_tcga_cor$cor) <- rownames(cl_tcga_cor$cor) %>% str_remove(" Cancer")

column_ha = HeatmapAnnotation(
  cell_line_diseas = rownames(cl_tcga_cor$cor),
  col = list(cell_line_diseas = disease_colors_combine2))
row_ha = rowAnnotation(TCGA_disease = colnames(cl_tcga_cor$cor),
                       col = list(TCGA_disease = tcga_colors))

pdf("results/Figure_3/tcga_cl_cor.pdf", width = 9, height = 8)
Heatmap(t(cl_tcga_cor$cor), show_column_names = T, show_row_names = T, 
        row_title = "TCGA cohorts", column_title = "Cell line diseases",
        width = nrow(cl_tcga_cor$cor)*unit(4.5, "mm"), 
        height = ncol(cl_tcga_cor$cor)*unit(4.5, "mm"),
        col = circlize::colorRamp2(c(0.68, 0.867), c("white", "#B2182B")),
        bottom_annotation = column_ha,
        cluster_rows = F, cluster_columns = F,
        # border = "white",
        right_annotation = row_ha,
        # col = colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B")),
        heatmap_legend_param = list(title = "Spearman's ρ", at = c(0.7, 0.75, 0.8, 0.85), direction = "horizontal"))
dev.off()



# Cell line TCGA gene classification alluvial aligned -------------------------------------------------------------

cl_tcga_aligned_alluvial <- list()
cl_tcga_aligned_alluvial$cl_classification <- 
  hpa_gene_classification(exp_landscape$combine_disease %>% filter(primary_disease %in% hpa_tcga$tcga_disease_match),
                          "ntpm", "primary_disease", "ensg_id", 4, 10)

identical(cl_tcga_aligned_alluvial$cl_classification$gene, exp_landscape$hpa_summary$Ensembl)

cl_tcga_aligned_alluvial$alluvial <- exp_landscape$hpa_summary %>%
  dplyr::select(Ensembl, cancer_specificity = 'RNA cancer specificity') %>%
  mutate(cancer_specificity = str_replace(cancer_specificity, pattern = "cancer|Cancer", replacement = "tissue") %>% str_to_sentence,
         cl_disease_tcga = cl_tcga_aligned_alluvial$cl_classification$spec_category %>% str_to_sentence)

table(cl_tcga_aligned_alluvial$alluvial$cancer_specificity)
table(cl_tcga_aligned_alluvial$alluvial$cancer_specificity) %>% sum
table(cl_tcga_aligned_alluvial$alluvial$cl_disease_tcga)
table(cl_tcga_aligned_alluvial$alluvial$cl_disease_tcga) %>% sum

cl_tcga_aligned_alluvial$alluvial <- cl_tcga_aligned_alluvial$alluvial[complete.cases(cl_tcga_aligned_alluvial$alluvial),]
cl_tcga_aligned_alluvial$alluvial_long <- 
  to_lodes_form(cl_tcga_aligned_alluvial$alluvial %>% dplyr::select(Ensembl, cl_disease_tcga, cancer_specificity), key = "category", axes = c(2,3))
cl_tcga_aligned_alluvial$alluvial_long$stratum <- 
  factor(cl_tcga_aligned_alluvial$alluvial_long$stratum, 
         levels = c("Tissue enriched", "Group enriched", "Tissue enhanced", "Low tissue specificity", "Not detected"))
ggplot(cl_tcga_aligned_alluvial$alluvial_long, aes(x = category, stratum = stratum, alluvium = alluvium, y = alluvium, fill = stratum)) +
  geom_lode() + geom_flow(curve_type = "cubic", alpha = 2/3, width = 0.45) + geom_stratum(alpha = 1, width = 0.45) + 
  theme_minimal() + scale_fill_manual(values = gene_category_pal) + guides(fill=guide_legend(ncol=3)) +
  theme(legend.position = "bottom", panel.grid = element_blank(), axis.text.x = element_text(colour = "black", size = 10),
        axis.text.y = element_blank(), axis.title.y = element_blank())
ggsave("results/Figure_S4/flow_specificity_cl_tcga_aligned.pdf", width = 6, height = 6, dpi = 600, bg = "white")



# Immune GO terms TCGA -------------------------------------------------------------

immune_go_genes <- list()
immune_go_genes$tcell <- go_specificity$cl_vs_cancer_cl@result$geneID[7] %>% str_split("/")
immune_go_genes$tcell <- immune_go_genes$tcell[[1]]
immune_go_genes$bcell <- go_specificity$cl_vs_cancer_cl@result$geneID[15] %>% str_split("/")
immune_go_genes$bcell <- immune_go_genes$bcell[[1]]

immune_go_genes$tcell_tcga <- readRDS("data/tcga_cancer_exp_avg.Rds") %>% 
  filter(Gene %in% immune_go_genes$tcell) %>% 
  mutate(fpkm = log2(fpkm + 1))
immune_go_genes$bcell_tcga <- readRDS("data/tcga_cancer_exp_avg.Rds") %>% 
  filter(Gene %in% immune_go_genes$bcell) %>% 
  mutate(fpkm = log2(fpkm + 1))
immune_go_genes$tcell_cl <- exp_landscape$combine_disease %>% 
  filter(ensg_id %in% immune_go_genes$tcell) %>% 
  mutate(ntpm = log2(ntpm + 1), primary_disease = str_remove(primary_disease, " Cancer"))
immune_go_genes$bcell_cl <- exp_landscape$combine_disease %>% 
  filter(ensg_id %in% immune_go_genes$bcell) %>% 
  mutate(ntpm = log2(ntpm + 1), primary_disease = str_remove(primary_disease, " Cancer"))

immune_go_genes$tcell_cl %>% 
  ggplot(aes(x=primary_disease, y=ntpm)) + 
  geom_boxplot(aes(fill = primary_disease), width = .6, outlier.shape = NA) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 11.5)) +
  ylab(expression(paste("Log"[2], "(nTPM+1)"))) +
  scale_fill_manual(values = disease_colors_combine2) +
  xlab("Cell line diseases") +
  ggtitle("T cell activation genes") +
  theme(panel.border = element_blank(), 
        panel.background = element_blank(),
        aspect.ratio = 0.3,
        # panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 15, colour = "black", face = "bold"), 
        axis.title = element_text(size = 15, colour = "black"), 
        axis.text.y = element_text(size = 13, colour = "black"), 
        axis.text.x = element_text(size = 13, colour = "black", angle = 90, hjust = 1, vjust = 0.5), 
        axis.line = element_line(), axis.ticks.length = unit(.1, "cm"),
        legend.position = "none", 
        legend.title = element_text(size = 11, colour = "black"),
        legend.text = element_text(size = 10, colour = "black"))
ggsave("results/Figure_S4/T_cell_genes_cl.pdf", height = 5, width = 7, dpi = 600, units = "in")

immune_go_genes$tcell_tcga %>% 
  ggplot(aes(x=Cancer, y=fpkm)) + 
  geom_boxplot(aes(fill = Cancer), width = .6, outlier.shape = NA) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 8.2)) +
  ylab(expression(paste("Log"[2], "(FPKM+1)"))) +
  scale_fill_manual(values = tcga_colors) +
  xlab("TCGA cohorts") +
  ggtitle("T cell activation genes") +
  theme(panel.border = element_blank(), 
        panel.background = element_blank(),
        aspect.ratio = 0.3,
        # panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 15, colour = "black", face = "bold"), 
        axis.title = element_text(size = 15, colour = "black"), 
        axis.text.y = element_text(size = 13, colour = "black"), 
        axis.text.x = element_text(size = 13, colour = "black", angle = 90, hjust = 1, vjust = 0.5), 
        axis.line = element_line(), axis.ticks.length = unit(.1, "cm"),
        legend.position = "none", 
        legend.title = element_text(size = 11, colour = "black"),
        legend.text = element_text(size = 10, colour = "black"))
ggsave("results/Figure_S4/T_cell_genes_tcga.pdf", height = 5, width = 7, dpi = 600, units = "in")





# tau plot -------------------------------------------------------------

tau <- list()
tau$cancer <-
  readRDS("data/tcga_cancer_exp_avg.Rds") %>%
  ungroup %>%
  dplyr::rename(gene = Gene, exp = fpkm) %>%
  mutate(exp = log10(exp+1)) %>% 
  calculate_tau

tau$cell_line_disease <-
  exp_landscape$combine_disease %>%
  ungroup %>%
  dplyr::rename(gene = ensg_id, exp = ntpm) %>%
  mutate(exp = log10(exp+1)) %>% 
  calculate_tau

tau$sc_type <-
  read_tsv("data/rna_single_cell_type.tsv") %>%
  dplyr::select(gene = Gene, cell_type = 'Cell type', nTPM) %>%
  group_by(cell_type, gene) %>%
  summarise(exp = mean(nTPM, na.rm = T)) %>%
  ungroup %>%
  mutate(exp = log10(exp+1)) %>% 
  calculate_tau


tau$plot <-
  tibble(gene = Reduce(intersect, list(names(tau$cancer), names(tau$cell_line_disease), names(tau$sc_type)))) %>%
  mutate(cell_line_disease_tau = tau$cell_line_disease[gene],
         cancer_tau = tau$cancer[gene],
         sc_type_tau = tau$sc_type[gene]) %>%
  left_join(exp_landscape$gene_disease %>% dplyr::select(gene, cell_line_spec = spec_category)) %>%
  left_join(exp_landscape$hpa_summary %>% dplyr::select(gene = Ensembl,
                                                        cancer_spec = `RNA cancer specificity`,
                                                        sc_spec = `RNA single cell type specificity`)) %>%
  mutate(cell_line_spec = str_replace(cell_line_spec, pattern = "tissue", replacement = "cell line disease") %>% str_to_sentence) %>%
  filter(cell_line_spec != "Not detected", cancer_spec != "Not detected", sc_spec != "Not detected") %>%
  mutate(cell_line_spec_short = ifelse(cell_line_spec %in% c("Low cell line disease specificity"), "not elevated", "elevated"),
         cancer_spec_short = ifelse(cancer_spec %in% c("Low cancer specificity"), "not elevated", "elevated"),
         sc_spec_short = ifelse(sc_spec %in% c("Low cell type specificity"), "not elevated", "elevated"),
         cl_vs_cancer = "Others", cl_vs_sc = "Others")

tau$plot$cl_vs_cancer[tau$plot$cell_line_spec_short == "elevated" & tau$plot$cancer_spec_short == "elevated"] = "Elevated in both"
tau$plot$cl_vs_cancer[tau$plot$cell_line_spec_short == "elevated" & tau$plot$cancer_spec_short != "elevated"] = "Elevated in cell line disease"
tau$plot$cl_vs_cancer[tau$plot$cell_line_spec_short != "elevated" & tau$plot$cancer_spec_short == "elevated"] = "Elevated in TCGA cancer"

tau$plot$cl_vs_sc[tau$plot$cell_line_spec_short == "elevated" & tau$plot$sc_spec_short == "elevated"] = "Elevated in both"
tau$plot$cl_vs_sc[tau$plot$cell_line_spec_short == "elevated" & tau$plot$sc_spec_short != "elevated"] = "Elevated in cell line disease"
tau$plot$cl_vs_sc[tau$plot$cell_line_spec_short != "elevated" & tau$plot$sc_spec_short == "elevated"] = "Elevated in single-cell type"



ggplot(tau$plot, aes(x = cell_line_disease_tau, y = cancer_tau)) +
  geom_point(aes(color = cl_vs_cancer), size = 1, alpha = 0.5, shape = 16) +
  theme_bw(base_size = 14) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 1.01)) + scale_y_continuous(expand = c(0,0), limits = c(0, 1.01)) +
  guides(color=guide_legend(ncol=2, override.aes = list(size = 3))) +
  scale_color_manual(values = c("Elevated in both" = "#DC0000",
                                "Elevated in cell line disease" = "#4DBBD5",
                                "Elevated in TCGA cancer" = "#00A087",
                                "Others" = "grey")) +
  theme(panel.background = element_rect(colour = "black",size = 1), aspect.ratio=1,
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 13, colour = "black"), axis.title.x = element_text(size = 16, colour = "black"),
        axis.title.y = element_text(size = 16, colour = "black"), axis.text.y = element_text(size = 13, colour = "black"),
        legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 13, colour = "black")) +
  xlab("Cell line disease tau") + ylab("TCGA cancer tau")
ggsave(file = "results/Figure_S4/CL_vs_Cancer_tau.pdf", dpi = 600, width = 5, height = 5)


ggplot(tau$plot, aes(x = cell_line_disease_tau, y = sc_type_tau)) +
  geom_point(aes(color = cl_vs_sc), size = 1, alpha = 0.5, shape = 16) +
  theme_bw(base_size = 14) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 1.01)) + scale_y_continuous(expand = c(0,0), limits = c(0, 1.01)) +
  guides(color=guide_legend(ncol=2, override.aes = list(size = 3))) +
  scale_color_manual(values = c("Elevated in both" = "#DC0000",
                                "Elevated in cell line disease" = "#4DBBD5",
                                "Elevated in single-cell type" = "#00A087",
                                "Others" = "grey")) +
  theme(panel.background = element_rect(colour = "black",size = 1), aspect.ratio=1,
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 13, colour = "black"), axis.title.x = element_text(size = 16, colour = "black"),
        axis.title.y = element_text(size = 16, colour = "black"), axis.text.y = element_text(size = 13, colour = "black"),
        legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 13, colour = "black")) +
  xlab("Cell line disease tau") + ylab("Single-cell type tau")
ggsave(file = "results/Figure_S4/CL_vs_SCtype_tau.pdf", dpi = 600, width = 5, height = 5)




ggplot(tau$plot, aes(x = cell_line_disease_tau, y = cell_line_spec)) +
  geom_violin(aes(fill = cell_line_spec)) +
  theme_bw() + xlab("Tau score") +
  scale_x_continuous(expand = c(0,0), limits = c(0, 1.05), breaks = c(0, 0.5, 1)) +
  guides(color=guide_legend(ncol=2, override.aes = list(size = 3))) +
  scale_y_discrete(limits=c("Cell line disease enriched", "Group enriched", "Cell line disease enhanced", "Low cell line disease specificity")) +
  scale_fill_manual(values = c("Cell line disease enriched" = "#E41A1C",
                               "Group enriched" = "#FF9D00",
                               "Cell line disease enhanced" = "#984EA3",
                               "Low cell line disease specificity" = "grey40")) +
  theme(panel.border = element_blank(), aspect.ratio = 0.4,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 15, colour = "black"),
        panel.background = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 13, colour = "black"),
        axis.line = element_line(), axis.ticks.length = unit(.1, "cm"),
        legend.position = "none")
ggsave(file = "results/Figure_S4/tau_cell_line_tau.pdf", dpi = 600, width = 4, height = 2.5)



ggplot(tau$plot, aes(x = cancer_tau, y = cancer_spec)) +
  geom_violin(aes(fill = cancer_spec)) +
  theme_bw() + xlab("Tau score") +
  scale_x_continuous(expand = c(0,0), limits = c(0, 1.05), breaks = c(0, 0.5, 1)) +
  guides(color=guide_legend(ncol=2, override.aes = list(size = 3))) +
  scale_y_discrete(limits=c("Cancer enriched", "Group enriched", "Cancer enhanced", "Low cancer specificity")) +
  scale_fill_manual(values = c("Cancer enriched" = "#E41A1C",
                               "Group enriched" = "#FF9D00",
                               "Cancer enhanced" = "#984EA3",
                               "Low cancer specificity" = "grey40")) +
  theme(panel.border = element_blank(), aspect.ratio = 0.4,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 15, colour = "black"),
        panel.background = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 13, colour = "black"),
        axis.line = element_line(), axis.ticks.length = unit(.1, "cm"),
        legend.position = "none")
ggsave(file = "results/Figure_S4/tau_cancer_tau.pdf", dpi = 600, width = 4, height = 2.5)


ggplot(tau$plot, aes(x = sc_type_tau, y = sc_spec)) +
  geom_violin(aes(fill = sc_spec)) +
  theme_bw() + xlab("Tau score") +
  scale_x_continuous(expand = c(0,0), limits = c(0, 1.05), breaks = c(0, 0.5, 1)) +
  scale_y_discrete(limits=c("Cell type enriched", "Group enriched", "Cell type enhanced", "Low cell type specificity")) +
  guides(color=guide_legend(ncol=2, override.aes = list(size = 3))) +
  scale_fill_manual(values = c("Cell type enriched" = "#E41A1C",
                               "Group enriched" = "#FF9D00",
                               "Cell type enhanced" = "#984EA3",
                               "Low cell type specificity" = "grey40")) +
  theme(panel.border = element_blank(), aspect.ratio = 0.4,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 15, colour = "black"),
        panel.background = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 13, colour = "black"),
        axis.line = element_line(), axis.ticks.length = unit(.1, "cm"),
        legend.position = "none")
ggsave(file = "results/Figure_S4/tau_sc_type_tau.pdf", dpi = 600, width = 4, height = 2.5)


