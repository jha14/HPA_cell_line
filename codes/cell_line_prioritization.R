pacman::p_load(tidyverse, ggpubr, xlsx, fgsea, ggrepel, ggdist, ggalluvial, foreach, doMC, org.Hs.eg.db)


# data preparation --------------------------------------------------------

hpa_tcga <- list()

hpa_tcga$tcga_disease_match <- c(BLCA = "Bladder Cancer",
                                 BRCA = "Breast Cancer",
                                 CHOL = "Bile Duct Cancer",
                                 CESC = "Cervical Cancer",
                                 COAD = "Colon/Colorectal Cancer",
                                 ESCA = "Esophageal Cancer",
                                 GBM  = "Brain Cancer",
                                 HNSC = "Head and Neck Cancer",
                                 KICH = "Kidney Cancer",
                                 KIRC = "Kidney Cancer",
                                 KIRP = "Kidney Cancer",
                                 LAML = "Leukemia",
                                 LGG  = "Brain Cancer",
                                 LIHC = "Liver Cancer",
                                 LUAD = "Lung Cancer",
                                 LUSC = "Lung Cancer",
                                 OV   = "Ovarian Cancer",
                                 PAAD = "Pancreatic Cancer",
                                 PRAD = "Prostate Cancer",
                                 READ = "Colon/Colorectal Cancer",
                                 SARC = "Sarcoma",
                                 SKCM = "Skin Cancer",
                                 STAD = "Gastric Cancer",
                                 TGCT = "Testis Cancer",
                                 THCA = "Thyroid Cancer",
                                 UCEC = "Endometrial/Uterine Cancer")

hpa_tcga$cell_line_annotation <- exp_landscape$combined_annotation_disease
hpa_tcga$cell_line_annotation$primary_or_metastasis[is.na(hpa_tcga$cell_line_annotation$primary_or_metastasis)] <- "NA"
hpa_tcga$cell_line_annotation$primary_or_metastasis <- factor(hpa_tcga$cell_line_annotation$primary_or_metastasis,
                                                              levels = c("Primary", "Metastasis", "NA"))

hpa_tcga$tcga_exp <-
  readRDS("../data_publish/tcga_exp.Rds") %>%
  dplyr::select(ensg_id, sample, ntpm)

hpa_tcga$cell_line_exp <-
  exp_landscape$combine %>%
  filter(ensg_id %in% unique(hpa_tcga$tcga_exp$ensg_id)) %>%
  pivot_wider(names_from = sample, values_from = ntpm) %>%
  as.data.frame %>%
  column_to_rownames("ensg_id")


hpa_tcga$tcga_clinical <-
  readRDS("../data_publish/tcga_clinical.Rds")[["final_clin"]] %>%
  dplyr::select(project_id, barcode,
                pathologic_stage = ajcc_pathologic_stage_simple,
                subtype = subtype_selected) %>%
  mutate(project_id = str_remove(project_id, "TCGA-"))



# cell line cohort prioritization data
hpa_tcga$cl_priori_cohort <- list()

hpa_tcga$cl_priori_cohort$tcga_cancer_exp_avg <-
  hpa_tcga$tcga_exp %>%
  left_join(hpa_tcga$tcga_clinical %>%
              dplyr::select(project_id, sample = barcode)) %>%
  group_by(project_id, ensg_id) %>%
  summarise(ntpm = mean(ntpm)) %>%
  pivot_wider(names_from = project_id, values_from = ntpm) %>%
  as.data.frame %>%
  column_to_rownames("ensg_id")
dim(hpa_tcga$cl_priori_cohort$tcga_cancer_exp_avg)



# cell line stage prioritization data
hpa_tcga$cl_priori_stage <- list()

hpa_tcga$cl_priori_stage$tcga_stage_exp_avg <-
  hpa_tcga$tcga_exp %>%
  left_join(hpa_tcga$tcga_clinical %>%
              dplyr::select(project_id, sample = barcode, pathologic_stage)) %>%
  filter(!is.na(pathologic_stage)) %>%
  group_by(project_id, ensg_id, pathologic_stage) %>%
  summarise(ntpm = mean(ntpm)) %>%
  mutate(cohort_stage = paste(project_id, pathologic_stage)) %>%
  ungroup() %>%
  dplyr::select(ensg_id, cohort_stage, ntpm) %>%
  pivot_wider(names_from = cohort_stage, values_from = ntpm) %>%
  as.data.frame %>%
  column_to_rownames("ensg_id")
dim(hpa_tcga$cl_priori_stage$tcga_stage_exp_avg)



# cell line stage prioritization data
hpa_tcga$cl_priori_subtype <- list()

hpa_tcga$cl_priori_subtype$tcga_subtype_exp_avg <-
  hpa_tcga$tcga_exp %>%
  left_join(hpa_tcga$tcga_clinical %>%
              dplyr::select(project_id, sample = barcode, subtype)) %>%
  filter(!is.na(subtype)) %>%
  group_by(project_id, ensg_id, subtype) %>%
  summarise(ntpm = mean(ntpm)) %>%
  mutate(cohort_subtype = paste(project_id, subtype)) %>%
  ungroup() %>%
  dplyr::select(ensg_id, cohort_subtype, ntpm) %>%
  pivot_wider(names_from = cohort_subtype, values_from = ntpm) %>%
  as.data.frame %>%
  column_to_rownames("ensg_id")
dim(hpa_tcga$cl_priori_subtype$tcga_subtype_exp_avg)


hpa_tcga$tcga_exp <- NULL


# align TCGA and cell line genes
identical(rownames(hpa_tcga$cl_priori_cohort$tcga_cancer_exp_avg), rownames(hpa_tcga$cell_line_exp))
identical(rownames(hpa_tcga$cl_priori_stage$tcga_stage_exp_avg), rownames(hpa_tcga$cell_line_exp))
identical(rownames(hpa_tcga$cl_priori_subtype$tcga_subtype_exp_avg), rownames(hpa_tcga$cell_line_exp))





# TCGA cohort gene classification -----------------------------------------

hpa_tcga$tcga_sig_barplot <- list()

hpa_tcga$tcga_sig_barplot$gene_cf <- 
  exp_landscape$gene_tcga_disease

hpa_tcga$tcga_sig_barplot$gene_short <- 
  hpa_tcga$tcga_sig_barplot$gene_cf %>% 
  dplyr::select(spec_category, enriched_diseases) %>% 
  filter(spec_category %in% c("enriched", "enhanced", "group enriched")) %>% 
  mutate(enriched_diseases = str_remove_all(enriched_diseases, "TCGA-"),
         spec_category = str_to_sentence(spec_category))

hpa_tcga$tcga_sig_barplot$gene_barplot <- 
  data.frame(matrix(0, length(hpa_tcga$tcga_disease_match), 3), 
             row.names = names(hpa_tcga$tcga_disease_match))

colnames(hpa_tcga$tcga_sig_barplot$gene_barplot) <- c("Group enriched", "Enhanced", "Enriched")

for(i in 1:nrow(hpa_tcga$tcga_sig_barplot$gene_short)){
  temp <- 
    str_split(hpa_tcga$tcga_sig_barplot$gene_short$enriched_diseases[i], pattern = ";", simplify = T)
  
  hpa_tcga$tcga_sig_barplot$gene_barplot[temp,hpa_tcga$tcga_sig_barplot$gene_short$spec_category[i]] <- 
    hpa_tcga$tcga_sig_barplot$gene_barplot[temp,hpa_tcga$tcga_sig_barplot$gene_short$spec_category[i]] + 1
}

hpa_tcga$tcga_sig_barplot$barplot_long <- 
  hpa_tcga$tcga_sig_barplot$gene_barplot %>% 
  rownames_to_column("cohort") %>% 
  pivot_longer(!cohort, names_to = "Specificity", values_to = "value") %>% 
  mutate(Specificity = factor(Specificity, levels = c("Enhanced", "Group enriched", "Enriched")))


hpa_tcga$tcga_sig_barplot$barplot_long %>% 
  ggplot(aes(fill=Specificity, x=cohort, y=value)) + 
  geom_bar(stat="identity", width = 0.75) + 
  xlab("TCGA cohort") + ylab("Number of genes") +
  theme_bw() + scale_fill_manual(values = gene_category_pal[1:3]) +
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
  scale_y_continuous(expand = c(0,0), limits = c(0, 1250), breaks = c(0, 300, 600, 900, 1200)) + ylab("Number of genes")
ggsave("results/Figure_S5/TCGA_cohort_bar.pdf", height = 4, width = 10, dpi = 600, units = "in")
dev.off()



# calculate TCGA cell line correlation ------------------------------------

# TGCT Teratoma
cor(hpa_tcga$cl_priori_cohort$tcga_cancer_exp_avg[,"TGCT"], hpa_tcga$cell_line_exp[,"CVCL_L280"], method = "spearman")

hpa_tcga$cl_priori_cohort$corr_results <- list()
for(i in names(hpa_tcga$tcga_disease_match)){
    print(i)
    temp_samples <- hpa_tcga$cell_line_annotation %>% 
      filter(primary_disease == hpa_tcga$tcga_disease_match[i]) %>% 
      dplyr::select(RRID) %>% pull
    cor_temp <- cor(hpa_tcga$cl_priori_cohort$tcga_cancer_exp_avg[,i], 
                    hpa_tcga$cell_line_exp[, temp_samples, drop = F], method = "spearman")
    hpa_tcga$cl_priori_cohort$corr_results[[i]] <- 
      data.frame(correlation = cor_temp[1,], RRID = colnames(cor_temp)) %>% 
      arrange(-correlation) %>% 
      mutate(rank = 1:length(cor_temp)) %>% 
      left_join(hpa_tcga$cell_line_annotation, by = "RRID")
}


temp_low <- c()
temp_high <- c()
hpa_tcga$cl_priori_cohort$corr_results$BLCA
for(i in names(hpa_tcga$cl_priori_cohort$corr_results)){
  temp_low <- c(temp_low, range(hpa_tcga$cl_priori_cohort$corr_results[[i]]$correlation)[1])
  temp_high <- c(temp_high, range(hpa_tcga$cl_priori_cohort$corr_results[[i]]$correlation)[2])
}
min(temp_low)
max(temp_high)
sapply(hpa_tcga$cl_priori_cohort$corr_results, nrow)



# calculate TCGA signature GSEA ------------------------------------

# load tcga cancer gene classification
hpa_tcga$cl_priori_cohort$tcga_cancer_gene_classification <- 
  exp_landscape$gene_tcga_disease %>% 
  filter(spec_category %in% c("enriched", "group enriched", "enhanced"))

hpa_tcga$cl_priori_cohort$tcga_disease_signature <- list()

for(i in 1:nrow(hpa_tcga$cl_priori_cohort$tcga_cancer_gene_classification)){
  cohort_temp <- 
    str_split(hpa_tcga$cl_priori_cohort$tcga_cancer_gene_classification$enriched_diseases[i], ";")[[1]] %>% 
    str_remove("TCGA-")
  for(j in cohort_temp)
    hpa_tcga$cl_priori_cohort$tcga_disease_signature[[j]] <- 
      c(hpa_tcga$cl_priori_cohort$tcga_disease_signature[[j]], 
        hpa_tcga$cl_priori_cohort$tcga_cancer_gene_classification$gene[i])
}

sapply(hpa_tcga$cl_priori_cohort$tcga_disease_signature, length)



#prepare logfc of each cell line relative to disease mean exp

hpa_tcga$cl_priori_cohort$cl_exp_gsea <- 
  get_gsea_exp(hpa_tcga$cell_line_exp, 
               exp_landscape$combine)


#calculate GSEA
hpa_tcga$cl_priori_cohort$gsea_res <- list()

registerDoMC(8)

for(i in names(hpa_tcga$tcga_disease_match)){
  print(i)
  disease_enriched_genes <- hpa_tcga$cl_priori_cohort$tcga_disease_signature[[i]]
  
  print(disease_enriched_genes[!(disease_enriched_genes %in% rownames(hpa_tcga$cell_line_exp))])
  temp_samples <- hpa_tcga$cell_line_annotation %>% 
    filter(primary_disease == hpa_tcga$tcga_disease_match[i]) %>% 
    dplyr::select(RRID) %>% pull
  
  temp_cl_exp <- hpa_tcga$cl_priori_cohort$cl_exp_gsea[temp_samples, , drop = F]
  
  disease_signature <- list()
  disease_signature[[i]] <- disease_enriched_genes
  
  hpa_tcga$cl_priori_cohort$gsea_res[[i]] <- 
    foreach(j = rownames(temp_cl_exp), .combine = rbind, .inorder = T) %dopar%
    cell_line_gsea(temp_cl_exp[j,], disease_signature)
    
  
  hpa_tcga$cl_priori_cohort$gsea_res[[i]] <- 
    hpa_tcga$cl_priori_cohort$gsea_res[[i]] %>% 
    mutate(RRID = temp_samples, 
           padj = p.adjust(padj, method = "BH")) %>% 
    arrange(-NES) %>% 
    mutate(rank = 1:nrow(.)) %>% 
    relocate(RRID, rank, .before = pathway) %>% 
    left_join(hpa_tcga$cell_line_annotation, by = "RRID")
  
  gc()
}
















# combine correlation and GSEA results ------------------------------------

hpa_tcga$cl_priori_cohort$combined_rank <- list()

for(i in names(hpa_tcga$tcga_disease_match)){
  hpa_tcga$cl_priori_cohort$combined_rank[[i]] <- 
    hpa_tcga$cl_priori_cohort$corr_results[[i]] %>% 
    dplyr::select(uniq_id, RRID, cell_line_name, correlation, correlation_rank = rank) %>% 
    left_join(hpa_tcga$cl_priori_cohort$gsea_res[[i]] %>% 
                dplyr::select(uniq_id, gsea_rank = rank, NES, pval, padj, primary_disease, 
                              primary_or_metastasis, sample_collection_site, data_source), by = "uniq_id") %>% 
    mutate(overall_rank = correlation_rank + gsea_rank) %>% 
    mutate(overall_rank = rank(overall_rank, ties.method = "min")) %>% 
    arrange(overall_rank) %>% 
    relocate(overall_rank, .after = cell_line_name) %>% 
    relocate(gsea_rank, .after = padj)
}

View(hpa_tcga$cl_priori_cohort$combined_rank$BLCA)

if(file.exists("results/Figure_4/cell_line_prioritization_cohort.xlsx"))
  unlink("results/Figure_4/cell_line_prioritization_cohort.xlsx")
for(i in names(hpa_tcga$cl_priori_cohort$combined_rank))
  write.xlsx(hpa_tcga$cl_priori_cohort$combined_rank[[i]], 
             file = "results/Figure_4/cell_line_prioritization_cohort.xlsx", 
             sheetName = i, row.names = F, append = T)



# visualization dotplot ---------------------------------------------------

correlation_range <- c()
nes_range <- c()
for(i in names(hpa_tcga$tcga_disease_match)){
  correlation_range <-
    c(correlation_range, hpa_tcga$cl_priori_cohort$combined_rank[[i]] %>%
        filter(overall_rank <= 10) %>%
        dplyr::select(correlation) %>%
        pull)
  nes_range <- c(nes_range,
                 hpa_tcga$cl_priori_cohort$combined_rank[[i]] %>%
                   filter(overall_rank <= 10) %>%
                   dplyr::select(NES) %>%
                   pull)
}
range(correlation_range)
range(nes_range)


plot_cell_line <- function(cohort, xbreaks, ybreaks){
  alpha_dot <- ifelse(hpa_tcga$cl_priori_cohort$combined_rank[[cohort]]$overall_rank < 6, 1, 0.4)
  p <- hpa_tcga$cl_priori_cohort$combined_rank[[cohort]] %>% 
    mutate(text_show = ifelse(overall_rank < 6, cell_line_name, NA)) %>%
    ggplot(aes(x = NES, y = correlation)) +
    geom_point(aes(colour = overall_rank < 6), alpha = alpha_dot, size = 3, shape = 16) +
    geom_text_repel(aes(label = text_show), max.overlaps = 100, size = 5) +
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

range(hpa_tcga$cl_priori_cohort$combined_rank[["UCEC"]]$NES)
range(hpa_tcga$cl_priori_cohort$combined_rank[["UCEC"]]$correlation)

hpa_tcga$cl_priori_cohort$cell_line_plot <- list()
hpa_tcga$cl_priori_cohort$cell_line_plot[["BLCA"]] <- plot_cell_line("BLCA", seq(0.6, 2.2, 0.4), seq(0.74, 0.84, 0.02))
hpa_tcga$cl_priori_cohort$cell_line_plot[["BRCA"]] <- plot_cell_line("BRCA", seq(0.7, 2.2, 0.5), seq(0.71, 0.80, 0.03))
hpa_tcga$cl_priori_cohort$cell_line_plot[["CHOL"]] <- plot_cell_line("CHOL", seq(1.1, 1.9, 0.2), seq(0.73, 0.77, 0.02))
hpa_tcga$cl_priori_cohort$cell_line_plot[["CESC"]] <- plot_cell_line("CESC", seq(1.6, 2.2, 0.2), seq(0.75, 0.79, 0.01))
hpa_tcga$cl_priori_cohort$cell_line_plot[["COAD"]] <- plot_cell_line("COAD", seq(0.7, 2.3, 0.4), seq(0.73, 0.82, 0.03))
hpa_tcga$cl_priori_cohort$cell_line_plot[["ESCA"]] <- plot_cell_line("ESCA", seq(1.1, 2.3, 0.4), seq(0.74, 0.79, 0.01))
hpa_tcga$cl_priori_cohort$cell_line_plot[["GBM"]]  <- plot_cell_line("GBM",  seq(0.6, 2.0, 0.7), seq(0.71, 0.80, 0.03))
hpa_tcga$cl_priori_cohort$cell_line_plot[["HNSC"]] <- plot_cell_line("HNSC", seq(1.0, 2.4, 0.7), seq(0.75, 0.81, 0.02))
hpa_tcga$cl_priori_cohort$cell_line_plot[["KICH"]] <- plot_cell_line("KICH", seq(0.8, 1.4, 0.2), seq(0.70, 0.76, 0.02))
hpa_tcga$cl_priori_cohort$cell_line_plot[["KIRC"]] <- plot_cell_line("KIRC", seq(1.1, 2.0, 0.3), seq(0.68, 0.76, 0.02))
hpa_tcga$cl_priori_cohort$cell_line_plot[["KIRP"]] <- plot_cell_line("KIRP", seq(1.0, 2.2, 0.4), seq(0.70, 0.80, 0.02))
hpa_tcga$cl_priori_cohort$cell_line_plot[["LAML"]] <- plot_cell_line("LAML", seq(0.8, 2.0, 0.4), seq(0.70, 0.82, 0.03))
hpa_tcga$cl_priori_cohort$cell_line_plot[["LGG"]]  <- plot_cell_line("LGG",  seq(0.6, 1.8, 0.4), seq(0.66, 0.76, 0.02))
hpa_tcga$cl_priori_cohort$cell_line_plot[["LIHC"]] <- plot_cell_line("LIHC", seq(0.7, 2.2, 0.5), seq(0.67, 0.77, 0.02))
hpa_tcga$cl_priori_cohort$cell_line_plot[["LUAD"]] <- plot_cell_line("LUAD", seq(0.5, 2.3, 0.6), seq(0.65, 0.77, 0.03))
hpa_tcga$cl_priori_cohort$cell_line_plot[["LUSC"]] <- plot_cell_line("LUSC", seq(0.6, 2.1, 0.5), seq(0.67, 0.79, 0.03))
hpa_tcga$cl_priori_cohort$cell_line_plot[["OV"]]   <- plot_cell_line("OV",   seq(0.8, 2.3, 0.5), seq(0.73, 0.79, 0.02))
hpa_tcga$cl_priori_cohort$cell_line_plot[["PAAD"]] <- plot_cell_line("PAAD", seq(0.7, 1.9, 0.3), seq(0.69, 0.77, 0.02))
hpa_tcga$cl_priori_cohort$cell_line_plot[["PRAD"]] <- plot_cell_line("PRAD", seq(0.7, 2.3, 0.4), seq(0.74, 0.78, 0.01))
hpa_tcga$cl_priori_cohort$cell_line_plot[["READ"]] <- plot_cell_line("READ", seq(0.8, 2.4, 0.4), seq(0.73, 0.83, 0.02))
hpa_tcga$cl_priori_cohort$cell_line_plot[["SARC"]] <- plot_cell_line("SARC", seq(1.0, 2.0, 0.2), seq(0.74, 0.79, 0.01))
hpa_tcga$cl_priori_cohort$cell_line_plot[["SKCM"]] <- plot_cell_line("SKCM", seq(0.8, 2.3, 0.5), seq(0.75, 0.81, 0.02))
hpa_tcga$cl_priori_cohort$cell_line_plot[["STAD"]] <- plot_cell_line("STAD", seq(0.6, 2.1, 0.5), seq(0.70, 0.78, 0.02))
hpa_tcga$cl_priori_cohort$cell_line_plot[["THCA"]] <- plot_cell_line("THCA", seq(0.9, 2.1, 0.4), seq(0.70, 0.78, 0.02))
hpa_tcga$cl_priori_cohort$cell_line_plot[["UCEC"]] <- plot_cell_line("UCEC", seq(0.8, 2.0, 0.3), seq(0.70, 0.79, 0.03))


ggarrange(plotlist = hpa_tcga$cl_priori_cohort$cell_line_plot, nrow = 5, ncol = 5)
ggsave("results/Figure_4/cell_line_prioritization.pdf", width = 15, height = 15)
dev.off()



# cell line prioritization stage ----------------------------------------

hpa_tcga$cl_priori_stage$cohorts <- 
  hpa_tcga$cl_priori_stage$tcga_stage_exp_avg %>% 
  colnames %>% 
  str_split(., " ", simplify = T) %>% 
  as.data.frame %>% 
  dplyr::select(V1) %>% 
  filter(V1 != "TGCT") %>% 
  pull %>% 
  unique

hpa_tcga$cl_priori_stage$cohort_match_samples <- list()
hpa_tcga$cl_priori_stage$disease_signature <- list()
hpa_tcga$cl_priori_stage$stage_corr <- list()
hpa_tcga$cl_priori_stage$cl_exp_gsea <- list()
hpa_tcga$cl_priori_stage$gsea_res <- list()

registerDoMC(8)
for(i in hpa_tcga$cl_priori_stage$cohorts){
  print(i)
  
  hpa_tcga$cl_priori_stage$stage_corr[[i]] <- list()
  
  hpa_tcga$cl_priori_stage$cohort_match_samples[[i]] <- 
    hpa_tcga$cell_line_annotation %>% 
    filter(primary_disease == hpa_tcga$tcga_disease_match[i]) %>% 
    dplyr::select(RRID) %>% pull
  
  temp_tcga <- 
    hpa_tcga$cl_priori_stage$tcga_stage_exp_avg %>% 
    dplyr::select(grep(i, colnames(.), value = T))
  
  colnames(temp_tcga) <- 
    str_replace_all(colnames(temp_tcga), pattern = " ", replacement = "_")
  
  hpa_tcga$cl_priori_stage$disease_signature[[i]] <- get_disease_signature(temp_tcga)
  hpa_tcga$cl_priori_stage$disease_signature[[i]][["all"]] <- 
    hpa_tcga$cl_priori_stage$disease_signature[[i]] %>% 
    unlist(use.names = F) %>% 
    unique
  
  print(sapply(hpa_tcga$cl_priori_stage$disease_signature[[i]], length))
  
  for(j in colnames(temp_tcga)){
    cor_temp <- cor(temp_tcga[, j], 
                    hpa_tcga$cell_line_exp[, hpa_tcga$cl_priori_stage$cohort_match_samples[[i]], drop = F], method = "spearman")
    
    hpa_tcga$cl_priori_stage$stage_corr[[i]][[j]] <- 
      data.frame(correlation = cor_temp[1,], 
                 RRID = colnames(cor_temp)) %>% 
      arrange(-correlation) %>% 
      mutate(rank = 1:length(cor_temp)) %>% 
      left_join(hpa_tcga$cell_line_annotation, by = "RRID")
  }


  hpa_tcga$cl_priori_stage$cl_exp_gsea[[i]] <- 
    get_gsea_exp(hpa_tcga$cell_line_exp[, hpa_tcga$cl_priori_stage$cohort_match_samples[[i]], drop = F], 
                 exp_landscape$combine)
  
  
  temp_samples <- rownames(hpa_tcga$cl_priori_stage$cl_exp_gsea[[i]])
  hpa_tcga$cl_priori_stage$gsea_res[[i]] <- list()
  
  for(j in names(hpa_tcga$cl_priori_stage$disease_signature[[i]])){
    print(paste("GSEA:", j))
    disease_signature <- list()
    disease_signature[[j]] <- hpa_tcga$cl_priori_stage$disease_signature[[i]][[j]]
    hpa_tcga$cl_priori_stage$gsea_res[[i]][[j]] <- 
      foreach(k = temp_samples, .combine = rbind, .inorder = T) %dopar%
      cell_line_gsea(hpa_tcga$cl_priori_stage$cl_exp_gsea[[i]][k,], disease_signature)
    
    hpa_tcga$cl_priori_stage$gsea_res[[i]][[j]] <- 
      hpa_tcga$cl_priori_stage$gsea_res[[i]][[j]] %>% 
      mutate(RRID = temp_samples, 
             padj = p.adjust(padj, method = "BH")) %>% 
      arrange(-NES) %>% 
      mutate(rank = 1:nrow(.)) %>% 
      relocate(RRID, rank, .before = pathway) %>% 
      left_join(hpa_tcga$cell_line_annotation, by = "RRID")
    
    gc()
  }
}



hpa_tcga$cl_priori_stage$combined_rank <- list()

for(i in hpa_tcga$cl_priori_stage$cohorts){
  
  hpa_tcga$cl_priori_stage$combined_rank[[i]] <- list()
  
  for(j in names(hpa_tcga$cl_priori_stage$stage_corr[[i]])){
    
    hpa_tcga$cl_priori_stage$combined_rank[[i]][[j]] <- 
      hpa_tcga$cl_priori_stage$stage_corr[[i]][[j]] %>% 
      dplyr::select(uniq_id, RRID, cell_line_name, correlation, correlation_rank = rank) %>% 
      left_join(hpa_tcga$cl_priori_stage$gsea_res[[i]][[j]] %>% 
                  dplyr::select(uniq_id, gsea_rank = rank, stage = pathway, NES, pval, padj, primary_disease, 
                                primary_or_metastasis, sample_collection_site, data_source), by = "uniq_id") %>% 
      mutate(overall_rank = correlation_rank + gsea_rank) %>% 
      mutate(overall_rank = rank(overall_rank, ties.method = "min")) %>% 
      mutate(stage = str_replace_all(stage, "_", " ")) %>% 
      arrange(overall_rank) %>% 
      relocate(stage, overall_rank, .after = cell_line_name) %>% 
      relocate(gsea_rank, .after = padj)
  }
  
  hpa_tcga$cl_priori_stage$combined_rank[[i]] <- 
    bind_rows(hpa_tcga$cl_priori_stage$combined_rank[[i]])
}


if(file.exists("results/Figure_4/cell_line_prioritization_stage.xlsx"))
  unlink("results/Figure_4/cell_line_prioritization_stage.xlsx")
for(i in names(hpa_tcga$cl_priori_stage$combined_rank))
  write.xlsx(hpa_tcga$cl_priori_stage$combined_rank[[i]], 
             file = "results/Figure_4/cell_line_prioritization_stage.xlsx", 
             sheetName = i, row.names = F, append = T)



# cell line prioritization subtype ----------------------------------------

hpa_tcga$cl_priori_subtype$cohorts <- 
  hpa_tcga$cl_priori_subtype$tcga_subtype_exp_avg %>% 
  colnames %>% 
  str_split(., " ", simplify = T) %>% 
  as.data.frame %>% 
  dplyr::select(V1) %>% 
  filter(V1 != "TGCT") %>% 
  pull %>% 
  unique

hpa_tcga$cl_priori_subtype$cohort_match_samples <- list()
hpa_tcga$cl_priori_subtype$disease_signature <- list()
hpa_tcga$cl_priori_subtype$subtype_corr <- list()
hpa_tcga$cl_priori_subtype$cl_exp_gsea <- list()
hpa_tcga$cl_priori_subtype$gsea_res <- list()

registerDoMC(8)
for(i in hpa_tcga$cl_priori_subtype$cohorts){
  print(i)
  
  hpa_tcga$cl_priori_subtype$subtype_corr[[i]] <- list()
  
  hpa_tcga$cl_priori_subtype$cohort_match_samples[[i]] <- 
    hpa_tcga$cell_line_annotation %>% 
    filter(primary_disease == hpa_tcga$tcga_disease_match[i]) %>% 
    dplyr::select(RRID) %>% pull
  
  temp_tcga <- 
    hpa_tcga$cl_priori_subtype$tcga_subtype_exp_avg %>% 
    dplyr::select(grep(paste0("^", i), colnames(.), value = T))
  
  colnames(temp_tcga) <- 
    str_replace_all(colnames(temp_tcga), pattern = " ", replacement = "_")

  hpa_tcga$cl_priori_subtype$disease_signature[[i]] <- get_disease_signature(temp_tcga)
  hpa_tcga$cl_priori_subtype$disease_signature[[i]][["all"]] <- 
    hpa_tcga$cl_priori_subtype$disease_signature[[i]] %>% 
    unlist(use.names = F) %>% 
    unique
  
  for(j in colnames(temp_tcga)){
    cor_temp <- cor(temp_tcga[, j], 
                    hpa_tcga$cell_line_exp[, hpa_tcga$cl_priori_subtype$cohort_match_samples[[i]], drop = F], method = "spearman")
    
    hpa_tcga$cl_priori_subtype$subtype_corr[[i]][[j]] <- 
      data.frame(correlation = cor_temp[1,], 
                 RRID = colnames(cor_temp)) %>% 
      arrange(-correlation) %>% 
      mutate(rank = 1:length(cor_temp)) %>% 
      left_join(hpa_tcga$cell_line_annotation, by = "RRID")
  }
  
  print(sapply(hpa_tcga$cl_priori_subtype$disease_signature[[i]], length))
  
  hpa_tcga$cl_priori_subtype$cl_exp_gsea[[i]] <- 
    get_gsea_exp(hpa_tcga$cell_line_exp[, hpa_tcga$cl_priori_subtype$cohort_match_samples[[i]], drop = F], 
                 exp_landscape$combine)
  
  
  temp_samples <- rownames(hpa_tcga$cl_priori_subtype$cl_exp_gsea[[i]])
  hpa_tcga$cl_priori_subtype$gsea_res[[i]] <- list()
  
  for(j in names(hpa_tcga$cl_priori_subtype$disease_signature[[i]])){
    print(paste("GSEA:", j))
    disease_signature <- list()
    disease_signature[[j]] <- hpa_tcga$cl_priori_subtype$disease_signature[[i]][[j]]
    hpa_tcga$cl_priori_subtype$gsea_res[[i]][[j]] <- 
      foreach(k = temp_samples, .combine = rbind, .inorder = T) %dopar%
      cell_line_gsea(hpa_tcga$cl_priori_subtype$cl_exp_gsea[[i]][k,], disease_signature)
    
    hpa_tcga$cl_priori_subtype$gsea_res[[i]][[j]] <- 
      hpa_tcga$cl_priori_subtype$gsea_res[[i]][[j]] %>% 
      mutate(RRID = temp_samples, 
             padj = p.adjust(padj, method = "BH")) %>% 
      arrange(-NES) %>% 
      mutate(rank = 1:nrow(.)) %>% 
      relocate(RRID, rank, .before = pathway) %>% 
      left_join(hpa_tcga$cell_line_annotation, by = "RRID")
    
    gc()
  }
}



hpa_tcga$cl_priori_subtype$combined_rank <- list()

for(i in hpa_tcga$cl_priori_subtype$cohorts){
  
  hpa_tcga$cl_priori_subtype$combined_rank[[i]] <- list()
  
  for(j in names(hpa_tcga$cl_priori_subtype$subtype_corr[[i]])){
    
    hpa_tcga$cl_priori_subtype$combined_rank[[i]][[j]] <- 
      hpa_tcga$cl_priori_subtype$subtype_corr[[i]][[j]] %>% 
      dplyr::select(uniq_id, RRID, cell_line_name, correlation, correlation_rank = rank) %>% 
      left_join(hpa_tcga$cl_priori_subtype$gsea_res[[i]][[j]] %>% 
                  dplyr::select(uniq_id, gsea_rank = rank, subtype = pathway, NES, pval, padj, primary_disease, 
                                primary_or_metastasis, sample_collection_site, data_source), by = "uniq_id") %>% 
      mutate(overall_rank = correlation_rank + gsea_rank) %>% 
      mutate(overall_rank = rank(overall_rank, ties.method = "min")) %>% 
      mutate(subtype = str_replace_all(subtype, "_", " ")) %>% 
      arrange(overall_rank) %>% 
      relocate(subtype, overall_rank, .after = cell_line_name) %>% 
      relocate(gsea_rank, .after = padj)
  }
  
  hpa_tcga$cl_priori_subtype$combined_rank[[i]] <- 
    bind_rows(hpa_tcga$cl_priori_subtype$combined_rank[[i]])
}


if(file.exists("results/Figure_4/cell_line_prioritization_subtype.xlsx"))
  unlink("results/Figure_4/cell_line_prioritization_subtype.xlsx")
for(i in names(hpa_tcga$cl_priori_subtype$combined_rank))
  write.xlsx(hpa_tcga$cl_priori_subtype$combined_rank[[i]],
             file = "results/Figure_4/cell_line_prioritization_subtype.xlsx",
             sheetName = i, row.names = F, append = T)



savehistory()
save.image()




# TCGA cohort gene signature in cell line -----------------------------------------

lihc_sig_plot <- list()
lihc_sig_plot$lihc_signature <- 
  exp_landscape$gene_tcga_disease %>% 
  dplyr::select(gene, spec_category, enriched_diseases) %>% 
  filter(spec_category %in% c("group enriched", "enhanced", "enriched"),
         grepl("LIHC", enriched_diseases))


lihc_sig_plot$lihc_signature_dot <- 
  log2(hpa_tcga$cell_line_exp[lihc_sig_plot$lihc_signature$gene, hpa_tcga$cl_priori_cohort$combined_rank$LIHC$RRID] + 1)


lihc_sig_plot$lihc_signature_dot %>% 
  rownames_to_column("ensg_id") %>% 
  pivot_longer(!ensg_id, names_to = "lihc_cell_line", values_to = "value") %>% 
  left_join(lihc_sig_plot$lihc_signature, by = c("ensg_id" = "gene")) %>% 
  mutate(spec_category = factor(spec_category %>% str_to_sentence, 
                                levels = c("Enhanced", "Group enriched", "Enriched"))) %>% 
  left_join(hpa_tcga$cl_priori_cohort$combined_rank$LIHC %>% 
              dplyr::select(RRID, cell_line_name), by = c("lihc_cell_line" = "RRID")) %>% 
  mutate(cell_line_name = factor(cell_line_name, 
                                 levels = hpa_tcga$cl_priori_cohort$combined_rank$LIHC$cell_line_name)) %>% 
  ggplot(aes(x=cell_line_name, y=value)) + 
  # stat_halfeye(adjust = .5, width = .7, .width = 0, justification = -.2, point_colour = NA) +
  geom_boxplot(width = .5, outlier.shape = NA, fill = "#D1CBE5") +
  geom_jitter(color = "#D1CBE5", size = 0.4, alpha = 0.5, shape = 16) +
  # geom_jitter(aes(color = spec_category), size = 0.3, alpha = 0.2, shape = 16) +
  scale_color_manual(values = gene_category_pal[1:3]) +
  guides(color = guide_legend(override.aes = list(size = 2))) + 
  scale_y_continuous(expand = c(0,0), limits = c(-0.00000000001, 15.5)) +
  ylab(expression(paste("Log"[2], "(nTPM+1)"))) +
  xlab("Liver cancer cell lines") +
  ggtitle("TCGA-LIHC signature") +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 15, colour = "black", face = "bold"), 
        axis.title = element_text(size = 15, colour = "black"), 
        panel.background = element_blank(), axis.text.y = element_text(size = 13, colour = "black"), 
        axis.text.x = element_text(size = 13, colour = "black", angle = 90, hjust = 1, vjust = 0.5), 
        axis.line = element_line(), axis.ticks.length = unit(.1, "cm"),
        legend.position = "right", 
        legend.title = element_text(size = 11, colour = "black"),
        legend.text = element_text(size = 10, colour = "black"))
ggsave("results/Figure_5/LIHC_cell_line.pdf", height = 4, width = 8.3, dpi = 600, units = "in")
dev.off()



lihc_sig_plot$signature_gsoa <- 
  lihc_sig_plot$lihc_signature$gene %>% 
  clusterProfiler::enrichGO(OrgDb = org.Hs.eg.db, 
                            keyType = "ENSEMBL", 
                            universe = exp_landscape$gene_tcga_disease$gene,
                            ont = "BP", pAdjustMethod = "BH")


x <- lihc_sig_plot$signature_gsoa@result[c(1:10),]
gr <- do.call(rbind, strsplit(x$GeneRatio, '/'))
range(as.numeric(gr[,1])/as.numeric(gr[,2]))

plot_GO_dot(lihc_sig_plot$signature_gsoa@result[c(1:10),],
            "LIHC signature", 
            c(0.0001, 0.15), c(0, 60), c(0.06, 0.15), c(0.06, 0.09, 0.12, 0.15))
ggsave(file = "results/Figure_5/LIHC_signature_gsoa.pdf", width = 8, height = 3.5, bg = "white")
dev.off()


exp_landscape$hpa_summary$Gene[exp_landscape$hpa_summary$Ensembl %in% 
                                 Reduce(intersect, list(hyper_testing$cell_line_disease %>% filter(disease == "Liver Cancer") %>% dplyr::select(Ensembl) %>% pull,
                                                        hyper_testing$single_cell_type %>% filter(cell_type == "Hepatocytes") %>% dplyr::select(Ensembl) %>% pull,
                                                        hyper_testing$tissue %>% filter(tissue == "Liver") %>% dplyr::select(Ensembl) %>% pull,
                                                        hyper_testing$tcga_cancer %>% filter(cohort == "TCGA-LIHC") %>% dplyr::select(Ensembl) %>% pull))] %>% sort



# ALB in liver cancer cell line -----------------------------------------

lihc_sig_plot$alb <- 
  hpa_tcga$cell_line_exp["ENSG00000163631", hpa_tcga$cl_priori_cohort$combined_rank$LIHC$RRID] %>% 
  t %>% as.data.frame %>% 
  mutate(ENSG00000163631 = log2(ENSG00000163631 + 1)) %>% 
  rownames_to_column("cell_line") %>% 
  left_join(hpa_tcga$cl_priori_cohort$combined_rank$LIHC %>% 
              dplyr::select(RRID, cell_line_name), by = c("cell_line" = "RRID")) %>% 
  mutate(cell_line_name = factor(cell_line_name, 
                                 levels = hpa_tcga$cl_priori_cohort$combined_rank$LIHC$cell_line_name))

ggplot(lihc_sig_plot$alb, aes(y=ENSG00000163631, x=cell_line_name)) + 
  geom_bar(position = "dodge", stat = "identity", fill = "#D1CBE5", width = 0.7) + theme_bw() + 
  ylab(expression(paste("Log"[2], "(nTPM+1) ALB"))) +
  ggtitle("ALB expression") +
  scale_y_continuous(expand = c(0,0), limits = c(0,15.5)) +
  xlab("Liver cancer cell lines") +
  theme(panel.border = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 15, colour = "black", face = "bold"), 
        axis.title = element_text(size = 15, colour = "black"), 
        panel.background = element_blank(), axis.text.y = element_text(size = 13, colour = "black"), 
        axis.text.x = element_text(size = 13, colour = "black", angle = 90, hjust = 1, vjust = 0.5), 
        axis.line = element_line(), axis.ticks.length = unit(.1, "cm"),
        legend.position = "right", 
        legend.title = element_text(size = 11, colour = "black"),
        legend.text = element_text(size = 10, colour = "black"))
ggsave("results/Figure_5/LIHC_cell_line_ALB.pdf", height = 4, width = 8.3, dpi = 600, units = "in")
dev.off()


# TCGA cohort liver cancer cell line correlation -----------------------------------------

lihc_sig_plot$liver_cell_line_exp <- 
  hpa_tcga$cell_line_exp[,hpa_tcga$cell_line_annotation %>% 
                           filter(primary_disease == "Liver Cancer") %>% 
                           dplyr::select(RRID) %>% pull]

lihc_sig_plot$liver_cell_line_exp <- 
  log2(lihc_sig_plot$liver_cell_line_exp + 1) %>% 
  rownames_to_column("Gene")

lihc_sig_plot$lihc_corr_plot <- 
  hpa_tcga$cl_priori_cohort$tcga_cancer_exp_avg %>% 
  rownames_to_column("Gene") %>% 
  dplyr::select(Gene, LIHC) %>% 
  mutate(LIHC = log2(LIHC+1)) %>% 
  left_join(lihc_sig_plot$liver_cell_line_exp) %>% 
  left_join(lihc_sig_plot$lihc_signature %>% dplyr::select(Gene = gene, spec_category)) %>% 
  mutate(spec_category = ifelse(is.na(spec_category), "not elevated", spec_category)) %>% 
  mutate(spec_category = factor(spec_category %>% str_to_sentence, 
                                levels = c("Enhanced", "Group enriched", "Enriched", "Not elevated"))) %>% 
  mutate(size = ifelse(spec_category == "Not elevated", 1, 2)) %>% 
  mutate(transparent = ifelse(spec_category == "Not elevated", 0.3, 0.75))


range(lihc_sig_plot$lihc_corr_plot$CVCL_0336)
ggplot(lihc_sig_plot$lihc_corr_plot, aes(x = LIHC, y = CVCL_0336)) +
  geom_point(aes(color = spec_category), size = lihc_sig_plot$lihc_corr_plot$size,
             alpha = lihc_sig_plot$lihc_corr_plot$transparent, shape = 16) +
  theme_bw(base_size = 14) +
  scale_color_manual(values = gene_category_pal1) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 16.8)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 16.8)) +
  theme(panel.background = element_rect(colour = "black",size = 1),aspect.ratio=1,
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 16, colour = "black"), axis.title.x = element_text(size = 20, colour = "black"), 
        axis.title.y = element_text(size = 20, colour = "black"), axis.text.y = element_text(size = 16, colour = "black"),
        legend.position = "none") + 
  xlab(expression(paste("Log"[2], "(nTPM+1) LIHC"))) + ylab(expression(paste("Log"[2], "(nTPM+1) Huh-7")))
ggsave(file = "results/Figure_5/LIHC_cell_line_huh7.pdf", dpi = 600, width = 4, height = 4)
dev.off()


range(lihc_sig_plot$lihc_corr_plot$CVCL_0027)
ggplot(lihc_sig_plot$lihc_corr_plot, aes(x = LIHC, y = CVCL_0027)) + 
  geom_point(aes(color = spec_category), size = lihc_sig_plot$lihc_corr_plot$size,
             alpha = lihc_sig_plot$lihc_corr_plot$transparent, shape = 16) +
  theme_bw(base_size = 14) +
  scale_color_manual(values = gene_category_pal1) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 16.8)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 16.8)) +
  theme(panel.background = element_rect(colour = "black",size = 1),aspect.ratio=1,
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 16, colour = "black"), axis.title.x = element_text(size = 20, colour = "black"), 
        axis.title.y = element_text(size = 20, colour = "black"), axis.text.y = element_text(size = 16, colour = "black"),
        legend.position = "none") + 
  xlab(expression(paste("Log"[2], "(nTPM+1) LIHC"))) + ylab(expression(paste("Log"[2], "(nTPM+1) Hep-G2")))
ggsave(file = "results/Figure_5/LIHC_cell_line_hepg2.pdf", dpi = 600, width = 4, height = 4)
dev.off()


range(lihc_sig_plot$lihc_corr_plot$CVCL_0077)
ggplot(lihc_sig_plot$lihc_corr_plot, aes(x = LIHC, y = CVCL_0077)) + 
  geom_point(aes(color = spec_category), size = lihc_sig_plot$lihc_corr_plot$size,
             alpha = lihc_sig_plot$lihc_corr_plot$transparent, shape = 16) +
  theme_bw(base_size = 14) +
  scale_color_manual(values = gene_category_pal1) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 16.8)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 16.8)) +
  theme(panel.background = element_rect(colour = "black",size = 1),aspect.ratio=1,
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 16, colour = "black"), axis.title.x = element_text(size = 20, colour = "black"), 
        axis.title.y = element_text(size = 20, colour = "black"), axis.text.y = element_text(size = 16, colour = "black"),
        legend.position = "none") + 
  xlab(expression(paste("Log"[2], "(nTPM+1) LIHC"))) + ylab(expression(paste("Log"[2], "(nTPM+1) SNU-398")))
ggsave(file = "results/Figure_5/LIHC_cell_line_SNU398.pdf", dpi = 600, width = 4, height = 4)
dev.off()


range(lihc_sig_plot$lihc_corr_plot$CVCL_0077)
ggplot(lihc_sig_plot$lihc_corr_plot, aes(x = LIHC, y = CVCL_0077)) + 
  geom_point(aes(color = spec_category), size = lihc_sig_plot$lihc_corr_plot$size,
             alpha = lihc_sig_plot$lihc_corr_plot$transparent, shape = 16) +
  theme_bw(base_size = 14) +
  guides(color=guide_legend(ncol=4, override.aes = list(size = 3))) +
  scale_color_manual(values = gene_category_pal1) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 16.8)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 16.8)) +
  theme(panel.background = element_rect(colour = "black",size = 1),aspect.ratio=1,
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 16, colour = "black"), axis.title.x = element_text(size = 20, colour = "black"), 
        axis.title.y = element_text(size = 20, colour = "black"), axis.text.y = element_text(size = 16, colour = "black"),
        legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 15, colour = "black")) + 
  xlab(expression(paste("Log"[2], "(nTPM+1) LIHC"))) + ylab(expression(paste("Log"[2], "(nTPM+1) SNU-398")))
ggsave(file = "results/Figure_5/LIHC_cell_line_legend.pdf", dpi = 600, width = 4, height = 4)
dev.off()


