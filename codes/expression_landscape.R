pacman::p_load(tidyverse, readxl, ggpubr, ComplexHeatmap, scales, ggsci, ggrepel, ggdendro, ggalluvial, 
               RColorBrewer, venneuler, VennDiagram, org.Hs.eg.db, ggdist)
source("codes/hpa_functions.R")


# annotation table ----------------------------------------------------

exp_landscape <- list()

exp_landscape$hpa_annotation <- read_excel("data/HPA_annotation.xlsx") %>% 
  dplyr::select(HPA_ID, RRID, cell_line_name, hpa_name, primary_disease, primary_or_metastasis, sample_collection_site) %>% 
  dplyr::rename(uniq_id = HPA_ID) %>% mutate(data_source = "HPA")

exp_landscape$ccle_annotation <- read_excel("data/CCLE_2019_annotation.xlsx") %>% 
  dplyr::select(sra_Run, RRID, cell_line_name, primary_disease, primary_or_metastasis, sample_collection_site) %>% 
  dplyr::rename(uniq_id = sra_Run) %>% mutate(data_source = "CCLE")

add_hpa_cell_lines <- exp_landscape$hpa_annotation %>% 
  filter(!(primary_disease %in% c("Non-Cancerous", "Unknown")), !(RRID %in% exp_landscape$ccle_annotation$RRID))

# add_hpa_cell_lines <- exp_landscape$hpa_annotation %>% 
#   filter(!is.na(primary_disease), !(RRID %in% exp_landscape$ccle_annotation$RRID))

overlap_cell_lines <- exp_landscape$hpa_annotation %>% 
  filter(RRID %in% exp_landscape$ccle_annotation$RRID) %>% 
  left_join(exp_landscape$ccle_annotation %>% dplyr::select(CCLE_id = uniq_id, RRID)) %>% 
  mutate(combined_uniq_id = paste0(uniq_id, "+", CCLE_id))

exp_landscape$combined_annotation <- 
  bind_rows(add_hpa_cell_lines %>% dplyr::select(-hpa_name), exp_landscape$ccle_annotation) %>% 
  left_join(overlap_cell_lines %>% dplyr::select(combined_uniq_id, RRID), by = "RRID")

exp_landscape$combined_annotation$uniq_id[!is.na(exp_landscape$combined_annotation$combined_uniq_id)] <- 
  exp_landscape$combined_annotation$combined_uniq_id[!is.na(exp_landscape$combined_annotation$combined_uniq_id)]
exp_landscape$combined_annotation$data_source[exp_landscape$combined_annotation$RRID %in% overlap_cell_lines$RRID] <- "HPA+CCLE"
exp_landscape$combined_annotation <- exp_landscape$combined_annotation %>% dplyr::select(-combined_uniq_id)

exp_landscape$combined_annotation_disease <- 
  exp_landscape$combined_annotation %>% 
  filter(!(primary_disease %in% c("Non-Cancerous", "Unknown")))



# expression data ----------------------------------------------------

exp_landscape$hpa <- 
  read_tsv("data/hpa_celline_all_samples_gene_tpm_103.tsv.gz") %>% 
  group_by(sample, ensg_id) %>% 
  summarise(ntpm = mean(ntpm)) %>% 
  mutate(sample = exp_landscape$hpa_annotation$uniq_id[match(sample, exp_landscape$hpa_annotation$hpa_name)])

exp_landscape$hpa_add <- exp_landscape$hpa %>% 
  filter(sample %in% c(add_hpa_cell_lines$uniq_id, overlap_cell_lines$uniq_id)) %>% 
  mutate(sample = exp_landscape$hpa_annotation$RRID[match(sample, exp_landscape$hpa_annotation$uniq_id)])

exp_landscape$ccle <- read_tsv("data/CCLE_2019_exp_103.tsv.gz") %>% dplyr::select(sample, ensg_id, ntpm) %>% 
  mutate(sample = exp_landscape$ccle_annotation$RRID[match(sample, exp_landscape$ccle_annotation$uniq_id)])

exp_landscape$combine <- 
  bind_rows(exp_landscape$hpa_add, exp_landscape$ccle) %>% 
  filter(sample %in% exp_landscape$combined_annotation_disease$RRID)
table(table(exp_landscape$combine$sample))

exp_landscape$combine <- exp_landscape$combine %>% 
  group_by(sample, ensg_id) %>% 
  summarise(ntpm = mean(ntpm, na.rm = T))
table(table(exp_landscape$combine$sample))
dim(exp_landscape$combine)
sum(exp_landscape$combine$sample %in% exp_landscape$combined_annotation_disease$RRID)

exp_landscape$combine_disease <- exp_landscape$combine %>% 
  left_join(exp_landscape$combined_annotation_disease %>% 
              dplyr::select(RRID, primary_disease), by = c("sample" = "RRID")) %>% 
  group_by(primary_disease, ensg_id) %>% summarise(ntpm = mean(ntpm, na.rm = T))

exp_landscape$combine_primary_disease <- exp_landscape$combine %>% 
  left_join(exp_landscape$combined_annotation_disease %>% 
              dplyr::select(RRID, primary_disease, primary_or_metastasis), 
            by = c("sample" = "RRID")) %>% 
  filter(primary_or_metastasis == "Primary") %>% 
  group_by(primary_disease, ensg_id) %>% summarise(ntpm = mean(ntpm, na.rm = T))

exp_landscape$combine_metastasis_disease <- exp_landscape$combine %>% 
  left_join(exp_landscape$combined_annotation_disease %>% 
              dplyr::select(RRID, primary_disease, primary_or_metastasis), 
            by = c("sample" = "RRID")) %>% 
  filter(primary_or_metastasis == "Metastasis") %>% 
  group_by(primary_disease, ensg_id) %>% summarise(ntpm = mean(ntpm, na.rm = T))


# calculate distribution and specificity ----------------------------------------------------

exp_landscape$gene_hpa_sample <- hpa_gene_classification(exp_landscape$hpa, "ntpm", "sample", "ensg_id", 4, 10)
exp_landscape$gene_ccle_sample <- hpa_gene_classification(exp_landscape$ccle, "ntpm", "sample", "ensg_id", 4, 10)
exp_landscape$gene_disease <- hpa_gene_classification(exp_landscape$combine_disease, "ntpm", "primary_disease", "ensg_id", 4, 10)
exp_landscape$gene_primary_disease <- hpa_gene_classification(exp_landscape$combine_primary_disease, "ntpm", "primary_disease", "ensg_id", 4, 10)
exp_landscape$gene_metastasis_disease <- hpa_gene_classification(exp_landscape$combine_metastasis_disease, "ntpm", "primary_disease", "ensg_id", 4, 10)

# Based on HPA version 21
exp_landscape$hpa_summary <- read_tsv("data/proteinatlas.tsv")

# table(exp_landscape$gene_hpa_sample$dist_category)
# table(exp_landscape$hpa_summary$`RNA cell line distribution`)


# barplot distribution ----------------------------------------------------

gene_dist_barplot <- data.frame(num_genes = c(table(exp_landscape$gene_hpa_sample$dist_category),
                                              table(exp_landscape$gene_ccle_sample$dist_category)),
                                category = rep(names(table(exp_landscape$gene_hpa_sample$dist_category)), 2),
                                data = c(rep("HPA cell line v21", 5), rep("CCLE 2019", 5)))
gene_dist_barplot$category[gene_dist_barplot$category == "detected in all"] <- "Detected in all cell lines"
gene_dist_barplot$category[gene_dist_barplot$category == "detected in many"] <- "Detected in many cell lines"
gene_dist_barplot$category[gene_dist_barplot$category == "detected in some"] <- "Detected in some cell lines"
gene_dist_barplot$category[gene_dist_barplot$category == "detected in single"] <- "Detected in a single cell line"
gene_dist_barplot$category[gene_dist_barplot$category == "not detected"] <- "Not detected"
gene_dist_barplot$category = factor(gene_dist_barplot$category,
                                    levels = c("Detected in all cell lines", 
                                               "Detected in many cell lines", 
                                               "Detected in some cell lines",
                                               "Detected in a single cell line", 
                                               "Not detected"))

ggplot(gene_dist_barplot, aes(fill=data, y=num_genes, x=category)) + 
  geom_bar(position="dodge", stat="identity") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_text(size = 15, colour = "black"),
        panel.background = element_blank(), plot.title = element_text(size = 15, colour = "black", face = "bold"),
        axis.text = element_text(colour = "black", size = 14), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle=320, hjust=0, size = 15), legend.position = "right", legend.title = element_blank(), 
        legend.text = element_text(size = 14, colour = "black")) + scale_y_continuous(expand = c(0,0), limits = c(0,7900)) +
  ylab("Number of genes") + geom_text(aes(label = num_genes), vjust = -0.3, size = 5) +
  scale_fill_manual(values = c("HPA cell line v21" = "#4DBBD5", "CCLE 2019" = "#E64B35"))
ggsave("results/Figure_2/HPA_vs_CCLE_distribution.pdf", height = 5, width = 9, dpi = 600, units = "in")



# alluvial distribution ----------------------------------------------------

identical(exp_landscape$gene_disease$gene, exp_landscape$hpa_summary$Ensembl)
identical(exp_landscape$gene_primary_disease$gene, exp_landscape$hpa_summary$Ensembl)
identical(exp_landscape$gene_metastasis_disease$gene, exp_landscape$hpa_summary$Ensembl)

gene_specificity_alluvial <- exp_landscape$hpa_summary %>%
  dplyr::select(Ensembl,
                cancer_specificity = 'RNA cancer specificity',
                sc_specificity = 'RNA single cell type specificity') %>%
  mutate(cancer_specificity = str_replace(cancer_specificity, pattern = "cancer|Cancer", replacement = "tissue") %>% str_to_sentence,
         sc_specificity = str_replace(sc_specificity, pattern = "cell type|Cell type", replacement = "tissue") %>% str_to_sentence,
         gene_disease = exp_landscape$gene_disease$spec_category %>% str_to_sentence,
         gene_primary_disease = exp_landscape$gene_primary_disease$spec_category %>% str_to_sentence,
         gene_metastasis_disease = exp_landscape$gene_metastasis_disease$spec_category %>% str_to_sentence)

table(gene_specificity_alluvial$cancer_specificity)
table(gene_specificity_alluvial$cancer_specificity) %>% sum
table(gene_specificity_alluvial$sc_specificity)
table(gene_specificity_alluvial$sc_specificity) %>% sum
table(gene_specificity_alluvial$gene_disease)
table(gene_specificity_alluvial$gene_disease) %>% sum

gene_specificity_alluvial_cancer <- gene_specificity_alluvial[complete.cases(gene_specificity_alluvial),]
gene_specificity_alluvial_long <- 
  to_lodes_form(gene_specificity_alluvial_cancer %>% 
                  dplyr::select(Ensembl, gene_disease, cancer_specificity), key = "category", axes = c(2,3))
gene_specificity_alluvial_long$stratum <- 
  factor(gene_specificity_alluvial_long$stratum, 
         levels = c("Tissue enriched", "Group enriched", "Tissue enhanced", "Low tissue specificity", "Not detected"))
ggplot(gene_specificity_alluvial_long, aes(x = category, stratum = stratum, alluvium = alluvium, y = alluvium, fill = stratum)) +
  geom_lode() + geom_flow(curve_type = "cubic", alpha = 2/3, width = 0.5) + geom_stratum(alpha = 1, width = 0.5) + 
  theme_minimal() + scale_fill_manual(values = gene_category_pal) + guides(fill=guide_legend(ncol=1)) +
  theme(legend.position = "right", panel.grid = element_blank(), axis.text.x = element_text(colour = "black", size = 10),
        axis.text.y = element_blank(), axis.title.y = element_blank())
ggsave("results/Figure_2/flow_specificity_disease.pdf", width = 6.3, height = 6, dpi = 600, bg = "white")


gene_specificity_alluvial_long <- 
  to_lodes_form(gene_specificity_alluvial %>% 
                  dplyr::select(Ensembl, gene_disease, sc_specificity), key = "category", axes = c(2,3))
gene_specificity_alluvial_long$stratum <- factor(gene_specificity_alluvial_long$stratum,
                                                 levels = c("Tissue enriched", "Group enriched", "Tissue enhanced", "Low tissue specificity", "Not detected"))
ggplot(gene_specificity_alluvial_long, aes(x = category, stratum = stratum, alluvium = alluvium, y = alluvium, fill = stratum)) +
  geom_lode() + geom_flow(curve_type = "cubic", alpha = 2/3, width = 0.5) + geom_stratum(alpha = 1, width = 0.5) + 
  theme_minimal() + scale_fill_manual(values = gene_category_pal) + guides(fill=guide_legend(ncol=1)) +
  theme(legend.position = "right", panel.grid = element_blank(), axis.text.x = element_text(colour = "black", size = 10),
        axis.text.y = element_blank(), axis.title.y = element_blank())
ggsave("results/Figure_S4/flow_specificity_sc.pdf", width = 6.3, height = 6, dpi = 600, bg = "white")




# compare elevated genes ----------------------------------------------------

identical(exp_landscape$gene_disease$gene, exp_landscape$hpa_summary$Ensembl)

compare_specificity <- 
  exp_landscape$hpa_summary %>%
  dplyr::select(Ensembl, cancer_elevation = 'RNA cancer specificity', sc_elevation = 'RNA single cell type specificity') %>%
  mutate(cancer_elevation = ifelse(cancer_elevation %in% c("Low cancer specificity", "Not detected"), "No specificity", "High specificity"),
         sc_elevation = ifelse(sc_elevation %in% c("Low cell type specificity", "Not detected"), "No specificity", "High specificity"),
         cell_line_disease_elevation = ifelse(exp_landscape$gene_disease$spec_category %in% c("low tissue specificity", "not detected"), "No specificity", "High specificity"),
         cl_vs_cancer = "Others", cl_vs_sc = "Others")

compare_specificity$cl_vs_cancer[compare_specificity$cell_line_disease_elevation == "High specificity" & compare_specificity$cancer_elevation == "High specificity"] = "High specificity in both"
compare_specificity$cl_vs_cancer[compare_specificity$cell_line_disease_elevation == "High specificity" & compare_specificity$cancer_elevation != "High specificity"] = "High specificity in cell line disease"
compare_specificity$cl_vs_cancer[compare_specificity$cell_line_disease_elevation != "High specificity" & compare_specificity$cancer_elevation == "High specificity"] = "High specificity in TCGA cancer"

compare_specificity$cl_vs_sc[compare_specificity$cell_line_disease_elevation == "High specificity" & compare_specificity$sc_elevation == "High specificity"] = "High specificity in both"
compare_specificity$cl_vs_sc[compare_specificity$cell_line_disease_elevation == "High specificity" & compare_specificity$sc_elevation != "High specificity"] = "High specificity in cell line disease"
compare_specificity$cl_vs_sc[compare_specificity$cell_line_disease_elevation != "High specificity" & compare_specificity$sc_elevation == "High specificity"] = "High specificity in single-cell type"


myCol <- brewer.pal(4, "Set2")
venn <- venn.diagram(x = list(compare_specificity %>% filter(cancer_elevation == "High specificity") %>% dplyr::select(Ensembl) %>% pull,
                              compare_specificity %>% filter(cell_line_disease_elevation == "High specificity") %>% dplyr::select(Ensembl) %>% pull),
                     filename = NULL, lty = 'blank', main.fontface = "bold",
                     fill=c(myCol[1],myCol[2]), main.fontfamily = "sans", 
                     lwd = 1, fontfamily = "sans", cat.fontfamily = "sans", cex=0.55, cat.cex=0.5,
                     margin = c(0.1,0.1,0.1,0.1), cat.default.pos = "outer",
                     cat.pos = c(-20, 20), cat.dist = c(0.04, 0.04),
                     category.names = c("TCGA cancer", "Cell line disease"))
pdf(file="results/Figure_S4/cl_cancer_elevation.pdf", width = 2, height = 2)
grid.draw(venn)
dev.off()



# GO analysis -------------------------------------------------------------

go_specificity <- list()
go_specificity$cl_vs_cancer_common <- 
  compare_specificity %>% filter(cl_vs_cancer == "High specificity in both") %>% dplyr::select(Ensembl) %>% pull %>% 
  clusterProfiler::enrichGO(OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", universe = compare_specificity$Ensembl,
                            ont = "BP", pAdjustMethod = "BH")

go_specificity$cl_vs_cancer_cl <- 
  compare_specificity %>% filter(cl_vs_cancer == "High specificity in cell line disease") %>% dplyr::select(Ensembl) %>% pull %>% 
  clusterProfiler::enrichGO(OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", universe = compare_specificity$Ensembl,
                            ont = "BP", pAdjustMethod = "BH")

go_specificity$cl_vs_cancer_can <- 
  compare_specificity %>% filter(cl_vs_cancer == "High specificity in TCGA cancer") %>% dplyr::select(Ensembl) %>% pull %>% 
  clusterProfiler::enrichGO(OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", universe = compare_specificity$Ensembl,
                            ont = "BP", pAdjustMethod = "BH")


plot_GO_dot = function(x, title, size_limit, color_limit, xlim, xbreaks){
  gr = do.call(rbind, strsplit(x$GeneRatio, '/'))
  x$GeneRatio = as.numeric(gr[,1])/as.numeric(gr[,2])
  x$Description = factor(x$Description, levels = x$Description[order(x$GeneRatio)])
  x$log_p_value = -log10(x$p.adjust)
  p = ggplot(x, aes(x = GeneRatio, y = Description)) +
    geom_point(aes(size = GeneRatio, fill = log_p_value), color = "grey60", shape = 21) + theme_bw() + #add colors
    labs(size = "Gene Ratio", fill = expression(paste("-Log"[10], " adj. P-value"))) +
    scale_fill_gradient(low="white", high="red", limits = color_limit,
                        guide = guide_colourbar(direction = "horizontal", barwidth = 5, barheight = 0.6)) +
    scale_size(limits = size_limit) + ggtitle(title) +
    scale_x_continuous(limits = xlim, breaks = xbreaks) +
    xlab("Gene Ratio") +
    theme(panel.border = element_rect(colour = "black",size = 1),
          legend.title = element_text(size = 10),
          legend.position = "bottom",
          axis.title.x = element_text(size = 12), axis.title.y = element_blank(),
          aspect.ratio=1, axis.line.x = element_blank(), axis.line.y = element_blank(),
          plot.title = element_text(size = 13, face="bold", hjust = 0.5),
          axis.text.x = element_text(size = 11),
          axis.text.y = element_text(size = 12))
  return(p)
}


ggarrange(plotlist = list(plot_GO_dot(go_specificity$cl_vs_cancer_common@result[1:8,], "High specificity genes in both", 
                                      c(0.0001, 0.1), c(0,53), c(0.009, 0.061), c(0.02, 0.04, 0.06)),
                          plot_GO_dot(go_specificity$cl_vs_cancer_common@result[1:8,], "High specificity genes in both", 
                                      c(0.0001, 0.1), c(0,53), c(0.009, 0.061), c(0.02, 0.04, 0.06)),
                          plot_GO_dot(go_specificity$cl_vs_cancer_cl@result[c(1,5,6,7,11,15,16,19),], 
                                      "High specificity genes in cell line disease", 
                                      c(0.0001, 0.1), c(0,53), c(0.029, 0.091), c(0.03, 0.06, 0.09)),
                          plot_GO_dot(go_specificity$cl_vs_cancer_can@result[1:8,], "High specificity genes in cancer", 
                                      c(0.0001, 0.1), c(0,53), c(0.019, 0.065), c(0.02, 0.04, 0.06))),
          ncol = 2, nrow = 2, legend = "bottom", align = "hv", common.legend = T)
ggsave(file = "results/Figure_2/GO_combined.pdf", width = 12, height = 5.5, bg = "white")



# hierarchical clustering -------------------------------------------------

hierarchical_disease <- list()

hierarchical_disease$exp = exp_landscape$combine_disease %>% 
  pivot_wider(names_from = primary_disease, values_from = ntpm) %>% 
  as.data.frame %>% column_to_rownames("ensg_id")

hierarchical_disease$dendro <- hierarchical_disease$exp %>% 
  cor(method = "spearman") %>% {1 - .} %>% as.dist %>% hclust %>% dendro_data

ggplot() + 
  geom_segment(data = segment(hierarchical_disease$dendro), aes(x = x, y = y, xend = xend, yend = yend)) + 
  geom_point(data = label(hierarchical_disease$dendro), aes(x = x, y = y, color = label), shape = 15, size = 6.2) +
  geom_text(data = label(hierarchical_disease$dendro), aes(x = x, y = y, label = label), size = 4) +
  theme_bw() + scale_color_manual(values = disease_colors_combine) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 10, colour = "black"), axis.title.y = element_blank(), 
        panel.background = element_blank(), axis.text.x = element_text(size = 8, colour = "black"),
        axis.line = element_blank(), axis.ticks.length = unit(.1, "cm"), axis.ticks.y = element_blank(),
        axis.text.y = element_blank(), legend.position = "none") + ylab("1 - Spearman's ρ") +
  coord_flip() + 
  scale_y_reverse(expand = c(0.2, 0))
ggsave("results/Figure_2/hierarchical_tree_disease.pdf", height = 6, width = 2.5, dpi = 600, units = "in")


hierarchical_disease$spec <- exp_landscape$gene_disease %>% dplyr::select(spec_category, enriched_tissues) %>% 
  filter(spec_category %in% c("group enriched", "tissue enhanced", "tissue enriched"))
hierarchical_disease$spec_matrix <- data.frame(matrix(0, length(unique(exp_landscape$combine_disease$primary_disease)), 3), 
                                               row.names = unique(exp_landscape$combine_disease$primary_disease))
colnames(hierarchical_disease$spec_matrix) <- c("group enriched", "tissue enhanced", "tissue enriched")
for(i in 1:nrow(hierarchical_disease$spec)){
  temp <- str_split(hierarchical_disease$spec$enriched_tissues[i], pattern = ";", simplify = T)
  hierarchical_disease$spec_matrix[temp,hierarchical_disease$spec$spec_category[i]] <- hierarchical_disease$spec_matrix[temp,hierarchical_disease$spec$spec_category[i]] + 1
}
hierarchical_disease$spec_matrix_long <- hierarchical_disease$spec_matrix %>% rownames_to_column("disease") %>% 
  pivot_longer(!disease, names_to = "Specificity", values_to = "value") %>% 
  mutate(Specificity = factor(Specificity %>% str_to_sentence %>% str_replace("Tissue", "Disease"), 
                              levels = c("Disease enhanced", "Group enriched", "Disease enriched")),
         disease = factor(disease, levels = hierarchical_disease$dendro$labels$label))

hierarchical_disease$spec_matrix_long %>% 
  ggplot(aes(fill=Specificity, x=value, y=disease)) + 
  geom_bar(stat="identity") + 
  theme_bw() + scale_fill_manual(values = gene_category_pal1) +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 10, colour = "black"), axis.title.y = element_blank(), 
        panel.background = element_blank(), axis.text = element_text(size = 8, colour = "black"), 
        axis.line = element_blank(), axis.ticks.length = unit(.1, "cm"),
        legend.position = "right", 
        legend.title = element_text(size = 11, colour = "black"),
        legend.text = element_text(size = 10, colour = "black")) + 
  scale_x_continuous(expand = c(0,0)) + xlab("Number of genes")
ggsave("results/Figure_2/hierarchical_bar_disease.pdf", height = 6, width = 5, dpi = 600, units = "in")



# hypergeometric testing ----------------------------------------------------

hyper_testing <- list()
hyper_testing$tissue <- exp_landscape$hpa_summary %>% 
  dplyr::select(Ensembl, tissue_specificity = 'RNA tissue specificity', tissue = 'RNA tissue specific nTPM') %>% 
  filter(tissue_specificity == 'Tissue enriched')
hyper_testing$tissue$tissue <- str_split(hyper_testing$tissue$tissue, pattern = ":", simplify = T)[,1]

hyper_testing$cancer <- exp_landscape$hpa_summary %>% 
  dplyr::select(Ensembl, cancer_specificity = 'RNA cancer specificity', cancer = 'RNA cancer specific FPKM') %>% 
  filter(cancer_specificity == 'Cancer enriched')
hyper_testing$cancer$cancer <- str_split(hyper_testing$cancer$cancer, pattern = ":", simplify = T)[,1]

hyper_testing$single_cell_type <- exp_landscape$hpa_summary %>% 
  dplyr::select(Ensembl, sc_specificity = 'RNA single cell type specificity', cell_type = 'RNA single cell type specific nTPM') %>% 
  filter(sc_specificity == 'Cell type enriched')
hyper_testing$single_cell_type$cell_type <- str_split(hyper_testing$single_cell_type$cell_type, pattern = ":", simplify = T)[,1]

hyper_testing$cell_line_disease <- exp_landscape$gene_disease %>% 
  dplyr::select(Ensembl = gene, disease_specificity = spec_category, disease = enriched_tissues) %>% 
  filter(disease_specificity == "tissue enriched")

hyper_testing$cell_line_primary_disease <- exp_landscape$gene_primary_disease %>% 
  dplyr::select(Ensembl = gene, disease_specificity = spec_category, disease = enriched_tissues) %>% 
  filter(disease_specificity == "tissue enriched")

hyper_testing$cell_line_metastasis_disease <- exp_landscape$gene_metastasis_disease %>% 
  dplyr::select(Ensembl = gene, disease_specificity = spec_category, disease = enriched_tissues) %>% 
  filter(disease_specificity == "tissue enriched")


table(hyper_testing$tissue$tissue)
table(hyper_testing$cancer$cancer)
table(hyper_testing$single_cell_type$cell_type)
table(hyper_testing$cell_line_disease$disease)
table(hyper_testing$cell_line_primary_disease$disease)
table(hyper_testing$cell_line_metastasis_disease$disease)


# disease
hyper_testing_results <- list()
hyper_testing_results$cell_line_disease <- hyper_testing_aggr(hyper_testing$cell_line_disease)
hyper_testing_results$cell_line_primary_disease <- hyper_testing_aggr(hyper_testing$cell_line_primary_disease)
hyper_testing_results$cell_line_metastasis_disease <- hyper_testing_aggr(hyper_testing$cell_line_metastasis_disease)

dim(hyper_testing_results$cell_line_disease$combined_adj)
dim(hyper_testing_results$cell_line_primary_disease$combined_adj)
dim(hyper_testing_results$cell_line_metastasis_disease$combined_adj)

split_line <- data.frame(HPA_combined = "Bone Marrow", cell_line_disease = "12split", value = c(10))
hyper_primeta_plot_data <- bind_rows(hyper_testing_results$cell_line_primary_disease$plot_data %>% 
                                       mutate(cell_line_disease = paste("1", cell_line_disease)),
                                     hyper_testing_results$cell_line_metastasis_disease$plot_data %>% 
                                       mutate(cell_line_disease = paste("2", cell_line_disease)),
                                     split_line)
levels_combined <- c(sort(unique(c(rownames(hyper_testing_results$cell_line_primary_disease$combined_adj)[1:12],
                                   rownames(hyper_testing_results$cell_line_metastasis_disease$combined_adj)[1:16]))),
                     sort(unique(c(rownames(hyper_testing_results$cell_line_primary_disease$combined_adj)[13:24],
                                   rownames(hyper_testing_results$cell_line_metastasis_disease$combined_adj)[17:25]))),
                     sort(unique(c(rownames(hyper_testing_results$cell_line_primary_disease$combined_adj)[25:47],
                                   rownames(hyper_testing_results$cell_line_metastasis_disease$combined_adj)[26:44]))))

hyper_primeta_plot_data$HPA_combined <- factor(as.vector(hyper_primeta_plot_data$HPA_combined),
                                               levels = levels_combined %>% str_to_title)


hyper_comparison_plot <- function(plot_data, ylab, width, height, file){
  plot_data %>% 
    ggplot(aes(x = HPA_combined, y = cell_line_disease)) +
    geom_point(aes(size = value, fill = value), shape = 21, color = "grey50") + theme_bw() + #add colors
    labs(size = expression(paste("-Log"[10], "(adj. P-value)")), fill = expression(paste("-Log"[10], "(adj. P-value)"))) +
    scale_fill_gradient(low="white", high="orangered3") +
    guides(size = guide_legend(ncol = 2)) + ylab(ylab) + xlab("HPA combined") +
    theme(panel.border = element_rect(colour = "black",size = 1), 
          axis.text.x = element_text(size = 11, angle = 270, hjust = 0, vjust = 0.5),
          axis.line = element_blank(),
          aspect.ratio = length(unique(plot_data$cell_line_disease))/length(unique(plot_data$HPA_combined)),
          axis.title.y = element_text(size = 12), 
          axis.text.y = element_text(size = 11))
  ggsave(paste0("results/", file), height = height, width = width, dpi = 600, units = "in")
}

hyper_comparison_plot(hyper_testing_results$cell_line_disease$plot_data, "Cell line diseases", 12, 10, "Figure_3/CL_disease_HPA_comparison.pdf")
hyper_comparison_plot(hyper_primeta_plot_data, "Primary cell lines (by disease)", 13.5, 11, "Figure_S5/CL_primary_metastasis_disease_HPA_comparison.pdf")



# Venn Diagram Liver ------------------------------------------------------

table(hyper_testing$tissue$tissue)
table(hyper_testing$cancer$cancer)
table(hyper_testing$single_cell_type$cell_type)
table(hyper_testing$cell_line_disease$disease)
table(hyper_testing$cell_line_primary_disease$disease)
table(hyper_testing$cell_line_metastasis_disease$disease)


myCol <- brewer.pal(4, "Set2")
venn <- venn.diagram(x = list(hyper_testing$cell_line_disease %>% filter(disease == "Liver Cancer") %>% dplyr::select(Ensembl) %>% pull,
                              hyper_testing$single_cell_type %>% filter(cell_type == "Hepatocytes") %>% dplyr::select(Ensembl) %>% pull,
                              hyper_testing$tissue %>% filter(tissue == "liver") %>% dplyr::select(Ensembl) %>% pull,
                              hyper_testing$cancer %>% filter(cancer == "liver cancer") %>% dplyr::select(Ensembl) %>% pull),
                     filename = NULL, lty = 'blank', main.fontface = "bold",
                     # fill=c(myCol[4], "#2142AB", "#36C26B", "#BF0000"),
                     fill=myCol[c(2,3,1,4)],
                     main.fontfamily = "sans", 
                     lwd = 1, fontfamily = "sans", cat.fontfamily = "sans", cex=0.6, cat.cex=0.6,
                     margin = c(0.1,0.1,0.1,0.1), cat.default.pos = "outer",
                     cat.pos = c(-15, 15, -25, 25), cat.dist = c(0.22, 0.22, 0.12, 0.12),
                     category.names = c("Cell line liver cancer", "Hepatocytes", "Liver tissue", "TCGA liver cancer"))
pdf(file="results/Figure_3/liver_overlap.pdf", width = 2, height = 2)
grid.draw(venn)
dev.off()



# Correlation primary/metastatic cell lines ----------------------------------------------------

pm_diff <- list()
pm_diff$annotation <- exp_landscape$combined_annotation_disease %>% filter(!is.na(sample_collection_site))
pm_diff$exp <- exp_landscape$combine %>% spread(sample, ntpm)
pm_diff$primary <- pm_diff$annotation %>% filter(primary_or_metastasis == "Primary")
pm_diff$metastasis <- pm_diff$annotation %>% filter(primary_or_metastasis == "Metastasis")

pm_diff$primary_cor <- c()
pm_diff$meta_disease_cor <- c()
pm_diff$meta_site_cor <- c()

# primary cell line correlation
for(i in unique(pm_diff$primary$primary_disease)){
  if(sum(pm_diff$primary$primary_disease == i) > 1){
    print(i)
    cor_temp = pm_diff$exp %>% 
      dplyr::select(pm_diff$primary %>% filter(primary_disease == i) %>% dplyr::select(RRID) %>% pull) %>% 
      cor(method = "spearman")
    pm_diff$primary_cor <- c(pm_diff$primary_cor, cor_temp[upper.tri(cor_temp, diag = F)])
  }
}

# metastatic cell line correlation (by disease)
for(i in unique(pm_diff$metastasis$primary_disease)){
  if(sum(pm_diff$metastasis$primary_disease == i) > 1){
    print(i)
    cor_temp = pm_diff$exp %>% 
      dplyr::select(pm_diff$metastasis %>% filter(primary_disease == i) %>% dplyr::select(RRID) %>% pull) %>% 
      cor(method = "spearman")
    pm_diff$meta_disease_cor <- c(pm_diff$meta_disease_cor, cor_temp[upper.tri(cor_temp, diag = F)])
  }
}

# metastatic cell line correlation (by collection site)
for(i in unique(pm_diff$metastasis$sample_collection_site)){
  if(sum(pm_diff$metastasis$sample_collection_site == i) > 1){
    print(i)
    cor_temp = pm_diff$exp %>% 
      dplyr::select(pm_diff$metastasis %>% filter(sample_collection_site == i) %>% dplyr::select(RRID) %>% pull) %>% 
      cor(method = "spearman")
    pm_diff$meta_site_cor <- c(pm_diff$meta_site_cor, cor_temp[upper.tri(cor_temp, diag = F)])
  }
}


pm_diff$plot_box <- rbind(data.frame(correlation = pm_diff$primary_cor, type = "Primary"),
                          data.frame(correlation = pm_diff$meta_disease_cor, type = "Metastatic (same disease)"),
                          data.frame(correlation = pm_diff$meta_site_cor, type = "Metastatic (same site)"))
pm_diff$plot_box$type <- factor(pm_diff$plot_box$type, levels = c("Primary", "Metastatic (same disease)", "Metastatic (same site)")[3:1])

pm_diff$plot_box %>% 
  ggplot(aes(x = type, y = correlation, fill = type)) +
  stat_halfeye(adjust = .5, width = .7, .width = 0, justification = -.2, point_colour = NA) +
  geom_boxplot(width = .2, outlier.shape = NA) +
  # geom_jitter(width = .05, alpha = .3) +
  theme_bw() + scale_y_continuous(limits = c(0.73, 1.05))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_text(size = 17, colour = "black"),
        panel.background = element_blank(), 
        axis.text = element_text(size = 13, colour = "black"), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle=330, hjust=0),
        plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm"),
        legend.position = "none",
        legend.title = element_blank()) + ylab("Spearman's p") +
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     comparisons = list(c("Metastatic (same disease)", "Metastatic (same site)"), c("Primary", "Metastatic (same disease)"), 
                                        c("Primary", "Metastatic (same site)")), label.y = c(0.98, 1.0, 1.03)) +
  scale_fill_manual(values = c("Primary" = "#619CFF", "Metastatic (same disease)" = "#00A087", "Metastatic (same site)" = "#F8766D"))
ggsave("results/Figure_3/raincloud_plot.pdf", height = 4.5, width = 4, dpi = 600, units = "in")
dev.off()


# Essential genes comparison ------------------------------------------

essential_genes <- list()
essential_genes$ccle <- 
  exp_landscape$gene_ccle_sample %>% 
  filter(dist_category == "detected in all") %>% 
  dplyr::select(gene) %>% pull
essential_genes$ccle = exp_landscape$hpa_summary$Gene[exp_landscape$hpa_summary$Ensembl %in% essential_genes$ccle]


essential_genes$cell_paper <- xlsx::read.xlsx("data/essential_genes_cell.xlsx", sheetIndex = 1)
essential_genes$cell_paper$BF_hct116_core = ifelse(essential_genes$cell_paper$BF_hct116 > 1.57, T, F)
essential_genes$cell_paper$BF_hela_core = ifelse(essential_genes$cell_paper$BF_hela > 15.47, T, F)
essential_genes$cell_paper$BF_gbm_core = ifelse(essential_genes$cell_paper$BF_gbm > 3.2, T, F)
essential_genes$cell_paper$BF_rpe1_core = ifelse(essential_genes$cell_paper$BF_rpe1 > 6.84, T, F)
essential_genes$cell_paper$BF_dld1_core = ifelse(essential_genes$cell_paper$BF_dld1 > 3.57, T, F)

essential_genes$cell_paper_genes <- 
  essential_genes$cell_paper[apply(essential_genes$cell_paper[,10:14], 1, function(x) return(sum(x == T, na.rm = T))) > 2,] %>% 
  dplyr::select(Gene) %>% pull


essential_genes$science_paper = xlsx::read.xlsx("data/essential_genes_science.xlsx", sheetIndex = 1)
essential_genes$science_paper$KBM7.sig <- ifelse(essential_genes$science_paper$KBM7.CS < -0.1 & essential_genes$science_paper$KBM7.adjusted.p.value < 0.05, T, F)
essential_genes$science_paper$K562.sig <- ifelse(essential_genes$science_paper$K562.CS < -0.1 & essential_genes$science_paper$K562.adjusted.p.value < 0.05, T, F)
essential_genes$science_paper$Jiyoye.sig <- ifelse(essential_genes$science_paper$Jiyoye.CS < -0.1 & essential_genes$science_paper$Jiyoye.adjusted.p.value < 0.05, T, F)
essential_genes$science_paper$Raji.sig <- ifelse(essential_genes$science_paper$Raji.CS < -0.1 & essential_genes$science_paper$Raji.adjusted.p.value < 0.05, T, F)

essential_genes$science_paper_genes <- 
  essential_genes$science_paper[essential_genes$science_paper$KBM7.sig,] %>% 
  dplyr::select(Gene) %>% pull


myCol <- brewer.pal(4, "Set2")
venn <- venn.diagram(x = list(essential_genes$cell_paper_genes,
                              essential_genes$science_paper_genes,
                              essential_genes$ccle),
                     filename = NULL, lty = 'blank', main.fontface = "bold",
                     fill=c(myCol[1:3]), main.fontfamily = "sans", 
                     lwd = 1, fontfamily = "sans", cat.fontfamily = "sans", cex=0.55, cat.cex=0.5,
                     margin = c(0.1,0.1,0.1,0.1), cat.default.pos = "outer",
                     cat.pos = c(-27, 27, 180), cat.dist = c(0.055, 0.055, 0.055),
                     category.names = c("Fitness genes\nHart et al. 2015", "Essential genes\nWang et al. 2015", "CCLE 2019"))
pdf(file="results/Figure_2/essential_gene.pdf", width = 2, height = 2)
grid.draw(venn)
dev.off()


essential_genes$go_ccle_specific <- 
  essential_genes$ccle[!(essential_genes$ccle %in% c(essential_genes$cell_paper_genes, essential_genes$science_paper_genes))] %>% 
  clusterProfiler::enrichGO(OrgDb = org.Hs.eg.db, keyType = "SYMBOL", universe = exp_landscape$hpa_summary$Gene, ont = "BP", pAdjustMethod = "BH")

essential_genes$go_common <- 
  Reduce(intersect, list(essential_genes$ccle, essential_genes$cell_paper_genes, essential_genes$science_paper_genes)) %>% 
  clusterProfiler::enrichGO(OrgDb = org.Hs.eg.db, keyType = "SYMBOL", universe = exp_landscape$hpa_summary$Gene, ont = "BP", pAdjustMethod = "BH")


ggarrange(plotlist = list(plot_GO_dot(essential_genes$go_common@result[1:10,],
                                      "970 common genes", 
                                      c(0.0001, 0.2), c(0,110), c(0.059, 0.19), c(0.06, 0.12, 0.18)),
                          plot_GO_dot(essential_genes$go_ccle_specific@result[c(2:9,14,15),], 
                                      "3,578 CCLE specific genes", 
                                      c(0.0001, 0.2), c(0,110), c(0.018, 0.065), c(0.02, 0.04, 0.06))),
          ncol = 1, nrow = 2, legend = "bottom", align = "hv", common.legend = T)
ggsave(file = "results/Figure_S4/GO_essential.pdf", width = 8, height = 6.5, bg = "white")
dev.off()
