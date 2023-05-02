pacman::p_load(tidyverse)

# hpa_tcga$cl_priori_cohort$combined_rank$BLCA %>% View

progeny_cytosig_priori <- list()

tcga110cl <- read_tsv("../data_publish/tcga110.tsv")
tcga110cl$Disease[!(tcga110cl$Disease %in% names(hpa_tcga$tcga_disease_match))]
tcga110cl <- 
  tcga110cl %>% 
  filter(Disease != "DLBC", Disease != "MESO")

rank_cor_results <- list()
rank_cor_results$estimate <- c()
rank_cor_results$pvalue <- c()
progeny_cytosig_selected_cl <- list()
tcga_selected_cl <- list()

for(i in names(hpa_tcga$tcga_disease_match)[-24]){
  progeny_cytosig_priori[[i]] <- 
    progeny_comparison$cellline_tcga_dist_results %>% 
    filter(grepl(i, category)) %>% 
    dplyr::select(cell_line_name = ccle, progeny_mse = mse, category) %>% 
    arrange(progeny_mse) %>% 
    mutate(progeny_rank = 1:length(progeny_mse)) %>% 
    left_join(cytosig_comparison$cellline_tcga_dist_results %>% 
                filter(grepl(i, category)) %>% 
                dplyr::select(cell_line_name = ccle, cytosig_mse = mse) %>% 
                arrange(cytosig_mse) %>% 
                mutate(cytosig_rank = 1:length(cytosig_mse))) %>% 
    mutate(overall_rank = progeny_rank + cytosig_rank) %>% 
    arrange(overall_rank) %>% 
    mutate(overall_rank = 1:length(overall_rank)) %>% 
    left_join(hpa_tcga$cell_line_annotation %>% 
                dplyr::select(uniq_id, RRID, cell_line_name, primary_disease, primary_or_metastasis, sample_collection_site, data_source)) %>% 
    relocate(uniq_id, RRID, overall_rank, .before = cell_line_name) %>% 
    dplyr::select(-category)
  
  exp_rank <- progeny_cytosig_priori[[i]] %>% 
    left_join(hpa_tcga$cl_priori_cohort$combined_rank[[i]] %>% 
                dplyr::select(uniq_id, exp_rank = overall_rank))
  rank_cor_results$estimate <- c(rank_cor_results$estimate, cor.test(exp_rank$overall_rank, exp_rank$exp_rank, method = "kendall")$estimate)
  rank_cor_results$pvalue <- c(rank_cor_results$pvalue, cor.test(exp_rank$overall_rank, exp_rank$exp_rank, method = "kendall")$p.value)
  
  progeny_cytosig_selected_cl[[i]] <- exp_rank %>% filter(overall_rank < 6) %>% dplyr::select(cell_line_name) %>% pull
  tcga_selected_cl[[i]] <- exp_rank %>% filter(exp_rank < 6) %>% dplyr::select(cell_line_name) %>% pull
}

if(file.exists("results/cl_prioritization/cell_line_prioritization_progeny_cytosig.xlsx"))
  unlink("results/cl_prioritization/cell_line_prioritization_progeny_cytosig.xlsx")
for(i in names(progeny_cytosig_priori))
  write.xlsx(progeny_cytosig_priori[[i]] %>% relocate(cell_line_name, .after = RRID) %>% as.data.frame,
             file = "results/cl_prioritization/cell_line_prioritization_progeny_cytosig.xlsx", 
             sheetName = i, row.names = F, append = T)


names(rank_cor_results$estimate) <- names(hpa_tcga$tcga_disease_match)[-24]
names(rank_cor_results$pvalue) <- names(hpa_tcga$tcga_disease_match)[-24]

rank_cor_results$estimate
rank_cor_results$pvalue

progeny_cytosig_selected_cl$all <- unlist(progeny_cytosig_selected_cl, use.names = F) %>% unique
tcga_selected_cl$all <- unlist(tcga_selected_cl, use.names = F) %>% unique

length(progeny_cytosig_selected_cl$all)
length(tcga_selected_cl$all)
sum(progeny_cytosig_selected_cl$all %in% tcga_selected_cl$all)
sum(tcga_selected_cl$all %in% tcga110cl$Cell.Line)



hyper_test <- function(list1, list2, all = 985)
  return(phyper(sum(list1 %in% list2)-1, length(list1), all-length(list1), length(list2), lower.tail=F))

hyper_test(tcga_selected_cl$all, tcga110cl$Cell.Line)
hyper_test(progeny_cytosig_selected_cl$all, tcga_selected_cl$all)




