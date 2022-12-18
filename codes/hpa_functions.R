

do_pca <- 
  function(wide_data, npcs = NULL) {
    suppressMessages(require(tidyverse))
    suppressMessages(require(pcaMethods))
    
    if(is.null(npcs)) {
      npcs <- min(dim(wide_data))
    }
    
    wide_data %>% 
      t() %>% 
      pca(nPcs = npcs) 
  }

get_pca_scores <- 
  function(pca_res, 
           use_R2cum_PCselection = F,
           R2cum_lim = 0.8,
           use_sDev_PCselection = F) {
    suppressMessages(require(tidyverse))
    
    pc_lim <- ncol(pca_res@scores)
    
    if(use_R2cum_PCselection) {
      pc_lim <- 
        which(pca_res@R2cum > R2cum_lim)[1]
    } else if(use_sDev_PCselection) {
      pc_lim <- 
        rev(which(pca_res@sDev >= 1))[1]
    }
    
    pca_res %>% 
      scores() %>% 
      {.[,1:pc_lim]} 
  }

do_umap <- 
  function(wide_data, 
           seed = 42, 
           n_neighbors = 15,
           n_components = 2, 
           ...) {
    suppressMessages(require(tidyverse))
    suppressMessages(require(uwot))
    
    set.seed(seed)
    
    umap_res <- 
      umap(wide_data, 
           n_neighbors = n_neighbors,
           n_components = n_components,
           ...) %>% 
      as.data.frame()
      
    rownames(umap_res) <- rownames(wide_data)
    colnames(umap_res) <- paste0("UMAP", 1:ncol(umap_res))
    
    umap_res
  }


hpa_gene_classification <- 
  function(data, expression_col, tissue_col, gene_col, enr_fold, max_group_n, det_lim = 1) {
    data_ <- 
      data %>% 
      dplyr::select(gene = gene_col,
                    expression = expression_col,
                    tissue = tissue_col) %>% 
      mutate(expression = round(expression, 4)) 
    
    if(any(is.na(data_$expression))) stop("NAs in expression column")
    if(any(is.na(data_$gene))) stop("NAs in gene column")
    if(any(is.na(data_$tissue))) stop("NAs in tissue column")
    
    n_groups <- length(unique(data_$tissue))
    
    gene_class_info <- 
      data_ %>%
      group_by(gene) %>%
      summarise(
        
        # Gene expression distribution metrics
        mean_exp = mean(expression, na.rm = T),
        min_exp = min(expression, na.rm = T),
        max_exp = max(expression, na.rm = T), 
        max_2nd = sort(expression)[length(expression)-1],
        
        # Expression frequency metrics
        n_exp = length(which(expression >= det_lim)),
        frac_exp = n_exp/length(expression[!is.na(expression)])*100,
        
        # Limit of enhancement metrics
        lim = max_exp/enr_fold, 
        
        exps_over_lim = list(expression[which(expression >= lim & expression >= det_lim)]),
        n_over = length(exps_over_lim[[1]]), 
        mean_over = mean(exps_over_lim[[1]]),
        min_over = ifelse(n_over == 0, NA,
                          min(exps_over_lim[[1]])),
        
        max_under_lim = max(expression[which(expression < min_over)], det_lim*0.1),
        
        
        exps_enhanced = list(which(expression/mean_exp >= enr_fold & expression >= det_lim)),
        
        
        
        
        # Expression patterns
        enrichment_group = paste(sort(tissue[which(expression >= lim & expression >= det_lim)]), collapse=";"),
        
        n_enriched = length(tissue[which(expression >= lim & expression >= det_lim)]),
        n_enhanced = length(exps_enhanced[[1]]), 
        enhanced_in = paste(sort(tissue[exps_enhanced[[1]]]), collapse=";"),
        n_na = n_groups - length(expression),
        max_2nd_or_lim = max(max_2nd, det_lim*0.1),
        tissues_not_detected = paste(sort(tissue[which(expression < det_lim)]), collapse=";"),
        tissues_detected = paste(sort(tissue[which(expression >= det_lim)]), collapse=";")) 
    
    
    gene_categories <- 
      gene_class_info %>%
      
      mutate(
        spec_category = case_when(n_exp == 0 ~ "not detected", 
                                  
                                  # Genes with expression fold times more than anything else are tissue enriched
                                  max_exp/max_2nd_or_lim >= enr_fold ~ "tissue enriched", 
                                  
                                  # Genes with expression fold times more than other tissues in groups of max group_n - 1 are group enriched
                                  max_exp >= lim &
                                    n_over <= max_group_n & n_over > 1 &
                                    mean_over/max_under_lim >= enr_fold ~ "group enriched", 
                                  
                                  # Genes with expression in tissues fold times more than the mean are tissue enhance
                                  n_enhanced > 0 ~ "tissue enhanced", 
                                  
                                  # Genes expressed with low tissue specificity
                                  T ~ "low tissue specificity"), 
        
        
        dist_category = case_when(frac_exp == 100 ~ "detected in all",
                                  frac_exp >= 33 ~ "detected in many",
                                  n_exp > 1 ~ "detected in some",
                                  n_exp == 1 ~ "detected in single",
                                  n_exp == 0 ~ "not detected"),
        
        spec_score = case_when(spec_category == "tissue enriched" ~ max_exp/max_2nd_or_lim,
                               spec_category == "group enriched" ~ mean_over/max_under_lim, 
                               spec_category == "tissue enhanced" ~ max_exp/mean_exp)) 
    
    
    
    
    ##### Rename and format
    gene_categories %>%
      mutate(enriched_tissues = case_when(spec_category %in% c("tissue enriched", "group enriched") ~ enrichment_group,
                                          spec_category == "tissue enhanced" ~ enhanced_in),
             n_enriched = case_when(spec_category %in% c("tissue enriched", "group enriched") ~ n_enriched,
                                    spec_category == "tissue enhanced" ~ n_enhanced)) %>%
      dplyr::select(gene, 
                    spec_category, 
                    dist_category, 
                    spec_score,
                    n_expressed = n_exp, 
                    fraction_expressed = frac_exp,
                    max_exp = max_exp,
                    enriched_tissues,
                    n_enriched,
                    n_na = n_na,
                    tissues_not_detected,
                    tissues_detected) 
    
    
    
  }


calculate_tau <- function(x){
  genes <- x %>%
    dplyr::select(gene) %>%
    pull %>%
    unique

  tau <- c()
  for(i in genes){
    exps <- x %>%
      ungroup %>%
      filter(gene == i) %>%
      dplyr::select(exp) %>%
      pull
    exps <- exps/max(exps)
    exps <- 1-exps
    tau <- c(tau, sum(exps)/(length(exps)-1))
  }
  names(tau) <- genes
  tau <- tau[!is.na(tau)]
  return(tau)
}


hyper_test <- function(list_1, list_2, all = 20090)
  return(phyper(sum(list_1 %in% list_2) - 1, length(list_1), all - length(list_1), length(list_2), lower.tail = F))



hyper_testing_comparison <- function(list1, list2){
  list1_unique <- list1[,3] %>% pull %>% unique
  list2_unique <- list2[,3] %>% pull %>% unique
  result_matrix <- data.frame(matrix(0, length(list1_unique), length(list2_unique)))
  rownames(result_matrix) <- list1_unique
  colnames(result_matrix) <- list2_unique
  for(i in list1_unique)
    for(j in list2_unique)
      result_matrix[i,j] <- hyper_test(list1$Ensembl[pull(list1[,3]) == i], list2$Ensembl[pull(list2[,3]) == j])
  return(result_matrix)
}



hyper_testing_aggr <- function(list1){
  result_list <- list()
  result_list$tissue <- hyper_testing_comparison(hyper_testing$tissue, list1) %>% arrange(rownames(.))
  result_list$cancer <- hyper_testing_comparison(hyper_testing$cancer, list1) %>% arrange(rownames(.))
  result_list$single_cell_type <- hyper_testing_comparison(hyper_testing$single_cell_type, list1) %>% arrange(rownames(.))
  
  result_list$combined <- rbind(result_list$tissue,
                                result_list$cancer,
                                result_list$single_cell_type)
  result_list$combined_adj <- p.adjust(as.matrix(result_list$combined), method = "BH")
  dim(result_list$combined_adj) <- dim(result_list$combined)
  colnames(result_list$combined_adj) <- colnames(result_list$combined)
  rownames(result_list$combined_adj) <- rownames(result_list$combined)
  
  filter_rows <- apply(result_list$combined_adj, 1, function(x) all(x>0.05))
  result_list$combined_adj <- result_list$combined_adj[!filter_rows,]
  filter_cols <- apply(result_list$combined_adj, 2, function(x) all(x>0.05))
  result_list$combined_adj <- result_list$combined_adj[,!filter_cols]
  
  result_list$combined_adj <- -log10(result_list$combined_adj)
  result_list$combined_adj[result_list$combined_adj< -log10(0.05)] = 0
  
  result_list$plot_data <- result_list$combined_adj %>% 
    as.data.frame() %>% rownames_to_column("HPA_combined") %>% 
    pivot_longer(cols = -HPA_combined, names_to = "cell_line_disease") %>% filter(value != 0)
  
  result_list$plot_data$HPA_combined <- 
    factor(str_to_title(result_list$plot_data$HPA_combined), 
           levels = unique(str_to_title(result_list$plot_data$HPA_combined)))
  return(result_list)
}


