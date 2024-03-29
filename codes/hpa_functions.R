

remove_genes <- 
  function(wide_dat,
           rm_NA_genes = F,
           rm_not_expressed = F,
           not_expressed_lim = 1,
           rm_no_variance = F) {
    suppressMessages(require(tidyverse))
    
    
    if(rm_NA_genes) {
      
      wide_dat <- 
        wide_dat %>% 
        {.[complete.cases(.), ]}
    }
    
    if(rm_not_expressed) {
      
      wide_dat <- 
        wide_dat[apply(wide_dat,
                       MARGIN = 1, 
                       max) >= not_expressed_lim,] 
    }
    
    if(rm_not_expressed) {
      
      wide_dat <- 
        wide_dat[apply(wide_dat,
                       MARGIN = 1, 
                       max) >= not_expressed_lim,] 
    }
    
    if(rm_no_variance) {
      
      wide_dat <- 
        wide_dat[apply(wide_dat,
                       MARGIN = 1, 
                       sd) > 0,] 
    }
    
    wide_dat
  }

scale_data <- 
  function(wide_data, 
           logp1_scale = F,
           zscore_scale = F, 
           max_scale = F) {
    suppressMessages(require(tidyverse))
    
    if(zscore_scale & max_scale) stop("Scaling: Choose either max scaling or z-score scaling")
    
    if(logp1_scale) {
      wide_data <- 
        log10(wide_data + 1)
    }
    if(zscore_scale) {
      wide_data <- 
        wide_data %>% 
        t() %>% 
        scale() %>% 
        t()
    }
    
    if(max_scale) {
      gene_max <- 
        apply(wide_data, 
            MARGIN = 1,
            max)
      
      wide_data <-
        sweep(wide_data, MARGIN=1, gene_max, `/`)
    }
    
    wide_data
  }

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



extract_variable_features <- 
  function(wide_data, 
           log1p_data = F,
           maxvars = 5000,
           ...) {
    
    # Whether to log scale data for variable extraction or not. 
    # This does not affect output data
    if(log1p_data) {
      wide_data <- 
        log1p(wide_data)
    }
    
    # Calculate mean and standar deviation
    vardata <- 
      wide_data %>% 
      {tibble(var = rownames(wide_data), 
              mean = apply(., 
                           MARGIN = 1, 
                           mean),
              sd = apply(., 
                         MARGIN = 1, 
                         sd))}
    
    # Make loess model
    vardata_loess <- 
      loess(sqrt(sd) ~ mean, data = vardata,
            ...)
    
    vardata_res <- 
      vardata %>%
      mutate(pred_sd = predict(vardata_loess, newdata = .),
             high = sqrt(sd) > pred_sd,
             diff = sqrt(sd) - pred_sd,
             top = rank(-diff) < maxvars) 
    
    # vardata_res %T>%
    #   {group_by(., high) %>%
    #       count %>%
    #       print} %>%
    #   mutate(diff = sqrt(sd) - pred_sd,
    #          top = rank(-diff) < maxvars) %>% 
    #   ggplot(aes(mean, sqrt(sd), color  = top)) +
    #   geom_point() +
    #   stat_function(fun = function(x) predict(vardata_loess, newdata = tibble(mean = x)),
    #                 color = "black")
    
    highvar_vars <- 
      vardata_res %>% 
      filter(top) %>% 
      pull(var)
    
    wide_data[highvar_vars,]
    
  }





HPA_sample_UMAP <- 
  function(wide_data) {
    wide_data %>% 
      # Remove genes that have max expression < 1 in dataset, 
      # and genes without variance in expression:
      remove_genes(rm_NA_genes = T,
                   rm_not_expressed = T, 
                   not_expressed_lim = 1,
                   rm_no_variance = T) %>% 
      
      # Extract highly variable variables: 
      # extract_variable_features(log1p_data = T) %>%
      
      # Scale data genewise, first log10(exp + 1), 
      # and then z-score:
      scale_data(logp1_scale = T, 
                 zscore_scale = T) %>% 
      
      # Perform PCA:
      do_pca() %>% 
      
      # Get PCA scores, selecting PCs that constitute 80% cumulative r2:
      get_pca_scores(use_R2cum_PCselection = T,
                     R2cum_lim = 0.8) %>% 
      
      # Perform UMAP with default settings:
      do_umap(n_neighbors = round(sqrt(dim(.)[1]))) %>% 
      
      # Add rownames as a column for saving:
      as_tibble(rownames = "sample")
  }



hpa_gene_classification <- 
  function(data, expression_col, disease_col, gene_col, enr_fold, max_group_n, det_lim = 1) {
    data_ <- 
      data %>% 
      dplyr::select(gene = all_of(gene_col),
                    expression = all_of(expression_col),
                    disease = all_of(disease_col)) %>% 
      mutate(expression = round(expression, 4)) 
    
    if(any(is.na(data_$expression))) stop("NAs in expression column")
    if(any(is.na(data_$gene))) stop("NAs in gene column")
    if(any(is.na(data_$disease))) stop("NAs in disease column")
    
    n_groups <- length(unique(data_$disease))
    
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
        enrichment_group = paste(sort(disease[which(expression >= lim & expression >= det_lim)]), collapse=";"),
        
        n_enriched = length(disease[which(expression >= lim & expression >= det_lim)]),
        n_enhanced = length(exps_enhanced[[1]]), 
        enhanced_in = paste(sort(disease[exps_enhanced[[1]]]), collapse=";"),
        n_na = n_groups - length(expression),
        max_2nd_or_lim = max(max_2nd, det_lim*0.1),
        diseases_not_detected = paste(sort(disease[which(expression < det_lim)]), collapse=";"),
        diseases_detected = paste(sort(disease[which(expression >= det_lim)]), collapse=";")) 
    
    
    gene_categories <- 
      gene_class_info %>%
      
      mutate(
        spec_category = case_when(n_exp == 0 ~ "not detected", 
                                  
                                  # Genes with expression fold times more than anything else are enriched
                                  max_exp/max_2nd_or_lim >= enr_fold ~ "enriched", 
                                  
                                  # Genes with expression fold times more than other diseases in groups of max group_n - 1 are group enriched
                                  max_exp >= lim &
                                    n_over <= max_group_n & n_over > 1 &
                                    mean_over/max_under_lim >= enr_fold ~ "group enriched", 
                                  
                                  # Genes with expression in diseases fold times more than the mean are enhance
                                  n_enhanced > 0 ~ "enhanced", 
                                  
                                  # Genes expressed with low specificity
                                  T ~ "low specificity"), 
        
        
        dist_category = case_when(frac_exp == 100 ~ "detected in all",
                                  frac_exp >= 33 ~ "detected in many",
                                  n_exp > 1 ~ "detected in some",
                                  n_exp == 1 ~ "detected in single",
                                  n_exp == 0 ~ "not detected"),
        
        spec_score = case_when(spec_category == "enriched" ~ max_exp/max_2nd_or_lim,
                               spec_category == "group enriched" ~ mean_over/max_under_lim, 
                               spec_category == "enhanced" ~ max_exp/mean_exp)) 
    
    
    
    
    ##### Rename and format
    gene_categories %>%
      mutate(enriched_diseases = case_when(spec_category %in% c("enriched", "group enriched") ~ enrichment_group,
                                          spec_category == "enhanced" ~ enhanced_in),
             n_enriched = case_when(spec_category %in% c("enriched", "group enriched") ~ n_enriched,
                                    spec_category == "enhanced" ~ n_enhanced)) %>%
      dplyr::select(gene, 
                    spec_category, 
                    dist_category, 
                    spec_score,
                    n_expressed = n_exp, 
                    fraction_expressed = frac_exp,
                    max_exp = max_exp,
                    enriched_diseases,
                    n_enriched,
                    n_na = n_na,
                    diseases_not_detected,
                    diseases_detected) 
    
    
    
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
  result_list$tcga_cancer <- hyper_testing_comparison(hyper_testing$tcga_cancer, list1) %>% arrange(rownames(.))
  result_list$tissue <- hyper_testing_comparison(hyper_testing$tissue, list1) %>% arrange(rownames(.))
  result_list$single_cell_type <- hyper_testing_comparison(hyper_testing$single_cell_type, list1) %>% arrange(rownames(.))
  
  result_list$combined <- rbind(result_list$tcga_cancer,
                                result_list$tissue,
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
    pivot_longer(cols = -HPA_combined, names_to = "cell_line_disease") %>% 
    filter(value != 0) %>% 
    mutate(HPA_combined = factor(HPA_combined, levels = unique(HPA_combined)))
  
  return(result_list)
}


plot_GO_dot = function(x, title, size_limit, color_limit, xlim, xbreaks){
  gr <- do.call(rbind, strsplit(x$GeneRatio, '/'))
  x$GeneRatio <- as.numeric(gr[,1])/as.numeric(gr[,2])
  x$Description <- factor(x$Description, levels = x$Description[order(x$GeneRatio)])
  x$log_p_value <- -log10(x$p.adjust)
  p <- ggplot(x, aes(x = GeneRatio, y = Description)) +
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


cell_line_gsea <- 
  function(cl_exp_gsea, disease_signature){
    gene_rank <- sort(cl_exp_gsea)
    results <- fgsea::fgsea(pathways = disease_signature, 
                                  stats = gene_rank, 
                                  eps = 1e-300, 
                                  minSize = 1, 
                                  maxSize = 10000)
}


get_disease_signature <- 
  function(temp_tcga){
    disease_signature <- list()
    for(i in colnames(temp_tcga)){
      temp_max <- 
        temp_tcga %>% 
        rownames_to_column("ensg_id") %>% 
        dplyr::select(-all_of(i)) %>% 
        pivot_longer(!ensg_id, names_to = "stage", values_to = "ntpm") %>% 
        group_by(ensg_id) %>% 
        summarise(max_ntpm = 4*max(ntpm)) %>% 
        suppressMessages
      
      disease_signature[[i]] <- 
        temp_tcga %>% 
        rownames_to_column("ensg_id") %>% 
        as_tibble() %>% 
        dplyr::select(ensg_id, target = all_of(i)) %>% 
        left_join(temp_max) %>% 
        filter(target > max_ntpm) %>% 
        dplyr::select(ensg_id) %>% 
        pull %>% 
        suppressMessages
    }
    
    return(disease_signature)
  }


get_gsea_exp <- 
  function(cell_line_exp, all_exp){
    cl_exp_gsea <- 
      cell_line_exp
    
    sd_genes <- apply(cl_exp_gsea, MARGIN = 1, sd)
    sum(sd_genes == 0)
    cl_exp_gsea <- 
      cl_exp_gsea[sd_genes != 0,]
    
    zero_exp <- apply(cl_exp_gsea, MARGIN = 1, function(x) sum(x==0))
    sum(zero_exp == ncol(cl_exp_gsea))
    
    cl_exp_gsea <- 
      cl_exp_gsea[zero_exp != ncol(cl_exp_gsea),]
    
    cl_mean_exp <- 
      all_exp %>% 
      filter(sample %in% colnames(cl_exp_gsea)) %>% 
      left_join(hpa_tcga$cell_line_annotation, 
                by = c("sample" = "RRID")) %>% 
      ungroup %>% 
      dplyr::select(ensg_id, ntpm, primary_disease) %>% 
      group_by(ensg_id, primary_disease) %>% 
      summarise(ntpm = mean(ntpm, na.rm = T)) %>% 
      dplyr::select(-primary_disease) %>%
      group_by(ensg_id) %>% 
      summarise(ntpm = mean(ntpm, na.rm = T)) %>% 
      filter(ensg_id %in% rownames(cl_exp_gsea)) %>% 
      suppressMessages()
    
    identical(rownames(cl_exp_gsea), cl_mean_exp$ensg_id)
    sum(cl_mean_exp$ntpm == 0)
    
    cl_exp_gsea <-
      log2(cl_exp_gsea/cl_mean_exp$ntpm + 1) %>% t
  }


