dotplot_gg <- function(rwhn, seed_names = NULL, n_terms = 20, remove_common = F){
  ggdf <- lapply(1:length(rwhn), function(i){
    rwhn[[i]] %>% 
      mutate(seed = ifelse(is.null(seed_names),
                           i,
                           seed_names[i])) %>%  
      mutate(rank = 1:nrow(.))
  }) %>% do.call(rbind, .)
  
  if(remove_common){
    ggdf <- ggdf %>% 
      group_by(name) %>% 
      mutate(rank_dif = (rank - mean(rank))) %>% 
      filter(rank_dif > 0) %>% 
      group_by(seed) %>%
      mutate(rank_flt = 1:n()) %>% 
      filter(rank_flt <= n_terms) %>% 
      ungroup() %>% 
      mutate(V1 = signif(V1, digits = 2))
  }else{
    ggdf <- ggdf %>% 
      group_by(name) %>% 
      mutate(rank_dif = (rank - mean(rank))) %>% 
      group_by(seed) %>%
      mutate(rank_flt = 1:n()) %>% 
      filter(rank_flt <= n_terms) %>% 
      ungroup() %>% 
      mutate(V1 = signif(V1, digits = 2)) %>% 
      mutate(name = factor(name, unique(name)))
  }
  
  dot <- ggplot(ggdf, aes(y = as.factor(seed), x = name)) +
    geom_count(aes(color = rank_flt)) +
    theme_bw() +
    scale_y_discrete(name = "Cluster") +
    scale_color_distiller(name = "RWHN Rank", palette = "OrRd") +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          plot.margin = margin(10,10,10,160)
    ) +
    guides(size = F)
  
  return(list(dot, ggdf))
}
