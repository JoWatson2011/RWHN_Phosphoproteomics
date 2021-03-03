overrepresentationAnalysis <- function(clustering, 
                                       colours,
                                       RWHN_sig = NULL,
                                       simplify = T,
                                       ylab = "GOBP Term", 
                                       database = "GO_Biological_Process_2018"){

  
  enrichedTerms <- lapply(unique(clustering)[order(unique(clustering))], function(i){
    ids <- names(clustering[clustering == i ] )
    cl_prots <- unique(gsub("_.*", "", ids))
    cl_prots <- cl_prots[cl_prots != ""]
    
    enriched <- enrichr(cl_prots, databases = database) %>% 
      .[[1]] %>%
      filter(Adjusted.P.value <= 0.05)  %>% 
      mutate(cluster = i)
    return(enriched)
  }) %>% 
    do.call(rbind, .)
  
  if(grepl("GO_", database)){
    enrichedTerms <-  enrichedTerms %>%
      separate(Term,
               into = c("Term", "GOID"),
               sep = " \\(",
               extra = "drop") %>%
      mutate(GOID = sub("\\)",
                        "",
                        GOID))
  }
  
  
  if(grepl("GO_Biological_Process", database) & simplify){
    simpleGO <- simplifyGOReqData()
    
    enrichedTerms_flt <- lapply(unique(clustering)[order(unique(clustering))], function(i){
      df <- enrichedTerms[enrichedTerms$cluster == i,]
      
      keepID <- simplifyGO(GOID = df$GOID, simplifyData = simpleGO)
      
      df_flt <- df %>% 
        filter(GOID %in% keepID) %>% 
        arrange(desc(Adjusted.P.value)) 
      if(nrow(df_flt) > 0){
        df_flt <- mutate(df_flt, rank = 1:n())
      }
      return(df_flt)
    }) %>% 
      do.call(rbind, .) %>% 
      mutate(V1 = signif(Adjusted.P.value, digits = 2),
             name = factor(Term, unique(Term)))

  } else {
    
    enrichedTerms_flt <- lapply(unique(clustering)[order(unique(clustering))], function(i){
      df <- enrichedTerms[enrichedTerms$cluster ==i, c("Term", "Adjusted.P.value", "cluster")] %>% 
        arrange(desc(Adjusted.P.value))
      
      if(nrow(df) > 0){
        df <- mutate(df, rank = 1:n())
      }
      
      return(df)
    }) %>% 
      do.call(rbind, .) %>% 
      mutate(V1 = signif(Adjusted.P.value, digits = 2),
             name = factor(Term, unique(Term)))
    
  }
  
  if(!is.null(RWHN_sig)){
    enrichedTerms_flt$Term <- tolower(enrichedTerms_flt$Term)
    enrichedTerms_flt <- merge(enrichedTerms_flt,RWHN_sig$data[,c("name", "seed")],
                               by.x = "Term", by.y = "name", all.x = T)
    enrichedTerms_flt$rwhn <- ifelse(enrichedTerms_flt$cluster == enrichedTerms_flt$seed, T, NA)
    enrichedTerms_flt <- dplyr::select(enrichedTerms_flt, -seed) %>%  unique()
  }
  
  
  gg <- ggplot(enrichedTerms_flt, aes(y = name, x = as.factor(cluster))) +
    geom_tile(aes(fill = Adjusted.P.value), color = "white") +
    scale_fill_gradient(breaks = seq(0, 0.05, 0.01), limits = c(0, 0.05),
                        low = colours[1], high = colours[2]) +
    guides(color = FALSE,
           fill = guide_colourbar(title="FDR", reverse = T)) +
    theme(legend.key.size = unit(.25, "cm"),
          legend.title = element_text(size = 5),
          legend.text = element_text(size = 5), 
          axis.text.y = element_text(size = 5),
          axis.title = element_text(size = 6),
          panel.background = element_rect(fill = "black"), 
          panel.grid = element_blank(),
          legend.margin = margin(0,0,0,0, "cm")
    ) +
    scale_x_discrete(position = "top") +
    xlab("Cluster") +
    ylab(ylab)
  
  if(!is.null(RWHN_sig)){
    gg <- gg + 
      geom_point(aes(shape = rwhn), show.legend = F) 
  }
  
  return(gg)
}
