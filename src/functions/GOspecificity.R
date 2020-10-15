GOspecific_vis <- function(ORA_terms, RWHN_terms){
  
  GOratio <- function(t){
    library(GO.db)

    GOancestors <- as.list(GOBPANCESTOR)
    GOoffspring <- as.list(GOBPOFFSPRING)

    tp <- length(GOancestors[[t]]) # 'parents' of query term
    tc <- length(GOoffspring[[t]]) # 'children' of query term

    spec <- tp / tc

    return(spec)
  }

  ORA_ratio <- lapply(1:max(ORA_terms$cluster), function(i){
    tmp <- sapply(unique(ORA_terms[ORA_terms$cluster == i,]$GOID), GOratio)
    if(length(tmp) > 0){
      df <- data.frame(Specificity = tmp) 
      df$Method <- "ORA"
      df$cluster <-  i
    }else{
      df <- data.frame()
    }
    return(df)
  }) %>% 
    do.call(rbind,.)
  RWHN_ratio <- lapply(1:max(RWHN_terms$seed), function(i){
    tmp <- sapply(unique(RWHN_terms[RWHN_terms$seed == i,]$GOID), GOratio)
    if(length(tmp) > 0){
      df <- data.frame(Specificity = tmp) 
      df$Method <- "RWHN"
      df$cluster <-  i
    }else{
      df <- data.frame()
    }
    return(df)
  }) %>% 
    do.call(rbind,.)
  
  gg_ratio <- rbind(ORA_ratio, RWHN_ratio) %>% 
    ggplot(aes(x = Method, y = Specificity, fill = Method)) +
    geom_jitter(color="black", size=0.4, alpha=0.9) +
    geom_boxplot() +
    facet_wrap(~cluster) +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ylab("Ratio Parent Terms:Child Terms")
  
  ##
  GOancestors <- as.list(GOBPANCESTOR)
  GOoffspring <- as.list(GOBPOFFSPRING)
  
  ORA_offspring <- lapply(1:max(ORA_terms$cluster), function(x){
    tmp <- sapply(unique(ORA_terms[ORA_terms$cluster == x,]$GOID), function(i) length(GOoffspring[[i]]))
    if(length(tmp) > 0){
      df <- data.frame(Specificity = tmp) 
      df$Method <- "ORA"
      df$cluster <-  x
    }else{
      df <- data.frame()
    }
    return(df)
  }) %>% do.call(rbind, .)
  
  RWHN_offspring <- lapply(1:max(RWHN_terms$seed), function(x){
    tmp <- sapply(unique(RWHN_terms[RWHN_terms$seed == x,]$GOID), function(i) length(GOoffspring[[i]]))
    if(length(tmp) > 0){
      df <- data.frame(Specificity = tmp) 
      df$Method <- "RWHN"
      df$cluster <-  x
    }else{
      df <- data.frame()
    }
    return(df)
  }) %>% do.call(rbind, .)
  
  gg_offspring <- rbind(ORA_offspring, RWHN_offspring) %>%
    ggplot(aes(x = Method, y = log10(Specificity), fill = Method)) +
    geom_jitter(color="black", size=0.4, alpha=0.9) +
    geom_boxplot() +
    geom_hline(yintercept = 1.399422) + # log10(mean(sapply(GOoffspring, length)))
    facet_wrap(~cluster) +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ylab("Log10 Number of Offspring per term")
  
  ##
  ORA_parentsofRWHN <- lapply(1:max(ORA_terms$cluster), function(x){
    tmp <- sapply(ORA_terms[ORA_terms$cluster == x,]$GOID, function(i){
      tc <- GOoffspring[[i]]
      
      tc_rwhn <- tc[tc %in% RWHN_terms]
      
      return(length(tc_rwhn))
  })
    if(length(tmp) > 0){
      df <- data.frame(ID = names(tmp),
        Specificity = tmp) 
      df$Method <- "ORA"
      df$cluster <-  x
    }else{
      df <- data.frame()
    }
    return(df)
  }) %>% do.call(rbind, .)
  
  RWHN_parentsofORA <- lapply(1:max(RWHN_terms$seed), function(x){
    tmp <- sapply(unique(RWHN_terms[RWHN_terms$seed == x,]$GOID), function(i){
      tc <- GOoffspring[[i]]
      
      tc_rwhn <- tc[tc %in% RWHN_terms[RWHN_terms$seed == x,]$GOID]
      
      return(length(tc_rwhn))
    })
    if(length(tmp) > 0){
      df <- data.frame(ID = names(tmp),
                       Specificity = tmp) 
      df$Method <- "RWHN"
      df$cluster <-  x
    }else{
      df <- data.frame()
    }
    return(df)
  }) %>% do.call(rbind, .)
    

  
  gg_parents  <- rbind(ORA_parentsofRWHN, RWHN_parentsofORA) %>% 
    filter(Specificity > 0) %>% 
    ggplot(aes(x = Method, y = Specificity, fill = ID)) + 
    geom_bar(position="stack", stat="identity") +
    facet_wrap(~cluster) +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ylab("Number of offspring terms found in other method output")
  
  return(list(gg_ratio, 
              gg_offspring,
              gg_parents))
}