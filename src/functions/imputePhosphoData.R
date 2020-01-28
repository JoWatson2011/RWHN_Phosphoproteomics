imputeTruncNorm <- function(v){         
  dis <- msm::rtnorm(length(v), 
                     lower = -(sd(v, na.rm = T)),
                     upper = (-(sd(v, na.rm = T)) + 1.0)
  )
  v[is.na(v)] <- sample(dis, length(v[is.na(v)]), replace = F)
  return(v)
}

filter_missing <- function(sty, allowed = 0){
  keep <- sty %>% 
    group_by(id) %>% 
    dplyr::select(id, grep("ratio", colnames(.))) %>% 
    pivot_longer(cols = -id, names_to = "experiment", values_to = "value") %>% 
    filter(sum(is.na(value)) <= allowed) %>% 
    dplyr::select(id) %>% 
    unique()
  
  flt <- filter(sty, id %in% keep$id)
  
  return(flt)
}