library(dplyr)
library(igraph)

model <- readr::read_csv("data/CFR_STY_2016.csv") %>% 
  filter(`Regulated by EGF` == "+" )
ptmsig <- readr::read_csv("data/data_PTMsigDB_all_sites_v1.9.0.csv")# %>% 
  #filter(grepl("PATH", category))

ptmsig_flt <- apply(model, 1, function(i){
  tmp <- paste0(i["Gene names"], "_", i["Swiss-Prot phosphosite"])
  
  df <- filter(ptmsig, grepl(tmp, site.annotation))
}) %>% do.call(rbind, .) %>% 
  group_split(signature) %>% 
  lapply(function(i){
    sites <- gsub(":.*", "", i$site.annotation)
    expand.grid(sites, sites, 
                stringsAsFactors = F)
  }) %>% 
  do.call(rbind, .) %>% 
  unique()

nw <- graph_from_data_frame(ptmsig_flt)


