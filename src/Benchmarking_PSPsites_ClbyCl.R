# From CRAN
library(dplyr)
library(tidyr)
library(igraph)
library(enrichR)
library(patchwork)
library(ggplot2)
# From Bioconductor
library(GOSemSim)
library(AnnotationDbi)
library(GO.db)
# From src/
source("src/functions/overrepresentationAnalysis.R")
source("src/functions/heatmap_RWHNsig.R")
source("src/functions/simplifyGO.R")
source("src/functions/simplifyGOReqData.R")
source("src/functions/constructHetNet.R")
source("src/functions/calculateRWHN.R")

set.seed(123)

data <- readr::read_tsv("data/Regulatory_sites.txt", skip = 2) %>% 
  filter(!is.na(ON_PROCESS),
         ORGANISM == "human",
         grepl("-p", MOD_RSD)) %>% 
  separate_rows(ON_PROCESS, sep = "; ") %>% 
  mutate(id_site = paste0(GENE,
                          "_",
                          gsub("-.*",
                               "",
                               MOD_RSD))) %>% 
  unique() 

# clustering <- data$ON_PROCESS
# names(clustering) <- data$id_site
# 
# saveRDS(hetNet, "results/data/PSP_benchmark_all_hetNet.rds")
# 
# seed_l <- lapply(unique(clustering)[order(unique(clustering))], function(x){
#   c(names(clustering[clustering == x]))
# })
rwhn <- lapply(unique(data$ON_PROCESS), function(s){
  print(s)
  clustering <- rep(1, nrow(data[data$ON_PROCESS == s,]))
  names(clustering) <- data[data$ON_PROCESS == s,]$id_site
  
  hetNet <- constructHetNet(clustering = clustering)
  
  rwhn <- calculateRWHN(edgelists = hetNet$edgelists,
                verti = hetNet$v,
                seeds = names(clustering),
                transitionProb = 0.7,
                restart = 0.7,
                weight_xy = 0.3,
                weight_yz = 0.7,
                eps = 1/10^12) %>%
    filter(name %in% hetNet$v[hetNet$v$layer=="func",]$v)
})
saveRDS(rwhn, "results/data/rwhn_all_benchmarks.rds")

heatmap_RWHN(rwhn, "GO")
