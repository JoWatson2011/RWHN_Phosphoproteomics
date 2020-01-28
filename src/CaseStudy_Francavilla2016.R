# From CRAN
library(dplyr)
library(tidyr)
library(igraph)
library(e1071)
library(enrichR)
# From Bioconductor
library(Mfuzz)
library(GOSemSim)
library(AnnotationDbi)
library(GO.db)
# From src/
source("src/imputePhosphoData.R")
source("src/mfuzz-ggplot.R")
source("src/simplifyGO.R")
source("src/simplifyGOReqData.R")
source("src/constructHetNet.R")
source("src/calculateRWHN.R")
source("src/dotplot_gg.R")

# Import data
cfr_sty <- data.table::fread(input = "data/CFR_STY_2016.csv",
                             select = c(1:9, 26:33)
) %>% 
  mutate(id = paste0(`Gene names`, "_", `Swiss-Prot phosphosite`)) %>% 
  filter(!duplicated(id))
colnames(cfr_sty) <- gsub(" ", "_", colnames(cfr_sty))



# Filter data with > 2 missing values and
# impute phosphorylated sites regulated by EGF or TGF-a
egf <- cfr_sty %>% 
  dplyr::select(id, grep("EGF", colnames(cfr_sty))) %>% 
  filter_missing(allowed = 2) %>% 
  mutate_at(vars(matches("ratio")), imputeTruncNorm)

tgf <- cfr_sty %>% 
  dplyr::select(id, grep("TGF", colnames(cfr_sty))) %>% 
  filter_missing(allowed = 2) %>% 
  mutate_at(vars(matches("ratio")), imputeTruncNorm)

# Cluster based on dynamics of phosphorylated sites
egf_fcm <- fuzzyC(egf, 6)
tgf_fcm <- fuzzyC(tgf, 6)

# Construct heterogeneous network
egf_mlnw <- constructHetNet(cfr_sty, egf, egf_fcm)
tgf_mlnw <- constructHetNet(cfr_sty, tgf, tgf_fcm)

## Run RWHN algorithm
# Recommend to run overnight or on HPC
# EGF:
seed_l <- lapply(1:max(egf_fcm$clustering), function(i){
  c(names(egf_fcm$clustering[egf_fcm$clustering == i]))
})
rwhn_egf <- lapply(seed_l, function(s){
  calculateRWHN(edgelists = edgelists_egf,
                #verti = v_egf[v_egf$v != "negative regulation of transcription by RNA polymerase II",],
                verti = v_egf,
                seeds = s,
                transitionProb = 0.7,
                restart = 0.7,
                weight_xy = 0.3,
                weight_yz = 0.7) %>%
    filter(name %in% v_egf[v_egf$layer=="func",]$v)
})
saveRDS(rwhn_egf, "results/data/rwhn_egf_clusters.rds")

#TGFa
seed_l <- lapply(1:max(tgf_fcm$clustering), function(i){
  c(names(tgf_fcm$clustering[tgf_fcm$clustering == i]))
})
rwhn_tgf <- lapply(seed_l, function(s){
  calculateRWHN(edgelists = edgelists_tgf,
                verti = v_tgf,
                seeds = s,
                transitionProb = 0.7,
                restart = 0.7,
                weight_xy = 0.3,
                weight_yz = 0.7) %>%
    filter(name %in% v_tgf[v_tgf$layer=="func",]$v)
})
saveRDS(rwhn_tgf, "results/data/rwhn_tgf_clusters.rds")

# visualise results with dot plot
dot_egf <- dotplot_gg(rwhn_egf_cl)
ggsave("results/figs/rwhn_egf_clustersAsSeeds.tif", dot_egf)
dot_tgf <- dotplot_gg(rwhn_tgf_cl)
ggsave("results/figs/rwhn_tgf_clustersAsSeeds.tif", dot_tgf)

# To determine the frequency of common sites
egf_com <- dot_egf[[2]] %>% filter(rank_dif == 0) %>% dplyr::select(egf_rank = rank, name) %>% unique()
tgf_com <- dot_tgf[[2]] %>% filter(rank_dif == 0) %>% dplyr::select(tgf_rank = rank, name) %>% unique()
com <- merge(egf_com, tgf_com, by = "name", all = T) %>% arrange(egf_rank, tgf_rank)
simpl <- simplifyGOReqData()
freq <- sapply(simpl$GO2Gene,length)
freq_p <- freq / sum(freq) * 100
termid <- filter(simpl$GOterms, TERM %in% com$name)
write.csv(termid, "Table4_commonTermsCFRdata.csv")

