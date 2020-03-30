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
source("src/functions/imputePhosphoData.R")
source("src/functions/mfuzz-ggplot.R")
source("src/functions/simplifyGO.R")
source("src/functions/simplifyGOReqData.R")
source("src/functions/constructHetNet.R")
source("src/functions/calculateRWHN.R")
source("src/functions/dotplot_gg.R")

# Import data
ruprecht_sty <- data.table::fread(input = "data/Ruprecht_STY_2017.csv",
                             select = c(13, 30, 37,
                                        39:46,
                                        93:94,
                                        100
                                        ),
                             skip = 1
) %>% 
  mutate(id = paste0(`Gene names`, "_", `Amino acid` ,`Position`)) %>% 
  filter(!duplicated(id)) %>% 
  mutate(`Gene names` = gsub(";.*", "", `Gene names`)
  )
colnames(ruprecht_sty) <- gsub(" ", "_", colnames(ruprecht_sty))



# Filter data with > 2 NA in experimental condition
par_lap <- ruprecht_sty %>% 
  dplyr::select(id, grep("M/L", colnames(.))) %>%
  filter(`t-test_SignificantM/L_(FDR_<_1%)` == "+") %>% 
  filter_missing(allowed = 1,  colnms = "^M/L")

res_lap <- ruprecht_sty %>% 
  dplyr::select(id, grep("H/L", colnames(.))) %>%
  filter(`t-test_SignificantH/L_(FDR_<_1%)` == "+") %>% 
  filter_missing(allowed = 1,  colnms = "^H/L") %>% 
  mutate_at


# Cluster based on dynamics of phosphorylated sites
par_lap_fcm <- fuzzyC(par_lap, 6)
res_lap_fcm <- fuzzyC(res_lap, 6)

# Construct heterogeneous network
egf_mlnw <- constructHetNet(ruprecht_sty, par_lap, par_lap_fcm)
tgf_mlnw <- constructHetNet(ruprecht_sty, res_lap, res_lap_fcm)

## Run RWHN algorithm
# Recommend to run overnight or on HPC
# EGF:
seed_l <- lapply(1:max(par_lap_fcm$clustering), function(i){
  c(names(par_lap_fcm$clustering[par_lap_fcm$clustering == i]))
})
rwhn_par <- lapply(seed_l, function(s){
  calculateRWHN(edgelists = egf_mlnw$edgelists,
                verti = egf_mlnw$v[egf_mlnw$v$v != "negative regulation of transcription by RNA polymerase II",],
                #verti = egf_mlnw$v,
                seeds = s,
                transitionProb = 0.7,
                restart = 0.7,
                weight_xy = 0.3,
                weight_yz = 0.7) %>%
    filter(name %in% egf_mlnw$v[egf_mlnw$v$layer=="func",]$v)
})
saveRDS(rwhn_par, "results/data/rwhn_egf_clusters.rds")

#TGFa
seed_l <- lapply(1:max(res_lap_fcm$clustering), function(i){
  c(names(res_lap_fcm$clustering[res_lap_fcm$clustering == i]))
})
rwhn_res <- lapply(seed_l, function(s){
  calculateRWHN(edgelists = tgf_mlnw$edgelists,
                verti = tgf_mlnw$v,
                seeds = s,
                transitionProb = 0.7,
                restart = 0.7,
                weight_xy = 0.3,
                weight_yz = 0.7) %>%
    filter(name %in% tgf_mlnw$v[tgf_mlnw$v$layer=="func",]$v)
})
saveRDS(rwhn_res, "results/data/rwhn_tgf_clusters.rds")

# visualise results with dot plot
dot_par <- dotplot_gg(rwhn_egf_cl)
ggsave("results/figs/rwhn_egf_clustersAsSeeds.tif", dot_par)
dot_res <- dotplot_gg(rwhn_tgf_cl)
ggsave("results/figs/rwhn_tgf_clustersAsSeeds.tif", dot_res)

# To determine the frequency of common sites
par_com <- dot_par[[2]] %>% filter(rank_dif == 0) %>% dplyr::select(egf_rank = rank, name) %>% unique()
res_com <- dot_res[[2]] %>% filter(rank_dif == 0) %>% dplyr::select(tgf_rank = rank, name) %>% unique()
com <- merge(par_com, res_com, by = "name", all = T) %>% arrange(egf_rank, tgf_rank)
simpl <- simplifyGOReqData()
freq <- sapply(simpl$GO2Gene,length)
freq_p <- freq / sum(freq) * 100
termid <- filter(simpl$GOterms, TERM %in% com$name)
write.csv(termid, "Table4_commonTermsCFRdata.csv")

