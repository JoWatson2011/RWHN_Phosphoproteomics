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
  mutate(`Gene names` = gsub(";.*", "", `Gene names`),
         id = paste0(`Gene names`, "_", `Amino acid` ,`Position`)) %>% 
  filter(!duplicated(id)) %>% 
  rename(`Gene names` = "gene.symbol")
colnames(ruprecht_sty) <- gsub(" ", "_", colnames(ruprecht_sty))


# PARENTAL + LAP M/L: Filter NAs and non significant
par_lap <- ruprecht_sty %>%
  dplyr::select(id, grep("M/L", colnames(.))) %>%
  filter(`t-test_SignificantM/L_(FDR_<_1%)` == "+" ) %>%
  na.omit() %>%
  dplyr::select(-`t-test_SignificantM/L_(FDR_<_1%)`
  ) %>%
  pivot_longer(-id,
               names_to = "rep",
               values_to = "ratio") %>%
  mutate(rep = gsub("_R[1234]$", "", rep)) %>%
  group_by(id, rep) %>%
  summarise(ratio = median(ratio, na_rm = T)) %>%
  pivot_wider(id_cols = id,
              names_from = rep,
              values_from = ratio)


# RESISTANT  + LAP H/L: Filter NAs
par_lap <- ruprecht_sty %>%
  dplyr::select(id, grep("H/L", colnames(.))) %>%
  filter(`t-test_SignificantH/L_(FDR_<_1%)` == "+" ) %>%
  na.omit() %>%
  dplyr::select(-`t-test_SignificantH/L_(FDR_<_1%)`
  ) %>%
  pivot_longer(-id,
               names_to = "rep",
               values_to = "ratio") %>%
  mutate(rep = gsub("_R[1234]$", "", rep)) %>%
  group_by(id, rep) %>%
  summarise(ratio = median(ratio, na_rm = T)) %>%
  pivot_wider(id_cols = id,
              names_from = rep,
              values_from = ratio)

# ALL COND. : Filter data with NA in an experimental condition
tot_lap <- ruprecht_sty %>%
  dplyr::select(id, grep("H[M]/L", colnames(.))) %>%
filter(`t-test_SignificantM/L_(FDR_<_1%)` == "+" |
        `t-test_SignificantH/L_(FDR_<_1%)` == "+") %>%
  filter(`t-test_SignificantM/L_(FDR_<_1%)` == "+" ) %>%
#  filter_missing(allowed = 1,  colnms = "^[M]/L") %>%
  na.omit() %>%
  # dplyr::select(-c(`t-test_SignificantM/L_(FDR_<_1%)`,
  #                  `t-test_SignificantH/L_(FDR_<_1%)`)
  #               ) %>%
  dplyr::select(-`t-test_SignificantM/L_(FDR_<_1%)`
  ) %>%
  pivot_longer(-id,
               names_to = "rep",
               values_to = "ratio") %>%
  mutate(rep = gsub("_R[1234]$", "", rep)) %>%
  group_by(id, rep) %>%
  summarise(ratio = median(ratio, na_rm = T)) %>%
  pivot_wider(id_cols = id,
              names_from = rep,
              values_from = ratio)

# Cluster based on dynamics of phosphorylated sites

# par_lap_fcm <- fuzzyC(par_lap, 3)
# 
# par_lap_hclust <- hclust(dist(par_lap[,2:3]))
# plot(par_lap_hclust)
# par_lap_cl <- cutree(par_lap_hclust, k = 5)

par_lap_cl <- kmeans(par_lap[,2], 4)$cluster
par_lap_cl <- kmeans(par_lap[,2:3], 4)$cluster
names(par_lap_cl) <- par_lap$id

#res_lap_fcm <- fuzzyC(res_lap, 6)

# Construct heterogeneous network
par_mlnw <- constructHetNet(ruprecht_sty, par_lap, par_lap_cl, modules = T)
#tgf_mlnw <- constructHetNet(ruprecht_sty, res_lap, res_lap_fcm)

## Run RWHN algorithm
# Recommend to run overnight or on HPC

par_mlnw$edgelists <- lapply(par_mlnw$edgelists, function(i)
  i %>% filter_all(any_vars(. != "negative regulation of transcription by RNA polymerase II"))
)
par_mlnw$v <- par_mlnw$v[par_mlnw$v$v != "negative regulation of transcription by RNA polymerase II",]

seed_l <- lapply(1:max(par_lap_cl), function(i){
  c(names(par_lap_cl[par_lap_cl == i]))
})

start <- Sys.time()
rwhn_par <- lapply(seed_l, function(s){
  calculateRWHN(edgelists = par_mlnw$edgelists,
                #verti = par_mlnw$v[par_mlnw$v$v != "negative regulation of transcription by RNA polymerase II",],
                verti = par_mlnw$v,
                seeds = s,
                transitionProb = 0.7,
                restart = 0.7,
                weight_xy = 0.3,
                weight_yz = 0.7) %>%
    filter(name %in% par_mlnw$v[par_mlnw$v$layer=="func",]$v)
})
end <- Sys.time()
  saveRDS(rwhn_par, "results/data/rwhn_ruprecht_clusters.rds")



# visualise results with dot plot
dot_par <- dotplot_gg(rwhn_par, n_terms = 27, remove_common = F)
ggsave("results/figs/rwhn_ruprecht_clustersAsSeeds.tiff", dot_par[[1]])

# To determine the frequency of common sites
par_com <- dot_par[[2]] %>% filter(rank_dif == 0) %>% dplyr::select(egf_rank = rank, name) %>% unique()
res_com <- dot_res[[2]] %>% filter(rank_dif == 0) %>% dplyr::select(tgf_rank = rank, name) %>% unique()
com <- merge(par_com, res_com, by = "name", all = T) %>% arrange(egf_rank, tgf_rank)
simpl <- simplifyGOReqData()
freq <- sapply(simpl$GO2Gene,length)
freq_p <- freq / sum(freq) * 100
termid <- filter(simpl$GOterms, TERM %in% com$name)
write.csv(termid, "Table4_commonTermsCFRdata.csv")

