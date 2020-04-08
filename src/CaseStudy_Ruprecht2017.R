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
  filter(!duplicated(id)) 
colnames(ruprecht_sty) <- gsub(" ", ".", colnames(ruprecht_sty))
colnames(ruprecht_sty)[grep("FDR.<.1%",
                            colnames(ruprecht_sty))] <- substr(grep("FDR.<.1%", colnames(ruprecht_sty),
                                                                      value = T),
                                                                 8, 
                                                                 21)
colnames(ruprecht_sty) <- gsub("/", "", colnames(ruprecht_sty))


# PARENTAL + LAP M/L: Filter NAs and non significant
par_lap <- ruprecht_sty %>%
  dplyr::select(id, grep("ML", colnames(.))) %>%
  filter(SignificantML == "+" ) %>%
  na.omit() %>%
  dplyr::select(-SignificantML
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
res_lap <- ruprecht_sty %>%
  dplyr::select(id, grep("HL", colnames(.))) %>%
  filter(SignificantHL == "+" ) %>%
  na.omit() %>%
  dplyr::select(-SignificantHL
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
  dplyr::select(id, grep("[HM]L", colnames(.))) %>%
  filter(SignificantML == "+" |
           SignificantHL == "+") %>%
  filter_missing(allowed = 1,  colnms = "^[HM]L_") %>%
  mutate_if(is.numeric, imputeTruncNorm) %>% 
  # dplyr::select(-c(SignificantML,
  #                  SignificantHL)
  # ) %>%
  pivot_longer(-c(id, SignificantML, SignificantHL),
               names_to = "rep",
               values_to = "ratio") %>%
  mutate(rep = gsub("_R[1234]$", "", rep)) %>%
  group_by(id, rep) %>%
  mutate(ratio = median(ratio, na_rm = T)) %>%
  unique() %>% 
  pivot_wider(#id_cols = id,
              names_from = rep,
              values_from = ratio)

# Cluster based on dynamics of phosphorylated sites
wss <- data.frame(type = c(rep("All significant", 15),
                           rep("Just * parental", 15), 
                           rep("Just * resistant", 15)),
                  cl = rep(1:15, 3),
                  wss = c(sapply(1:15, 
                                 function(k){kmeans(tot_lap[,4:5], 
                                                    k, nstart=50,
                                                    iter.max = 15 )$tot.withinss}),
                          sapply(1:15, 
                                 function(k){kmeans(par_lap[,2], 
                                                    k, nstart=50,
                                                    iter.max = 15 )$tot.withinss}),
                          sapply(1:15, 
                                 function(k){kmeans(res_lap[,2], 
                                                    k, nstart=50,
                                                    iter.max = 15 )$tot.withinss})),
                  stringsAsFactors = F)
ggplot(wss, aes(x = cl, y = wss, color = type)) +
  geom_point() + 
  geom_line() +
  xlab("Number of clusters K") +
  ylab("Total within-clusters sum of squares") +
  scale_color_discrete()


par_lap_cl <- kmeans(par_lap[,2], 2)$cluster
names(par_lap_cl) <- par_lap$id

res_lap_cl <- kmeans(res_lap[,2], 4)$cluster
names(res_lap_cl) <- res_lap$id

tot_lap_cl <- kmeans(tot_lap[,4:5], 5)$cluster
names(tot_lap_cl) <- tot_lap$id


# Construct heterogeneous network
par_mlnw <- constructHetNet(ruprecht_sty, par_lap, par_lap_cl, modules = T)
res_mlnw <- constructHetNet(ruprecht_sty, res_lap, res_lap_cl, modules = T)
tot_mlnw <- constructHetNet(ruprecht_sty, tot_lap, tot_lap_cl, modules = T)

## Run RWHN algorithm
# Recommend to run overnight or on HPC

# par_mlnw$edgelists <- lapply(par_mlnw$edgelists, function(i)
#   i %>% filter_all(any_vars(. != "negative regulation of transcription by RNA polymerase II"))
# )
# par_mlnw$v <- par_mlnw$v[par_mlnw$v$v != "negative regulation of transcription by RNA polymerase II",]

seed_par <- lapply(1:max(par_lap_cl), function(i){
  c(names(par_lap_cl[par_lap_cl == i]))
})

seed_res <- lapply(1:max(res_lap_cl), function(i){
  c(names(res_lap_cl[res_lap_cl == i]))
})

seed_tot <- lapply(1:max(tot_lap_cl), function(i){
  c(names(tot_lap_cl[tot_lap_cl == i]))
})

paste("start par", Sys.time())
rwhn_par <- lapply(seed_par, function(s){
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

paste("start res", Sys.time())
rwhn_res <- lapply(seed_res, function(s){
  calculateRWHN(edgelists = res_mlnw$edgelists,
                #verti = res_mlnw$v[res_mlnw$v$v != "negative regulation of transcription by RNA polymerase II",],
                verti = res_mlnw$v,
                seeds = s,
                transitionProb = 0.7,
                restart = 0.7,
                weight_xy = 0.3,
                weight_yz = 0.7) %>%
    filter(name %in% res_mlnw$v[res_mlnw$v$layer=="func",]$v)
})

paste("start tot", Sys.time())
rwhn_tot <- lapply(seed_tot, function(s){
  calculateRWHN(edgelists = tot_mlnw$edgelists,
                #verti = tot_mlnw$v[tot_mlnw$v$v != "negative regulation of transcription by RNA polymerase II",],
                verti = tot_mlnw$v,
                seeds = s,
                transitionProb = 0.7,
                restart = 0.7,
                weight_xy = 0.3,
                weight_yz = 0.7) %>%
    filter(name %in% tot_mlnw$v[tot_mlnw$v$layer=="func",]$v)
})

end <- Sys.time()
#saveRDS(rwhn_par, "results/data/rwhn_ruprecht_clusters.rds")



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

