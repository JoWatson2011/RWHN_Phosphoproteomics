# From CRAN
library(dplyr)
library(tidyr)
library(igraph)
library(e1071)
library(enrichR)
library(ggplot2)
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
colnames(ruprecht_sty) <- c("amino.acid", "localisation.prob", "position",
                            "HL_R1", "HL_R2", "HL_R3", "HL_R4",
                            "ML_R1", "ML_R2", "ML_R3", "ML_R4",
                            "SignificantML", "SignificantHL",
                            "gene.symbol", "id")

# PARENTAL + LAP M/L: Filter NAs and non significant
par_lap <- ruprecht_sty %>%
  dplyr::select(id, grep("ML", colnames(.))) %>%
  filter(SignificantML == "+" ) %>%
  dplyr::select(-SignificantML
  ) %>%
  pivot_longer(-id,
               names_to = "rep",
               values_to = "ratio") %>%
  mutate(ratio = imputeLCMD::impute.QRILC(as.matrix(ratio))[[1]],
         rep = gsub("_R[1234]$", "", rep)) %>% 
  group_by(id, rep) %>%
  summarise(ratio = median(ratio, na_rm = T)) %>%
  pivot_wider(id_cols = id,
              names_from = rep,
              values_from = ratio)


# RESISTANT  + LAP H/L
res_lap <- ruprecht_sty %>%
  dplyr::select(id, grep("HL", colnames(.))) %>%
  filter(SignificantHL == "+" )  %>% 
  dplyr::select(-SignificantHL
  ) %>%
  pivot_longer(-id,
               names_to = "rep",
               values_to = "ratio") %>%
  mutate(ratio = imputeLCMD::impute.QRILC(as.matrix(ratio))[[1]],
         rep = gsub("_R[1234]$", "", rep)) %>%
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
  # dplyr::select(-c(SignificantML,
  #                  SignificantHL)
  # ) %>%
  pivot_longer(-c(id, SignificantML, SignificantHL),
               names_to = "rep",
               values_to = "ratio") %>%
  mutate(ratio = imputeLCMD::impute.QRILC(as.matrix(ratio))[[1]],
         rep = gsub("_R[1234]$", "", rep)) %>%
  group_by(id, rep) %>%
  mutate(ratio = median(ratio, na_rm = T)) %>%
  unique() %>% 
  pivot_wider(#id_cols = id,
              names_from = rep,
              values_from = ratio)

# Cluster based on dynamics of phosphorylated sites
set.seed(1)
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


sites <- data.frame(
  gene = c("HNRNPU", "SRRM2", "DDX23", "NCBP1", "SRRM2", 
           "PCBP1", "PRPF4B", "SRRM2", "CDK1", "CDK1", 
           "SRC", "IGR1R", "LDAH", "ALDOA", "ALDOA", 
           "PFKP", "GAPDH", "ENO1", "PKM", "PCAM", "PYGB",
           "PGM1", "PGM2", "PFKB2", "HSP90AB1", "BAD", "FOXO3", 
           "SRC"),
  sites = c("S59", "S1132", "S107", "S22", "S1987", "S264", 
            "S431", "S970", "T161", "NA", "NA", "NA", "Y16", 
            "S39", "S46", "S386", "S83", "Y44", "Y175", "Y92", 
            "T59", "S117", "S165", "S466", "S255", "S99", "T32",
            "S17"),
  paperDescription = c("Splicesome", "Splicesome", "Splicesome", 
                       "Splicesome", "Splicesome", "Splicesome", 
                       "Splicesome", "Splicesome", "CC. strongest observed activation of any kinase", 
                       NA, NA, NA, "Glycotic enzymes. Highly incresed in resistance",
                       "Glycotic enzymes. Highly incresed in resistance", "Glycotic enzymes. 
                       Highly incresed in resistance", "Glycotic enzymes. Highly incresed in resistance", 
                       "Glycotic enzymes. Highly incresed in resistance", 
                       "Glycotic enzymes. Highly incresed in resistance", 
                       "Glycotic enzymes. Highly incresed in resistance", 
                       "Glycotic enzymes. Highly incresed in resistance", 
                       "Glycogen catabolism", "Glycogen catabolism", "Glycogen catabolism",
                       NA, NA, NA, NA, NA),
  stringsAsFactors = F
) %>% 
  mutate(id = paste0(gene, "_", sites)) %>% 
  merge(., data.frame(id = names(tot_lap_cl), tot_cl = tot_lap_cl), by = "id",
        all.x = T) %>% 
  merge(., data.frame(id = names(res_lap_cl), res_cl = res_lap_cl), by = "id",
        all.x = T)
write.csv(sites, "results/data/ruprecht_clusters.csv")


# Construct heterogeneous network
par_mlnw <- constructHetNet(stytxt = ruprecht_sty, 
                            phosphoData =  par_lap, 
                            clustering = par_lap_cl,
                            modules = T,
                            enrichrLib =  "KEGG_2019_Human",
                            pval = 0.05)

res_mlnw <- constructHetNet(stytxt = ruprecht_sty,
                            phosphoData = res_lap,
                            clustering =  res_lap_cl,
                            modules = T,
                            enrichrLib =  "KEGG_2019_Human",
                            pval = 0.05)

tot_mlnw <- constructHetNet(stytxt = ruprecht_sty,
                            phosphoData =  tot_lap[,-c(2:3)],
                            clustering =  tot_lap_cl,
                            modules = T,
                            enrichrLib =  "KEGG_2019_Human",
                            pval = 0.05)

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
                verti = tot_mlnw$v,
                seeds = s,
                transitionProb = 0.7,
                restart = 0.7,
                weight_xy = 0.3,
                weight_yz = 0.7) %>%
    filter(name %in% tot_mlnw$v[tot_mlnw$v$layer=="func",]$v)
})

end <- Sys.time()

saveRDS(rwhn_par, "results/data/rwhn_ruprecht_par.rds")
saveRDS(rwhn_res, "results/data/rwhn_ruprecht_res.rds")
saveRDS(rwhn_tot, "results/data/rwhn_ruprecht_tot.rds")

rwhn_par <- readRDS("results/data/rwhn_ruprecht_par.rds")
rwhn_res <- readRDS("results/data/rwhn_ruprecht_res.rds")
rwhn_tot <- readRDS("results/data/rwhn_ruprecht_tot.rds")


# visualise results with dot plot
dot_par_all <- dotplot_gg(rwhn_par, n_terms = 30, remove_common = F)
dot_res_all <- dotplot_gg(rwhn_res, n_terms = 30, remove_common = F)
dot_tot_all <- dotplot_gg(rwhn_tot, n_terms = 30, remove_common = F)
dot_par_diff <- dotplot_gg(rwhn_par, n_terms = 20, remove_common = T)
dot_res_diff <- dotplot_gg(rwhn_res, n_terms = 30, remove_common = T)
dot_tot_diff <- dotplot_gg(rwhn_tot, n_terms =30, remove_common = T)


dot_res_diff[[1]]  <- dot_res_diff[[1]] + 
    ggtitle("RWHN results from \`lapatanib-resistant\` network")
dot_tot_diff[[1]]  <- dot_tot_diff[[1]] +
  ggtitle("RWHN results from \'Total\' network")

gg <- dot_tot_diff[[1]]  / dot_res_diff[[1]] + plot_layout(guides = "collect")
gg <- gg + plot_annotation(tag_levels = "A")

ggsave("results/figs/rwhn_ruprecht_kegg_patchwork.tiff", width = 11.7, height = 10, units= "in")

ggsave("results/figs/rwhn_ruprecht_parental.tiff", dot_par_all[[1]], width = 11.7, height = 7, units= "in")
ggsave("results/figs/rwhn_ruprecht_resistant.tiff", dot_res_diff[[1]], width = 11.7, height = 7, units= "in")
ggsave("results/figs/rwhn_ruprecht_tot.tiff", dot_tot_diff[[1]], width = 11.7, height = 7, units= "in")

######
par_mlnw_GO <- constructHetNet(stytxt = ruprecht_sty,
                            phosphoData =  par_lap, 
                            clustering = par_lap_cl,
                            modules = T,
                            pval = 0.05)

res_mlnw_GO <- constructHetNet(stytxt = ruprecht_sty,
                            phosphoData = res_lap,
                            clustering =  res_lap_cl,
                            modules = T,
                            pval = 0.05)

tot_mlnw_GO <- constructHetNet(stytxt = ruprecht_sty,
                            phosphoData =  tot_lap[,-c(2:3)],
                            clustering =  tot_lap_cl,
                            modules = T,
                            pval = 0.05)



paste("start res", Sys.time())
rwhn_res_GO <- lapply(seed_res, function(s){
  calculateRWHN(edgelists = res_mlnw_GO$edgelists,
                verti = res_mlnw_GO$v,
                seeds = s,
                transitionProb = 0.7,
                restart = 0.7,
                weight_xy = 0.3,
                weight_yz = 0.7) %>%
    filter(name %in% res_mlnw_GO$v[res_mlnw_GO$v$layer=="func",]$v)
})

paste("start tot", Sys.time())
rwhn_tot_GO <- lapply(seed_tot, function(s){
  calculateRWHN(edgelists = tot_mlnw_GO$edgelists,
                verti = tot_mlnw_GO$v,
                seeds = s,
                transitionProb = 0.7,
                restart = 0.7,
                weight_xy = 0.3,
                weight_yz = 0.7) %>%
    filter(name %in% tot_mlnw_GO$v[tot_mlnw_GO$v$layer=="func",]$v)
})
paste("end", Sys.time())

saveRDS(rwhn_res_GO, "results/data/rwhn_res_GO.rds")
saveRDS(rwhn_tot_GO, "results/data/rwhn_tot_GO.rds")


dot_res_diff <- dotplot_gg(rwhn_res_GO, n_terms = 30, remove_common = T, size = 2)
dot_res_diff[[1]]  <- dot_res_diff[[1]] + 
  theme(axis.text.x = element_text(size = 7)) +
  ggtitle("RWHN results from \`lapatanib-resistant\` network")
dot_tot_diff <- dotplot_gg(rwhn_tot_GO, n_terms =30, remove_common = T, size = 2)
dot_tot_diff[[1]]  <- dot_tot_diff[[1]] + theme(axis.text.x = element_text(size = 7))  +
  ggtitle("RWHN results from \'Total\' network")


gg <- dot_tot_diff[[1]]  / dot_res_diff[[1]] + plot_layout(guides = "collect")
gg <- gg + plot_annotation(tag_levels = "A")

ggsave("results/figs/rwhn_ruprecht_patchwork.tiff", width = 11.7, height = 10, units= "in")
# 
# ggsave("results/figs/rwhn_ruprecht_resistant_GO.tiff", dot_res_diff[[1]], width = 11.7, height = 7, units= "in")
# ggsave("results/figs/rwhn_ruprecht_tot_GO.tiff", dot_tot_diff[[1]], width = 11.7, height = 7, units= "in")
