# From CRAN
library(dplyr)
library(tidyr)
library(igraph)
library(e1071)
library(enrichR)
library(patchwork)
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
cfr_sty <- data.table::fread(input = "data/CFR_STY_2016.csv",
                             select = c(1:9, 26:33)
) %>% 
  mutate(id = paste0(`Gene names`, "_", `Swiss-Prot phosphosite`)) %>% 
  filter(!duplicated(id)) %>% 
  rename(`Gene names` = "gene.symbol")
colnames(cfr_sty) <- gsub(" ", "_", colnames(cfr_sty))



# Filter data with > 2 missing values and
# impute phosphorylated sites regulated by EGF or TGF-a
# egf <- cfr_sty %>% 
#   dplyr::select(id, grep("EGF", colnames(cfr_sty))) %>% 
#   filter_missing(allowed = 2) %>% 
#   mutate_at(vars(matches("ratio")), imputeTruncNorm)

impute_fun <- function(i){
  data <- imputeLCMD::impute.QRILC(as.matrix(i))[[1]]
  return(data)
}

egf <- cfr_sty %>% 
  dplyr::select(id, grep("EGF", colnames(cfr_sty))) %>% 
 # filter_at(vars(matches("ratio")), ~ . > 0.5 | . < -1) %>% 
  filter_missing(allowed = 2, colnms = "ratio") %>%
  mutate_at(vars(matches("ratio")), impute_fun)
  
tgf <- cfr_sty %>% 
  dplyr::select(id, grep("TGF", colnames(cfr_sty))) %>% 
#  filter_at(vars(matches("ratio")), ~ . > 0.5 | . < -1) %>%
  filter_missing(allowed = 2, colnms = "ratio") %>%  
  mutate_at(vars(matches("ratio")), impute_fun)


# Cluster based on dynamics of phosphorylated sites
egf_fcm <- fuzzyC(egf, 6)
tgf_fcm <- fuzzyC(tgf, 6)

# Construct heterogeneous network
egf_mlnw <- constructHetNet(cfr_sty, egf, egf_fcm$clustering)
tgf_mlnw <- constructHetNet(cfr_sty, tgf, tgf_fcm$clustering)

## Run RWHN algorithm
# Recommend to run overnight or on HPC
# EGF:
seed_l <- lapply(1:max(egf_fcm$clustering), function(i){
  c(names(egf_fcm$clustering[egf_fcm$clustering == i]))
})
rwhn_egf <- lapply(seed_l, function(s){
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
saveRDS(rwhn_egf, "results/data/rwhn_egf_clusters.rds")

#TGFa
seed_l <- lapply(1:max(tgf_fcm$clustering), function(i){
  c(names(tgf_fcm$clustering[tgf_fcm$clustering == i]))
})
rwhn_tgf <- lapply(seed_l, function(s){
  calculateRWHN(edgelists = tgf_mlnw$edgelists,
                verti = tgf_mlnw$v[tgf_mlnw$v$v != "negative regulation of transcription by RNA polymerase II",],
                seeds = s,
                transitionProb = 0.7,
                restart = 0.7,
                weight_xy = 0.3,
                weight_yz = 0.7) %>%
    filter(name %in% tgf_mlnw$v[tgf_mlnw$v$layer=="func",]$v)
})
saveRDS(rwhn_tgf, "results/data/rwhn_tgf_clusters.rds")

# visualise results with dot plot
rwhn_egf <- readRDS("results/data/rwhn_egf_clusters.rds")
rwhn_tgf <- readRDS("results/data/rwhn_tgf_clusters.rds")

dot_egf <- dotplot_gg(rwhn_egf,remove_common = T, size = 2, n_terms = 20)
dot_egf[[1]] <- dot_egf[[1]] + 
  theme(axis.text.x = element_text(size = 5.5),
        legend.key.size = unit(.3, "cm"),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7)) +
  ggtitle("RWHN ranks from EGF network")

dot_tgf <- dotplot_gg(rwhn_tgf,remove_common = T, size = 2, n_terms = 20) 
dot_tgf[[1]] <- dot_tgf[[1]] + 
  theme(axis.text.x = element_text(size = 5.5),
        legend.key.size = unit(.3, "cm"),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7)) +
  ggtitle("RWHN ranks from TGF-a network")

dots <- dot_egf[[1]] / dot_tgf[[1]] + plot_layout(guides = "collect")
ggsave("results/figs/rwhn_francavilla.tiff", dots, width = 8.3, height = 8, units = "in")


data <- rbind(tgf_fcm$g$data, egf_fcm$g$data) %>% 
  mutate(stim = c(rep("TGF-a", nrow(tgf_fcm$g$data)),
                  rep("EGF", nrow(egf_fcm$g$data)))
  )

raf1 <- ggplot() +
  geom_line(aes(x = Time, y = expression, group = Identifier),
            #data = subset(data, clusterindex %in% c(3, 4, 5)),]
            data = data,
            colour = alpha("grey", 0.7)) + 
  geom_line(aes(x = Time, y = expression, 
                group = Identifier,color = Identifier),
            data = filter(data,
                          grepl("RAF1", Identifier)),
            size = 1) +
  facet_grid(stim ~ clusterindex) +
  scale_x_discrete(limits = c("1", "8", "40", "90")) +
  geom_hline(yintercept = 0, 
             colour="#666666", alpha = 0.5) + 
  theme_minimal() +
  ylab("Phosphorylation changes from baseline") + labs(fill= "Membership")


ggsave("results/figs/rwhn_francavilla_raf1.tiff", raf1, width = 8.3, height = 3.7, units = "in")

# To determine the frequency of common sites
egf_com <- dot_egf[[2]] %>% filter(rank_dif == 0) %>% dplyr::select(egf_rank = rank, name) %>% unique()
tgf_com <- dot_tgf[[2]] %>% filter(rank_dif == 0) %>% dplyr::select(tgf_rank = rank, name) %>% unique()
com <- merge(egf_com, tgf_com, by = "name", all = T) %>% arrange(egf_rank, tgf_rank)
simpl <- simplifyGOReqData()
freq <- sapply(simpl$GO2Gene,length)
freq_p <- freq / sum(freq) * 100
termid <- filter(simpl$GOterms, TERM %in% com$name)
write.csv(termid, "results/data/Table4_commonTermsCFRdata.csv")

