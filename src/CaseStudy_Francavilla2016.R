# From CRAN
library(dplyr)
library(tidyr)
library(igraph)
library(e1071)
library(enrichR)
library(patchwork)
library(ggplot2)
# From Bioconductor
library(Mfuzz)
library(GOSemSim)
library(AnnotationDbi)
library(GO.db)
# From src/
source("src/functions/overrepresentationAnalysis.R")
source("src/functions/heatmap_RWHNsig.R")
source("src/functions/imputePhosphoData.R")
source("src/functions/mfuzz-ggplot.R")
source("src/functions/simplifyGO.R")
source("src/functions/simplifyGOReqData.R")
source("src/functions/constructHetNet.R")
source("src/functions/calculateRWHN.R")

set.seed(123)

# Import data
cfr_sty <- data.table::fread(input = "data/CFR_STY_2016.csv",
                             select = c(1:9, 26:35)
) %>% 
  mutate(id = paste0(`Gene names`, "_", `Swiss-Prot phosphosite`)) %>% 
  filter(!duplicated(id)) %>% 
  rename("Gene names" = "gene.symbol")
colnames(cfr_sty) <- gsub(" ", "_", colnames(cfr_sty))



# Filter data with > 2 missing values and
# impute phosphorylated sites regulated by EGF or TGF-a

impute_fun <- function(i){
  data <- imputeLCMD::impute.QRILC(as.matrix(i))[[1]]
  return(data)
}

egf <- cfr_sty %>% 
  dplyr::select(id, grep("EGF", colnames(cfr_sty))) %>%
  filter(Regulated_by_EGF == "+") %>% 
 # filter_at(vars(matches("ratio")), ~ . > 0.5 | . < -1) %>% 
  dplyr::select(-Regulated_by_EGF) %>% 
  filter_missing(allowed = 2, colnms = "ratio") %>%
  mutate_at(vars(matches("ratio")), impute_fun)
  
tgf <- cfr_sty %>% 
  dplyr::select(id, grep("TGF", colnames(cfr_sty))) %>% 
  filter(Regulated_by_TGFalfa == "+") %>% 
#  filter_at(vars(matches("ratio")), ~ . > 0.5 | . < -1) %>%
  dplyr::select(-Regulated_by_TGFalfa) %>% 
  filter_missing(allowed = 2, colnms = "ratio") %>%  
  mutate_at(vars(matches("ratio")), impute_fun)


# Cluster based on dynamics of phosphorylated sites
egf_fcm <- fuzzyC(egf, 6)
tgf_fcm <- fuzzyC(tgf, 6)


data <- rbind(tgf_fcm$g$data, egf_fcm$g$data) %>% 
  mutate(stim = c(rep("TGF-a", nrow(tgf_fcm$g$data)),
                  rep("EGF", nrow(egf_fcm$g$data)))
  )

ggcl <- ggplot(data, aes(x = Time, y = expression, group = Identifier, color = maxMembership)) +
  geom_line() +
  facet_grid(stim ~ clusterindex) +
  scale_x_discrete(limits = c("1", "8", "40", "90")) +
  geom_hline(yintercept = 0, 
             colour="#666666", alpha = 0.5) + 
  theme_minimal() +
  ylab("Phosphorylation changes from baseline") + labs(fill= "Membership") +
  scale_colour_gradient(low = "white", high = "red")

ggsave("results/figs/francavilla_clusters.tiff", ggcl,
       width = 8.3, height = 3.7, units = "in", dpi = 300)

# Construct heterogeneous network
egf_mlnw <- constructHetNet(egf, egf_fcm$clustering)
tgf_mlnw <- constructHetNet(tgf, tgf_fcm$clustering)


saveRDS(egf_mlnw, "results/data/egf_mlnw.rds")
saveRDS(tgf_mlnw, "results/data/tgf_mlnw.rds")

## Run RWHN algorithm
# Recommend to run overnight or on HPC
# EGF:
seed_l <- lapply(1:max(egf_fcm$clustering), function(i){
  c(names(egf_fcm$clustering[egf_fcm$clustering == i]))
})
rwhn_egf <- lapply(seed_l, function(s){
  calculateRWHN(edgelists = egf_mlnw$edgelists,
                verti = egf_mlnw$v,
                seeds = s,
                transitionProb = 0.7,
                restart = 0.7,
                weight_xy = 0.3,
                weight_yz = 0.7,
                eps = 1/10^12) %>%
    filter(name %in% egf_mlnw$v[egf_mlnw$v$layer=="func",]$v)
})
saveRDS(rwhn_egf, "results/data/rwhn_egf_clusters.rds")

#TGFa
seed_l <- lapply(1:max(tgf_fcm$clustering), function(i){
  c(names(tgf_fcm$clustering[tgf_fcm$clustering == i]))
})
rwhn_tgf <- lapply(seed_l, function(s){
  calculateRWHN(edgelists = tgf_mlnw$edgelists,
                verti = tgf_mlnw$v,
                seeds = s,
                transitionProb = 0.7,
                restart = 0.7,
                weight_xy = 0.3,
                weight_yz = 0.7,
                eps = 1/10^12) %>%
    filter(name %in% tgf_mlnw$v[tgf_mlnw$v$layer=="func",]$v)
})
saveRDS(rwhn_tgf, "results/data/rwhn_tgf_clusters.rds")



######
# Export, Visualise, etc.
######
rwhn_egf <- readRDS("results/data/rwhn_egf_clusters.rds")
rwhn_tgf <- readRDS("results/data/rwhn_tgf_clusters.rds")

rwhn_egf_df <- lapply(1:length(rwhn_egf), function(i)
  mutate(rwhn_egf[[i]], 
         seed = i,
         rank = 1:nrow(rwhn_egf[[i]])
  )
) %>% 
  do.call(cbind, .) 
readr::write_csv(rwhn_egf_df, "results/data/CFR_egf_rwhn_results.csv")

rwhn_tgf_df <- lapply(1:length(rwhn_tgf), function(i)
  mutate(rwhn_tgf[[i]], 
         seed = i,
         rank = 1:nrow(rwhn_tgf[[i]])
  )
) %>% 
  do.call(cbind, .) 
readr::write_csv(rwhn_tgf_df, "results/data/CFR_tgf_rwhn_results.csv")


## EGF Filter top 5%
sighm_egf <- heatmap_RWHN(rwhn_output = rwhn_egf, ylab = "GOBP")



## tgf Filter top 5%

sighm_tgf <- heatmap_RWHN(rwhn_tgf, ylab = "GOBP")


######################
# Standard ORA analysis
######################

## EGF
enrichedTerms_egf <- overrepresentationAnalysis(clustering = egf_fcm$clustering,
                                                RWHN_sig = sighm_egf,
                                                colours = c("#effff6","#168d49"))


## TGF
enrichedTerms_tgf <- overrepresentationAnalysis(clustering = tgf_fcm$clustering,
                                                RWHN_sig = sighm_tgf, 
                                                colours = c("#effff6","#168d49"))


##################

ggsave("results/figs/Francavilla_RWHN_ORA.pdf",
       sighm_egf + enrichedTerms_egf + sighm_tgf + enrichedTerms_tgf +
         plot_layout(ncol = 2,
                     widths= c(1.8, 1),
                     heights = c(1.3, 1),
                     guides=  "collect"),
       width = 182,
       height = 150,
       units = "mm",
       dpi = "print"
)

#######
# Final visualisation for paper
#######


ora_egf <- enrichedTerms_egf$data %>% 
  dplyr::select(term = Term, value = Adjusted.P.value, cluster) %>% 
  mutate(method = "ORA",
         term = tolower(term))
rwhn_egf <- sighm_egf$data %>% 
  dplyr::select(term = name, value = rank, cluster = seed) %>% 
  mutate(method = "RWHN",
         term = tolower(term))

rwhn_egf_gg <- rwhn_egf %>% 
  rbind(ora_egf) %>% 
  mutate(cluster = ifelse(method == "RWHN", cluster, NA),
         value = ifelse(method == "RWHN", value, NA)
  ) %>% 
  arrange(term) %>% 
  ggplot(aes(y = term,
             x = as.factor(cluster)
  )
  ) +
  geom_tile(aes(fill = value), color = "white") +
  scale_fill_gradient(low = "#5bd670", high = "#def6e2",
                      limits = c(1, 40)) +
  guides(color = FALSE,
         fill = guide_colourbar(title="Rank", reverse = T)) +
  theme_minimal() +
  theme(legend.key.size = unit(.25, "cm"),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 6),
        legend.margin = margin(0,0,0,0, "cm"),
        strip.background = element_blank(),
        strip.text.y = element_blank()
  ) +
  scale_x_discrete(breaks = factor(c(1:5)), limits = factor(c(1:5)), 
                   position = "top") +
  xlab("Cluster") +
  ylab("GOBP") 

ora_egf_gg <- ora_egf %>% 
  rbind(rwhn_egf) %>% 
  mutate(cluster = ifelse(method == "ORA", cluster, NA),
         value = ifelse(method == "ORA", value, NA)
  ) %>% 
  arrange(term) %>% 
  ggplot(aes(y = term,
             x = as.factor(cluster)
  )
  ) +
  geom_tile(aes(fill = value), color = "white") +
  scale_fill_gradient(low = "#845bd6", high = "#e6def6",
                      limits = c(0, 0.05)) +
  guides(color = FALSE,
         fill = guide_colourbar(title="Rank", reverse = T)) +
  theme_minimal() +
  theme(legend.key.size = unit(.25, "cm"),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 5),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 6),
        legend.margin = margin(0,0,0,0, "cm")
  ) +
  scale_x_discrete(breaks = factor(c(1:5)), limits = factor(c(1:5)), 
                   position = "top") +
  xlab("Cluster") +
  ylab("GOBP") 
##
# TGF
ora_tgf <- enrichedTerms_tgf$data %>% 
  dplyr::select(term = Term,  value = Adjusted.P.value, cluster, ) %>% 
  mutate(method = "ORA",
         term = tolower(term))
rwhn_tgf <- sighm_tgf$data %>% 
  dplyr::select(term = name, value = rank, cluster = seed) %>% 
  mutate(method = "RWHN",
         term = tolower(term))

rwhn_tgf_gg <- rwhn_tgf %>% 
  rbind(ora_tgf) %>% 
  mutate(cluster = ifelse(method == "RWHN", cluster, NA),
         value = ifelse(method == "RWHN", value, NA)
  ) %>% 
  arrange(term) %>% 
  ggplot(aes(y = term,
             x = as.factor(cluster)
  )
  ) +
  geom_tile(aes(fill = value), color = "white") +
  scale_fill_gradient(low = "#5bd670", high = "#def6e2",
                      limits= c(1,40)) +
  guides(color = FALSE,
         fill = guide_colourbar(title="Rank", reverse = T)) +
  theme_minimal() +
  theme(legend.key.size = unit(.25, "cm"),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        axis.text = element_text(size = 5),
        axis.title.x = element_text(size = 6),
        axis.title.y = element_blank(),
        legend.margin = margin(0,0,0,0, "cm"),
        strip.background = element_blank(),
        strip.text.y = element_blank()
  ) +
  scale_x_discrete(breaks = factor(c(1:5)), limits = factor(c(1:5)), 
                   position = "top") +
  xlab("Cluster") +
  ylab("GOBP") 

ora_tgf_gg <- ora_tgf %>% 
  rbind(rwhn_tgf) %>% 
  mutate(cluster = ifelse(method == "ORA", cluster, NA),
         value = ifelse(method == "ORA", value, NA)
  ) %>% 
  arrange(term) %>% 
  ggplot(aes(y = term,
             x = as.factor(cluster)
  )
  ) +
  geom_tile(aes(fill = value), color = "white") +
  scale_fill_gradient(low = "#845bd6", high = "#e6def6",
                      limits = c(0, 0.05)) +
  guides(color = FALSE,
         fill = guide_colourbar(title="Rank", reverse = T)) +
  theme_minimal() +
  theme(legend.key.size = unit(.25, "cm"),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 5),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 6),
        legend.margin = margin(0,0,0,0, "cm")
  ) +
  scale_x_discrete(breaks = factor(c(1:5)), limits = factor(c(1:5)), 
                   position = "top") +
  xlab("Cluster") +
  ylab("GOBP") 

ggsave(
  "results/figs/Francavilla_RWHN_ORA_v2.tiff",
  rwhn_egf_gg +
    ora_egf_gg + 
    rwhn_tgf_gg + 
    ora_tgf_gg +
    plot_layout(nrow = 1, guides= "collect"),
  width = 7,
  height = 9,
  units = "in",
  dpi = 300
)
ggsave(
  "results/figs/Francavilla_RWHN_ORA_v2.pdf",
  rwhn_egf_gg +
    ora_egf_gg + 
    rwhn_tgf_gg + 
    ora_tgf_gg +
    plot_layout(nrow = 1, guides= "collect"),
  width = 7,
  height = 9,
  units = "in",
  dpi = 300
)
