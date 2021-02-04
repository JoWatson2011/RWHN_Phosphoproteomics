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
  filter(!duplicated(id))
colnames(cfr_sty) <- gsub(" ", "_", colnames(cfr_sty))

psp <- readr::read_tsv("data/Regulatory_sites.txt", skip = 2) %>% 
  filter(!is.na(ON_PROCESS),
         ORGANISM == "human",
         grepl("-p", MOD_RSD)) %>% 
  separate_rows(ON_PROCESS, sep = "; ") %>% 
  mutate(id_site = paste0(GENE,
                          "_",
                          gsub("-.*",
                               "",
                               MOD_RSD))) %>% 
  filter(ON_PROCESS == "cell growth, inhibited") %>%
  dplyr::select(id_site) %>% 
  distinct()

psp_vec <- rep(7, nrow(psp)) 
names(psp_vec) <- psp$id_site
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

egf_fcm$clustering <- c(egf_fcm$clustering,
                        psp_vec[!names(psp_vec)
                                                 %in% 
                                                   names(egf_fcm$clustering)]
)
egf_wNew <- rbind(egf,
             data.frame(id = names(egf_fcm$clustering)[!names(egf_fcm$clustering) %in% egf$id],
                        EGF_1_min_average_ratio = sample(seq(0.9, 1.1, by = 0.01),
                                                         length(egf_fcm$clustering[egf_fcm$clustering == 7]),
                                                         replace = T),
                        EGF_8_min_average_ratio = sample(seq(0.9, 1.1, by = 0.01),
                                                         length(egf_fcm$clustering[egf_fcm$clustering == 7]),
                                                         replace = T),
                        EGF_40_min_average_ratio = sample(seq(0.9, 1.1, by = 0.01),
                                                          length(egf_fcm$clustering[egf_fcm$clustering == 7]),
                                                          replace = T),
                        EGF_90__min_average_ratio = sample(seq(0.9, 1.1, by = 0.01),
                                                           length(egf_fcm$clustering[egf_fcm$clustering == 7]),
                                                           replace = T)
                        )
             )
                        
tgf_fcm <- fuzzyC(tgf, 6)
tgf_fcm$clustering <- c(tgf_fcm$clustering,
                        psp_vec[!names(psp_vec)
                                                 %in% 
                                                   names(tgf_fcm$clustering)]
)
tgf_wNew <- rbind(tgf,
             data.frame(id = names(tgf_fcm$clustering)[!names(tgf_fcm$clustering) %in% tgf$id],
                        TGFalfa_1_min_average_ratio = sample(seq(0.9, 1.1, by = 0.01),
                                                             length(tgf_fcm$clustering[tgf_fcm$clustering == 7]),
                                                             replace = T),
                        TGFalfa_8_min_average_ratio = sample(seq(0.9, 1.1, by = 0.01),
                                                             length(tgf_fcm$clustering[tgf_fcm$clustering == 7]),
                                                             replace = T),
                        TGFalfa_40_min_average_ratio = sample(seq(0.9, 1.1, by = 0.01),
                                                              length(tgf_fcm$clustering[tgf_fcm$clustering == 7]),
                                                              replace = T),
                        TGFalfa_90_min_average_ratio = sample(seq(0.9, 1.1, by = 0.01),
                                                              length(tgf_fcm$clustering[tgf_fcm$clustering == 7]),
                                                              replace = T)
             )
)

# Construct heterogeneous network
egf_mlnw <- constructHetNet(phosphoData = egf_wNew, clustering = egf_fcm$clustering)
tgf_mlnw <- constructHetNet(phosphoData = tgf_wNew, clustering = tgf_fcm$clustering)


saveRDS(egf_mlnw, "results/data/egf_mlnw_spikeIn.rds")
saveRDS(tgf_mlnw, "results/data/tgf_mlnw_spikeIn.rds")

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
saveRDS(rwhn_egf, "results/data/rwhn_egf_spikeIn.rds")

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
saveRDS(rwhn_tgf, "results/data/rwhn_tgf_spikeIn.rds")



######
# Export, Visualise, etc.
######
# rwhn_egf <- readRDS("results/data/rwhn_egf_spikeIn_Growth.rds")
# rwhn_tgf <- readRDS("results/data/rwhn_tgf_spikeIn_Growth.rds")
rwhn_egf <- readRDS("results/data/rwhn_egf_spikeIn_Growth.rds")
rwhn_tgf <- readRDS("results/data/rwhn_tgf_spikeIn_Growth.rds")


## EGF Filter top 5%
sighm_egf <- heatmap_RWHN(rwhn_egf, ylab = "GOBP Term")

## tgf Filter top 5%

sighm_tgf <- heatmap_RWHN(rwhn_tgf, "GOBP Term")

######################
# Standard ORA analysis
######################

## EGF
enrichedTerms_egf <- overrepresentationAnalysis(clustering = egf_fcm$clustering[egf_fcm$clustering == 7],
                                                RWHN_sig = sighm_egf,
                                                colours = c("#effff6","#168d49"))


## TGF
enrichedTerms_tgf <- overrepresentationAnalysis(clustering = tgf_fcm$clustering[tgf_fcm$clustering ==7],
                                                RWHN_sig = sighm_tgf,
                                                colours = c("#effff6","#168d49"))

# Visualise results of Cl 7 in all analysis types.


ora_vis <- list(EGF = enrichedTerms_egf$data[enrichedTerms_egf$data$cluster==7,],
                `TGF` = enrichedTerms_tgf$data[enrichedTerms_tgf$data$cluster==7,])
rwhn_vis <- list(EGF = sighm_egf$data[sighm_egf$data$seed == 7,],
                 TGF = sighm_tgf$data[sighm_tgf$data$seed == 7,])

vis <- rbind(
  do.call(rbind, lapply(names(ora_vis), function(i){
    ora_vis[[i]] %>% 
      arrange(Adjusted.P.value) %>% 
      mutate(rank = 1:nrow(.)) %>% 
      dplyr::select(Term, rank) %>% 
      mutate(Data = i,
             Analysis = "ORA")
  })
  ),
  do.call(rbind, lapply(names(rwhn_vis), function(i){
    rwhn_vis[[i]] %>%
      dplyr::select(Term = name, rank) %>% 
      mutate(Data = i,
             Analysis = "RWHN")
  })
  )
)
ggplot(vis, aes(x = as.factor(data), y = reorder(Term, rank))) +
  geom_tile(aes(fill = rank), colour = "white") +
  scale_fill_gradient(breaks = seq(from = 1,  to = max(vis$rank), by = 5),
                      low = "#b8dfeb", high = "#06779a") +
  
  scale_x_discrete(position = "top") +
  guides(color = FALSE,
         fill = guide_colourbar(title="Rank", reverse = T)) +
  xlab("Analysis type") +
  ylab("GO") +
  facet_wrap(~analysis) +
  theme(legend.key.size = unit(.25, "cm"),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5), 
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 6),
        panel.background = element_rect(fill = "black"), 
        panel.grid = element_blank(),
        legend.margin = margin(0,0,0,0, "cm"),
        panel.spacing.x=unit(0.1, "lines")
  )

num <- vis %>% 
  group_by(Data, Analysis) %>% 
  summarise(`No. terms` = n(), .groups = "keep")  %>% 
  arrange(Analysis)

# For clarity, excluded terms that appear in <3 of the scenarios.
# Really... I'm only interested in the 
vis <- rbind(
  do.call(rbind, lapply(names(ora_vis), function(i){
    ora_vis[[i]] %>% 
      arrange(Adjusted.P.value) %>% 
      mutate(rank = 1:nrow(.)) %>% 
      dplyr::select(Term, rank) %>% 
      mutate(data = i,
             analysis = "ORA")
  })
  ),
  do.call(rbind, lapply(names(rwhn_vis), function(i){
    rwhn_vis[[i]] %>%
      dplyr::select(Term = name, rank) %>% 
      mutate(data = i,
             analysis = "RWHN")
  })
  )
) %>% 
  group_by(Term) %>% 
  filter(n() > 2)

gg_vis <- ggplot(vis, aes(x = as.factor(data), y = reorder(Term, rank))) +
  geom_tile(aes(fill = rank), colour = "white") +
  scale_fill_gradient(breaks = seq(from = 1,  to = max(vis$rank), by = 5),
                      labels = c("1", rep("", length(seq(from = 1,  to = max(vis$rank), by = 5))-1)),
                      low = "#b8dfeb", high = "#06779a") +
  guides(color = FALSE,
         fill = guide_colourbar(title="Rank", reverse = T)) +
  xlab("Data subset from  Francavilla et al. (2016)") +
  ylab("GOBP") +
  facet_wrap(~analysis)+
  theme(legend.key.size = unit(.25, "cm"),
        strip.text = element_text(size = 7.5),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5), 
        axis.text.y = element_text(size = 5),
        axis.text.x = element_text(size = 5, angle = 45),
        axis.title = element_text(size = 5),
        panel.background = element_rect(fill = "black"), 
        panel.grid = element_blank(),
        plot.margin = margin(0,0,0,0, "cm"),
        legend.margin = margin(0,0,0,0, "cm")
  ) +
  scale_x_discrete(position = "top") 

gg_vis / gridExtra::tableGrob(num, row = NULL) + plot_annotation(tag_levels = "A")

ggsave("results/figs/CFR_SPIKEIN_GROWTH.tiff",
       gg_vis /
         gridExtra::tableGrob(num, 
                              row = NULL,
                              theme = gridExtra::ttheme_default(
                                base_size = 8
                              )) + 
         plot_annotation(tag_levels = "A")  & 
         theme(plot.tag = element_text(size = 8)),
       width = 3.33, height = 4, units = "in", dpi = 300)
