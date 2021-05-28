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

set.seed(1)

sty_sig <- readr::read_csv("data/Santra_etal_2016_Sig.csv") %>% 
  rowwise() %>% 
  filter(sum(across(WT:M1)) > 0) %>% 
  ungroup() %>% 
  #dplyr::select(`Gene names`, `Positions within proteins`) %>% 
  mutate(across(c(`Gene names`, `Positions within proteins`), 
                gsub, pattern = ";.*", replacement = "")) %>% 
  mutate(id = paste0(`Gene names`, "_", `Positions within proteins`)) %>% 
  unique() #%>%                       ### In case I decide to cluster by compartment
  # dplyr::select(id, WT:M1) %>% 
  # pivot_longer(-id) %>%  
  # filter(value == 0) %>% 
  # group_by(id) %>% 
  # summarise(name = paste(name, collapse = " "))

sty <- readr::read_csv("data/Santra_etal_2016.csv") %>% 
  mutate(across(c(`Gene names`, `Positions within proteins`), 
                gsub, pattern = ";.*", replacement = "")) %>% 
  mutate(id = paste0(`Gene names`, "_", `Positions within proteins`)) %>% 
  filter(id %in% sty_sig$id)  %>% 
  dplyr::select(id, starts_with("Intensity ")) %>% 
  mutate(across(where(is.numeric), ~ ifelse(. == 0, NA, .))) %>% 
  mutate(across(where(is.numeric), log))

sty[,grep("Intensity ", colnames(sty))] <- limma::normalizeBetweenArrays(sty[,grep("Intensity ", colnames(sty))],
                                                                          method = "quantile")
sty_fin <- sty %>% 
  pivot_longer(-id) %>%
  mutate(name =  gsub("_[1-3]$", "", name)) %>% 
  mutate(value = imputeQRLIC.mod(as.matrix(value))) %>% 
  group_by(id) %>% 
  mutate(value = (value - mean(value, na.rm = T))/
           sd(value, na.rm = T)) %>% 
  unique() %>% 
  pivot_wider(values_fn = median) %>% 
  na.omit() %>% 
  ungroup()


set.seed(1)
wss <- data.frame(
  cl = 1:15,
  wss = sapply(1:15,
           function(k) {
             kmeans(sty_fin[, 2:7],
                    k, nstart = 50,
                    iter.max = 15)$tot.withinss
           })
)
elbow <- ggplot(wss, aes(x = cl, y = wss)) +
  geom_point() +
  geom_line() +
  xlab("Number of clusters (k)") +
  ylab("Total within-clusters sum of squares") +
  theme_minimal() +
  theme(legend.justification = c(1, 1),
        legend.position = c(1, 1),
        legend.box.background = element_rect(color = "black", fill = "white"),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        axis.title = element_text(size = 5),
        axis.text = element_text(size = 5)
  )

ggsave("results/figs/Santra_elbow.pdf", elbow, width = 4)

  set.seed(1)
  cl <- kmeans(sty_fin[,2:7], 6)$cluster
  names(cl) <- sty_fin$id
  
  sty_fin$cluster <- cl
  gg_cl <- sty_fin %>%
    pivot_longer(-c(id, cluster)) %>%
    mutate(name = gsub("Intensity ", "", name)) %>% 
    arrange(cluster) %>%
    ggplot(aes(x =  name, y = forcats::fct_reorder(id, cluster), fill = value)) +
    geom_tile() +
    scale_x_discrete(limits = c("Con", "WT", "KDEL", "M1", "CD8", "LCK"),
                     labels = c("Endogenous", "Unlocalised",
                                "Endoplasmic Reticulum",
                                "Golgi apparatus", "Plasma membrane", 
                                "Lipid Rafts")) +
    scale_fill_gradient2(low = "red", mid = "black", high = "green", 
                         name = "Z-score\nIntensity") +
    facet_wrap(~cluster, ncol = 1, scales = "free_y", strip.position = "right") +
    ylab("Regulated sites") +
    xlab("HRASV12 localisation") +
    theme(
      legend.key.size = unit(.25, "cm"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = 5),
      legend.text = element_text(size = 5), 
      axis.title = element_text(size = 6),
      legend.margin = margin(0,0,0,0, "cm"),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 5, angle = 45, vjust = 1, hjust = 1),
      panel.spacing.y=unit(0.1, "lines"),
        axis.ticks.y = element_blank(),
      strip.text = element_text(size = 5)
    )


ggsave("results/figs/Santra_clusters.tiff",
       elbow + gg_cl + plot_annotation(tag_levels = "A"),
       width = 7,
       dpi = 300
)

mlnw <- constructHetNet(phosphoData = dplyr::select(sty_fin, -cluster),
                        clustering =  cl,
                        modules = T,
                        stringPath = "data/STRINGexpmtgene_highconf.rds",
                        pval = 0.05)
saveRDS(mlnw, "results/data/Santra_mlnw.rds")

seed <- lapply(1:max(cl), function(i){
  names(cl[cl==i])
})

rwhn <- lapply(seed, function(s){
  calculateRWHN(edgelists = mlnw$edgelists,
                verti = mlnw$v,
                seeds = s,
                transitionProb = 0.7,
                restart = 0.7,
                weight_xy = 0.3,
                weight_yz = 0.7,
  ) %>%
    filter(name %in% mlnw$v[mlnw$v$layer=="func",]$v)
})

saveRDS(rwhn, "results/data/rwhn_santra.rds")

#####
rwhn <- readRDS("results/data/rwhn_santra.rds")
names(rwhn) <- as.character(1:length(rwhn))
sighm <- heatmap_RWHN(rwhn_output = rwhn[c("1", "4", "6")],
                      ylab = "GOBP Term")


ora <- overrepresentationAnalysis(clustering = cl[cl %in% c(1,4,6)],
                                  RWHN_sig = sighm,
                                  simplify= F,
                                  colours = c("#effff6","#168d49")
)

ggsave("results/figs/Santra_Clusters.pdf",
  gg_cl +theme(plot.margin = unit(c(0,0,0,0), "cm")),
  width = 2,
  height = 3.3,
  units = "in",
  dpi = "print"
)
ggsave("results/figs/Santra_RWHN.pdf",
       sighm +theme(plot.margin = unit(c(0,0,0,0), "cm"),
                    axis.text.y = element_text(size = 4.5)),
       width = 3.5,
       height = 7,
       units = "in",
       dpi = "print"
)
ggsave("results/figs/Santra_ORA.pdf",
       ora + 
         scale_y_discrete(labels = gsub(" \\(go:.*$", "", unique(ora$data$Term))) +
         theme(axis.text.y = element_text(size = 4.5),
                   plot.margin = unit(c(0,0,0,0), "cm")),
       width = 3.5,
       height = 7,
       units = "in",
       dpi = "print"
)

##########
# Final visualisation for paper
##########
ora <- ora$data %>% 
  dplyr::select(term = Term,  value = Adjusted.P.value, cluster, ) %>% 
  mutate(method = "ORA")
rwhn <- sighm$data %>% 
  dplyr::select(term = name, value = rank, cluster = seed) %>% 
  mutate(method = "RWHN")

rwhn_gg <- rwhn %>% 
  rbind(ora) %>% 
  mutate(cluster = ifelse(method == "RWHN", cluster, NA),
         value = ifelse(method == "RWHN", value, NA)
  ) %>% 
  arrange(term) %>% 
  ggplot(aes(y = term,
             x = as.factor(cluster)
  )
  ) +
  geom_tile(aes(fill = value), color = "white") +
  scale_fill_gradient(low = "#5bd670", high = "#def6e2") +
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
  scale_x_discrete(#breaks = factor(c(1:5)), limits = factor(c(1:5)), 
                   position = "top", na.translate = F) +
  xlab("Cluster") +
  ylab("GOBP") 

ora_gg <- ora %>% 
  rbind(rwhn) %>% 
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
  scale_x_discrete(#breaks = factor(c(1:5)), limits = factor(c(1:5)), 
                   position = "top", na.translate = F) +
  xlab("Cluster") +
  ylab("GOBP") 

ggsave(
  "results/figs/Santra_RWHN_ORA_v2.tiff",
    rwhn_gg + 
    ora_gg +
    plot_layout(nrow = 1, guides= "collect"),
  width = 7,
  height = 9,
  units = "in",
  dpi = 300
)

ggsave(
  "results/figs/Santra_RWHN_ORA_v2.pdf",
  rwhn_gg + 
    ora_gg +
    plot_layout(nrow = 1, guides= "collect"),
  width = 7,
  height = 9,
  units = "in",
  dpi = 300
)
