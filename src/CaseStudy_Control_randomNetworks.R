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
source("src/functions/mfuzz-ggplot.R")
source("src/functions/simplifyGO.R")
source("src/functions/simplifyGOReqData.R")
source("src/functions/constructHetNet.R")
source("src/functions/calculateRWHN.R")

# Import multilayer heterogenous network for model data
mlnw_model <- readRDS("results/data/mlnw_model.rds")

fcm <- mlnw_model$fcm
v <- mlnw_model$v
edgelists <- mlnw_model$edgelists

seed_l <- lapply(1:max(fcm$clustering), function(i){
  c(names(fcm$clustering[fcm$clustering == i]))
})

rwhn_random <- lapply(seed_l, function(s){
  list <- lapply(1:100, function(i){
    ranks <- calculateRWHN(edgelists = edgelists,
                           verti = v,
                           seeds = s,
                           transitionProb = 0.7,
                           restart = 0.7,
                           weight_xy = 0.3,
                           weight_yz = 0.7,
                           random = T) %>%
      filter(name %in% v[v$layer=="func",]$v)
    print(i)
    return(ranks)
  })
  return(list)
})

# Import 'actual' result
rwhn <- readRDS("results/data/rwhn_model.rds")

# Wrangle data to long formats and combine
rwhn_random_gg <- lapply(1:length(rwhn), function(x){
  df <- lapply(1:length(rwhn_random[[x]]), function(i){
    df1 <- rwhn_random[[x]][[i]] %>%
      mutate(iter = i,
             rank = 1:nrow(.),
             actual = F)
  }) %>% do.call(rbind, .) %>%
    mutate(seed = x)
})

rwhn_gg <- lapply(1:length(rwhn), function(x){
  df <- rwhn[[x]] %>% 
    mutate(iter = 0,
           rank = 1:nrow(.),
           actual = T,
           seed = x)
})

gg <- lapply(1:length(rwhn_gg), function(x) rbind(rwhn_gg[[x]], rwhn_random_gg[[x]]) %>% 
               mutate(name = factor(x = name, levels = rwhn_gg[[x]]$name, order = T)))

# Visualise with ggplot
gs <- lapply(gg, function(i){
  ggplot(data = subset(i, actual == F), aes(x = rank, y = name, group =1)) +
    stat_density2d(geom = "tile",
                   aes(fill = ..density..),
                   contour = F) +
    geom_line(data = subset(i, actual == T), 
              aes(group = 1), alpha = 0.3) +
    geom_point(data = subset(i, actual == T),
               aes(group = 1,
                   color = rank),
               #show.legend = T,
               size = 0.3) +
    scale_color_gradient(name = "True Rank", low = "red4", high = "white") +
    scale_x_continuous(name = "RWHN rank",
                       breaks = seq(from = 50, to = max(i$rank), by = 50),
                       expand = c(0,0)) +
    scale_y_discrete(name = "GO Term") +
    scale_fill_gradient(name = "KDE of random\nnetworks",
                        limits = c(0, max(sapply(gg, function(i){ max(MASS::kde2d(as.numeric(i$name), i$rank)$z)}))),
                        low = "slategray1", high = "steelblue4",
                        na.value = "white"
    ) +
    theme(axis.title = element_text(size = 5),
         # legend.text = element_text(size = 5),
         # legend.title= element_text(size = 5),
         axis.text.x = element_text(size = 5),
         axis.text.y = element_blank(),
         axis.ticks.y = element_blank(),
         axis.line = element_line(colour = "black"),
         # legend.key.height = unit(2, "cm")
    ) +
    ggtitle(i$seed)
})


#Export

library(patchwork)
layout <- "
AABBCC
AABBCC
#DDEE#
#DDEE#
"

ggsave("results/figs/controls/model_random_nw_patchwork.tiff",
  gs[[1]] + 
    gs[[2]] + 
    gs[[3]] +
    gs[[4]] +
    gs[[5]] +
    plot_layout(nrow = 2, guides = "collect"),
  width = 115,
  height = 65,
  units = "mm",
  dpi = "print"
)


ggsave("results/figs/controls/model_random_nw_cl_legend.tiff", plot = gs[[1]], dpi = "print")



