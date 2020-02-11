library(parallel)
library(dplyr)
library(igraph)
library(ggplot2)
library(ggridges)
source("src/functions/calculateRWHN.R")

# Parameter ranges
parameters <- expand.grid(trans = seq(0.2, 0.8, 0.1),
                          restart = seq(0.2, 0.8, 0.1),
                          eta_xy = seq(0.2, 0.8, 0.1),
                          eta_yz = seq(0.2, 0.8, 0.1))

# Import multilayer network
mlnw <- readRDS("results/data/mlnw_model.rds")

# Speed up processing with library Parallel
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores, type = "FORK")

# Calculate RWHN with toy network, each cluster set to seed
# and each set of parameters
rwhn <- parApply(cl, parameters, 1, function(i){
  lapply(1:5, function(s){
  calculateRWHN(edgelists = mlnw$edgelists,
                verti = mlnw$v,
                seeds = names(mlnw$fcm$clustering[mlnw$fcm$clustering == s]),
                transitionProb = i[1],
                restart = i[2],
                weight_xy = i[3],
                weight_yz = i[4]) %>%
    filter(name %in% mlnw$v[mlnw$v$layer=="func",]$v)
  })
})
stopCluster(cl)

# Aggrgate and visualise
forhm <- lapply(1:length(rwhn), function(i){
  lapply(1:5, function(x){
    rwhn[[i]][[x]] %>% 
      mutate(l = i,
             seed = x,
             rank = 1:nrow(.))
  }) %>% do.call(rbind, .)
}) %>%  do.call(rbind, .)

g <- ggplot(forhm, aes(x = rank, y = name, fill = name)) +
  geom_density_ridges(color = "grey") +
  facet_wrap(~ seed) +
  ylab("GO Term") +
  xlab("RWHN Rank") +
  theme_ridges() + 
  theme(legend.position = "none",
        axis.text.y = element_blank())

ggsave(filename = "results/figs/controls/parameterChanges.tiff", g, width = 8.3, height = 11.7, units = "in")
