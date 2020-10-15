library(parallel)
library(dplyr)
library(tidyr)
library(igraph)
library(ggplot2)
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
cl <- makeCluster(no_cores)

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
      mutate(l = paste(parameters[i,], collapse = " "),
             seed = x,
             rank = 1:nrow(.)) %>% 
      select(-V1) 
  }) %>% do.call(rbind, .)
}) %>%  do.call(rbind, .)


hm <- ggplot(tmp, aes(x = rank, 
                      y = reorder(name, rank), 
                      fill = term_at_pos
                      )) +
  geom_tile() +
  facet_wrap(~ seed, ncol = 1) +
  scale_fill_gradient("Frequency\nat rank", 
                      low = "#FFFFFF", 
                      high = "#132B43") +
  scale_y_discrete("GO Term") +
  scale_x_continuous("Rank") +
  geom_segment(aes(x = 0, y = 10, xend = 10, yend = 10),
               color = "red", size = 0.2) +
  geom_segment(aes(x = 10, y = 0, xend = 10, yend = 10), 
               color = "red", size = 0.2) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        axis.text.y = element_text(size = 2),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        legend.position = "left",
        legend.text = element_text(size = 4),
        legend.title = element_text(size = 4))
bar <- tmp %>% 
  group_by(name) %>% 
  top_n(1) %>% 
  ggplot(aes(x =  reorder(name, rank),
             y = term_at_pos, 
             fill = term_at_pos
            )) +
  geom_col() + 
  coord_flip() +
  scale_fill_gradient(low = "#FFFFFF", high = "#132B43") +
  facet_wrap(~ seed, ncol = 1) +
  scale_y_continuous("Times at most frequent rank") +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1)
        )

library(patchwork)
ggsave("results/figs/controls/parameterChanges_all.pdf", hm + bar, 
       width = 8.3, height = 11.7+2*(11.7/3), units = "in")

######


# Parameter ranges
parameters <- data.frame(trans = c(seq(0.2, 0.8, 0.1), rep(0.5, 21)),
                         restart = c(rep(0.5, 7), seq(0.2, 0.8, 0.1), rep(0.5, 14)),
                         eta_xy = c( rep(0.5, 14), seq(0.2, 0.8, 0.1), rep(0.5, 7)),
                         eta_yz = c(rep(0.5, 21), seq(0.2, 0.8, 0.1)), 
                         param_altered = c(rep("Transition Probability", 7),
                                           rep("Restart Probability", 7),
                                           rep("Protein Node Weighting", 7),
                                           rep("Function Node Weighting.", 7)),
                         stringsAsFactors = F
) 

no_cores <- detectCores() - 1
cl <- makeCluster(no_cores, type = "FORK")

# Calculate RWHN with toy network, each cluster set to seed
# and each set of parameters
rwhn <- parApply(cl, parameters[,1:4], 1, function(i){
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

forhm <- lapply(1:length(rwhn), function(i){
  lapply(1:5, function(x){
    rwhn[[i]][[x]] %>% 
      mutate(l = paste(parameters[i,1:4], collapse = " "),
             param_altered = parameters[i,5],
             seed = x,
             rank = 1:nrow(.)) %>% 
      select(-V1) 
  }) %>% do.call(rbind, .)
}) %>%  do.call(rbind, .) %>% 
  group_by(name, rank, seed, param_altered) %>% 
  summarise(n = n()) %>% 
  ungroup %>% 
  mutate(n = 100/7*n,
         seed = paste("Cluster", seed))


g <- ggplot(forhm, aes(x = rank, y = name, fill = n)) +
  facet_grid(seed ~ param_altered) 
  geom_tile() +
  scale_fill_gradient("", low = "red", high = "blue") +
  ylab("GO Term") +
  xlab("RWHN Rank") +
  theme_ridges() + 
  theme(legend.text = element_text(size = 8),
        plot.title = element_text(size = 8),
        axis.text.y = element_blank(), 
        axis.text.x = element_text(size = 8),
        axis.title = element_text(size = 8), 
        axis.title.x = element_text(hjust = 0.5),
        axis.title.y = element_text(hjust = 0.5),
        strip.text = element_text(size = 8)
  )

ggsave(filename = paste0("results/figs/controls/parameterChanges_fct.tiff"),
       plot = g + theme(legend.position = "none"),
       width = 16.2, height = 15, unit = "cm")
ggsave(filename = "results/figs/controls/parameterChanges_lgd.tiff",
       plot = g,
       width = 16.2, height = 15, unit = "cm")

