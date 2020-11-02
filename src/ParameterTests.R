library(dplyr)
library(tidyr)
library(igraph)
library(ggplot2)
source("src/functions/calculateRWHN.R")

# Import multilayer network
mlnw <- readRDS("results/data/mlnw_model.rds")



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

# Calculate RWHN with toy network, each cluster set to seed
# and each set of parameters
rwhn <- apply(parameters[,1:4], 1, function(i){
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
  facet_grid(seed ~ param_altered)  +
  geom_tile(aes(fill = n)) +
  scale_fill_gradient("% times\nfound\nat rank", 
                      low = "red", high = "blue", 
                      na.value = "grey") +
  ylab("GO Term") +
  xlab("RWHN Rank") +
  ggridges::theme_ridges() +
  theme(legend.text = element_text(size = 8),
        plot.title = element_text(size = 8),
        axis.text.y = element_blank(), 
        axis.text.x = element_text(size = 8),
        axis.title = element_text(size = 8), 
        axis.title.x = element_text(hjust = 0.5),
        axis.title.y = element_text(hjust = 0.5),
        strip.text = element_text(size = 8)
  )
g

ggsave(filename = paste0("results/figs/controls/parameterChanges_fct.tiff"),
       plot = g,
       width = 18.2, height = 15, unit = "cm", dpi = "print")
ggsave(filename = paste0("results/figs/controls/parameterChanges_fct.pdf"),
       plot = g,
       width = 18.2, height = 15, unit = "cm", dpi = "print")

