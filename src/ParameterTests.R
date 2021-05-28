library(dplyr)
library(tidyr)
library(igraph)
library(ggplot2)
library(reticulate)

sp <- reticulate::import("scipy.stats")
source("src/functions/calculateRWHN.R")

# Import multilayer network
mlnw <- readRDS("results/data/mlnw_model.rds")



# Parameter ranges
parameters <- data.frame(trans = c(seq(0.1, 0.9, 0.1), rep(0.5, 27)),
                         restart = c(rep(0.5, 9), seq(0.1, 0.9, 0.1), rep(0.5, 18)),
                         eta_xy = c( rep(0.5, 18), seq(0.1, 0.9, 0.1), rep(0.5, 9)),
                         eta_yz = c(rep(0.5, 27), seq(0.1, 0.9, 0.1)), 
                         param_altered = c(rep("Transition Probability", 9),
                                           rep("Restart Probability", 9),
                                           rep("Protein Weighting", 9),
                                           rep("Function Weighting", 9)),
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

wt <- lapply(unique(parameters$param_altered), function(i){
  
  index <- which(grepl(i, parameters$param_altered))
  col <- ifelse(i == "Transition Probability", "trans", 
                ifelse(i == "Restart Probability", "restart",
                       ifelse(i == "Protein Weighting", "eta_xy",
                              ifelse(i == "Function Weighting", "eta_yz", NA))))


  dat <- expand.grid(res1 = index,
                     res2 = index)
  
  wt <- lapply(1:5, function(cl){
    apply(dat, 1, function(res){
    
    res1 <- res["res1"]
    res2 <- res["res2"]
    
    x <- as.numeric(as.factor(rwhn[[res1]][[cl]]$name))
    y <- as.numeric(as.factor(rwhn[[res2]][[cl]]$name))
    wtr <- sp$weightedtau(x, y)
    
    fin <- data.frame(
      res1 = parameters[res1,col],
      res2 = parameters[res2,col],
      wt = wtr$correlation,
      cl = paste("Cluster", cl)
    )
    return(fin)
  })
})   %>% 
    bind_rows() %>% 
    mutate(param = i)
}) %>%
  bind_rows() %>% 
  ggplot(aes(x = as.character(res1), y = as.character(res2), fill = wt)) +
  geom_tile() +
  facet_grid(cl ~ param) +
  xlab("Parameter value") +
  ylab("Parameter value") +
  # scale_x_discrete(labels = seq(0.1, 0.9, 0.1)) +
  # scale_y_discrete(labels = seq(0.1, 0.9, 0.1)) + 
  scale_fill_viridis_c("Weighted\nTau\nCorrelation")+#, limits = c(0.5, 1)) +
  theme(plot.background = element_blank(),
        panel.grid = element_blank(),
        legend.key.size = unit(.2, "cm"),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        axis.text.x = element_text(size = 5, angle = 45),
        axis.title = element_blank(),
        legend.margin = margin(0,0,0,0, "cm"),
        strip.text = element_text(size = 5),
        plot.margin = margin(0,0,0,0,"mm")
  )

wt
ggsave("results/figs/controls/parameterChanges_wt_fct.pdf",
       wt,
       height = 4, width = 4, dpi = 300)


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

