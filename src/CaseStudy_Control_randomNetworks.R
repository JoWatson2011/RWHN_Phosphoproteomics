library(phosphoRWHN)
library(patchwork)
library(tidyverse)

mlnw_model <- readRDS("results/data/mlnw_model.rds")

fcm <- mlnw_model$fcm
v <- mlnw_model$v
edgelists <- mlnw_model$edgelists

seed_l <- lapply(1:max(fcm$cluster), function(i){
  c(names(fcm$clustering[fcm$clustering == i]))
})

rwhn_random <- lapply(seed_l, function(s){
  list <- lapply(1:100, function(i){
    ranks <- calculateRWHN(hetNet = mlnw_model,
                           seeds = s,
                           transitionProb = 0.7,
                           restart = 0.7,
                           eta_xy = 0.3,
                           eta_yz = 0.7,
                           random = T,
                           filterFunctions = F)
    return(ranks)
  })
  return(list)
})

saveRDS(rwhn_random, "results/data/rwhn_model_randomPermutation.rds")

# Import
#rwhn_random <- readRDS("results/data/rwhn_model_randomPermutation.rds")
# Import 'actual' result
rwhn <- readRDS("results/data/rwhn_model.rds")

sighm <- heatmap_RWHN(rwhn,#, ylab = "GOBP Term",
                      colours = c(low = "#ffe6e8", high = "#f74451"),
                      removeCommon = T)

# Wrangle data to long formats and combine for density plot visualisation
rwhn_random_gg <- lapply(1:length(rwhn), function(x){
  df <- lapply(1:length(rwhn_random[[x]]), function(i){
    df1 <- rwhn_random[[x]][[i]] %>%
      filter(name %in% v$v[v$layer=="func"]) %>% 
      mutate(iter = i,
             rank = 1:nrow(.),
             actual = F)
  }) %>% do.call(rbind, .) %>%
    mutate(seed = x)
})
rwhn_gg <- lapply(1:length(rwhn), function(x){
  df <- rwhn[[x]] %>% 
    filter(name %in% v$v[v$layer=="func"]) %>% 
    mutate(iter = 0,
           rank = 1:nrow(.),
           actual = T,
           seed = x)
})
gg <- lapply(1:length(rwhn_gg), function(x) 
  rbind(rwhn_gg[[x]], rwhn_random_gg[[x]]) %>% 
    mutate(name = factor(
      x = name, 
      levels = rwhn_gg[[x]]$name, 
      order = T))
  )

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
    scale_fill_viridis_c(
      
      name = "KDE of random\nnetworks",
      limits = c(0, max(sapply(gg, function(i){ max(MASS::kde2d(as.numeric(i$name), i$rank)$z)}))),
      na.value = "grey"
    ) +
    theme(axis.title = element_text(size = 5),
        legend.text = element_text(size = 5),
         legend.title= element_text(size = 5),
         axis.text.x = element_text(size = 5),
         axis.text.y = element_blank(),
         axis.ticks.y = element_blank(),
         axis.line = element_line(colour = "black"),
         legend.key.height = unit(0.25, "cm"),
         legend.key.width = unit(0.25, "cm"),
         title = element_text(size = 5)
    ) +
    ggtitle(paste("Cluster", i$seed))
})

collate <- gs[[1]] + 
  gs[[2]] + 
  gs[[3]] +
  gs[[4]] +
  gs[[5]] +
  plot_layout(nrow = 1, guides = "collect")


ggsave("results/figs/controls/model_random_nw_patchwork.pdf",
  collate,
  width = 7,
  height = 2,
  units = "in",
  dpi = "print"
)
ggsave("results/figs/controls/model_random_nw_patchwork.tiff",
       collate,
       width = 7,
       height = 2,
       units = "in",
       dpi = "print"
)

#####
# Calculate divergence from random
# of true ranks
#####


sighm_1 <- heatmap_RWHN(rwhn,#, ylab = "GOBP Term",
                        colours = c(low = "#ffe6e8", high = "#f74451"),
                        removeCommon = T, pct_cutoff = 0.01) %>% 
  .$data %>% 
  mutate(pct_cutoff = 0.01)
sighm_5 <- heatmap_RWHN(rwhn,#, ylab = "GOBP Term",
                      colours = c(low = "#ffe6e8", high = "#f74451"),
                      removeCommon = T, pct_cutoff = 0.05) %>% 
  .$data %>% 
  mutate(pct_cutoff = 0.05)
sighm_10 <- heatmap_RWHN(rwhn,#, ylab = "GOBP Term",
                        colours = c(low = "#ffe6e8", high = "#f74451"),
                        removeCommon = T, pct_cutoff = 0.1) %>% 
  .$data %>% 
  mutate(pct_cutoff = 0.1)
sighm_15 <- heatmap_RWHN(rwhn,#, ylab = "GOBP Term",
                         colours = c(low = "#ffe6e8", high = "#f74451"),
                         removeCommon = T, pct_cutoff = 0.15) %>% 
  .$data %>% 
  mutate(pct_cutoff = 0.15)


ptests <- 
  lapply(1:length(seed_l),function(i){
    lapply(v$v[v$layer=="func"], function(term){
      
      ranks <-  sapply(rwhn_random[[i]],
                       function(x) ifelse(term %in% x$name, x$rank[x$name == term], NA))
      
      if(all(is.na(ranks))){
        return(data.frame())
      }
      
      true_rank <- rwhn[[i]]$rank[rwhn[[i]]$name == term]
      t_p <- t.test(ranks, mu = true_rank)$p.value
      t_padj <- p.adjust(t_p, "fdr")
      mann_p <- wilcox.test( ranks, mu = true_rank)$p.value
      mann_padj <- p.adjust(mann_p, "fdr")
  
      res <- data.frame(
        term,
        true_rank,
        true_prob =  rwhn[[i]]$V1[rwhn[[i]]$name == term],
        t_p,
        t_padj,
        mann_p,
        mann_padj
      )
      return(res)
    }) %>% 
      bind_rows() %>% 
      mutate(sig_t = ifelse(t_padj < 0.05, T, F),
             sig_mann = ifelse(mann_padj < 0.05, T, F),
             insighm_1 = ifelse(term %in% sighm_1[sighm_1$seed == i,]$name, T, F),
             insighm_5 = ifelse(term %in% sighm_5[sighm_5$seed == i,]$name, T, F),
             insighm_10 = ifelse(term %in% sighm_10[sighm_10$seed ==i,]$name, T, F),
             insighm_15 = ifelse(term %in% sighm_15[sighm_15$seed == i,]$name, T, F),
             cluster = i) %>% 
      arrange(true_rank)
  }) %>% 
  bind_rows()

limits <- c(
  "sig_mann",
  "insighm_1",
  "insighm_5",
  "insighm_10",
  "insighm_15"
)

names <- c(
  "Mann Whitney-U\n< 0.05",
  "1%",
  "5%",
  "10%",
  "15%"
)
  
mannwhitney_hm <- ptests %>% 
  pivot_longer(c(sig_mann, insighm_1, insighm_5,
                 insighm_10, insighm_15), 
               names_to = "test", values_to = "p") %>% 
  ggplot(aes(x = factor(test, levels = limits), y = term, fill = p)) +
  geom_tile() + 
  scale_x_discrete(name = "", 
                   labels = names
  ) +
  xlab("GOBP Term") +
  scale_fill_manual(name = "", values = c("TRUE" = "#1ecbe1", "FALSE" = "#E1341E")) +
  facet_wrap(~cluster, nrow = 1) +
  theme(panel.background = element_blank(),
        plot.background = element_blank(),
        axis.title = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        axis.text.x = element_text(size = 5, angle = 45, hjust = 1),
        legend.key.height = unit(0.25, "cm"),
        legend.key.width = unit(0.25, "cm"),
        legend.text = element_text(size = 5),
        plot.margin = margin(0,0,0,0, "cm"),
        legend.margin = margin(0,0,0,0,"cm")
  ) 
ggsave("results/figs/controls/mannWhitney.pdf",
       mannwhitney_hm,
       width = 7,
       height = 5,
       units = "in",
       dpi = "print")

