## Using packaged code for convenience
library(phosphoRWHN)
library(dplyr)
library(tidyr)
library(Mfuzz)
library(igraph)
library(ggplot2)
library(reticulate)
source("src/functions/imputePhosphoData.R")
source("src/functions/mfuzz-ggplot.R")

sp <- reticulate::import("scipy.stats")


##############
# Process data
##############
# Import
# Import data
cfr_sty <- data.table::fread(input = "data/CFR_STY_2016.csv",
                             select = c(1:9, 26:35)
) %>% 
  mutate(id = paste0(`Gene names`, "_", `Swiss-Prot phosphosite`)) %>% 
  filter(!duplicated(id))
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

set.seed(123)
egf_fcm <- fuzzyC(egf, 6)

# Run RWHN with heterogeneous etworks constructed
# with different STRING / PPI confidence scores
stringTests_hetNet <- lapply(seq(0.1, 0.9, 0.1), function(x){
  
  hetNet <- constructHetNet(clustering = egf_fcm$clustering,
                            phosphoData = egf,
                            stringConf = x)
})

seed_l <- lapply(1:max(egf_fcm$clustering), function(i){
  c(names(egf_fcm$clustering[egf_fcm$clustering == i]))
})

stringTests_RWHN <- lapply(stringTests_hetNet, function(hetNet){
  rwhn <- lapply(seed_l, function(s){
    calculateRWHN(hetNet = hetNet,
                  seeds = s, 
                  transitionProb = 0.7, 
                  restart = 0.7,
                  eta_xy =  0.3,
                  eta_yz = 0.7,
                  eps = 1/10^12
    )
  })
  return(rwhn)
})
#Save
# saveRDS(stringTests_RWHN, "results/data/STRINGtests_EGF.rds")

##################################
# Correlation between RWHN outputs
# when STRING confidence changed
##################################

#stringTests_RWHN <- readRDS("results/data/STRINGtests_EGF.rds")
dat <- expand.grid(
  data.frame(res1 = 1:length(stringTests_RWHN),
             res2 = 1:length(stringTests_RWHN)))

wt <- lapply(1:max(egf_fcm$clustering), function(clus){
  apply(dat, 1, function(res){
    
    res1 <- res["res1"]
    res2 <- res["res2"]
    
    x <- as.numeric(as.factor(stringTests_RWHN[[res1]][[clus]]$name))
    y <- as.numeric(as.factor(stringTests_RWHN[[res2]][[clus]]$name))
    wtr <- sp$weightedtau(x, y)
    
    fin <- data.frame(
      res1 = res1,
      res2 = res2,
      cl = paste("Cluster", clus),
      wt = wtr$correlation
    )
    return(fin)
  }) %>% 
    bind_rows()
}) %>% 
  bind_rows() %>% 
  ggplot(aes(x = as.character(res1), y = as.character(res2), fill = wt)) +
  geom_tile(color = "grey") +
  geom_text(aes(label = round(wt, digits = 2)), angle = 45, size = 1.5) +
  facet_wrap(~cl) +
  scale_x_discrete(labels = seq(0.1, 0.9, 0.1)) +
  scale_y_discrete(labels = seq(0.1, 0.9, 0.1)) + 
  scale_fill_gradient("Weighted Tau\nCorrelation", 
                    high = "#9cd64f", low = "#e2f311",
                    #limits = c(0, 1)
                    ) +
  xlab("STRING experimental confidence cut-off") +
  ylab("STRING experimental confidence cut-off") +
  theme(plot.background = element_blank(),
        panel.grid = element_blank(),
        legend.key.size = unit(.25, "cm"),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        axis.text.x = element_text(size = 5, angle = 45),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        legend.margin = margin(0,0,0,0, "cm"),
        strip.text = element_text(size = 5)
  ) 


# ggsave(
#   "results/figs/controls/STRINGconfidence_cor_legend.pdf",
#   wt ,
#   width = 4,
#   height = 4,
#   dpi = 300
# )
##################################
# Inter-connectedness of protein layer
# given different STRING confidence thresholds
##################################


#PPIs <- 
calcs <- data.frame(
  stringConf = seq(0.1, 0.9, 0.1),
  Components = sapply(stringTests_hetNet, function(x){
    nw <- components(graph_from_data_frame(x$edgelists$y, directed = F))$no
  }),
  # The adhesion of a graph is the minimum number of edges needed to remove to 
  # obtain a graph which is not strongly connected. 
  # This is the same as the edge connectivity of the graph.
  GraphAdhesion = sapply(stringTests_hetNet, function(x){
    nw <- edge_connectivity(graph_from_data_frame(x$edgelists$y, directed = F))
  }),
  # Transitivity measures the probability that the adjacent vertices of a vertex are connected.
  # This is sometimes also called the clustering coefficient.
  # https://en.wikipedia.org/wiki/Clustering_coefficient
  AvgLocalClusteringCoefficient = sapply(stringTests_hetNet, function(x){
    nw <-  transitivity(graph_from_data_frame(x$edgelists$y, directed = F), 
                        type = "global")
  }),
  # The vertex connectivity of a graph is the minimum vertex connectivity of all 
  # (ordered) pairs of vertices in the graph. In other words this is the minimum 
  # number of vertices needed to remove to make the graph not strongly connected. 
  # (If the graph is not strongly connected then this is zero.)
  GraphCohesion = sapply(stringTests_hetNet, function(x){
    nw <-  vertex_connectivity(graph_from_data_frame(x$edgelists$y, directed = F))
  })
) #%>% 
 # tidyr::pivot_longer(cols = -stringConf)
tidyr::pivot_longer(calcs, cols = -stringConf) %>% 
  ggplot(aes(x = as.factor(stringConf), y = value, group = name)) +
  geom_point() +
  geom_line() +
  xlab("STRING confidence score threshold") +
  ylab("n") +
  theme_minimal() +
  facet_wrap( ~ name, scales = "free_y")

# data.frame(
#   stringConf = seq(0.1, 0.9, 0.1),
lapply(1:length(stringTests_hetNet), function(x){
  
  nw <- graph_from_data_frame(
    stringTests_hetNet[[x]]$edgelists$y,
    directed = F)
  
  data.frame(
    stringConf  = seq(0.1, 0.9, 0.1)[x],
    components = components(nw)$csize,
    total = length(V(nw))
  )
  }) %>% 
  bind_rows() %>% 
  tidyr::pivot_longer(cols = -stringConf) %>% 
  ggplot(aes(x = as.factor(stringConf), y = value, group = name)) +
  geom_point() +
  facet_wrap( ~ name, scales = "free_y") +
  ylab("log10(number vertices)")



randos <- lapply(1:length(stringTests_hetNet), function(i){
  orig <-  graph_from_data_frame(
    stringTests_hetNet[[i]]$edgelists$y,
    directed = F)
  
  randos <- lapply(1:100, function(x){
    erdos.renyi.game(n = length(V(orig)),
                     p.or.m = length(E(orig)),
                     type = "gnm"
    )
  })
  
  res <- data.frame(
    stringConf = seq(0.1, 0.9, 0.1)[i],
    Components = sapply(randos, function(x){
      components(x)$no
    }),
    AvgLocalClusteringCoefficient =  sapply(randos, function(x){
      transitivity(x, type = "localaverage")
    })
  )
  return(res)
}) %>% 
  bind_rows() %>% 
  mutate(random = T)

dat <- calcs %>% 
  dplyr::select(stringConf, Components, AvgLocalClusteringCoefficient) %>% 
  mutate(random = F) %>% 
  rbind(randos)


comps <- ggplot(dat[dat$random == T,], aes(x = Components)) +
  geom_histogram() +
  geom_vline(data = dat[dat$random == F,], 
             mapping = aes(xintercept = Components), color = "red",
             size = 0.25) +
  xlab("No. components") +
  facet_wrap(~ stringConf, ncol = 2) +
  scale_x_continuous(breaks = c(1:8)) +
  theme(plot.background = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size = 5),
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        legend.margin = margin(0,0,0,0, "cm"),
        strip.text = element_text(size = 5)
  ) 


# ggsave(
#   "results/figs/controls/STRINGconfidence_components.pdf",
#   comps,
#   width = 2.8,
#   height = 7.7,
#   dpi = 300
# )
