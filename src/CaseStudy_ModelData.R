library(phosphoRWHN)
library(Mfuzz)
library(tidyverse)
library(GO.db)

model <- readr::read_csv("data/model_phospho.csv") %>% 
  mutate(id = paste0(prot, "_", site)) %>% 
  dplyr::select(id, t1:t5)

set.seed(1)
inputdata <- model %>% 
  tibble::remove_rownames() %>%
  tibble::column_to_rownames(var = "id") 
phospep <- new("ExpressionSet", 
               exprs = as.matrix(inputdata))
phospep.z <- standardise(phospep)      #Avg expression of each peptide = 0, sd = 1
optimalM <- mestimate(phospep.z)   #set fuzzifier
cl <- mfuzz(phospep.z, c = 5, m = optimalM)


# hetNet <- constructHetNet(clustering = cl$cluster, modules = F)
# saveRDS(list(edgelists = hetNet$edgelists,
#              v = hetNet$v,
#              fcm = cl), 
#         "results/data/mlnw_model.rds")

# 2. Run RWHN, using sites in each cluster as seeds
hetNet <- readRDS("results/data/mlnw_model.rds")
seed_l <- lapply(1:max(cl$cluster), function(i){
  c(names(cl$cluster[cl$cluster == i]))
})
set.seed(1)
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

saveRDS(rwhn, "results/data/rwhn_model.rds")
#######
# Export, Visualise, etc.
#######
# rwhn <- readRDS("results/data/rwhn_model.rds")
# 
# rwhn_df <- lapply(1:length(rwhn), function(i)
#   mutate(rwhn[[i]],
#          seed = i,
#          rank = 1:nrow(rwhn[[i]])
#          )
#   ) %>%
#   do.call(cbind, .)
# readr::write_csv(rwhn_df, "results/data/model_rwhn_results.csv")

sighm <-
  heatmap_RWHN(rwhn,
               database = "GOBP",
               colours = c(low = "#ffe6e8", high = "#f74451"))



######################
# Standard GO analysis
######################

enrichedTerms <-overrepresentationAnalysis(clustering = cl$cluster, vis = T)

#####################
# 
# ggsave("results/figs/model_RWHN_ORA.tiff",
#   sighm + enrichedTerms +
#     plot_layout(widths = c(2, 1)),
#   width = 180,
#   height = 80,
#   units = "mm", 
#   dpi = "print" 
# )
####
# Visualisation of results in GO tree form
# Added 30APR21
####

library(GO.db)
library(ggraph)
library(igraph)

RWHN_IDS <- AnnotationDbi::select(GO.db,
                                  unique(sighm$data$name),
                                  "GOID",
                                  "TERM") %>% 
  na.omit() %>% 
  merge(sighm$data, by.x= "TERM", by.y = "name") %>% 
  dplyr::select(TERM, GOID, seed)

GO_RWHN_ancestors <- as.list(GOBPANCESTOR)[unique(RWHN_IDS$GOID)]
# Is root node (all) in each?
#map(GO_RWHN_ancestors, function(i) "all" %in% i)

GO_RWHN_dirparents <- as.list(GOBPPARENTS)[unique(RWHN_IDS$GOID)] 
GO_RWHN_dirparents <- map(names(GO_RWHN_dirparents), function(i)
  data.frame(parent = GO_RWHN_dirparents[[i]], child = i)
) %>% 
  bind_rows()

GO_RWHN_allparents <- map(GO_RWHN_ancestors, function(i){
  listParents <- as.list(GOBPPARENTS)[i]
  
  dfParents <- map(na.omit(names(listParents)), function(x)
    data.frame(parent = listParents[[x]], child = x)
  ) %>% 
    bind_rows()
}) %>% 
  bind_rows()


GO_RWHN_comb_el <- rbind(GO_RWHN_allparents,
                         GO_RWHN_dirparents)
GO_RWHN_comb_v <- data.frame(name = unique(c(
  GO_RWHN_comb_el$parent, GO_RWHN_comb_el$child
))) %>%
  left_join(RWHN_IDS, by = c("name" = "GOID")) %>%
  group_by(name) %>%
  summarise(seed = paste(seed, collapse = ","),
            TERM = paste(unique(TERM), collapse = ",")) %>%
  mutate(
    seed = ifelse(seed == "NA", NA, seed),
    TERM = ifelse(TERM == "NA", NA, TERM),
    inRWHN = ifelse(seed == "", F, T),
    isAll = ifelse(name == "all", name, NA)
  )


## ORA
ORA_ID <- dplyr::select(enrichedTerms$data, TERM = Term, GOID, cluster)

GO_ORA_ancestors <- as.list(GOBPANCESTOR)[unique(ORA_ID$GOID)]
# Is root node (all) in each?
#map(GO_ORA_ancestors, function(i) "all" %in% i)

GO_ORA_dirparents <- as.list(GOBPPARENTS)[unique(ORA_ID$GOID)] 
GO_ORA_dirparents <- map(names(GO_ORA_dirparents), function(i)
  data.frame(parent = GO_ORA_dirparents[[i]], child = i)
) %>% 
  bind_rows()

GO_ORA_allparents <- map(GO_ORA_ancestors, function(i){
  listParents <- as.list(GOBPPARENTS)[i]
  
  dfParents <- map(na.omit(names(listParents)), function(x)
    data.frame(parent = listParents[[x]], child = x)
  ) %>% 
    bind_rows()
}) %>% 
  bind_rows()


GO_ORA_comb_el <- rbind(GO_ORA_allparents,
                        GO_ORA_dirparents)
GO_ORA_comb_v <- data.frame(name = unique(c(
  GO_ORA_comb_el$parent, GO_ORA_comb_el$child
))) %>%
  left_join(ORA_ID, by = c("name" = "GOID")) %>%
  group_by(name) %>%
  summarise(cluster = paste(cluster, collapse = ","),
            TERM = paste(unique(TERM), collapse = ",")) %>%
  mutate(
    cluster = ifelse(cluster == "NA", NA, cluster),
    TERM = ifelse(TERM == "NA", NA, TERM),
    inORA = ifelse(cluster == "", F, T),
    label = ifelse(inORA == T, TERM, NA)
  )


# Combine
GO_comb_el <- rbind(
  GO_ORA_comb_el,
  GO_RWHN_comb_el 
) %>% distinct()
GO_comb_v <- rbind(
  (dplyr::select(GO_RWHN_comb_v, name, inRes = inRWHN, cluster = seed) %>% 
     mutate(inRes = ifelse(inRes == T, "RWHN", NA))
  ),
  (dplyr::select(GO_ORA_comb_v, name, inRes = inORA, cluster) %>% 
     mutate(inRes = ifelse(inRes == T, "ORA", NA))
  )
) %>% 
  group_by(name) %>% 
  summarise(inRes = paste(na.omit(inRes), collapse = "&")) %>%  
  mutate(inRes = ifelse(inRes == "", "NA", inRes),
         size = ifelse(inRes != "NA", 2, 1))

nw <- simplify(
  graph_from_data_frame(GO_comb_el, vertices = GO_comb_v)
)

gg <- ggraph(nw) +
  geom_edge_bend(alpha = 0.25, color = "grey", edge_width = 0.2) +
  geom_node_point(aes(color = inRes, size = inRes), alpha = 0.7) +
  scale_color_manual(values = c("ORA" = "#845bd6", 
                                "RWHN&ORA" = "#5bc7d6",
                                "RWHN" = "#5bd670",
                                "NA" = "grey"
  )) +
  scale_size_manual(values = c("ORA" = 3, 
                               "RWHN&ORA" = 3,
                               "RWHN" = 3,
                               "NA" = 2
  )) +
  theme(legend.key.size = unit(.25, "cm"),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.position = "top",
        panel.background = element_blank(),
        plot.margin= margin(0,0,0,0, "pt"),
        plot.background = element_blank())

for(i in unique(gg$data$y)[unique(gg$data$y) != min(unique(gg$data$y))]){
  gg <-  gg +
    geom_hline(yintercept = i - 0.5, alpha = 0.2, colour = "#cba3cd") 
}

gg <- gg +
  geom_text(data = data.frame(
    x = -10, 
    y = unique(gg$data$y)[order(unique(gg$data$y))][1:11],
    label = unique(gg$data$y)[order(unique(gg$data$y))][1:11]
  ), aes(x = x, y = y, label = label), colour = "black"
  )


ora_with_nw <- left_join(enrichedTerms$data,
                         gg$data[,c("y", "name")],
                         by = c("GOID" = "name")) %>% 
  dplyr::select(term = Term, GOID , value = Adjusted.P.value, cluster, y) %>% 
  mutate(method = "ORA")
rwhn_with_nw <- left_join(sighm$data, 
                          AnnotationDbi::select(GO.db,
                                                sighm$data$name,
                                                "GOID",
                                                "TERM"),
                          by = c("name" = "TERM")) %>%
  left_join(gg$data[,c("y", "name")], by = c("GOID"= "name")) %>% 
  dplyr::select(term = name, GOID, value = rank, cluster = seed, y) %>% 
  mutate(method = "RWHN")


rwhn_with_nw_gg <- rwhn_with_nw %>% 
  filter(!is.na(cluster)) %>%
  rbind(ora_with_nw) %>% 
  mutate(cluster = ifelse(method == "RWHN", cluster, NA),
         value = ifelse(method == "RWHN", value, NA)
  ) %>% 
  filter(!is.na(y)) %>% 
  arrange(desc(y)) %>% 
  mutate(y = factor(y, levels = as.character(c(10:1)))) %>% 
  ggplot(aes(y = term,
             x = as.factor(cluster)
  )
  ) +
  geom_tile(aes(fill = value), color = "white") +
  scale_fill_gradient(low = "#5bd670", high = "#def6e2") +
  guides(color = FALSE,
         fill = guide_colourbar(title="Rank", reverse = T)) +
  theme_minimal() +
  theme(panel.grid = element_line(size = 0.5),
        legend.key.size = unit(.25, "cm"),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5, angle = 45),
        legend.position = "top",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 5),
        legend.margin = margin(0,0,0,0, "cm"),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        panel.spacing = unit(0.2, "mm"), 
        panel.border = element_rect(color = "grey", fill = "transparent")
  ) +
  scale_x_discrete(breaks = factor(c(1:5)), limits = factor(c(1:5)), 
                   position = "top") +
  xlab("Cluster") +
  ylab("GOBP") +
  facet_grid(y ~., scales = "free", space = "free")


ora_with_nw_gg <-  ora_with_nw %>% 
  filter(!is.na(cluster) | !is.na(term)) %>%
  rbind(rwhn_with_nw) %>% 
  mutate(cluster = ifelse(method == "ORA", cluster, NA),
         value = ifelse(method == "ORA", value, NA)
  ) %>% 
  filter(!is.na(y)) %>% 
  arrange(desc(y)) %>% 
  mutate(y = factor(y, levels = as.character(c(10:1)))) %>% 
  ggplot(aes(y = term,
             x = as.factor(cluster)
  )
  ) +
  geom_tile(aes(fill = value), color = "white") +
  scale_fill_gradient(breaks = seq(0, 0.05, 0.01), limits = c(0, 0.05),
                      low = "#845bd6", high = "#e6def6") +
  guides(color = FALSE,
         fill = guide_colourbar(title="FDR")) +
  theme_minimal() +
  theme(panel.grid = element_line(size = 0.5),
    legend.key.size = unit(.25, "cm"),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5, angle = 45, hjust = 1),
        legend.position = "top",
        axis.text.x = element_text(size = 5),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 5),
        #      panel.background = element_rect(fill = "white"),
        legend.margin = margin(0,0,0,0, "cm"),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        panel.spacing = unit(0.2, "mm"), 
        panel.border = element_rect(color = "grey", fill = "transparent"),
        plot.margin = margin(0,0,0,0,"mm")
        
  ) +
  scale_x_discrete(breaks = factor(c(1:5)), limits = factor(c(1:5)), 
                   position = "top") +
  xlab("Cluster") +
  ylab("GOBP") +
  facet_grid(as.factor(y) ~., scales = "free", space = "free") 

library(patchwork)
pw <- rwhn_with_nw_gg + 
  ora_with_nw_gg +
  gg +
  plot_layout(
    widths = c(1,1,3)
  )
pw

ggsave("results/figs/Model_with_nw.pdf",
       pw, dpi = 300, width = 7, height =5.8, units = "in")
