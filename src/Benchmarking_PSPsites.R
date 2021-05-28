# From CRAN
library(dplyr)
library(tidyr)
library(igraph)
library(enrichR)
library(patchwork)
library(ggplot2)
library(ggraph)
# From Bioconductor
library(GOSemSim)
library(AnnotationDbi)
library(GO.db)
# From src/
source("src/functions/overrepresentationAnalysis.R")
source("src/functions/heatmap_RWHNsig.R")
source("src/functions/simplifyGO.R")
source("src/functions/simplifyGOReqData.R")
source("src/functions/constructHetNet.R")
source("src/functions/calculateRWHN.R")

set.seed(123)

data <- readr::read_tsv("data/Regulatory_sites.txt", skip = 2) %>% 
  filter(!is.na(ON_PROCESS),
         ORGANISM == "human",
         grepl("-p", MOD_RSD)) %>% 
  separate_rows(ON_PROCESS, sep = "; ") %>% 
  mutate(id_site = paste0(GENE,
                          "_",
                          gsub("-.*",
                               "",
                               MOD_RSD))) %>% 
  unique() 

clustering <- data$ON_PROCESS
names(clustering) <- data$id_site

hetNet <- constructHetNet(clustering = clustering)

saveRDS(hetNet, "results/data/PSP_benchmark_all_hetNet.rds")

seed_l <- lapply(unique(clustering)[order(unique(clustering))], function(x){
  c(names(clustering[clustering == x]))
})
rwhn <- lapply(seed_l, function(s){
  calculateRWHN(edgelists = hetNet$edgelists,
                verti = hetNet$v,
                seeds = s,
                transitionProb = 0.7,
                restart = 0.7,
                weight_xy = 0.3,
                weight_yz = 0.7,
                eps = 1/10^12) %>%
    filter(name %in% hetNet$v[hetNet$v$layer=="func",]$v)
})
saveRDS(rwhn, "results/data/rwhn_all_benchmarks.rds")

#####
# Assessment
#####
rwhn <- readRDS("results/data/rwhn_all_benchmarks.rds")
sighm <- heatmap_RWHN(rwhn, ylab = "GOBP")

###SAVED BELOW
# ora <- overrepresentationAnalysis(clustering = clustering,
#                                   RWHN_sig = sighm,
#                                   colours = c("#effff6","#168d49"))
# saveRDS(ora, "data/PSP_allsites_ORA.rds")

ora <- readRDS("data/PSP_allsites_ORA.rds")


psp_mapped <- readr::read_tsv("data/Regulatory_sites_GOmapped.tsv") %>%
  dplyr::select(ON_PROCESS, GOID) %>% 
  unique() 

####################
# Visualize as tree
####################
library(GO.db)
library(ggraph)
library(igraph)

psp_tree_sep <- lapply(unique(psp_mapped$GOID), function(i){
  psp_offspring <- as.list(GOBPOFFSPRING)[[i]]
  
  dirchild <- data.frame(parent = i, 
                         child = GOBPCHILDREN[[i]])
  children <- lapply(psp_offspring, function(offspr){
    data.frame(parent = offspr,
               child = GOBPCHILDREN[[offspr]])
  }) %>% 
    bind_rows() %>% 
    rbind(dirchild)
  
  return(children)
})

psp_tree <- na.omit(bind_rows(psp_tree_sep))

RWHN_IDs <- AnnotationDbi::select(GO.db,
                                  unique(sighm$data$name),
                                  "GOID",
                                  "TERM") %>% 
  na.omit() %>% 
  merge(sighm$data, by.x= "TERM", by.y = "name") %>% 
  dplyr::select(TERM, GOID) %>% 
  unique()
ORA_IDS <- dplyr::select(ora$data, TERM = Term, GOID)

psp_tree_v <- data.frame(name = 
                           unique(
                             c(
                               psp_tree$parent,
                               psp_tree$child
                               )
                             )
) %>% 
  na.omit() %>% 
  mutate(inRes = ifelse(
    name %in% RWHN_IDs$GOID &
      name %in% ORA_IDS$GOID, "RWHN&ORA",
    ifelse(name %in% RWHN_IDs$GOID, "RWHN",
           ifelse(name %in% ORA_IDS$GOID, "ORA",
                  "")
    )
  )
  ) 

nw <- simplify(
  graph_from_data_frame(psp_tree, vertices = psp_tree_v)
)

gg <- ggraph(nw)

levels <- gg$data %>%
  group_by(y, inRes) %>% 
  summarise(n = n()) %>% 
  filter(inRes != "") %>%
  pivot_wider(names_from = inRes, 
              values_from = n) 
levels$y_reorder <- levels$y[order(levels$y, decreasing = T)]
levels_gg <- levels %>%
  pivot_longer(cols = -c(y, y_reorder),
               names_to = "inRes", 
               values_to = "n") %>% 
  ggplot(aes(x = as.factor(y_reorder), y = n, fill = inRes)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_y_continuous(name = "Number of GO terms",
                     n.breaks = 10) +
  xlab("Hierachy level (distance from top level term)") +
  #ylab("Number of GO terms") +
  scale_fill_viridis_d(name = "Results") +
  theme(#plot.background = element_blank(), 
      #  panel.background = element_blank(),
      #  panel.grid.major.x = element_blank(),
       # panel.grid.major.y = element_line(color = "grey"),
    legend.justification = c(1, 1),
    legend.position = c(1, 1),
    legend.key.size = unit(.25, "cm"),
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 5), 
    axis.text.y = element_text(size = 5),
    legend.margin = margin(0,0,0,0, "cm"),
    axis.text = element_text(size = 5),
    axis.title = element_text(size = 5)
  )
ggsave("results/figs/controls/PSP_Hierachy_distribution.pdf",
       levels_gg,
       width = 3.33,
       height = 2.5,
       dpi = 300
)
ggsave("results/figs/controls/PSP_Hierachy_distribution.tiff",
       levels_gg,
       width = 3.33,
       height = 2.5,
       dpi = 300
)

### Distribution of terms
psp_mapped <-
  readr::read_tsv("data/Regulatory_sites_GOmapped.tsv") %>%
  separate_rows(offspring, sep = ";") %>%
  pivot_longer(
    cols = -c(parentProc, ON_PROCESS),
    names_to = "parent",
    values_to = "GOID"
  ) %>%
  #dplyr::select(-name) %>%
  mutate(parent = ifelse(parent == "GOID", "parent", "offspring")) %>%
  unique()

psp_mapped$TERM <- AnnotationDbi::select(GO.db,
                                         keys = psp_mapped$GOID,
                                         columns = "TERM",
                                         
                                         keytype = "GOID")$TERM
# Import PSP site > annotation
psp_sites <- readr::read_tsv("data/Regulatory_sites.txt", skip = 2) %>% 
  filter(!is.na(ON_PROCESS),
         ORGANISM == "human",
         grepl("-p", MOD_RSD)) %>% 
  separate_rows(ON_PROCESS, sep = "; ") %>% 
  mutate(id_site = paste0(GENE,
                          "_",
                          gsub("-.*",
                               "",
                               MOD_RSD))) %>% 
  unique() 

# Load GOBP hierachy
offspring <- as.list(GOBPOFFSPRING)

# Extract child/offspring terms of mapped PSP terms
psp <-
  offspring[unique(psp_mapped[psp_mapped$parent == "parent",]$GOID)]

distr <- data.frame(
  term = AnnotationDbi::select(
    GO.db,
    keys = names(psp),
    columns = "TERM",
    keytype = "GOID"
  )$TERM,
  n = sapply(psp, length)
)
distr_gg <- ggplot(distr, aes(x = term, y = n, fill = term)) +
  geom_bar(stat = "identity") +
  guides(fill = "none") +
  coord_flip() +
  theme(plot.title.position = "plot") +
  ggtitle("Functional annotations in PSP",
          subtitle = paste("Total annotations of 4013 sites =", sum(distr$n)))


## 2) Which GO branch are predicted functional annotations found in...
# 2a) in the network?

# 2b) ...in the predicted annotations?

rwhn_tmp <- list(rwhn)

  funcs <-
    sighm$data %>%
    dplyr::select(name) %>%
    distinct()
  
  
  funcs$GOID <- AnnotationDbi::select(GO.db,
                                      keys = funcs$name,
                                      columns = "GOID",
                                      keytype = "TERM")$GOID
rwhn_go_branch <- data.frame(
    GOID = names(psp),
    inMLNW = sapply(psp, function(x)
      sum(funcs$GOID %in% x))
  ) %>% 
  mutate(TERM = AnnotationDbi::select(
    GO.db,
    keys = .$GOID,
    columns = "TERM",
    keytype = "GOID"
  )$TERM) %>%
  ggplot(aes(x = TERM, y = inMLNW, fill = TERM)) +
  geom_bar(stat = "identity")  +
  guides(fill = "none") +
  coord_flip() +
  ggtitle("RWHN output")
distr_gg + rwhn_go_branch + plot_layout(widths = c(1,5))

# 2c) ORA, PSP
ora <- readRDS("data/PSP_allsites_ORA.rds")$data

ora_go_branch <- data.frame(
  GOID = names(psp),
  inORA = sapply(psp, function(x)
    sum(ora$GOID %in% x))
) %>% 
  mutate(TERM = AnnotationDbi::select(
    GO.db,
    keys = .$GOID,
    columns = "TERM",
    keytype = "GOID"
  )$TERM) 

psp_distrs <- 
  data.frame(
    TERM = distr$term,
    n = distr$n,
    source = "Phosphosite Plus"
  ) %>% 
  ggplot(aes(x = TERM, y = n, fill= TERM )) +
  geom_bar(stat = "identity")  +
  guides(fill = "none") +
  facet_wrap(~  source, scales = "free_x") +
  coord_flip() +
  theme(plot.title = element_text(size = 5),
        strip.text = element_text(size = 5),
        legend.key.size = unit(.25, "cm"),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5), 
        axis.text.y = element_text(size = 5),
        axis.text.x = element_text(size = 4.5, angle = 45, hjust = 1, vjust = 1),
        axis.title = element_text(size = 5),
        legend.margin = margin(0,0,0,0, "cm")#,
        #plot.title.position = "plot"
  ) +
  ylab("No. sites annotated") +
  xlab("")


psp_distrs_results <- rbind(
  data.frame(
    TERM = ora_go_branch$TERM,
    n = ora_go_branch$inORA,
    source = "ORA results"
  ),
  data.frame(
    TERM = filter(rwhn_go_branch$data)$TERM,
    n = filter(rwhn_go_branch$data)$inMLNW,
    source = "RWHN results"
  )
) %>% 
  mutate(source = factor(source, levels = c("ORA results", 
                                            "RWHN results"))
  ) %>% 
  ggplot(aes(x = TERM, y = n, fill= TERM )) +
  geom_bar(stat = "identity")  +
  guides(fill = "none") +
  facet_wrap(~  source, scales = "free_x") +
  coord_flip() +
  theme(plot.title = element_text(size = 5),
        strip.text = element_text(size = 5),
        legend.key.size = unit(.25, "cm"),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5), 
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 4.5, angle = 45, hjust = 1, vjust = 1),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size = 5),
        legend.margin = margin(0,0,0,0, "cm")#,
        #plot.title.position = "plot"
  ) +
  ylab("No. GOBP Terms") +
  xlab("")


psp_distrs + 
  psp_distrs_results +
  plot_layout(widths = c(1,2)) +
  theme(plot.margin = margin(0,0,0,0, "cm"))

ggsave("results/figs/PSP_funcDistr.pdf",
       psp_distrs + 
         psp_distrs_results + 
         plot_layout(widths = c(1,2)) +
         plot_annotation(tag_levels = 'A') &
         theme(plot.margin = margin(0.05,0.05,0.05,0.05, "cm"),
               plot.tag = element_text(size = 6)),
       width = 3.33,
       height = 2.5,
       dpi = 300)