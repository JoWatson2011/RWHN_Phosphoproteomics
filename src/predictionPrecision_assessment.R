####
# Set up
####

#Import networks and RWHN ranks for prediction
egf_mlnw <- readRDS("results/data/egf_mlnw.rds")
tgf_mlnw <- readRDS("results/data/tgf_mlnw.rds")
rwhn_egf <- readRDS("results/data/rwhn_egf_clusters.rds")
rwhn_tgf <- readRDS("results/data/rwhn_tgf_clusters.rds")
psp_mlnw <- readRDS("results/data/PSP_benchmark_all_hetNet.rds")
rwhn_psp <- readRDS("results/data/rwhn_all_benchmarks.rds")
#Load libraries
library(GO.db)
library(tidyverse)
library(igraph)
library(ggraph)
library(patchwork)
source("src/functions/heatmap_RWHNsig.R")
#Import "True" annotations from PSP, mapped to GO terms
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

######
# Analysis
######

# What is the distribution of PSP functional annotations:
# which GO "branches" do they fall into ?
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

distr

## 2) Which GO branch are predicted functional annotations found in...
# 2a) in the network?
mlnw <- list(egf = egf_mlnw,
             tgf = tgf_mlnw,
             psp = psp_mlnw)

mlnw_go_branch <- sapply(1:length(mlnw), function(i) {
  funcs <- mlnw[[i]]$v[mlnw[[i]]$v$layer == "func",]
  funcs$GOID <- AnnotationDbi::select(GO.db,
                                      keys = funcs$v,
                                      columns = "GOID",
                                      keytype = "TERM")$GOID
  df <- data.frame(
    GOID = names(psp),
    inMLNW = sapply(psp, function(x)
      sum(funcs$GOID %in% x)),
    data = names(mlnw)[i]
  )
}, simplify = F, USE.NAMES = T) %>%
  do.call(rbind, .) %>%
  mutate(TERM = AnnotationDbi::select(
    GO.db,
    keys = .$GOID,
    columns = "TERM",
    keytype = "GOID"
  )$TERM) %>%
  ggplot(aes(x = TERM, y = inMLNW, fill = TERM)) +
  geom_bar(stat = "identity")  +
  guides(fill = "none") +
  facet_wrap(~ data, scales = "free") +
  coord_flip() +
  ggtitle("Function layer of network")

distr_gg + mlnw_go_branch + plot_layout(widths = c(1,5))


# 2b) ...in the predicted annotations?
rwhn <- list(egf = rwhn_egf,
             tgf = rwhn_tgf,
             psp = rwhn_psp)

rwhn_go_branch <- sapply(1:length(rwhn), function(i) {
  funcs <-
    heatmap_RWHN(rwhn_output = rwhn[[i]], ylab = "GOBP Term")$data %>%
    select(name) %>%
    distinct()
  
  
  funcs$GOID <- AnnotationDbi::select(GO.db,
                                      keys = funcs$name,
                                      columns = "GOID",
                                      keytype = "TERM")$GOID
  df <- data.frame(
    GOID = names(psp),
    inMLNW = sapply(psp, function(x)
      sum(funcs$GOID %in% x)),
    data = names(mlnw)[i]
  )
}, simplify = F, USE.NAMES = T) %>%
  do.call(rbind, .) %>%
  mutate(TERM = AnnotationDbi::select(
    GO.db,
    keys = .$GOID,
    columns = "TERM",
    keytype = "GOID"
  )$TERM) %>%
  ggplot(aes(x = TERM, y = inMLNW, fill = TERM)) +
  geom_bar(stat = "identity")  +
  guides(fill = "none") +
  facet_wrap(~ data) +
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
    TERM = filter(rwhn_go_branch$data, data == "psp")$TERM,
    n = filter(rwhn_go_branch$data, data == "psp")$inMLNW,
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


## 3) Is functional enrichment of terms derived from a particular part of the PPI?

networks <- lapply(1:length(mlnw), function(i) {
  nw <- rbind(mlnw[[i]]$edgelists$y,
              mlnw[[i]]$edgelists$yz,
              mlnw[[i]]$edgelists$z) %>%
    graph_from_data_frame(vertices = mlnw[[i]]$v[mlnw[[i]]$v$layer %in% c("prot", "func"),],
                          directed = F)
  
  
  tmp <- lapply(unique(psp_mapped$parentProc), function(i) {
    vids <-
      V(nw)[V(nw)$name %in% psp_mapped[psp_mapped$parentProc == i,]$TERM]
    
    if (length(vids) > 0) {
      nws <- make_ego_graph(nw, nodes = vids)
      
      vids <-
        do.call(rbind, lapply(nws, as_data_frame, what = "vertices")) %>%
        mutate(psp = i)
      
      return(vids)
    } else{
      return(data.frame())
    }
  }) %>%
    do.call(rbind, .)
  
  #tmp %>% group_by(name, psp) %>%  summarise(n = n()) %>%  arrange(desc(n))
  
  newnw <- nw
  V(newnw)$nbrs <- ifelse(V(newnw)$name %in% tmp$name, T, F)
  
  ggraph(induced_subgraph(newnw, v = V(newnw)$nbrs), layout = "nicely") +
    geom_edge_link() +
    geom_node_point(aes(color = layer), size = 2) +
    ggtitle(
      paste(
        "Induced subgraph of protein - function edges in",
        names(mlnw)[i],
        "network"
      ),
      subtitle = paste(
        length(V(newnw)[V(newnw)$nbrs == T & V(nw)$layer == "prot"]),
        "/" ,
        length(V(newnw)[V(newnw)$layer == "prot"]),
        "proteins.",
        length(V(newnw)[V(newnw)$nbrs == T &
                          V(newnw)$layer == "func"]),
        "functions."
      )
    )
})

# 4) What is the overlap between "True" and "predicted" annotations

true <- psp_mapped %>% 
  filter(parent == "parent") %>% 
  select(ON_PROCESS, GOID) %>% 
  unique() %>%
  merge(psp_sites[,c("id_site", "ON_PROCESS")], by = "ON_PROCESS")

forUpset <- lapply(1:length(rwhn), function(i){ 
  funcs <-
    heatmap_RWHN(rwhn_output = rwhn[[i]], ylab = "GOBP Term")$data %>% 
    distinct()
  
  funcs$GOID <- AnnotationDbi::select(GO.db,
                                      keys = funcs$name,
                                      columns = "GOID",
                                      keytype = "TERM")$GOID
  return(funcs)
})
