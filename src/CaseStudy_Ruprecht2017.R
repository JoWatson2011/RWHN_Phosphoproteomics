# From CRAN
library(dplyr)
library(tidyr)
library(igraph)
library(enrichR)
library(ggplot2)
library(patchwork)
# From Bioconductor
library(GOSemSim)
library(AnnotationDbi)
library(GO.db)
# From src/
source("src/functions/overrepresentationAnalysis.R")
source("src/functions/heatmap_RWHNsig.R")
source("src/functions/imputePhosphoData.R")
source("src/functions/simplifyGO.R")
source("src/functions/simplifyGOReqData.R")
source("src/functions/constructHetNet.R")
source("src/functions/calculateRWHN.R")

# Import data
ruprecht_sty <-
  data.table::fread(
    input = "data/Ruprecht_STY_2017.csv",
    select = c(13, 30, 37,
               39:46,
               93:94,
               100),
    skip = 1
  ) %>%
  mutate(
    `Gene names` = gsub(";.*", "", `Gene names`),
    id = paste0(`Gene names`, "_", `Amino acid` , `Position`)
  ) %>%
  filter(!duplicated(id))
colnames(ruprecht_sty) <-
  c(
    "amino.acid",
    "localisation.prob",
    "position",
    "HL_R1",
    "HL_R2",
    "HL_R3",
    "HL_R4",
    "ML_R1",
    "ML_R2",
    "ML_R3",
    "ML_R4",
    "SignificantML",
    "SignificantHL",
    "gene.symbol",
    "id"
  )

# RESISTANT  + LAP H/L
set.seed(1)
res_lap <- ruprecht_sty %>%
  dplyr::select(id, grep("HL", colnames(.))) %>%
  filter(SignificantHL == "+")  %>%
  filter_missing(allowed = 1,  colnms = "^HL_") %>%
  #  dplyr::select(-SignificantHL) %>%
  pivot_longer(-c(id, SignificantHL),
               names_to = "rep",
               values_to = "ratio") %>%
  mutate(ratio = imputeQRLIC.mod(as.matrix(ratio)),
         rep = gsub("_R[1234]$", "", rep)) %>%
  group_by(id, rep) %>%
  mutate(ratio = median(ratio, na_rm = T)) %>%
  unique() %>%
  pivot_wider(id_cols = id,
              names_from = rep,
              values_from = ratio)

# ALL COND. : Filter data with NA in an experimental condition
set.seed(1)
tot_lap <- ruprecht_sty %>%
  dplyr::select(id, grep("[HM]L", colnames(.))) %>%
  filter(SignificantML == "+" |
           SignificantHL == "+") %>%
  filter_missing(allowed = 1,  colnms = "^[HM]L_") %>%
  pivot_longer(-c(id, SignificantML, SignificantHL),
               names_to = "rep",
               values_to = "ratio") %>%
  mutate(ratio = imputeQRLIC.mod(as.matrix(ratio)),
         rep = gsub("_R[1234]$", "", rep)) %>%
  group_by(id, rep) %>%
  mutate(ratio = median(ratio, na_rm = T)) %>%
  unique() %>%
  pivot_wider(#id_cols = id,
    names_from = rep,
    values_from = ratio)

#Cluster based on dynamics of phosphorylated sites
set.seed(1)
wss <- data.frame(
  type = c(
    rep("Parental+lap. or Resistant+lap. cells", 15),
    rep("Resistant+lap. cells", 15)
  ),
  col = c(rep("#00798c",15),
  rep("#d1495b",15)),
  cl = rep(1:15, 2),
  wss = c(
    sapply(1:15,
           function(k) {
             kmeans(tot_lap[, 4:5],
                    k, nstart = 50,
                    iter.max = 15)$tot.withinss
           }),
    sapply(1:15,
           function(k) {
             kmeans(res_lap[, 2],
                    k, nstart = 50,
                    iter.max = 15)$tot.withinss
           })
  ),
  stringsAsFactors = F
)
elbow <- ggplot(wss, aes(x = cl, y = wss, color = type)) +
  geom_point() +
  geom_line() +
  scale_color_discrete("Sites significantly changing in...") +
  xlab("Number of clusters (k)") +
  ylab("Total within-clusters sum of squares") +
  theme_minimal() +
  theme(legend.justification = c(1, 1),
        legend.position = c(1, 1),
        legend.box.background = element_rect(color = "black", fill = "white"),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        axis.title = element_text(size = 5),
        axis.text = element_text(size = 5)
  )


# Cluster based on dynamics of phosphorylated sites
set.seed(1)
res_lap_cl <- kmeans(res_lap[, 2], 4)$cluster
names(res_lap_cl) <- res_lap$id


set.seed(1)
tot_lap_cl <- kmeans(tot_lap[, 4:5], 5)$cluster
names(tot_lap_cl) <- tot_lap$id



ggResCl <- lapply(1:max(res_lap_cl), function(i){
  ruprecht_sty %>% 
    filter(id %in% names(res_lap_cl[res_lap_cl == i])) %>% 
    filter(SignificantHL == "+") %>% 
    dplyr::select(grep("HL_", colnames(.))) %>% 
    pivot_longer(cols = everything()) %>% 
    summarise(mean = mean(value, na.rm = T),
              sd = sd(value, na.rm = T),
              .groups = "keep") %>% 
    mutate(cl = i,
           exp = "res")
}) %>% do.call(rbind, .) %>% 
  ggplot(aes(x = as.factor(cl), color = as.factor(cl))) +
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), size = 1, width = .2) +
  geom_point(aes(y = mean), size = 2) +
  theme_minimal() +
  theme(strip.text = element_blank(),
              panel.border = element_rect(fill = NA, color = "black"), 
              legend.position = "none",
        axis.title = element_text(size = 5),
        axis.text.y = element_text(size = 5)
        ) +
  ylab("Average change from Parental (Untreated)") +
  xlab("Cluster")

ggTotCl <- lapply(1:max(tot_lap_cl), function(i){
  ruprecht_sty %>% 
    filter(id %in% names(tot_lap_cl[tot_lap_cl == i])) %>% 
    filter(SignificantHL == "+" | SignificantML == "+") %>% 
    dplyr::select(grep("[HM]L_", colnames(.))) %>% 
    pivot_longer(cols = everything()) %>% 
    mutate(name = gsub("_R[1234]", "", name)) %>% 
    group_by(name) %>% 
    summarise(mean = mean(value, na.rm = T),
              sd = sd(value, na.rm =T),
              .groups = "keep") %>% 
    mutate(namecl = paste0(i, "_", name),
           cl = i)
}) %>% do.call(rbind, .) %>% 
  ggplot(aes(x = factor(name, levels = c("ML", "HL")), y = mean, group = cl, color = as.factor(cl))) +
  geom_line(size = 1) +
  geom_point() +
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width = .2, size = 1) +
  facet_wrap(~cl, nrow = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.border = element_rect(fill = NA, color = "black"), 
        legend.position = "none",
        axis.title = element_text(size = 5),
        axis.text = element_text(size = 5)
        ) +
  scale_x_discrete(name = "",
                   labels = c("Parental +lap.", "Resistant + lap.")) +
  ylab("Average change from Parental (Untreated)")

gg_cl <- elbow + (ggTotCl / ggResCl) + 
  plot_layout(widths = c(1,2)) +
  plot_annotation(tag_levels = 'A')
#ggsave("results/figs/Ruprecht_Clusters.tiff", gg_cl, width = 18.2, height = 10, units = "cm", dpi = "print")

# Construct heterogeneous network
res_mlnw <- constructHetNet(phosphoData = res_lap,
                            clustering =  res_lap_cl,
                            modules = T,
                            enrichrLib =  "KEGG_2019_Human",
                            stringPath = "data/STRINGexpmtgene_highconf.rds",
                            pval = 0.05)

tot_mlnw <- constructHetNet(phosphoData =  tot_lap[,-c(2:3)],
                            clustering =  tot_lap_cl,
                            modules = T,
                            enrichrLib =  "KEGG_2019_Human",
                            stringPath = "data/STRINGexpmtgene_highconf.rds",
                            pval = 0.05)

## Run RWHN algorithm
# Recommend to run overnight or on HPC

seed_res <- lapply(1:max(res_lap_cl), function(i){
  c(names(res_lap_cl[res_lap_cl == i]))
})

seed_tot <- lapply(1:max(tot_lap_cl), function(i){
  c(names(tot_lap_cl[tot_lap_cl == i]))
})


paste("start res", Sys.time())
rwhn_res <- lapply(seed_res, function(s){
  calculateRWHN(edgelists = res_mlnw$edgelists,
                verti = res_mlnw$v,
                seeds = s,
                transitionProb = 0.7,
                restart = 0.7,
                weight_xy = 0.3,
                weight_yz = 0.7,
                ) %>%
    filter(name %in% res_mlnw$v[res_mlnw$v$layer=="func",]$v)
})

paste("start tot", Sys.time())
rwhn_tot <- lapply(seed_tot, function(s){
  calculateRWHN(edgelists = tot_mlnw$edgelists,
                verti = tot_mlnw$v,
                seeds = s,
                transitionProb = 0.7,
                restart = 0.7,
                weight_xy = 0.3,
                weight_yz = 0.7) %>%
    filter(name %in% tot_mlnw$v[tot_mlnw$v$layer=="func",]$v)
})

end <- Sys.time()

saveRDS(rwhn_res, "results/data/rwhn_ruprecht_res.rds")
saveRDS(rwhn_tot, "results/data/rwhn_ruprecht_tot.rds")


#########
# Visualisation, etc. 
#########
rwhn_res <- readRDS("results/data/rwhn_ruprecht_res.rds")
rwhn_tot <- readRDS("results/data/rwhn_ruprecht_tot.rds")


rwhn_res_df <- lapply(1:length(rwhn_res), function(i)
  mutate(rwhn_res[[i]], 
         seed = i,
         rank = 1:nrow(rwhn_res[[i]])
  )
) %>% 
  do.call(cbind, .) 
readr::write_csv(rwhn_res_df, "results/data/Ruprecht_res_rwhn_results.csv")

rwhn_tot_df <- lapply(1:length(rwhn_tot), function(i)
  mutate(rwhn_tot[[i]], 
         seed = i,
         rank = 1:nrow(rwhn_tot[[i]])
  )
) %>% 
  do.call(cbind, .) 
readr::write_csv(rwhn_tot_df, "results/data/Ruprecht_tot_rwhn_results.csv")


# Visualise Results
sighm_res <- heatmap_RWHN(rwhn_output = rwhn_res,
                          ylab = "KEGG Pathway", 
                          colours = c("#8fc0d1","#209dc9"))


## tot Filter top 5%
sighm_tot <- heatmap_RWHN(rwhn_tot,
                          ylab = "KEGG Pathway",
                          colours = c("#8fc0d1","#209dc9"))


######################
# Standard ORA analysis
######################

## Resistant
enrichedTerms_res <- overrepresentationAnalysis(clustering = res_lap_cl,
                                                RWHN_sig = sighm_res,
                                                ylab = "KEGG Pathway",
                                                colours = c("#effff6","#168d49"),
                                                database = "KEGG_2019_Human")


## Total
enrichedTerms_tot <- overrepresentationAnalysis(tot_lap_cl, 
                                                RWHN_sig = sighm_tot,
                                                ylab = "KEGG Pathway",
                                                colours = c("#d5d2e2","#168d49"),
                                                database = "KEGG_2019_Human")

ggsave("results/figs/ruprecht_RWHN_ORA.pdf", 
       sighm_tot + sighm_res +
         enrichedTerms_tot + enrichedTerms_res + 
         plot_layout(ncol = 2, 
                     #widths= c(2, 1),
                     heights = c(1,2)),
       width = 180,
       height = 150,
       units = "mm",
       dpi = "print"
)



paste("RWHN:",
      length(unique(sighm_res$data$name)),
      "ORA:",
      length(unique(enrichedTerms_res$data$Term))
)

paste("RWHN:",
      length(unique(sighm_tot$data$name)),
      "ORA:",
      length(unique(enrichedTerms_tot$data$Term))
)

lapply(1:max(res_lap_cl), function(i){
  paste("RWHN:",
        nrow(sighm_res$data[sighm_res$data$seed == i,]),
        "ORA:",
        nrow(enrichedTerms_res$data[enrichedTerms_res$data$cluster == i,])
  )
        })


#########################
# Final visualisation for paper
#########################
ora_res <- enrichedTerms_res$data %>% 
  dplyr::select(term = Term, value = Adjusted.P.value, cluster) %>% 
  mutate(method = "ORA")
rwhn_res <- sighm_res$data %>% 
  dplyr::select(term = name, value = rank, cluster = seed) %>% 
  mutate(method = "RWHN")

rwhn_res_gg <- rwhn_res %>% 
  rbind(ora_res) %>% 
  mutate(cluster = ifelse(method == "RWHN", cluster, NA),
         value = ifelse(method == "RWHN", value, NA)
  ) %>% 
  arrange(term) %>% 
  ggplot(aes(y = term,
             x = as.factor(cluster)
  )
  ) +
  geom_tile(aes(fill = value), color = "white") +
  scale_fill_gradient(low = "#5bd670", high = "#def6e2",
                      limits = c(1, 15)) +
  guides(color = FALSE,
         fill = guide_colourbar(title="Rank", reverse = T)) +
  theme_minimal() +
  theme(legend.key.size = unit(.25, "cm"),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 6),
        legend.margin = margin(0,0,0,0, "cm"),
        strip.background = element_blank(),
        strip.text.y = element_blank()
  ) +
  scale_x_discrete(#breaks = factor(c(1:5)), limits = factor(c(1:5)), 
                   position = "top", na.translate = FALSE) +
  xlab("Cluster") +
  ylab("GOBP") 

ora_res_gg <- ora_res %>% 
  rbind(rwhn_res) %>% 
  mutate(cluster = ifelse(method == "ORA", cluster, NA),
         value = ifelse(method == "ORA", value, NA)
  ) %>% 
  arrange(term) %>% 
  ggplot(aes(y = term,
             x = as.factor(cluster)
  )
  ) +
  geom_tile(aes(fill = value), color = "white") +
  scale_fill_gradient(low = "#845bd6", high = "#e6def6",
                      limits = c(0, 0.05)) +
  guides(color = FALSE,
         fill = guide_colourbar(title="Rank", reverse = T)) +
  theme_minimal() +
  theme(legend.key.size = unit(.25, "cm"),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 5),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 6),
        legend.margin = margin(0,0,0,0, "cm")
  ) +
  scale_x_discrete(#breaks = factor(c(1:5)), limits = factor(c(1:5)), 
                   position = "top", na.translate = FALSE) +
  xlab("Cluster") +
  ylab("GOBP") 
##
# TGF
ora_tot <- enrichedTerms_tot$data %>% 
  dplyr::select(term = Term,  value = Adjusted.P.value, cluster, ) %>% 
  mutate(method = "ORA")
rwhn_tot <- sighm_tot$data %>% 
  dplyr::select(term = name, value = rank, cluster = seed) %>% 
  mutate(method = "RWHN")

rwhn_tot_gg <- rwhn_tot %>% 
  rbind(ora_tot) %>% 
  mutate(cluster = ifelse(method == "RWHN", cluster, NA),
         value = ifelse(method == "RWHN", value, NA)
  ) %>% 
  arrange(term) %>% 
  ggplot(aes(y = term,
             x = as.factor(cluster)
  )
  ) +
  geom_tile(aes(fill = value), color = "white") +
  scale_fill_gradient(low = "#5bd670", high = "#def6e2",
                      limits= c(1,15)) +
  guides(color = FALSE,
         fill = guide_colourbar(title="Rank", reverse = T)) +
  theme_minimal() +
  theme(legend.key.size = unit(.25, "cm"),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 6),
        legend.margin = margin(0,0,0,0, "cm"),
        strip.background = element_blank(),
        strip.text.y = element_blank()
  ) +
  scale_x_discrete(#breaks = factor(c(1:5)), limits = factor(c(1:5)), 
                   position = "top", na.translate = FALSE) +
  xlab("Cluster") +
  ylab("GOBP") 

ora_tot_gg <- ora_tot %>% 
  rbind(rwhn_tot) %>% 
  mutate(cluster = ifelse(method == "ORA", cluster, NA),
         value = ifelse(method == "ORA", value, NA)
  ) %>% 
  arrange(term) %>% 
  ggplot(aes(y = term,
             x = as.factor(cluster)
  )
  ) +
  geom_tile(aes(fill = value), color = "white") +
  scale_fill_gradient(low = "#845bd6", high = "#e6def6",
                      limits = c(0, 0.05)) +
  guides(color = FALSE,
         fill = guide_colourbar(title="Rank", reverse = T)) +
  theme_minimal() +
  theme(legend.key.size = unit(.25, "cm"),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 5),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 6),
        legend.margin = margin(0,0,0,0, "cm")
  ) +
  scale_x_discrete(#breaks = factor(c(1:5)), limits = factor(c(1:5)), 
                   position = "top", na.translate = FALSE) +
  xlab("Cluster") +
  ylab("GOBP") 

ggsave(
  "results/figs/Ruprecht_RWHN_ORA_v2.tiff",
  rwhn_tot_gg +
    ora_tot_gg + 
    rwhn_res_gg + 
    ora_res_gg +
    plot_layout(nrow = 1, guides= "collect"),
  width = 7,
  height = 6,
  units = "in",
  dpi = 300
)

ggsave(
  "results/figs/Ruprecht_RWHN_ORA_v2.pdf",
  rwhn_tot_gg +
    ora_tot_gg + 
    rwhn_res_gg + 
    ora_res_gg +
    plot_layout(nrow = 1, guides= "collect"),
  width = 7,
  height = 6,
  units = "in",
  dpi = 300
)
