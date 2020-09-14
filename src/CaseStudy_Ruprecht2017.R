# From CRAN
library(dplyr)
library(tidyr)
library(igraph)
library(e1071)
library(enrichR)
library(ggplot2)
library(patchwork)
# From Bioconductor
library(Mfuzz)
library(GOSemSim)
library(AnnotationDbi)
library(GO.db)
# From src/
source("src/functions/imputePhosphoData.R")
source("src/functions/mfuzz-ggplot.R")
source("src/functions/simplifyGO.R")
source("src/functions/simplifyGOReqData.R")
source("src/functions/constructHetNet.R")
source("src/functions/calculateRWHN.R")
source("src/functions/dotplot_gg.R")

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

# PARENTAL + LAP M/L: Filter NAs and non significant
par_lap <- ruprecht_sty %>%
  dplyr::select(id, grep("ML", colnames(.))) %>%
  filter(SignificantML == "+") %>%
  dplyr::select(-SignificantML) %>%
  pivot_longer(-id,
               names_to = "rep",
               values_to = "ratio") %>%
  mutate(ratio = imputeLCMD::impute.QRILC(as.matrix(ratio))[[1]],
         rep = gsub("_R[1234]$", "", rep)) %>%
  group_by(id, rep) %>%
  summarise(ratio = median(ratio, na_rm = T)) %>%
  pivot_wider(id_cols = id,
              names_from = rep,
              values_from = ratio)


# RESISTANT  + LAP H/L
res_lap <- ruprecht_sty %>%
  dplyr::select(id, grep("HL", colnames(.))) %>%
  filter(SignificantHL == "+")  %>%
  dplyr::select(-SignificantHL) %>%
  pivot_longer(-id,
               names_to = "rep",
               values_to = "ratio") %>%
  mutate(ratio = imputeLCMD::impute.QRILC(as.matrix(ratio))[[1]],
         rep = gsub("_R[1234]$", "", rep)) %>%
  group_by(id, rep) %>%
  summarise(ratio = median(ratio, na_rm = T)) %>%
  pivot_wider(id_cols = id,
              names_from = rep,
              values_from = ratio)

# ALL COND. : Filter data with NA in an experimental condition
tot_lap <- ruprecht_sty %>%
  dplyr::select(id, grep("[HM]L", colnames(.))) %>%
  filter(SignificantML == "+" |
           SignificantHL == "+") %>%
  filter_missing(allowed = 1,  colnms = "^[HM]L_") %>%
  # dplyr::select(-c(SignificantML,
  #                  SignificantHL)
  # ) %>%
  pivot_longer(-c(id, SignificantML, SignificantHL),
               names_to = "rep",
               values_to = "ratio") %>%
  mutate(ratio = imputeLCMD::impute.QRILC(as.matrix(ratio))[[1]],
         rep = gsub("_R[1234]$", "", rep)) %>%
  group_by(id, rep) %>%
  mutate(ratio = median(ratio, na_rm = T)) %>%
  unique() %>%
  pivot_wider(#id_cols = id,
    names_from = rep,
    values_from = ratio)

# Cluster based on dynamics of phosphorylated sites
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



res_lap_cl <- kmeans(res_lap[, 2], 4)$cluster
names(res_lap_cl) <- res_lap$id

tot_lap_cl <- kmeans(tot_lap[, 4:5], 5)$cluster
names(tot_lap_cl) <- tot_lap$id


ggResCl <- lapply(1:max(res_lap_cl), function(i){
  ruprecht_sty %>% 
    filter(id %in% names(res_lap_cl[res_lap_cl == i])) %>% 
    filter(SignificantHL == "+") %>% 
    dplyr::select(grep("HL_", colnames(.))) %>% 
    pivot_longer(cols = everything()) %>% 
    summarise(mean = mean(value, na.rm = T),
              sd = sd(value, na.rm = T)) %>% 
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
              sd = sd(value, na.rm =T)) %>% 
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
ggsave("results/figs/Ruprecht_Clusters.tiff", gg_cl, width = 18.2, height = 10, units = "cm")

# Construct heterogeneous network
par_mlnw <- constructHetNet(stytxt = ruprecht_sty, 
                            phosphoData =  par_lap, 
                            clustering = par_lap_cl,
                            modules = T,
                            enrichrLib =  "KEGG_2019_Human",
                            pval = 0.05)

res_mlnw <- constructHetNet(stytxt = ruprecht_sty,
                            phosphoData = res_lap,
                            clustering =  res_lap_cl,
                            modules = T,
                            enrichrLib =  "KEGG_2019_Human",
                            pval = 0.05)

tot_mlnw <- constructHetNet(stytxt = ruprecht_sty,
                            phosphoData =  tot_lap[,-c(2:3)],
                            clustering =  tot_lap_cl,
                            modules = T,
                            enrichrLib =  "KEGG_2019_Human",
                            pval = 0.05)

## Run RWHN algorithm
# Recommend to run overnight or on HPC

# par_mlnw$edgelists <- lapply(par_mlnw$edgelists, function(i)
#   i %>% filter_all(any_vars(. != "negative regulation of transcription by RNA polymerase II"))
# )
# par_mlnw$v <- par_mlnw$v[par_mlnw$v$v != "negative regulation of transcription by RNA polymerase II",]

seed_par <- lapply(1:max(par_lap_cl), function(i){
  c(names(par_lap_cl[par_lap_cl == i]))
})

seed_res <- lapply(1:max(res_lap_cl), function(i){
  c(names(res_lap_cl[res_lap_cl == i]))
})

seed_tot <- lapply(1:max(tot_lap_cl), function(i){
  c(names(tot_lap_cl[tot_lap_cl == i]))
})

paste("start par", Sys.time())
rwhn_par <- lapply(seed_par, function(s){
  calculateRWHN(edgelists = par_mlnw$edgelists,
                verti = par_mlnw$v,
                seeds = s,
                transitionProb = 0.7,
                restart = 0.7,
                weight_xy = 0.3,
                weight_yz = 0.7) %>%
    filter(name %in% par_mlnw$v[par_mlnw$v$layer=="func",]$v)
})

paste("start res", Sys.time())
rwhn_res <- lapply(seed_res, function(s){
  calculateRWHN(edgelists = res_mlnw$edgelists,
                verti = res_mlnw$v,
                seeds = s,
                transitionProb = 0.7,
                restart = 0.7,
                weight_xy = 0.3,
                weight_yz = 0.7) %>%
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

saveRDS(rwhn_par, "results/data/rwhn_ruprecht_par.rds")
saveRDS(rwhn_res, "results/data/rwhn_ruprecht_res.rds")
saveRDS(rwhn_tot, "results/data/rwhn_ruprecht_tot.rds")

rwhn_par <- readRDS("results/data/rwhn_ruprecht_par.rds")
rwhn_res <- readRDS("results/data/rwhn_ruprecht_res.rds")
rwhn_tot <- readRDS("results/data/rwhn_ruprecht_tot.rds")


## res Filter top 5%
rwhn_res_flt <-  lapply(1:length(rwhn_res), function(i){
  rwhn_res[[i]] %>% 
    mutate(seed = i)
}) %>% do.call(rbind, .) %>% 
  group_by(name) %>% 
  mutate(rank_dif = (rank - mean(rank)),
         color = ifelse(rank_dif == 0, T, NA),
         V1 = signif(V1, digits = 2)) %>% 
  filter(rank_dif > 0) %>%                    # Filter terms that appear in the same position in all conditions
  ungroup() %>%
  group_split(seed) %>% 
  lapply(., function(i){
    pct <- i[1,]$V1                             # Filter top 5% of terms ( with a messy loop!!)
    df <- i[1,]
    
    for(x in 2:nrow(i)){
      if(pct < 0.05){
        df <- rbind(df, i[x,])
        pct <- pct + i[x,]$V1
      }else{
        break
      }
    }
    
    return(df)
  }) %>% 
  do.call(rbind,.)

sighm_res <- ggplot(rwhn_res_flt, aes(x = as.factor(seed), y = name)) +
  geom_tile(aes(fill = V1, color = color)) +
  scale_color_manual(values = c("red", NA), name = "Probability") +
  guides(color = FALSE) +
  xlab("Seed nodes") +
  ylab("GOBP Term") +
  theme(legend.key.size = unit(.5, "cm"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8), 
        title = element_text(size = 8),
        axis.text.y = element_text(size = 4.5),
        panel.background = element_rect(fill = "black"), 
        panel.grid = element_blank()) +
  scale_x_discrete(position = "top")  

ggsave(filename = "results/figs/rwhn_sig_Ruprecht_res.tiff",
       plot = sighm_res,
       width = 182,
       height = 79,
       units = "mm")  

## tot Filter top 5%
rwhn_tot_flt <- lapply(1:length(rwhn_tot), function(i){
  rwhn_tot[[i]] %>% 
    mutate(seed = i)
}) %>% do.call(rbind, .) %>% 
  group_by(name) %>% 
  mutate(rank_dif = (rank - mean(rank)),
         color = ifelse(rank_dif == 0, T, NA),
         V1 = signif(V1, digits = 2)) %>% 
  filter(rank_dif > 0) %>%                    # Filter terms that appear in the same position in all conditions
  ungroup() %>% 
  group_split(seed) %>% 
  lapply(., function(i){
    pct <- i[1,]$V1                             # Filter top 5% of terms ( with a messy loop!!)
    df <- i[1,]
    
    for(x in 2:nrow(i)){
      if(pct < 0.05){
        df <- rbind(df, i[x,])
        pct <- pct + i[x,]$V1
      }else{
        break
      }
    }
    
    return(df)
  }) %>% 
  do.call(rbind,.)

sighm_tot <- ggplot(rwhn_tot_flt, aes(x = as.factor(seed), y = name)) +
  geom_tile(aes(fill = V1, color = color)) +
  scale_color_manual(values = c("red", NA), name = "Probability") +
  guides(color = FALSE) +
  xlab("Seed nodes") +
  ylab("GOBP Term") +
  theme(legend.key.size = unit(.5, "cm"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8), 
        title = element_text(size = 8),
        axis.text.y = element_text(size = 4.5),
        panel.background = element_rect(fill = "black"), 
        panel.grid = element_blank()) +
  scale_x_discrete(position = "top")  

ggsave(filename = "results/figs/rwhn_sig_Ruprecht_tot.tiff",
       plot = sighm_tot,
       width = 182,
       height = 79,
       units = "mm")  

# visualise results with dot plot

dot_res_diff <- dotplot_gg(rwhn_res, n_terms = 30, remove_common = T, col = "PuBu", size = 3)
dot_tot_diff <- dotplot_gg(rwhn_tot, n_terms =30, remove_common = T, col = "PuBu", size = 3)


dot_res_diff[[1]]  <- dot_res_diff[[1]] + 
  theme(axis.text.x = element_text(size = 5),
        legend.key.size = unit(.5, "cm"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8), 
        title = element_text(size = 8),
        plot.margin = margin(10,10,10,50)
  ) +
  xlab("KEGG pathways") +
  ggtitle("RWHN results from \`lapatanib-resistant\` network")
dot_tot_diff[[1]]  <- dot_tot_diff[[1]] +
  theme(axis.text.x = element_text(size = 5),
        legend.key.size = unit(.5, "cm"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8), 
        title = element_text(size = 8),
        plot.margin = margin(10,10,10,50)
  ) +
  xlab("KEGG pathway")
  ggtitle("RWHN results from \'Total\' network")

gg <- dot_tot_diff[[1]]  / dot_res_diff[[1]] + plot_layout(guides = "collect")
gg <- gg + plot_annotation(tag_levels = "A")

ggsave("results/figs/rwhn_ruprecht_kegg_patchwork.tiff", gg, width = 18.2, height = 15, units = "cm")

######################
# Standard ORA analysis
######################

## Resistant
enrichedTerms_res <- lapply(1:max(res_lap_cl), function(i){
  ids <- names(res_lap_cl[res_lap_cl == i ] )
  cl_prots <- unique(gsub("_.*", "", ids))
  
  enriched <- enrichr(cl_prots, databases = "KEGG_2019_Human") %>% 
    .[[1]] %>%
    filter(Adjusted.P.value <= 0.05)  %>% 
    mutate(cluster = i)
  return(enriched)
}) %>% 
  do.call(rbind, .)

enrichedTerms_flt_res <- lapply(1:max(res_lap_cl), function(i){
  df <- enrichedTerms_res[enrichedTerms_res$cluster ==i, ] %>% 
    arrange(desc(Adjusted.P.value)) %>% 
    slice_min(Adjusted.P.value, n = 10) %>% 
    mutate(rwhn = apply(., 1, function(x) {
      if(x["Term"] %in% dot_res_diff[[2]]$name){
        ifelse(i %in% dot_res_diff[[2]][dot_res_diff[[2]]$name == x["Term"],]$seed, T, NA)
      }else{
        NA
      }
    }
    )
    )
  if(nrow(df) > 0){
    df <- mutate(df, rank = 1:n())
  }
}) %>% 
  do.call(rbind, .) %>% 
  mutate(V1 = signif(Adjusted.P.value, digits = 2),
         name = factor(Term, unique(Term))) %>% 
  ggplot(aes(y = name, x = as.factor(cluster))) +
  geom_tile(aes(fill = as.factor(cluster))) +
  theme_bw() +
  geom_point(aes(shape = rwhn)) +
  theme(legend.key.size = unit(.5, "cm"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8), 
        title = element_text(size = 8),
        #axis.text.x = element_text(size = 8, angle = 30, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 4.5),
        legend.position = "none", 
        panel.background = element_rect(fill = "black"), 
        panel.grid = element_blank()
        
  ) +
  scale_x_discrete(position = "top") +
  ylab("Pathway") +
  xlab("") +
  ggtitle("Resitant")

ggsave(filename = "results/figs/standardORA_Ruprecht_Res.tiff",
       plot = enrichedTerms_flt_res,
       width = 55,
       height = 80,
       units = "mm")

## Total
enrichedTerms_tot <- lapply(1:max(tot_lap_cl), function(i){
  ids <- names(tot_lap_cl[tot_lap_cl == i ] )
  cl_prots <- unique(gsub("_.*", "", ids))
  
  enriched <- enrichr(cl_prots, databases = "KEGG_2019_Human") %>% 
    .[[1]] %>%
    filter(Adjusted.P.value <= 0.05)  %>% 
    mutate(cluster = i)
  return(enriched)
}) %>% 
  do.call(rbind, .)

enrichedTerms_flt_tot <- lapply(1:max(tot_lap_cl), function(i){
  df <- enrichedTerms_tot[enrichedTerms_tot$cluster ==i, c("Term", "Adjusted.P.value", "cluster")] %>% 
    arrange(desc(Adjusted.P.value)) %>% 
    slice_min(Adjusted.P.value, n = 10) %>% 
    mutate(rwhn = apply(., 1, function(x) {
      if(x["Term"] %in% dot_tot_diff[[2]]$name){
        ifelse(i %in% dot_tot_diff[[2]][dot_tot_diff[[2]]$name == x["Term"],]$seed, T, NA)
      }else{
        NA
      }
    }
    )
    )

  return(df)
}) %>% 
  do.call(rbind, .) %>% 
  mutate(name = factor(Term, unique(Term))) %>% 
  ggplot(aes(y = name, x = as.factor(cluster))) +
  geom_tile(aes(fill = as.factor(cluster))) +
  theme_bw() +
 # geom_point(aes(shape = rwhn)) +
  theme(legend.key.size = unit(.5, "cm"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8), 
        title = element_text(size = 8),
        axis.text.y = element_text(size = 4.5),
        legend.position = "none", 
        panel.background = element_rect(fill = "black"), 
        panel.grid = element_blank()
        
  ) +
  scale_x_discrete(position = "top") +
  ylab("Pathway") +
  xlab("") +
  ggtitle("Total")

ggsave(filename = "results/figs/standardORA_Ruprecht_Tot.tiff",
       plot = enrichedTerms_flt_tot,
       width = 55,
       height = 80,
       units = "mm")




######
par_mlnw_GO <- constructHetNet(stytxt = ruprecht_sty,
                            phosphoData =  par_lap, 
                            clustering = par_lap_cl,
                            modules = T,
                            pval = 0.05)

res_mlnw_GO <- constructHetNet(stytxt = ruprecht_sty,
                            phosphoData = res_lap,
                            clustering =  res_lap_cl,
                            modules = T,
                            pval = 0.05)

tot_mlnw_GO <- constructHetNet(stytxt = ruprecht_sty,
                            phosphoData =  tot_lap[,-c(2:3)],
                            clustering =  tot_lap_cl,
                            modules = T,
                            pval = 0.05)



paste("start res", Sys.time())
rwhn_res_GO <- lapply(seed_res, function(s){
  calculateRWHN(edgelists = res_mlnw_GO$edgelists,
                verti = res_mlnw_GO$v,
                seeds = s,
                transitionProb = 0.7,
                restart = 0.7,
                weight_xy = 0.3,
                weight_yz = 0.7) %>%
    filter(name %in% res_mlnw_GO$v[res_mlnw_GO$v$layer=="func",]$v)
})

paste("start tot", Sys.time())
rwhn_tot_GO <- lapply(seed_tot, function(s){
  calculateRWHN(edgelists = tot_mlnw_GO$edgelists,
                verti = tot_mlnw_GO$v,
                seeds = s,
                transitionProb = 0.7,
                restart = 0.7,
                weight_xy = 0.3,
                weight_yz = 0.7) %>%
    filter(name %in% tot_mlnw_GO$v[tot_mlnw_GO$v$layer=="func",]$v)
})
paste("end", Sys.time())

saveRDS(rwhn_res_GO, "results/data/rwhn_ruprecht_res_GO.rds")
saveRDS(rwhn_tot_GO, "results/data/rwhn_ruprecht_tot_GO.rds")

rwhn_res_GO <- readRDS("results/data/rwhn_ruprecht_res_GO.rds")
rwhn_tot_GO <- readRDS("results/data/rwhn_ruprecht_tot_GO.rds")

dot_res_diff <- dotplot_gg(rwhn_res_GO, n_terms = 20, remove_common = T, size = 2)
dot_tot_diff <- dotplot_gg(rwhn_tot_GO, n_terms = 20, remove_common = T, size = 2)

dot_res_diff[[1]]  <- dot_res_diff[[1]] + 
  theme(axis.text.x = element_text(size = 4),
        legend.key.size = unit(.5, "cm"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8), 
        title = element_text(size = 8),
        plot.margin = margin(10,10,10,20)
  ) +
  xlab("GOBP Terms") +
  ggtitle("RWHN results from \`lapatanib-resistant\` network")
dot_tot_diff[[1]]  <- dot_tot_diff[[1]] +
  theme(axis.text.x = element_text(size = 4),
        legend.key.size = unit(.5, "cm"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8), 
        title = element_text(size = 8),
        plot.margin = margin(10,10,10,20)
) +
  xlab("GOBP Terms") +
  ggtitle("RWHN results from \'Total\' network")


gg <- dot_tot_diff[[1]]  / dot_res_diff[[1]] + plot_layout(guides = "collect")
gg <- gg + plot_annotation(tag_levels = "A")

ggsave("results/figs/rwhn_ruprecht_GO_patchwork.tiff", gg, width = 18.2, height = 15, units = "cm")
