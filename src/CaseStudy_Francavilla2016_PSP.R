# From CRAN
library(dplyr)
library(tidyr)
library(igraph)
library(e1071)
library(enrichR)
library(patchwork)
library(ggplot2)
# From Bioconductor
library(Mfuzz)
library(GOSemSim)
library(AnnotationDbi)
library(GO.db)
# From src/
source("src/functions/overrepresentationAnalysis.R")
source("src/functions/heatmap_RWHNsig.R")
source("src/functions/imputePhosphoData.R")
source("src/functions/mfuzz-ggplot.R")
source("src/functions/simplifyGO.R")
source("src/functions/simplifyGOReqData.R")
source("src/functions/constructHetNet.R")
source("src/functions/calculateRWHN.R")

set.seed(123)

# Import data
psp <- readr::read_tsv("data/Regulatory_sites.txt", skip = 2) %>% 
  filter(!is.na(ON_PROCESS),
         ORGANISM == "human",
         grepl("-p", MOD_RSD)) %>% 
  separate_rows(ON_PROCESS, sep = "; ") %>% 
  mutate(id_site = paste0(PROTEIN,
                          "_",
                          gsub("-.*",
                               "",
                               MOD_RSD))) %>% 
  unique() 


cfr_sty <- data.table::fread(input = "data/CFR_STY_2016.csv",
                             select = c(1:9, 26:35)
) %>% 
  mutate(id = paste0(`Gene names`, "_", `Swiss-Prot phosphosite`)) %>% 
  filter(!duplicated(id),
         id %in% psp$id_site) %>% 
  rename(`Gene names` = "gene.symbol")
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
  
tgf <- cfr_sty %>% 
  dplyr::select(id, grep("TGF", colnames(cfr_sty))) %>% 
  filter(Regulated_by_TGFalfa == "+") %>% 
#  filter_at(vars(matches("ratio")), ~ . > 0.5 | . < -1) %>%
  dplyr::select(-Regulated_by_TGFalfa) %>% 
  filter_missing(allowed = 2, colnms = "ratio") %>%  
  mutate_at(vars(matches("ratio")), impute_fun)


# Cluster based on dynamics of phosphorylated sites
egf_fcm <- fuzzyC(egf, 6)
tgf_fcm <- fuzzyC(tgf, 6)

data <- rbind(tgf_fcm$g$data, egf_fcm$g$data) %>% 
  mutate(stim = c(rep("TGF-a", nrow(tgf_fcm$g$data)),
                  rep("EGF", nrow(egf_fcm$g$data)))
  )

ggcl <- ggplot(data, aes(x = Time, y = expression, group = Identifier, color = maxMembership)) +
  geom_line() +
  facet_grid(stim ~ clusterindex) +
  scale_x_discrete(limits = c("1", "8", "40", "90")) +
  geom_hline(yintercept = 0, 
             colour="#666666", alpha = 0.5) + 
  theme_minimal() +
  ylab("Phosphorylation changes from baseline") + labs(fill= "Membership") +
  scale_colour_gradient(low = "white", high = "red")

# ggsave("results/figs/francavilla_clusters.tiff", ggcl,
#        width = 8.3, height = 3.7, units = "in", dpi = 300)

# Construct heterogeneous network
egf_mlnw <- constructHetNet(phosphoData = egf, clustering = egf_fcm$clustering, simplify = F)
tgf_mlnw <- constructHetNet(tgf, tgf_fcm$clustering, simplify = F)

## Run RWHN algorithm
# Recommend to run overnight or on HPC
# EGF:
seed_l <- lapply(1:max(egf_fcm$clustering), function(i){
  c(names(egf_fcm$clustering[egf_fcm$clustering == i]))
})
rwhn_egf <- lapply(seed_l, function(s){
  calculateRWHN(edgelists = egf_mlnw$edgelists,
                verti = egf_mlnw$v,
                seeds = s,
                transitionProb = 0.7,
                restart = 0.7,
                weight_xy = 0.3,
                weight_yz = 0.7,
                eps = 1/10^12) %>%
    filter(name %in% egf_mlnw$v[egf_mlnw$v$layer=="func",]$v)
})
#saveRDS(rwhn_egf, "results/data/rwhn_egf_clusters.rds")

#TGFa
seed_l <- lapply(1:max(tgf_fcm$clustering), function(i){
  c(names(tgf_fcm$clustering[tgf_fcm$clustering == i]))
})
rwhn_tgf <- lapply(seed_l, function(s){
  calculateRWHN(edgelists = tgf_mlnw$edgelists,
                verti = tgf_mlnw$v,
                seeds = s,
                transitionProb = 0.7,
                restart = 0.7,
                weight_xy = 0.3,
                weight_yz = 0.7,
                eps = 1/10^12) %>%
    filter(name %in% tgf_mlnw$v[tgf_mlnw$v$layer=="func",]$v)
})
#saveRDS(rwhn_tgf, "results/data/rwhn_tgf_clusters.rds")



# ######
# # Export, Visualise, etc.
# ######
# rwhn_egf <- readRDS("results/data/rwhn_egf_clusters.rds")
# rwhn_tgf <- readRDS("results/data/rwhn_tgf_clusters.rds")
# 
# rwhn_egf_df <- lapply(1:length(rwhn_egf), function(i)
#   mutate(rwhn_egf[[i]], 
#          seed = i,
#          rank = 1:nrow(rwhn_egf[[i]])
#   )
# ) %>% 
#   do.call(cbind, .) 
# readr::write_csv(rwhn_egf_df, "results/data/CFR_egf_rwhn_results.csv")
# 
# rwhn_tgf_df <- lapply(1:length(rwhn_tgf), function(i)
#   mutate(rwhn_tgf[[i]], 
#          seed = i,
#          rank = 1:nrow(rwhn_tgf[[i]])
#   )
# ) %>% 
#   do.call(cbind, .) 
# readr::write_csv(rwhn_tgf_df, "results/data/CFR_tgf_rwhn_results.csv")
# 
# 
# ## EGF Filter top 5%
sighm_egf <- heatmap_RWHN(rwhn_output = rwhn_egf, ylab = "GOBP Term")
# 
# ggsave(filename = "results/figs/rwhn_sig_CFR_egf.tiff",
#        plot = sighm_egf,
#        width = 182,
#        height = 85,
#        units = "mm")  
# 
# ## tgf Filter top 5%
# 
sighm_tgf <- heatmap_RWHN(rwhn_tgf, "GOBP Term")
# 
# ggsave(filename = "results/figs/rwhn_sig_CFR_tgf.tiff",
#        plot = sighm_tgf,
#        width = 182,
#        height = 79,
#        units = "mm")  
# 
# 
# ######################
# # Standard ORA analysis
# ######################
# 
# ## EGF
enrichedTerms_egf <- overrepresentationAnalysis(clustering = egf_fcm$clustering,
                                                RWHN_sig = sighm_egf,
                                                colours = c("#effff6","#168d49"))
# 
# ggsave(filename = "results/figs/standardORA_Francavilla_EGF.tiff",
#        plot = enrichedTerms_egf,
#        width = 100,
#        height = 80,
#        units = "mm")
# 
# ## TGF
enrichedTerms_tgf <- overrepresentationAnalysis(clustering = tgf_fcm$clustering,
                                                RWHN_sig = sighm_tgf,
                                                colours = c("#effff6","#168d49"))
# 
# ggsave(filename = "results/figs/standardORA_Francavilla_TGF.tiff",
#        plot = enrichedTerms_tgf,
#        width = 100,
#        height = 80,
#        units = "mm")
# 
# 
# ggsave("results/figs/Francavilla_RWHN_ORA.pdf", 
#        sighm_egf + enrichedTerms_egf + sighm_tgf + enrichedTerms_tgf +
#          plot_layout(ncol = 2, 
#                      widths= c(1.8, 1),
#                      heights = c(1.3, 1),
#                      guides=  "collect"),
#        width = 182,
#        height = 150,
#        units = "mm",
#        dpi = "print"
# )
#        
# ####
# # F-SCORE
# ####
# 
# #RWHN
rwhn_f1 <- lapply(list(
  egf = list(cl = egf_fcm$clustering,
             rwhn_output = sighm_egf$data),
  tgf = list(cl = tgf_fcm$clustering,
             rwhn_output = sighm_tgf$data)
), function(i){
  true_mapped <- readr::read_tsv("data/Regulatory_sites_GOmapped.tsv") %>%
    arrange(ON_PROCESS)

  true <- readr::read_tsv("data/Regulatory_sites.txt", skip = 2) %>%
    mutate(id_site = paste0(PROTEIN, "_", gsub("-.*", "", MOD_RSD))) %>%
    filter(ORGANISM == "human" &
             grepl("-p", MOD_RSD),
           !is.na(ON_PROCESS)) %>%
    separate_rows(ON_PROCESS, sep ="; ") %>%
    dplyr::select(id_site, ON_PROCESS) %>%
    distinct() %>%
    merge(data.frame(id_site = names(i$cl), cl = i$cl), by = "id_site") %>%
    merge(true_mapped[,c("ON_PROCESS", "GOID", "offspring")], by = "ON_PROCESS")

  true_list <- true %>%
    group_by(cl) %>%
    group_split() %>%
    lapply(function(i){
      x <- i %>%
        separate_rows(offspring, sep = ";")
      unique(c(x$GOID, x$offspring))
    })

  library(GO.db)


  pred <- i$rwhn_output %>%
    arrange(seed)

  pred$GOID <- (AnnotationDbi::select(GO.db,
                                      pred$name,
                                      "GOID",
                                      keytype = "TERM"))$GOID

  pred_list <- pred %>%
    group_by(seed) %>%
    group_split() %>%
    lapply(function(x){
      ID <- x$GOID

      offspring <- do.call(c, as.list(GOBPOFFSPRING)[ID])

      return(c(ID, offspring))
    })
  sapply(1:length(pred_list), function(k){
    sapply(1:length(true_list), function(x){
      tp <- sum(pred_list[[k]] %in% true_list[[x]])
      fn <- sum(!true_list[[k]] %in% pred_list[[x]])
      fp <- sum(!pred_list[[k]] %in% true_list[[x]])

      # recall <- TruePositives / (TruePositives + FalseNegatives)
      recall <- tp / (tp + fn)
      # precision <- TruePositives / (TruePositives + FalsePositives)
      precision <- tp / (tp + fp)
      # F-Measure = (2 * Precision * Recall) / (Precision + Recall)
      F1 <- (2 * precision * recall) / (precision + recall)

      return(F1)
    })
  }) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Predicted") %>%
    pivot_longer(cols = -Predicted,
                 names_to = "True",
                 values_to = "F1") %>%
    ggplot(aes(x = Predicted, y = True, fill = F1)) +
    geom_tile() +
    scale_fill_viridis_c(limits = c(0, 1)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
})

###ORA
ora_f1 <- lapply(list(
  egf = list(cl = egf_fcm$clustering,
             ora_output = enrichedTerms_egf$data),
  tgf = list(cl = tgf_fcm$clustering,
             ora_output = enrichedTerms_tgf$data)
), function(i){

  true_mapped <- readr::read_tsv("data/Regulatory_sites_GOmapped.tsv") %>%
    arrange(ON_PROCESS)

  true <- readr::read_tsv("data/Regulatory_sites.txt", skip = 2) %>%
    mutate(id_site = paste0(PROTEIN, "_", gsub("-.*", "", MOD_RSD))) %>%
    filter(ORGANISM == "human" &
             grepl("-p", MOD_RSD),
           !is.na(ON_PROCESS)) %>%
    separate_rows(ON_PROCESS, sep ="; ") %>%
    dplyr::select(id_site, ON_PROCESS) %>%
    distinct() %>%
    merge(data.frame(id_site = names(i$cl), cl = i$cl), by = "id_site") %>%
    merge(true_mapped[,c("ON_PROCESS", "GOID", "offspring")], by = "ON_PROCESS")

  true_list <- true %>%
    group_by(cl) %>%
    group_split() %>%
    lapply(function(i){
      x <- i %>%
        separate_rows(offspring, sep = ";")
      unique(c(x$GOID, x$offspring))
    })

  library(GO.db)
  pred <- i$ora_output %>%
    arrange(cluster)

  pred_list <- pred %>%
    group_by(cluster) %>%
    group_split() %>%
    lapply(function(x){
      ID <- x$GOID

      offspring <- do.call(c, as.list(GOBPOFFSPRING)[ID])

      return(c(ID, offspring))
    })

  names(pred_list) <- unique(pred$cluster)

  true_list <- true_list[unique(pred$cluster)]
  names(true_list) <- unique(pred$cluster)
  sapply(names(pred_list), function(k){
    sapply(names(true_list), function(x){
      tp <- sum(pred_list[[k]] %in% true_list[[x]])
      fn <- sum(!true_list[[k]] %in% pred_list[[x]])
      fp <- sum(!pred_list[[k]] %in% true_list[[x]])

      # recall <- TruePositives / (TruePositives + FalseNegatives)
      recall <- tp / (tp + fn)
      # precision <- TruePositives / (TruePositives + FalsePositives)
      precision <- tp / (tp + fp)
      # F-Measure = (2 * Precision * Recall) / (Precision + Recall)
      F1 <- (2 * precision * recall) / (precision + recall)

      return(F1)
    })
  }) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Predicted") %>%
    pivot_longer(cols = -Predicted,
                 names_to = "True",
                 values_to = "F1") %>%
    ggplot(aes(x = Predicted, y = True, fill = F1)) +
    geom_tile() +
    scale_fill_viridis_c(limits = c(0, 1)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
})

