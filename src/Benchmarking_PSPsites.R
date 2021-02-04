# From CRAN
library(dplyr)
library(tidyr)
library(igraph)
library(enrichR)
library(patchwork)
library(ggplot2)
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

# This data contains proteins that are known to be involved in broadly similar processes,
# with sites that have functional annotations in PSP.
# We then cluster the sites of those proteins based on their PSP annotation.
# Then apply the alogirthm
# We would expect sites with the same functional annotation to rank similar GOBP terms higher.
# We can check this using the F-score, etc.
data <- readRDS("data/PSP_Benchmarks.rds")
hetNet <- lapply(data, function(i){
  dat <- i %>%  
    dplyr::select(-GO) %>% 
    mutate(id_site = paste0(SYMBOL,
                            "_",
                            gsub("-.*",
                                 "",
                                 MOD_RSD)
    )
    ) %>% 
    unique() %>% 
    group_by(id_site) %>% 
    filter(n() > 1)
  
  clustering <- dat$ON_PROCESS
  names(clustering) <- dat$id_site
  
  hetNet <- constructHetNet(clustering = clustering)
  
  return(list(
    hetNet = hetNet,
    clustering = clustering)
)
})
saveRDS(hetNet, "PSP_benchmarks_hetNet.rds")

rwhn <- lapply(hetNet, function(i){
  
  seed_l <- lapply(unique(i$clustering)[order(unique(i$clustering))], function(x){
  c(names(i$clustering[i$clustering == x]))
})
  rwhn <- lapply(seed_l, function(s){
    calculateRWHN(edgelists = i$hetNet$edgelists,
                  verti = i$hetNet$v,
                  seeds = s,
                  transitionProb = 0.7,
                  restart = 0.7,
                  weight_xy = 0.3,
                  weight_yz = 0.7,
                  eps = 1/10^12) %>%
      filter(name %in% i$hetNet$v[i$hetNet$v$layer=="func",]$v)
  })
})
saveRDS(rwhn, "results/data/rwhn_benchmarks.rds")


####
rwhn <- readRDS("results/data/rwhn_benchmarks.rds")
clustering <- lapply(readRDS("results/data/PSP_benchmarks_hetNet.rds"), 
                     function(i){
                       i$clustering
                     }
)

rwhn <- sapply(names(rwhn), function(i) {
  
  x <- rwhn[[i]]
  names(x) <- unique(clustering[[i]])[order(unique(clustering[[i]]))]
  
  return(x)
  
}, simplify = F, USE.NAMES = T)

# tmp <- sapply(names(rwhn), function(i){
# 
#   rwhn_s <- rwhn[[i]]
#   
#   rwhn_df <- lapply(1:length(rwhn_s), function(x)
#     mutate(rwhn_s[[x]], 
#            seed = unique(clustering[[i]])[order(unique(clustering[[i]]))][x],
#            rank = 1:nrow(rwhn_s[[x]])
#     )
#   ) %>% 
#     do.call(cbind, .) 
#   readr::write_csv(rwhn_df, paste0("results/data/PSP_Benchmark", i, "_rwhn_results.csv"))
#   
#   ## Filter top 5%
#   sighm <- heatmap_RWHN(rwhn_output = rwhn[[i]], ylab = "GOBP Term") + ggtitle(i) 
#   
#   ggsave(filename = paste0("results/figs/PSP_Benchmark", i, "_rwhn_results.tiff"),
#          plot = sighm,
#          width = 182,
#          height = 85,
#          units = "mm")  
#   
#   ## ORA
#   enrichedTerms <- overrepresentationAnalysis(clustering = clustering[[i]],
#                                                   RWHN_sig = sighm,
#                                                   colours = c("#effff6","#168d49"))
#   
#   ggsave(filename = paste0("results/figs/PSP_Benchmark", i, "_ORA_results.tiff"),
#          plot = enrichedTerms,
#          width = 100,
#          height = 80,
#          units = "mm")
#   
#   return(list(
#     ORA = enrichedTerms$data,
#     RWHN = sighm$data
#     )
#   )
# }, simplify = F, USE.NAMES = T)


tmp <- lapply(names(rwhn), function(i){
  sighm <- heatmap_RWHN(rwhn_output = rwhn[[i]], ylab = "GOBP Term") + ggtitle(i) 

  return(sighm$data)
})
names(tmp) <- names(rwhn)
  
true <- readr::read_tsv("data/Regulatory_sites_GOmapped.tsv") %>% 
  arrange(ON_PROCESS) 

true_list <- true %>% 
  group_by(ON_PROCESS) %>% 
  group_split() %>% 
  lapply(function(i){
    x <- i %>% 
      separate_rows(offspring, sep = ";")
    unique(c(x$GOID, x$offspring))
  })
names(true_list) <- true$ON_PROCESS

library(GO.db)


vis <- lapply(names(tmp), function(i){
  pred <- tmp[[i]] %>% 
    filter(seed %in% names(true_list)) %>% 
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
  
  names(pred_list) <- paste0(unique(pred$seed))
  
  sapply(names(pred_list), function(k){ 
    sapply(names(true_list)[names(true_list) %in% names(pred_list)], function(x){
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
  })
})

lapply(vis, function(mat){
  mat %>% 
        as.data.frame() %>%
        tibble::rownames_to_column("Predicted") %>%
        pivot_longer(cols = -Predicted,
                     names_to = "True",
                     values_to = "Number_of_terms") %>%
        ggplot(aes(x = Predicted, y = True, fill = Number_of_terms)) +
        geom_tile() +
        scale_fill_viridis_c(limits = c(0, 1)) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
})

lapply(1:length(vis), function(i){
  df <- vis[[i]] %>% 
    as.data.frame() %>%
    tibble::rownames_to_column("Predicted") %>%
    pivot_longer(cols = -Predicted,
                 names_to = "True",
                 values_to = "F1")
  
  frames <- data.frame(Predicted = rownames(vis[[i]]), True = colnames(vis[[i]]))
  
  g <- ggplot(df, aes(x = Predicted, y = True, fill = F1)) +
    geom_tile() +
    scale_fill_viridis_c(limits = c(0, 1)) +
    geom_tile(data=frames, aes(x = Predicted, y = True), size=1, fill=NA, colour="black") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle(names(data)[i])
  
  ggsave(paste0("results/figs/controls/benchmark_tests/", names(data)[i], ".tiff"),
         g)
  })


# #pred$GOID <- 
# 
# vis <-  lapply(names(tmp), function(i){
#   pred <- tmp[[i]] %>% 
#     filter(seed %in% names(true_list)) %>% 
#     arrange(seed)
#   
#   pred$GOID <- (AnnotationDbi::select(GO.db,
#                                       pred$name,
#                                       "GOID",
#                                       keytype = "TERM"))$GOID
#   pred_list <- pred %>% 
#     group_by(seed) %>% 
#     group_split() %>% 
#     lapply(function(x) x$GOID)
#   
#   names(pred_list) <- paste0(unique(pred$seed), "_PRED")
#   
#   #Filter true list for those terms found in clusters of benchmark data
#   mat <- sapply(true_list[true$ON_PROCESS[true$ON_PROCESS %in% pred$seed]], function(k){ 
#     sapply(pred_list, function(x){
#       sum(x %in% k)
#     })
#   })
#   
#   
#   g <- mat %>% 
#     as.data.frame() %>% 
#     tibble::rownames_to_column("Predicted") %>% 
#     pivot_longer(cols = -Predicted,
#                  names_to = "True",
#                  values_to = "Number_of_terms") %>% 
#     ggplot(aes(x = Predicted, y = True, fill = Number_of_terms)) +
#     geom_tile() +
#     scale_fill_viridis_b() +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
#     ggtitle(i)
#   
#   return(g)
# })
# 
# lapply(vis, function(i) 
#   ggsave(paste0("results/data/benchmark_tests/", i$labels$title, ".tiff"),
#          i)
# )
