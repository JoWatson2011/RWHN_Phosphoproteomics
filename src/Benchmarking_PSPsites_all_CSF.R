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
names(rwhn) <- unique(clustering)[order(unique(clustering))]


sighm <- heatmap_RWHN(rwhn_output = rwhn, ylab = "GOBP Term")
pred_rwhn <- sapply(unique(sighm$data$seed)[order(unique(sighm$data$seed))], 
                    function(i){
                      
                      goid <- AnnotationDbi::select(GO.db,
                                            sighm$data[sighm$data$seed == i,]$name,
                                            "GOID",
                                            "TERM"
                                            )$GOID
                      goid[!is.na(goid)]
                    } , simplify = F, USE.NAMES = T)

anc <- as.list(GOBPANCESTOR)
rwhn_ancestor <- sapply(pred_rwhn, function(i){
  unique(do.call(c, lapply(i, function(x) anc[[x]])))
})

offsp <-as.list(GOBPOFFSPRING)
rwhn_offspring <- sapply(pred_rwhn, function(i){
  unique(do.call(c, offsp[unique(i)]))
}, USE.NAMES = T)

###SAVED BELOW
# ora <- overrepresentationAnalysis(clustering = clustering,
#                                   RWHN_sig = sighm,
#                                   colours = c("#effff6","#168d49"))
# saveRDS(ora, "data/PSP_allsites_ORA.rds")

ora <- readRDS("data/PSP_allsites_ORA.rds")
pred_ora <- sapply(unique(ora$data$cluster)[order(unique(ora$data$cluster))], 
       function(i) {
         goid <- ora$data[ora$data$cluster == i,]$GOID
         
         unique(goid[!is.na(goid)])
         
         }, simplify = F, USE.NAMES = T)

ora_ancestor <- sapply(pred_ora, function(i){
  unique(do.call(c, lapply(i, function(x) anc[[x]])))
})
ora_offspring <- sapply(pred_ora, function(i){
  unique(do.call(c, offsp[unique(i)]))
}, USE.NAMES = T)


psp_mapped <- readr::read_tsv("data/Regulatory_sites_GOmapped.tsv") %>%
  dplyr::select(ON_PROCESS, GOID) %>% 
  unique() 

psp_sites <- readr::read_tsv("data/Regulatory_sites.txt", skip = 2) %>% 
  mutate(id_site = paste0(PROTEIN, "_", gsub("-.*", "", MOD_RSD))) %>% 
  filter(ORGANISM == "human" &
           grepl("-p", MOD_RSD),
         !is.na(ON_PROCESS)) %>% 
  separate_rows(ON_PROCESS, sep ="; ") %>% 
  dplyr::select(id_site, ON_PROCESS) %>% 
  distinct() %>% 
  merge(data.frame(id_site = names(clustering), cl = clustering), by = "id_site") %>% 
  merge(psp_mapped[,c("ON_PROCESS", "GOID")], by = "ON_PROCESS")

#TRUE annots
psp_exact <- sapply(unique(psp_sites$cl),
                    function(i) unique(psp_mapped[psp_mapped$ON_PROCESS==i,]$GOID), USE.NAMES = T)
psp_exact <- psp_exact[!(names(psp_exact) %in% c("carcinogenesis, altered",
                                                 "carcinogenesis, inhibited",
                                                 "carcinogenesis, induced",
                                                 "neural plasticity"))]
psp_exact <- psp_exact[order(names(psp_exact))]

psp_offspring <- sapply(psp_exact, function(i){
  goid <- unique(do.call(c, as.list(GOBPOFFSPRING)[unique(i)]))
  
  goid[!is.na(goid)]
}, USE.NAMES = T)




venn <- sapply(names(pred_ora), function(i){
  forVenn <- list(
    RWHN = unique(c(pred_rwhn[[i]], rwhn_ancestor[[i]], rwhn_offspring[[i]])),
    ORA = unique(c(pred_ora[[i]], ora_ancestor[[i]], ora_offspring[[i]])),
    PSP = psp_offspring[[i]]
  )
  
  # names(forVenn) <- paste(names(forVenn),":\n",
  #                         sapply(forVenn, length),
  #                         "terms"
  #                         
  # )
  ggVennDiagram::ggVennDiagram(forVenn, label = "count") + 
    ggtitle(i, "With offspring / ancestors of RWHN and ORA results") +
    # scale_fill_gradientn(name = "Enrichment (p-value)",
    #                      colors = c("white", "grey", "red"),
    #                      values = scales::rescale(c(10000, 200, 0)),
    #                      guide = "colorbar") +
    scale_fill_gradient(low="white",high = "white") +
    guides(fill = "none")  
}, simplify = F, USE.NAMES = T)

sapply(names(venn), function(i) ggsave(paste0("results/figs/PSP_tests/",i,".tiff"), venn[[i]]))

#phyper(
# overlap,
# set 1,
# set 2,
# not in set 1 / 2
#)

overlap<-sum(rwhn_ancestor[[1]] %in% psp_offspring[[1]])
set1 <- length(rwhn_ancestor[[1]])
set2 <- length(unique(psp_offspring[[1]]))
bg <- sum(length(unique(unname(do.call(c, psp_offspring)))), set1, set2)
phyper(
  overlap,
  set1,
  bg - set2,
  set1
  #lower.tail = F
)

## Does PSP_EXACT appear in ancestors of ORA/RWHN terms?

categories <- names(pred_rwhn)[names(pred_rwhn) %in% c("apoptosis, altered",
                                                       "apoptosis, induced",
                                                       "apoptosis, inhibited",
                                                       "cell cycle regulation",
                                                       "cell differentiation, altered",
                                                       "cell differentiation, induced",
                                                       "cell motility, altered",
                                                       "cell motility, induced",
                                                       "cell motility, inhibited",
                                                       "chromatin organization, altered",
                                                       "DNA repair, altered",
                                                       "DNA repair, induced",
                                                       "DNA repair, inhibited",
                                                       "signaling pathway regulation",
                                                       "transcription, altered",
                                                       "transcription, induced",
                                                       "transcription, inhibited" 
                                                       
)]

ora_ancs_gg <- sapply(psp_exact[categories], function(true) {
  sapply(pred_ora[categories],function(pred){ 
    sum(
      sapply(pred, function(GOID){
        true %in% anc[[GOID]]
      })
    )/length(pred) * 100
  })
})  %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("ORA") %>% 
  pivot_longer(cols = -ORA, names_to = "PSP") %>%  
  mutate(value = ifelse(value == 0, NA, value)) %>%
  ggplot(aes(x = ORA, y= PSP)) +
  geom_tile(aes(fill = value), color = "white") +
  scale_fill_gradient2(
    low ="white",
    high = "black",
    limit = c(0, 100),
    space = "Lab",
    name = "%",
    na.value = "white"
  ) +
  scale_x_discrete(drop = F) +
  scale_y_discrete(drop = F) +
  xlab("ORA") +
  ylab("PSP") +
  theme(legend.key.size = unit(.25, "cm"),
        legend.title = element_text(size = 4.5),
        legend.text = element_text(size = 4.5), 
        axis.text.y = element_text(size = 4.5),
        panel.grid = element_blank(),
        legend.margin = margin(0,0,0,0, "cm"),
        axis.text.x = element_text(size = 4.5, angle = 90, hjust = 1),
        axis.title.x = element_text(size = 6),
        axis.title.y = element_text(size = 6, angle = 90),
        axis.ticks = element_blank()
  ) +
  geom_tile(data = data.frame(PSP = names(psp_exact[categories]),
                              ORA = names(pred_ora[categories])),
            color = "#C19D9D", fill = NA) 

rwhn_ancs_gg <- sapply(psp_exact[categories], function(true) {
  sapply(pred_rwhn[categories],function(pred){ 
    sum(
      sapply(pred, function(GOID){
        true %in% anc[[GOID]]
      })
    )/length(pred) *100
  })
}) %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("RWHN") %>% 
  pivot_longer(cols = -RWHN, names_to = "PSP") %>%  
  mutate(value = ifelse(value == 0, NA, value)) %>%
  ggplot(aes(x = RWHN, y= PSP)) +
    geom_tile(aes(fill = value), color = "white") +
    scale_fill_gradient2(
      low ="white",
      high = "black",
      limit = c(0, 100),
      space = "Lab",
      name = "%",
      na.value = "white"
    ) +
  scale_x_discrete(drop = F) +
  scale_y_discrete(drop = F) +
  xlab("RWHN") +
  ylab("PSP") +
  theme(legend.key.size = unit(.25, "cm"),
        legend.title = element_text(size = 4.5),
        legend.text = element_text(size = 4.5), 
        axis.text.y = element_text(size = 4.5),
        panel.grid = element_blank(),
        legend.margin = margin(0,0,0,0, "cm"),
        axis.text.x = element_text(size = 4.5, angle = 90, hjust = 1),
        axis.title.x = element_text(size = 6),
        axis.title.y = element_text(size = 6, angle = 90),
        axis.ticks = element_blank()
  ) +
  geom_tile(data = data.frame(PSP = names(psp_exact[categories]),
    RWHN = names(pred_rwhn[categories])),
            color = "#C19D9D", fill = NA) 

ora_ancs_gg + rwhn_ancs_gg + 
  plot_layout(guides = 'collect')


ggsave("results/figs/PSP_anc_term.pdf", 
       ora_ancs_gg + rwhn_ancs_gg+ 
         plot_layout(guides = 'collect'),
       width = 7,
       height = 3.5,
       units = "in",
       dpi = "print"
)



rwhn_ora_ancs_gg <- sapply(pred_rwhn[categories], function(r) {
  sapply(categories,function(o){ 
    sum(
      sapply(r, function(GOID){
        any(GOID %in% ora_ancestor[[o]]| GOID %in% ora_offspring[[o]])
      })
    )/length(r) * 100
  })
})  %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("ORA") %>% 
  pivot_longer(cols = -ORA, names_to = "RWHN") %>%  
  mutate(value = ifelse(value == 0, NA, value)) %>%
  ggplot(aes(x = ORA, y= RWHN)) +
  geom_tile(aes(fill = value), color = "white") +
  scale_fill_gradient2(
    low ="white",
    high = "black",
    limit = c(0, 100),
    space = "Lab",
    name = "%",
    na.value = "white"
  ) +
  scale_x_discrete(drop = F) +
  scale_y_discrete(drop = F) +
  xlab("ORA") +
  ylab("RWHN") +
  theme(legend.key.size = unit(.25, "cm"),
        legend.title = element_text(size = 4.5),
        legend.text = element_text(size = 4.5), 
        axis.text.y = element_text(size = 4.5),
        panel.grid = element_blank(),
        legend.margin = margin(0,0,0,0, "cm"),
        axis.text.x = element_text(size = 4.5, angle = 90, hjust = 1),
        axis.title.x = element_text(size = 6),
        axis.title.y = element_text(size = 6, angle = 90),
        axis.ticks = element_blank()
  ) +
  geom_tile(data = data.frame(RWHN = names(pred_rwhn[categories]),
                              ORA = names(pred_ora[categories])),
            color = "#C19D9D", fill = NA) 



# 
# psp_mapped <- readr::read_tsv("data/Regulatory_sites_GOmapped.tsv") %>%
#   dplyr::select(ON_PROCESS, GOID) %>% 
#   unique() 
# 
# psp_sites <- readr::read_tsv("data/Regulatory_sites.txt", skip = 2) %>% 
#   mutate(id_site = paste0(PROTEIN, "_", gsub("-.*", "", MOD_RSD))) %>% 
#   filter(ORGANISM == "human" &
#            grepl("-p", MOD_RSD),
#          !is.na(ON_PROCESS)) %>% 
#   separate_rows(ON_PROCESS, sep ="; ") %>% 
#   dplyr::select(id_site, ON_PROCESS) %>% 
#   distinct() %>% 
#   merge(psp_mapped[,c("ON_PROCESS", "GOID")], by = "ON_PROCESS")
# 
# #TRUE annots
# psp_exact <- sapply(unique(psp_sites$ON_PROCESS),
#                     function(i) unique(psp_sites[psp_sites$ON_PROCESS==i,]$GOID), USE.NAMES = T)
# 
# psp_offspring <- sapply(psp_exact, function(i){
#   unique(do.call(c, as.list(GOBPOFFSPRING)[unique(i)]))
# }, USE.NAMES = T)
# 
# #PREDICTED annots
# rwhn_results <- sapply(unique(sighm$data$seed), function(x){
#   na.omit(AnnotationDbi::select(GO.db,
#                                 sighm$data[sighm$data$seed == x,]$name,
#                                 "GOID",
#                                 "TERM"
#   )$GOID)
# }, USE.NAMES = T)
# 
# rwhn_offspring <- sapply(rwhn_results, function(i){
#   unique(do.call(c, as.list(GOBPOFFSPRING)[unique(i)]))
# }, USE.NAMES = T)
# 
# rwhn_ancestors <- sapply(rwhn_results, function(i){
#   unique(do.call(c, as.list(GOBPANCESTOR)[unique(i)]))
# }, USE.NAMES = T)
# 
# rwhn_allGObranch <- sapply(1:length(rwhn_results), function(i){
#   c(rwhn_offspring[[i]], rwhn_ancestors[[i]])
# }, USE.NAMES = T)
# 
# #  names(psp_offspring) <- 
# #    if(is.null(names(psp_offspring))){
# #      paste0(1:length(psp_offspring), "_pred")
# # }else{
# #    paste0(names(psp_offspring), "_pred")
# # }
# 
# 
# sapply(psp_offspring, function(true){
#   sapply(rwhn_results, function(pred){
#     sum(true %in% pred)/length(pred)*100
#   })
# }) %>% 
#   as.data.frame() %>% 
#   tibble::rownames_to_column("true") %>%
#   pivot_longer(cols = -true,
#                names_to ="Pred_cl",
#                values_to = "pct") %>% 
#   mutate(Pred_cl = gsub("V", "", Pred_cl),
#          border = ifelse(true == Pred_cl, "black", NA))  %>% 
#   ggplot(aes(x = true, y = Pred_cl)) +
#   geom_tile(aes(fill = pct, color = border)) + guides(color = "none")
# 
# 
