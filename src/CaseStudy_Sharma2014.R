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
source("src/functions/mfuzz-ggplot.R")
source("src/functions/simplifyGO.R")
source("src/functions/imputePhosphoData.R")
source("src/functions/simplifyGOReqData.R")
source("src/functions/constructHetNet.R")
source("src/functions/calculateRWHN.R")

set.seed(123)

# Import data
sharma_sty <- data.table::fread(input = "data/Sharmaetal_2014.csv"
                                ) %>% 
  as.data.frame() %>% 
  filter(Class1 == "Yes") %>% 
  mutate(`Gene names` = gsub(";.*", "", `Gene names`),
    id = paste0(`Gene names`, "_", `Amino acid`, `Position`)) %>% 
  filter(!duplicated(id)) %>% 
  dplyr::select(-grep("pY_PV", colnames(.)))

colnames(sharma_sty)[grep("Intensity.*[CEN]", 
                          colnames(sharma_sty))] <- gsub("Intensity ", "",
                gsub("pY_", "",
                     gsub("[15]+_", "",
                          colnames(sharma_sty)[grep("Intensity.*[CEN]", 
                                                    colnames(sharma_sty))]
                     )
                )
)

expcols <- grep("^[CEN][1-6]", colnames(sharma_sty), value = T)

norm <- sharma_sty %>% 
  dplyr::select(id, all_of(expcols)) %>% 
  filter(rowSums(.[,expcols] > 0) >= 1) %>% 
  pivot_longer(-id) %>%
  group_by(id) %>% 
  #filter(sum(value == 0) > 1) %>% 
  mutate(value = ifelse(value == 0, NA, value)) %>% 
  mutate(value = log(value)) %>% 
  pivot_wider(values_fn = mean)
norm[,unique(expcols)] <-   limma::normalizeBetweenArrays(norm[,unique(expcols)], method = "quantile")
norm_med <- norm %>% 
  ungroup() %>% 
  pivot_longer(-id) %>%
  mutate(name =  gsub("[1-6]$", "", name)) %>% 
  mutate(value = imputeQRLIC.mod(as.matrix(value))) %>%
  group_by(id) %>%
  group_split() %>%
  lapply( function(i){

    df <- i
    dat <- i[,c("name", "value")]
    av <- aov(value ~ name, data = dat)
    p <- summary(av)[[1]][1, 5]
    df <- cbind(i, p)
    return(df)
  }) %>%
  do.call(rbind, .) %>%
  mutate(adjp = p.adjust(p, method = "fdr")) %>%
  filter(adjp < 0.05) %>%
  dplyr::select(name, id, value) %>% 
  group_by(id) %>% 
  mutate(value = (value - mean(value, na.rm = T))/
                          sd(value, na.rm = T)) %>% 
  unique() %>% 
  pivot_wider(values_fn = median) %>% 
  na.omit()
  
set.seed(1)
cl <- cutree(hclust(dist(norm_med)), k = 5)
names(cl) <- norm_med$id

norm_med$cl <- cl


gg_cl <- norm_med %>%
  pivot_longer(-c(id, cl)) %>%
  arrange(cl) %>%
  ggplot(aes(x =  name, y = forcats::fct_reorder(id, cl), fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "red", mid = "black", high = "green") +
  facet_wrap(~cl, ncol = 1, scales = "free_y", strip.position = "right") +
  theme(
    axis.text.y = element_blank(),
    panel.spacing.y=unit(0.1, "lines")
  )


mlnw <- constructHetNet(phosphoData = dplyr::select(norm_med, -cl),
                        clustering =  cl,
                        modules = F,
                        enrichrLib =  "GO_Biological_Process_2018",
                        stringPath = "data/STRINGexpmtgene_highconf.rds",
                        pval = 0.05)


seed <- lapply(1:max(cl), function(i){
  names(cl[cl==i])
})

rwhn <- lapply(seed, function(s){
  calculateRWHN(edgelists = mlnw$edgelists,
                verti = mlnw$v,
                seeds = s,
                transitionProb = 0.7,
                restart = 0.7,
                weight_xy = 0.3,
                weight_yz = 0.7,
  ) %>%
    filter(name %in% mlnw$v[mlnw$v$layer=="func",]$v)
})


saveRDS(rwhn, "results/data/rwhn_sharma.rds")

############

rwhn <- readRDS("results/data/rwhn_sharma.rds")

sighm <- heatmap_RWHN(rwhn_output = rwhn, ylab = "GOBP Term")

ora <- overrepresentationAnalysis(clustering = cl,
                                      RWHN_sig = sighm,
                                      colours = c("#effff6","#168d49")
) 

ora$data %>% 
  filter(rwhn == T) %>% 
  ggplot(aes(x = as.factor(cluster), y = Term, fill = Adjusted.P.value)) +
  geom_tile()

ora$data %>% 
  filter(grepl("mitosi[sc]", Term) | 
           grepl("cell cycle", Term)) %>%  
  group_by(Term) %>% 
  filter(n() < 5) %>% 
  summarise(cluster = paste(cluster, collapse = " "),
             .groups = "keep")




  pct <- rwhn[[2]][1,]$V1                             # Filter top 5% of terms ( with a messy loop!!)
  df <- rwhn[[2]][1,]
  for(x in 2:nrow(rwhn[[2]])){
    if(pct < 0.05){
      df <- rbind(df, rwhn[[2]][x,])
      pct <- pct + rwhn[[2]][x,]$V1
    }else{
      break
    }
  }
df$rank <- 1:nrow(df) 




