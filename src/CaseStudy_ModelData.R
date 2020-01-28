# From CRAN
library(dplyr)
library(tidyr)
library(igraph)
library(e1071)
library(enrichR)
# From Bioconductor
library(Mfuzz)
library(GOSemSim)
library(AnnotationDbi)
library(GO.db)
# From src/
source("src/mfuzz-ggplot.R")
source("src/simplifyGO.R")
source("src/simplifyGOReqData.R")
source("src/constructHetNet.R")
source("src/calculateRWHN.R")
source("src/dotplot_gg.R")

prots <- data.frame(node = c("RAF1", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3", "JUND", "DUSP6", "RPS6KA3"),
                    level = c(1, 2, 2, 3, 3, 4, 5, 5),
                    stringsAsFactors = F)
### 
# Phospho
###
model <- readr::read_csv("data/model_phospho.csv")

fuzzyC <- function(tidydata, clus){
  set.seed(1)
  inputdata <- tidydata %>% 
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(var = "id") # rownames must be id - see mfuzz.ggplot function
  phospep <- new("ExpressionSet", 
                 exprs = as.matrix(inputdata) )
  phospep.z <- standardise(phospep)      #Avg expression of each peptide = 0, sd = 1
  optimalM <- mestimate(phospep.z)   #set fuzzifier
  # tmp  <- Dmin(phospep.z,m=optimalM,crange=seq(4,20,1),repeats=3,visu=TRUE)    #How many clusters? Note: stochastic
  cl <- mfuzz(phospep.z, c = clus, m = optimalM)
  a <- mFuzz.ggplot(exprs(phospep.z), cl, returnData = T)
  a$g <- a$g + ggtitle(substr(a$exp$sample[1], 1, 6), clus) + theme_dark()
  
  return(list(data = phospep.z, exp = a$exp, g = a$g, clustering = cl$cluster))
}

dat <- model %>%
  mutate(id = paste0(prot, "_", site)) %>%
  dplyr::select(id, t1, t2, t3, t4, t5)

pp <- fuzzyC(dat, cl = 5)

phos <- lapply(1:max(pp$clustering), function(x){
  xpnd <- pp$clustering[pp$clustering == x]
  
  return(expand.grid(names(xpnd), names(xpnd)) %>% 
           mutate(cl = x))
}) %>% 
  do.call(rbind, .) %>% 
  graph_from_data_frame(directed = F) %>% 
  simplify() %>% 
  igraph::as_data_frame(what = "edges") %>% 
  dplyr::select(phos1 = 1, phos2 = 2)

###
# Phospho >> Prot
###
phos_prot <- graph_from_data_frame(phos, directed = F) %>%
  igraph::as_data_frame("vertices") %>% 
  mutate(prot = apply(., 1, function(x) strsplit(x[1], "_")[[1]][1])) %>% 
  dplyr::select(phos = name, prot)

###
# Prot
###
STRING.expmt.gene <- readRDS("data/STRINGexpmtgene_lowconf.rds")
prot <- STRING.expmt.gene[(STRING.expmt.gene$protein1 %in% prots$node) | (STRING.expmt.gene$protein2 %in% prots$node), c(2:3)] %>%
  na.omit() %>% 
  graph_from_data_frame(directed = F) %>% 
  igraph::simplify(remove.multiple = F,remove.loops = T) %>% 
  igraph::as_data_frame() %>% 
  dplyr::select(prot1 = 1, prot2 = 2)
prots <- c(prot$prot1, prot$prot2)

###
# Prot > Func
###
enrichedTerms <- enrichr(prots, databases = "GO_Biological_Process_2018") %>%
  .$GO_Biological_Process_2018 %>%
  filter(Adjusted.P.value < 0.01) %>%
  separate(Term,
           into = c("Term", "GOID"),
           sep = " \\(",
           extra = "drop") %>%
  mutate(GOID = sub("\\)",
                    "",
                    GOID)) %>%
  arrange(desc(Adjusted.P.value)) 

# More Terms = Weaker Analysis!!
keepID <- simplifyGO(GOID = enrichedTerms$GOID, simplifyData = simplifyGOReqData())

prot_func <- enrichedTerms %>% 
  filter(GOID %in% keepID) %>% 
  separate_rows(Genes, sep = ";") %>%
  dplyr::select(func = Term, prot = Genes, ID = GOID)

###
# Func
###
hsGO <- godata('org.Hs.eg.db', ont="BP")
semSim <- expand.grid(unique(prot_func$ID), unique(prot_func$ID), stringsAsFactors = F)
semSim$sim <- apply(semSim, 1, function(i) goSim(i[1], i[2], semData = hsGO, measure = "Wang"))
func <- semSim %>% filter(sim < 1 & sim > 0.7)
func$func1 <- (AnnotationDbi::select(GO.db,
                                     keys = func$Var1, columns = "TERM",
                                     keytype = "GOID"))$TERM
func$func2 <- (AnnotationDbi::select(GO.db,
                                     keys = func$Var2, columns = "TERM",
                                     keytype = "GOID"))$TERM
func <- dplyr::select(func, func1, func2) %>% na.omit %>% unique()

###
# Vertices
###
v <- rbind(
  data.frame(v = unique(phos_prot$phos), layer = "phos", stringsAsFactors = F),
  data.frame(v = unique(phos$phos1), layer = "phos", stringsAsFactors = F),
  data.frame(v = unique(phos$phos2), layer = "phos", stringsAsFactors = F),
  data.frame(v = unique(prot_func$prot), layer = "prot", stringsAsFactors = F),
  data.frame(v = unique(phos_prot$prot), layer = "prot", stringsAsFactors = F),
  data.frame(v = unique(prot$prot1), layer = "prot", stringsAsFactors = F),
  data.frame(v = unique(prot$prot2), layer = "prot", stringsAsFactors = F),
  data.frame(v = unique(prot_func$func), layer = "func", stringsAsFactors = F),
  data.frame(v = unique(func$func1), layer = "func",  stringsAsFactors = F),
  data.frame(v = unique(func$func2), layer = "func",  stringsAsFactors = F)
) %>% 
  unique()

edgelists <- list(x = phos, y = prot, z = func,
                  xy = phos_prot[,c("phos", "prot")], yx = phos_prot[,c("prot", "phos")],
                  yz = prot_func[,c("prot", "func")], zy = prot_func[,c("func", "prot")]
) %>% 
  lapply(dplyr::rename, "to" = 2, "from" = 1)

# Save heterogenous multilayer network
saveRDS(list(v = v,
             edgelists = edgelists,
             fcm = pp),
        "results/data/mlnw_model.rds")


seed_l <- lapply(1:max(pp$clustering), function(i){
  c(names(pp$clustering[pp$clustering == i]))
})

rwhn <- lapply(seed_l, function(s){
  calculateRWHN(edgelists = edgelists,
                verti = v,
                seeds = s,
                transitionProb = 0.7,
                restart = 0.7,
                weight_xy = 0.3,
                weight_yz = 0.7) %>%
    filter(name %in% v[v$layer=="func",]$v)
})

saveRDS(rwhn, "results/data/rwhn_model.rds")

dot <- dotplot_gg(rwhn)
ggsave(filename = "results/figs/rwhn_model_dotplot.tiff", 
       plot = dot,
       width = 209.804,
       height = 142,
       units = "mm")