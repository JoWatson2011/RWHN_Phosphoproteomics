library(tidyverse)
library(org.Hs.eg.db)
library(GO.db)

regsites <- read_tsv("data/Regulatory_sites.txt", skip = 2) %>% 
  separate_rows(ON_PROCESS, sep = "; ")

regsites_sum <- regsites %>% 
  mutate(parentProc = gsub(", .*", "", ON_PROCESS)) %>% 
  dplyr::select(ON_PROCESS, parentProc) %>% 
  unique()

parentProc_GOID <- tribble(
  ~GOID, ~parentProc,
  "GO:0006915", "apoptosis", # "Apoptotic process"
  "GO:0006914", "autophagy",
  # No term for carcinogenesis.
  "GO:0007155", "cell adhesion",
  "GO:0007049", "cell cycle regulation",  # "Cell cycle". Child term "regulation of cell cycle" may be more relevant.
  "GO:0030154", "cell differentiation",
  "GO:0016049", "cell growth",
  "GO:0048870", "cell motility",
  "GO:0006325", "chromatin organization",
  "GO:0007010", "cytoskeletal reorganization",
  "GO:0006281", "DNA repair",
  "GO:0006897", "endocytosis",
  "GO:0006887", "exocytosis",
  # No term for "neural plasticity",
  "GO:0008380", "RNA splicing",
  "GO:0043487", "RNA stability", # "Regulation of RNA stability"
  "GO:0007165", "signaling pathway regulation", # "signal transduction"
  "GO:0010467", "transcription", #"Gene expression". DNA templated/RNA templated in child terms
  "GO:0006412", "translation"
)

regsites_GO <- merge(regsites_sum, parentProc_GOID, by = "parentProc") %>% 
  merge(regsites, by = "ON_PROCESS", all.y = T) %>% 
  filter(!is.na(ON_PROCESS) &
           grepl("-p", MOD_RSD) &
           ORGANISM == "human"
  ) %>% 
  mutate(id_site = paste0(GENE, "_", gsub("-p", "", MOD_RSD)))


regsites_GO_offspring <- lapply(parentProc_GOID$GOID, function(i){
  offspring <- as.list(GOBPOFFSPRING)[[i]]
  
  # filtered <- geneEnrichment::simplifyGO(data.frame(GOID = offspring, p = 0.01),
  #                                      adjpcol = "p")
  
  df <- data.frame(GOID = i, 
                   # offspring = paste(offspring[offspring %in% filtered$GOID], 
                                     # collapse = ";")
                   offspring = paste(offspring, collapse = ";")
  )
  return(df)
}) %>% 
  do.call(rbind, .) %>% 
  merge(parentProc_GOID, by = "GOID") %>% 
  merge(regsites_sum, by = "parentProc") %>% 
  filter(!is.na(ON_PROCESS)
  )
  


write_tsv(regsites_GO_offspring, "data/Regulatory_sites_GOmapped.tsv")

# human sites in PSP with process annotation = 4013
# human genes = 1491

# of EGF data: 88 / 788 w associated function in PSP
# of TGF 79 / 727 w associated function in PSP


go_signaling <- as.list(GOBPOFFSPRING)[["GO:0007154"]]

signaling_prots <- AnnotationDbi::select(org.Hs.eg.db,
                      c("GO:0007154", go_signaling),
                      "SYMBOL", "GO") %>% 
  filter(SYMBOL %in% regsites_GO$GENE) %>% 
  merge(regsites_GO[,c("GENE", "MOD_RSD", "ON_PROCESS")],
        by.x = "SYMBOL", by.y = "GENE") %>% 
  dplyr::select(-c(EVIDENCE, ONTOLOGY)) %>% 
  unique()

go_cellCycle <- as.list(GOBPOFFSPRING)[["GO:0007049"]]
cellCycle_prots <- AnnotationDbi::select(org.Hs.eg.db,
                                         c("GO:0007049", go_cellCycle),
                                         "SYMBOL", "GO") %>% 
  filter(SYMBOL %in% regsites_GO$GENE) %>% 
  merge(regsites_GO[,c("GENE", "MOD_RSD", "ON_PROCESS")],
        by.x = "SYMBOL", by.y = "GENE") %>% 
  dplyr::select(-c(EVIDENCE, ONTOLOGY)) %>% 
  unique()

go_cellGrowth <- as.list(GOBPOFFSPRING)[["GO:0016049"]]
cellGrowth_prots <- AnnotationDbi::select(org.Hs.eg.db,
                                         c("GO:0016049", go_cellGrowth),
                                         "SYMBOL", "GO") %>% 
  filter(SYMBOL %in% regsites_GO$GENE) %>% 
  merge(regsites_GO[,c("GENE", "MOD_RSD", "ON_PROCESS")],
        by.x = "SYMBOL", by.y = "GENE") %>% 
  dplyr::select(-c(EVIDENCE, ONTOLOGY)) %>% 
  unique()


go_DNArepair <- as.list(GOBPOFFSPRING)[["GO:0006281"]]
DNArepair_prots <- AnnotationDbi::select(org.Hs.eg.db,
                                          c("GO:0006281", go_DNArepair),
                                          "SYMBOL", "GO") %>% 
  filter(SYMBOL %in% regsites_GO$GENE) %>% 
  merge(regsites_GO[,c("GENE", "MOD_RSD", "ON_PROCESS")],
        by.x = "SYMBOL", by.y = "GENE") %>% 
  dplyr::select(-c(EVIDENCE, ONTOLOGY))%>% 
  unique()

saveRDS(
  list(
    cellCycle = cellCycle_prots,
    cellGrowth = cellGrowth_prots,
    DNArepair = DNArepair_prots,
    signaling = signaling_prots
  ), 
  "data/PSP_Benchmarks.rds")


