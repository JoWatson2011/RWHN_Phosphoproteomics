#Read in STRING edges
STRING_flt <- function(path, conf){
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(tidyr))
  suppressPackageStartupMessages(library(biomaRt))
  
  STRING.expmt <- readr::read_table2(path) %>% 
  dplyr::select(1:2,7) %>%
  filter(experimental >= conf) %>% 
  tibble::rowid_to_column("ID")
  STRING.expmt$protein1 <- sub("9606.", "", STRING.expmt$protein1, fixed=T)
  STRING.expmt$protein2 <- sub("9606.", "", STRING.expmt$protein2, fixed=T)
  
  
  #ENSEMBL --> GENE name
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl" )
  genes1 <- unique(c(STRING.expmt$protein1, STRING.expmt$protein2))
  G_list <- getBM(filters = "ensembl_peptide_id", 
                  attributes= c("ensembl_peptide_id", "hgnc_symbol"),
                  values=genes1,
                  mart= ensembl)
  
  #How many duplicates? How many missing gene names? How many not matched
  G_list.duplicates <- G_list[duplicated(G_list$ensembl_peptide_id),]
  sum(G_list$ensembl_peptide_id %in% G_list.duplicates[,1]) # each dupl. occurs twice (most likely)
  G_list.unfound <- G_list[G_list$hgnc_symbol == "",1]
  G_list.unincluded <- genes1[!(genes1 %in% G_list$ensembl_peptide_id)]
  table(G_list.unfound %in% G_list.unincluded)  
  
  #match ensembl ids to gene symbols
  STRING.prot1 <- merge(STRING.expmt[,c(1:2,4)],G_list,by.x="protein1",by.y ="ensembl_peptide_id",all = T, sort = T)
  STRING.prot2 <- merge(STRING.expmt[,c(1,3:4)],G_list,by.x="protein2",by.y ="ensembl_peptide_id",all = T, sort = F)
  STRING.expmt.gene <- merge(STRING.prot1, STRING.prot2, by = "ID") %>% 
    filter(!duplicated(ID))
  
  #verify that it worked...
  nrow(STRING.expmt) == nrow(STRING.expmt.gene)
  identical(STRING.expmt$experimental, STRING.expmt.gene$experimental.x)
  identical(STRING.expmt$experimental, STRING.expmt.gene$experimental.y)
  
  #remove ensembl ids
  STRING.expmt.gene <- STRING.expmt.gene[,c(1,4,7,6)]
  colnames(STRING.expmt.gene) <- c("ID", "protein1","protein2","experimental")
  
  return(STRING.expmt.gene)
}
  
STRING_lowconf <- STRING_flt("data/9606.protein.links.detailed.v11.0.txt", 400)
saveRDS(STRING_lowconf, "data/STRINGexpmtgene_lowconf.rds")
STRING_highconf <- STRING_flt("data/9606.protein.links.detailed.v11.0.txt", 700)
saveRDS(STRING_highconf, "data/STRINGexpmtgene_highconf.rds")
