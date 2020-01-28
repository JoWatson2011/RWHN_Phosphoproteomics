constructHetNet <- function(stytxt, phosphoData, mfuzz_output, stringPath = "data/STRINGexpmtgene_lowconf.rds", modules = T, GOpval = 0.01){
  ## Phospho
  phos <- lapply(1:max(mfuzz_output$clustering), function(x){
    cl <- names(mfuzz_output$clustering[mfuzz_output$clustering == x])
    xpnd <- phosphoData %>% 
      filter(id %in% cl) %>% 
      tibble::column_to_rownames(var = "id") %>% 
      t() %>% 
      cor() %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column(var = "phos1") %>% 
      gather("phos2", "value", -phos1) %>% 
      filter(value < as.double(1.0), value > 0.99) %>% 
      unique
    return(xpnd) 
  }) %>% 
    do.call(rbind, .) %>% 
    graph_from_data_frame(directed = F) %>% 
    simplify() %>% 
    igraph::as_data_frame(what = "edges") %>% 
    dplyr::select(phos1 = from, phos2 = to)
  
  ## Phospho >> Prot
  phos_prot <- graph_from_data_frame(phos, directed = F) %>%
    igraph::as_data_frame("vertices") %>% 
    mutate(prot = apply(., 1, function(x) strsplit(x[1], "_")[[1]][1])) %>% 
    dplyr::select(phos = name, prot)
  
  ### INCLUDE EXTRA NODES OR DROP CONFIDENCE LEVEL???
  ## Prot
  prots <- unique(phos_prot$prot)
  STRING.expmt.gene <- readRDS(stringPath)
  prot <- STRING.expmt.gene[(STRING.expmt.gene$protein1 %in% prots) | (STRING.expmt.gene$protein2 %in% prots), c(2:3)] %>%
    na.omit %>% 
    graph_from_data_frame(directed = F) %>% 
    igraph::simplify(remove.multiple = F,remove.loops = T) %>% 
    igraph::as_data_frame() %>% 
    dplyr::select(prot1 = from, prot2 = to)
  
  ## Protein >> Func
  
  enrichedTerms <- lapply(1:max(mfuzz_output$clustering), function(i){
    ids <- names(mfuzz_output$clustering[mfuzz_output$clustering == i ] )
    cl_prots <- unique(stytxt[stytxt$id %in% ids,]$Gene_names)

    enrichedTerms <- if(modules){
      el <- STRING.expmt.gene %>% 
        filter(protein1 %in% cl_prots & protein2 %in% cl_prots) %>%
        dplyr::select(2:3) %>% 
        na.omit
      nw <- graph_from_data_frame(el, directed = F) %>% 
        igraph::simplify(remove.multiple = F,remove.loops = T) %>% 
        igraph::cluster_louvain()
      cluster_nw <- data.frame(gene.symbol = names(igraph::membership(nw)),
                               module = as.integer(igraph::membership(nw)),
                               stringsAsFactors = F)
      
      cluster_nw <- lapply(1:max(cluster_nw$module), function(x){
        sample <- cluster_nw %>%
          filter(module == x)
        
        enriched <- enrichr(sample$gene.symbol, databases = "GO_Biological_Process_2018") %>% 
          .$GO_Biological_Process_2018 %>%
          filter(Adjusted.P.value < GOpval) %>%
          separate(Term,
                   into = c("Term", "GOID"),
                   sep = " \\(",
                   extra = "drop") %>%
          mutate(GOID = sub("\\)",
                            "",
                            GOID))
        return(enriched)
      }) %>% 
        do.call(rbind, .)
    }else{
      enrichr(cl_prots, databases = "GO_Biological_Process_2018") %>% 
        .$GO_Biological_Process_2018 %>%
        filter(Adjusted.P.value < GOpval) %>%
        separate(Term,
                 into = c("Term", "GOID"),
                 sep = " \\(",
                 extra = "drop") %>%
        mutate(GOID = sub("\\)",
                          "",
                          GOID))
    }
  }) %>% do.call(rbind, .)
  
  keepID <- simplifyGO(GOID = enrichedTerms$GOID, simplifyData = simplifyGOReqData())
  
  prot_func <- enrichedTerms %>% 
    filter(GOID %in% keepID) %>% 
    separate_rows(Genes, sep = ";") %>%
    dplyr::select(func = Term, prot = Genes, ID = GOID)
  
  
  ## Func
  hsGO <- godata('org.Hs.eg.db', ont="BP")
  semSim <- expand.grid(prot_func$ID, prot_func$ID, stringsAsFactors = F)
  semSim$sim <- apply(semSim, 1, function(i) goSim(i[1], i[2], semData = hsGO, measure = "Wang"))
  func <- semSim %>% filter(sim < 1 & sim > 0.7)
  func$func1 <- (AnnotationDbi::select(GO.db,
                                       keys = func$Var1, columns = "TERM",
                                       keytype = "GOID"))$TERM
  func$func2 <- (AnnotationDbi::select(GO.db,
                                       keys = func$Var2, columns = "TERM",
                                       keytype = "GOID"))$TERM
  func <- dplyr::select(func, func1, func2) %>% unique() %>% na.omit()
  
  
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
  
  edgelists <- edgelists <- list(x = phos, y = prot, z = func,
                                 xy = phos_prot[,c("phos", "prot")], yx = phos_prot[,c("prot", "phos")],
                                 yz = prot_func[,c("prot", "func")], zy = prot_func[,c("func", "prot")]
  ) %>% 
    lapply(dplyr::rename, "to" = 2, "from" = 1)
  
  return(list(v = v,
             edgelists = edgelists)
  )
}
