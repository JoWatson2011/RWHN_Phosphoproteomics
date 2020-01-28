calculateRWHN <- function(edgelists, verti, seeds, transitionProb, 
                          restart, weight_xy, weight_yz, eps = 1/10^6,
                          random = F){
  
  heterogenousNetwork <- lapply(edgelists, function(i){
    nw <- i %>% 
      graph_from_data_frame(directed = F) %>% simplify
    bimap <- bipartite.mapping(nw)
    
    if(bimap[[1]] == T){
      if(random){
        incd <- matrix(data = sample.int(2,
                                        length(unique(i$from))*length(unique(i$to)),
                                        TRUE), 
                      nrow = length(unique(i$from)), 
                      ncol = length(unique(i$to)), 
                      dimnames = list(unique(i$from), unique(i$to))
        )
        randomnw <- graph_from_incidence_matrix(incd)
        bimap_random <- bipartite.mapping(randomnw)
        mat <- get.incidence(randomnw, types= bimap_random$type, attr=NULL, names=TRUE, sparse=T)
      }else{
        mat <- get.incidence(nw, types= bimap$type, attr=NULL, names=TRUE, sparse=T)
      }
    }else{
      if(random){
        names <- unique(c(i$from, i$to))
        adj <- matrix(data = sample.int(2,
                                         length(names)*length(names),
                                         TRUE), 
                       nrow = length(names), 
                       ncol = length(names), 
                       dimnames = list(names, names)
        )
        randomnw <- graph_from_adjacency_matrix(adj)
        mat <- get.adjacency(randomnw, sparse = T)
      }else{
        mat <- get.adjacency(nw, sparse = T, )
      }
    }
    return(as.matrix(mat)) 
  })
  
  verti$transition <- ifelse(verti$layer == "phos", "x", 
                             ifelse(verti$layer == "prot", "y",
                                    ifelse(verti$layer == "func", "z", NA)))
  
  transitionMat <- calculateTransitionMatrix(trans = transitionProb, 
                                             vert = verti,
                                             hetNet = heterogenousNetwork)
  
  rwhn <- getRWHN(transMat = transitionMat,
              gamma = restart,
              eta_xy = weight_xy,
              eta_yz = weight_yz,
              seeds = seeds,
              eps = eps) #%>% 
     #filter(name %in% verti[verti$layer == "func",]$v) %>% 
     #mutate(funcRank = 1:nrow(.))

  return(rwhn)
}

calculateTransitionMatrix <- function(trans, hetNet, vert){
  
  intra_tm <- function(intra, vCols, vRows, trans){
    
    # Which nodes have degree > 0 in bipartite network?
    present <- intra[rowSums(intra) > 0,]
    
    # eq1
    present_tp <- (trans * present) / rowSums(present)
    
    # rbind with nodes of degree = 0 in bipartite network
    tp <- rbind(present, intra[rowSums(intra)==0,])
    
    # return(tp)
    
    # add rows for nodes in inter but not intra networks
    rows <- vRows[!(vRows %in% rownames(intra))]
    cols <- vCols[!(vCols %in% colnames(intra))]
    
    fin_tp <- matrix(0,
                     nrow = (nrow(tp) + length(rows)),
                     ncol = (ncol(tp) + length(cols))
    )
    fin_tp[1:nrow(tp), 1:ncol(tp)] <- tp
    
    colnames(fin_tp) <- c(colnames(tp), cols)
    rownames(fin_tp) <- c(rownames(tp), rows)
    
    return(fin_tp)
  }
  
  xy_tp <- intra_tm(intra = hetNet[["xy"]], 
                    vCols =  vert[vert$transition == "y",]$v, 
                    vRows = vert[vert$transition == "x",]$v,
                    trans = trans)
  yx_tp <- t(xy_tp)
  yz_tp <- intra_tm(intra = hetNet[["yz"]], 
                    vCols =  vert[vert$transition == "z",]$v, 
                    vRows = vert[vert$transition == "y",]$v,
                    trans = trans)
  zy_tp <- t(yz_tp)
  
  inter_tm <- function(inter, intra_names, trans){
    ### Which nodes in inter have degree 0 in intra ?
    
    #notHetNames <- colnames(inter)[colnames(inter) %in% colnames(intra_tp[ , colSums(intra_tp) == 0])]
    #notHet <- inter[rownames(inter) %in% notHetNames, colnames(inter) %in% notHetNames]
    
    notHetNames <- colnames(inter)[!(colnames(inter) %in% intra_names)]
    notHet <- inter[rownames(inter) %in% notHetNames, colnames(inter) %in% notHetNames]
    
    # eq2 when (B)i,j = 0
    if(class(notHet) %in% c("data.frame", "matrix")){
      notHet_tp <- (notHet)/rowSums(notHet)
    }else{
      notHet_tp <- data.frame()
    }
    
    # Which nodes are in inter and intra?
    hetNames <- colnames(inter)[colnames(inter) %in% intra_names]
    het <- inter[rownames(inter) %in% hetNames, colnames(inter) %in% hetNames]
    
    # eq2 when (B)i,j > 0
    if(class(het) %in% c("data.frame", "matrix")){
      het_tp <- ((1-trans)*het)/rowSums(het)
    }else if(length(hetNames) != 0){
      het_tp <- data.frame()
      for(i in 1:length(hetNames)){
        het_tp[1,i] <- 1-trans
        colnames(het_tp)[i] <- rownames(het_tp)[i] <- hetNames[i]
      }
    } else {
      het_tp <- data.frame()
    }
    
    # which nodes aren't in inter but are in intra?
    onlyHetNames <- intra_names[!(intra_names %in% colnames(inter))]
    onlyHet <- matrix(rep.int(0, (length(onlyHetNames)^2)),
                      nrow = length(onlyHetNames),
                      ncol = (length(onlyHetNames)),
                      dimnames = list(onlyHetNames, onlyHetNames))
    
    # combine
    tp <- data.frame()
    
    if(ncol(notHet_tp) != 0 & ncol(het_tp) != 0){
      row <- nrow(notHet_tp)
      tp[1:row, 1:row] <- notHet_tp
      row <- row + 1
      
      tp[row:(row + nrow(het_tp) - 1), row:(row + nrow(het_tp) - 1)] <- het_tp
      row <- row + nrow(het_tp)
      
      if(ncol(onlyHet) != 0){
        tp[row:(row + nrow(onlyHet) - 1), row:(row + nrow(onlyHet) - 1)] <- onlyHet
        
        colnames(tp) <- c(colnames(notHet_tp), colnames(het_tp), colnames(onlyHet))
        rownames(tp) <- c(rownames(notHet_tp), rownames(het_tp), rownames(onlyHet))
        
      }else{
        colnames(tp) <- c(colnames(notHet_tp), colnames(het_tp))
        rownames(tp) <- c(rownames(notHet_tp), rownames(het_tp))
      }
    }else if(ncol(notHet_tp) == 0 & ncol(het_tp) != 0){
      row <- nrow(het_tp)
      
      tp[1:row, 1:row] <- het_tp
      row <- row + 1
      
      if(ncol(onlyHet) != 0){
        tp[row:(row + nrow(onlyHet) - 1), row:(row + nrow(onlyHet) - 1)] <- onlyHet
        colnames(tp) <- c(colnames(het_tp), colnames(onlyHet))
        rownames(tp) <- c(rownames(het_tp), rownames(onlyHet))
      } else {
        colnames(tp) <- colnames(het_tp)
        rownames(tp) <- rownames(het_tp)
      }
    }else if(ncol(notHet_tp) != 0 & ncol(het_tp) == 0){
      row <- nrow(notHet_tp)
      
      tp[1:row, 1:row] <- notHet_tp
      row <- row + 1
      
      if(ncol(onlyHet) != 0){
        tp[row:(row + nrow(onlyHet) - 1), row:(row + nrow(onlyHet) - 1)] <- onlyHet
        colnames(tp) <- c(colnames(notHet_tp), colnames(onlyHet))
        rownames(tp) <- c(rownames(notHet_tp), rownames(onlyHet))
      } else {
        colnames(tp) <- colnames(notHet_tp)
        rownames(tp) <- rownames(notHet_tp)
      }    
    }else{
      print("Situation not accounted for")
    }
    
    tp <- as.matrix(tp)
    tp[is.na(tp)] <- 0
    tp <- tp[order(match(rownames(tp), intra_names)), order(match(colnames(tp), intra_names))]
    
    return(tp)
  }
  
  x_tp <- inter_tm(inter = hetNet[["x"]], 
                   intra_names = colnames(hetNet[["yx"]]), 
                   trans = trans)
  
  y_tp <- inter_tm(inter = hetNet[["y"]], 
                   intra_names = unique(c(colnames(hetNet[["xy"]]), colnames(hetNet[["zy"]]))),
                   trans = trans)
  
  z_tp <- inter_tm(inter = hetNet[["z"]], 
                   intra_names = colnames(hetNet[["yz"]]),
                   trans = trans)
  
  # # create the full transition matrix
  # # |siteMatrix (x)   Site2Prot (yx)         0         |
  # # |Prot2Site (xy)   protMatrix (y)   prot2func (zy)  |
  # # |     0           func2prot (yz)   funcMatrix(z)   | 
  
  top0 <- matrix(0, 
                 nrow = nrow(x_tp),
                 ncol = ncol(z_tp),
                 dimnames = list(rownames(x_tp), colnames(z_tp)))
  btm0 <- matrix(0, 
                 nrow = nrow(z_tp),
                 ncol = ncol(x_tp),
                 dimnames = list(rownames(z_tp), colnames(x_tp)))
  
  # tmp1 <- cbind(x_tp, xy_tp[rownames(x_tp), colnames(y_tp), drop = F], top0)
  # tmp2 <- cbind(yx_tp[rownames(y_tp), colnames(x_tp), drop = F], y_tp, yz_tp[rownames(y_tp), colnames(z_tp), drop = F])
  # tmp3 <- cbind(btm0, zy_tp[rownames(z_tp), colnames(y_tp), drop = F], z_tp)
  tmp1 <- cbind(x_tp, xy_tp, top0)
  tmp2 <- cbind(yx_tp, y_tp, yz_tp)
  tmp3 <- cbind(btm0, zy_tp, z_tp)
  
  M1 <- rbind(tmp1, tmp2, tmp3)
  
  return(list(M1 = M1, x_tp = x_tp, y_tp = y_tp, z_tp = z_tp))
}


getRWHN <- function(transMat, gamma, eta_xy, eta_yz, seeds , eps){
  M1 <- transMat[["M1"]]
  x <- transMat[["x_tp"]]
  y <- transMat[["y_tp"]]
  z <- transMat[["z_tp"]]
  
  seedScores <- as.data.frame(rownames(x))
  seedScores$Scores <- ifelse(seedScores[,1] %in% seeds, 1/sum(seedScores[,1] %in% seeds), 0)
  
  probabilityVector <- c(seedScores$Score, rep.int((1/nrow(M1) * eta_xy), nrow(y)), rep.int((1/nrow(M1) * eta_yz), nrow(z)))
  #probabilityVector <- c(seedScores$Score, rep.int(0, nrow(y)), rep.int(0, nrow(z)))
  
  
  iter <- 0
  p0 <- probabilityVector
  p1 <- c(rep(0,length(probabilityVector)))
  
  while (sum(abs(p0 - p1 )) > eps){
    p0 <- p1
    p1 <- ((1- gamma) * t(M1)) %*% p0 + gamma * probabilityVector
    p1 <- p1 / sum(p1)
    
    if(!any(is.na(p1))){
      iter <- iter + 1
    }else{
      p1 <- p0
      break
    }
  } 
  
  pi <- data.frame(V1 = p1[,1],
                   name = c(rownames(x), rownames(y), rownames(z)),
                   stringsAsFactors = F)
  pi <- pi[order(pi[,1], decreasing = T),]
  # pi <- pi[pi$name %in% rownames(z),]
  pi$rank <- 1:nrow(pi)
  
  return(pi)
}
