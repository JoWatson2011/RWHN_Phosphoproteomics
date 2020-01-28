# Remove semantically similar GO terms.
# Based on code by Jamie Soul:
# https://github.com/soulj/SkeletalVis-Pipeline/blob/master/src/SYBIL%20Systems%20Biology/GOEnrichment.R

simplifyGO <- function(GOID, simplifyData){
  suppressPackageStartupMessages(library("dplyr"))
  suppressPackageStartupMessages(library("tidyr"))
  suppressPackageStartupMessages(library("GOSemSim"))
  
  
  semSim <- expand.grid(unique(GOID), unique(GOID), stringsAsFactors = F)
  semSim$sim <- apply(semSim, 1, function(i) goSim(i[1], i[2], semData = simplifyData$semDat, measure = "Wang"))
  
  semSim[is.na(semSim)] <- 0
  semSim <- semSim[!is.na(semSim$sim),]
  semSim <- semSim[ order(semSim$sim ,decreasing = T),]
  
  #get the simiar term pairs
  semSim <- semSim[semSim$Var1 != semSim$Var2,]
  semSim <- semSim[semSim$sim > 0.4,]
  
  
  #mark high fequency terms
  semSim$remove <- apply(semSim,1,function(x) {
    if (x[1] %in% simplifyData$highFreqTerms){
      return(x[1])
    }
    if (x[2] %in% simplifyData$highFreqTerms){
      return(x[2])
    } else {
      return(NA)
    }
  })
  remove <- na.omit(semSim$remove)
  semSim <- semSim[is.na(semSim$remove),]
  
  if(nrow(semSim) == 0){
    NA
  }else{
    for (i in 1:nrow(semSim)){
      Var1 <- semSim[i,"Var1"]
      Var2 <- semSim[i,"Var2"]
      if(Var2 %in% simplifyData$childTerms[[Var1]]){
        remove <- c(remove, Var2)
        next
      } else if (Var1 %in% simplifyData$childTerms[[Var2]])
        remove <- c(remove, Var1)
      next
    }
  }
    
    flt <- unique(GOID[!(GOID %in% remove)])
    
    return(flt)
}