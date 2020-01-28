# From Jamie Soul; adapted from jspaezp on GitHub
fuzzyC <- function(tidydata, clus){
  
  set.seed(1)
  inputdata <- tidydata %>% 
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(var = "id") # rownames must be id - see mfuzz.ggplot function
  phospep <- new("ExpressionSet", 
                 exprs = as.matrix(inputdata) )
  phospep.z <- standardise(phospep)      #Avg expression of each peptide = 0, sd = 1
  optimalM <- mestimate(phospep.z)   #set fuzzifier
  cl <- mfuzz(phospep.z, c = clus, m = optimalM)
  a <- mFuzz.ggplot(exprs(phospep.z), cl, returnData = T)
  a$g <- a$g + ggtitle(substr(a$exp$sample[1], 1, 6), clus)
  
  return(list(data = phospep.z, exp = a$exp, g = a$g, clustering = cl$cluster))
}


mFuzz.ggplot <- function(data, clustering,
                         centre = F,returnData=T,sort.columns=T) {
  
  suppressPackageStartupMessages(library("Biobase"))
  suppressPackageStartupMessages(library("Mfuzz"))
  suppressPackageStartupMessages(library("tidyr"))
  suppressPackageStartupMessages(library("tidyselect"))
  suppressPackageStartupMessages(library("dplyr"))
  suppressPackageStartupMessages(library("ggplot2"))
  
  
  clusterindex <- clustering$cluster
  
  # data frame with Membership values
  memship <- clustering$membership 
  colnames(memship) <- paste("membership", 
                             seq_along(memship[1,]), 
                             sep = ("")) 
  
  exp <-data
  
  # This chunk replaces col names by numbers if 
  # more than 1 is character only 
  # or when sort.columns is FALSE
  
  all.char.cols <- !grepl("\\d", colnames(exp))
  if ((sum(all.char.cols) > 1) | !sort.columns) {
    colnames(exp) <- seq_along(all.char.cols)    
  }
  
  exp <- data.frame(exp , 
                    Identifier = rownames(data),
                    clusterindex, memship,
                    stringsAsFactors = F) 
  
  # Transform data frame into a ggplot-compatible format
  exp <- exp %>% 
    gather(sample, 
           expression ,
           - Identifier,
           - clusterindex,
           - contains("membership")) %>% 
    mutate(Time = gsub("(\\w*\\D+(?=([0-9]+)))|((?<=\\d)\\D+$)", 
                       "", 
                       sample,
                       perl = TRUE)) %>%
    #  this regular expression deletes all characters and numbers prior to 
    #  the last number in the string z.b. AA00AA00__00 -> 00 else keeps the string
    mutate(Time = gsub("^\\D*$", # this needs to be fixed, bug when seveal character cols ...
                       "0", 
                       Time,
                       perl = TRUE)) %>%
    mutate(Time = factor(Time, ordered=T))
  
  exp$maxMembership <- exp %>%  
    dplyr::select(contains("membership")) %>%
    apply(., 1, max) 
  
  g <- ggplot(exp,
              aes(x = Time, y = expression)) +
    geom_line(
      aes(group = Identifier,  
          colour = maxMembership)
    ) + 
    scale_colour_gradient(low = "white", high = "red")
  
  # Center plotting when centre == TRUE
  if (centre) {
    centers <- clustering$centers %>% 
      data.frame(., 
                 clusterindex = rownames(.),
                 stringsAsFactors = F) %>% 
      gather(sample, 
             Centre,
             - clusterindex) %>% 
      mutate(Time = gsub("(\\w*\\D+(?=([0-9]+)))|((?<=\\d)\\D+$)", 
                         "", 
                         sample,
                         perl = TRUE)) %>%
      #  this regular expression deletes all characters and numbers prior to 
      #  the last number in the string z.b. AA00AA00__00 -> 00 else keeps the string
      mutate(Time = gsub("^\\D*$", # this needs to be fixed, bug when all character names
                         "0", 
                         Time,
                         perl = TRUE)) %>%
      mutate(Time = factor(Time, ordered=T))
    g <- g + geom_line(data = centers, aes(x = Time, y = Centre))
  }
  
  g <- g + facet_grid(. ~ clusterindex)
  g <- g + geom_hline(yintercept = 0, 
                      colour="#666666", alpha = 0.5) + 
    ylab("Fold change") + labs(fill= "Membership")
  
  return(list(g=g,exp=exp))
}