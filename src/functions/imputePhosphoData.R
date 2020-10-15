filter_missing <- function(sty, allowed = 0, colnms){
  keep <- sty %>% 
    group_by(id) %>% 
    dplyr::select(id, grep(colnms, colnames(.))) %>% 
    pivot_longer(cols = -id, names_to = "experiment", values_to = "value") %>% 
    filter(sum(is.na(value)) <= allowed) %>% 
    dplyr::select(id) %>% 
    unique()
  
  flt <- filter(sty, id %in% keep$id)
  
  return(flt)
}

imputeQRLIC.mod <- function (dataSet.mvs, tune.sigma = 1) 
{
  nFeatures = dim(dataSet.mvs)[1]
  nSamples = dim(dataSet.mvs)[2]
  dataSet.imputed = dataSet.mvs
  QR.obj = list()
  for (i in 1:nSamples) {
    curr.sample = dataSet.mvs[, i]
    pNAs = length(which(is.na(curr.sample)))/length(curr.sample)
    upper.q = 0.99
    q.normal = qnorm(seq((pNAs + 0.001), (upper.q + 0.001), 
                         (upper.q - pNAs)/(upper.q * 100)), mean = 0, sd = 1)
    q.curr.sample = quantile(curr.sample, probs = seq(0.001, 
                                                      (upper.q + 0.001), 0.01), na.rm = T)
    temp.QR = lm(q.curr.sample ~ q.normal)
    QR.obj[[i]] = temp.QR
    mean.CDD = temp.QR$coefficients[1]
    sd.CDD = as.numeric(temp.QR$coefficients[2])
    data.to.imp = tmvtnorm::rtmvnorm(n = nFeatures, mean = mean.CDD, 
                           sigma = sd.CDD * tune.sigma, upper = qnorm((pNAs + 
                                                                         0.001), mean = mean.CDD, sd = sd.CDD), algorithm = c("gibbs"))
    curr.sample.imputed = curr.sample
    curr.sample.imputed[which(is.na(curr.sample))] = data.to.imp[which(is.na(curr.sample))]
    dataSet.imputed[, i] = curr.sample.imputed
  }
  return(dataSet.imputed)
}
