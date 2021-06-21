### This script contains helper functions for model selection and summarizing
### rate estimates. 

# Old corSumm with standard error
corSumm <- function(results, rate.cat){
  # This function takes a list of results from a corHMM analysis and the number
  # of rate categories in the results. Creates a list of datatables. 
  # aic is a datatable of aic scores from each tree
  # aic.summ is a summary of aic scores with median, mean, and std error
  # par is datable of rate estimates and associated aic for each tree
  # par.melt is long version of par for easier plotting and summarising
  # par.summ is summary of pars with median, mean, and std error
  
  ### AIC scores ###
  aic <- unlist(sapply(results, "[", "AICc"))
  aic <- data.table(aic)
  # Summarise aic scores
  aic.summ <- aic[, .(med = median(aic), mean = mean(aic), 
                      std.err = sd(aic)/sqrt(.N))]
  
  ### Parameter esimates ###
  # Extract list of rate matrices, convert to vectors, and remove na's
  par <- sapply(results, "[", "solution")
  par <- lapply(par, as.vector)
  par <- lapply(par, na.omit)
  
  # Make vector of names for rates based on number of rate cats
  if(rate.cat == 1) 
  {par.names <- c("aic","Z.to.A", "A.to.Z")
  } else if(rate.cat == 2) 
  {par.names <- c("aic","ZS.to.AS", "AF.to.AS", 
                  "AS.to.ZS", "ZF.to.ZS", "AS.to.AF", 
                  "ZF.to.AF", "ZS.to.ZF", "AF.to.ZF")
  
  } else if(rate.cat == 3) 
  {par.names <- c("aic","ZS.to.AS", "AM.to.AS", "AS.to.ZS", "ZM.to.ZS", 
                  "AS.to.AM", "ZM.to.AM", "AF.to.AM", "ZS.to.ZM", 
                  "AM.to.ZM", "ZF.to.ZM", "AM.to.AF", "ZF.to.AF", 
                  "ZM.to.ZF", "AF.to.ZF")
  
  } else
  {
    par.names <- c("aic", "ZS.to.AS", "AMS.to.AS", "AS.to.ZS", "ZMS.to.ZS",
                   "AS.to.AMS", "ZMS.to.AMS", "AMF.to.AMS", "ZS.to.ZMS",
                   "AMS.to.ZMS", "ZMF.to.ZMS", "AMS.to.AMF", "ZMF.to.AMF",
                   "AF.to.AMF", "ZMS.to.ZMF", "AMF.to.ZMF", "ZF.to.ZMF",
                   "AMF.to.AF", "ZF.to.ZF", "ZMF.to.ZF", "AF.to.ZF")
  }
  # Turn list of rate vectors into list of 1 row datatables
  par <- lapply(par, function(x) setDT(as.list(x)))
  # Name each element of list by aic score
  names(par) <- aic[, aic]
  # Bind list into datatable with aic id col
  par <- rbindlist(par, idcol = TRUE)
  setnames(par, par.names)
  par[, aic := as.numeric(aic)]
  # Melt datatable for easier plotting and summary stats
  par.melt <- melt(par, id.vars = "aic", variable.name = "FromTo", value = "rate")
  # Summarise paramter estimates
  par.summ <- par.melt[, .(med = median(rate), mean = mean(rate), 
                           std.err = sd(rate)/sqrt(.N)), 
                       by = FromTo]
  
  # Create and return list of datatables
  results.list <- list(aic = aic, aic.summ = aic.summ, par = par, 
                       par.melt = par.melt, par.summ = par.summ)
  return(results.list)
}


bootMed <- function(x, n)
{
  bootmed <- apply(matrix(sample(x, rep=TRUE, n*length(x)), nrow=n), 1, median)
  bootmed <- data.table(median = bootmed)
  q <- quantile(bootmed[,median], c(0.025, 0.975))
  m <- median(x)
  # blist <- list(q)
  # return(blist)
  qm <- c(q[1], m, q[2])
  names(qm) <- c("2.5%", "median", "97.5%")
  return(qm)
}