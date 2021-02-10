# Load packages
library(data.table)
library(phytools)
library(phangorn)
library(here)

# Suppress scientific notation
options(scipen = 999)

# Set working directory
setwd(here("analysis/stree_asr"))

# Define functions-------------------------------------------------------------

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
  if(rate.cat == 1) {par.names <- c("aic","Z.to.A", "A.to.Z")
  } else if(rate.cat == 2) {par.names <- c("aic","ZS.to.AS", "AF.to.AS", 
                                           "AS.to.ZS", "ZF.to.ZS", "AS.to.AF", 
                                           "ZF.to.AF", "ZS.to.ZF", "AF.to.ZF")
  
  } else {par.names <- c("aic","ZS.to.AS", "AM.to.AS", "AS.to.ZS", "ZM.to.ZS", 
                         "AS.to.AM", "ZM.to.AM", "AF.to.AM", "ZS.to.ZM", 
                         "AM.to.ZM", "ZF.to.ZM", "AM.to.AF", "ZF.to.AF", 
                         "ZM.to.ZF", "AF.to.ZF")
  
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

treePaint <- function(tree, anc, rate.cat)
{
  if (rate.cat == 1)
  {
    as.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,1]>0.75)]
    as.paint <- unlist(Descendants(tree, node = as.edges, type = "children"))
    
    zs.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,2]>0.75)]
    zs.paint <- unlist(Descendants(tree, node = zs.edges, type = "children"))
    
    plist <- list(as.paint, zs.paint)
    nvec <- c("A", "Z")
    
    for (i in 1:length(plist))
    {
      tryCatch(
        {
          tree <- paintBranches(tree, edge = plist[[i]], state = nvec[i])
        }, error = function(e)
        {
          message("Some regimes could not be printed on tree.")
        })
    }
    
  }
  else if (rate.cat == 2)
  {
    as.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,1]>0.75)]
    as.paint <- unlist(Descendants(tree, node = as.edges, type = "children"))
    
    zs.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,2]>0.75)]
    zs.paint <- unlist(Descendants(tree, node = zs.edges, type = "children"))
    
    af.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,3]>0.75)]
    af.paint <- unlist(Descendants(tree, node = af.edges, type = "children"))
    
    zf.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,4]>0.75)]
    zf.paint <- unlist(Descendants(tree, node = zf.edges, type = "children"))
    
    plist <- list(as.paint, zs.paint, af.paint, zf.paint)
    nvec <- c("AS", "ZS", "AF", "ZF")
    
    for (i in 1:length(plist))
    {
      tryCatch(
        {
          tree <- paintBranches(tree, edge = plist[[i]], state = nvec[i])
        }, error = function(e)
        {
          message("Some regimes could not be printed on tree.")
        })
    }
  } else
  {
    # Designate edges that should be painted in each state
    as.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,1]>0.75)]
    as.paint <- unlist(Descendants(tree, node = as.edges, type = "children"))
    
    zs.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,2]>0.75)]
    zs.paint <- unlist(Descendants(tree, node = zs.edges, type = "children"))
    
    am.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,3]>0.75)]
    am.paint <- unlist(Descendants(tree, node = am.edges, type = "children"))
    
    zm.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,4]>0.75)]
    zm.paint <- unlist(Descendants(tree, node = zm.edges, type = "children"))
    
    af.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,5]>0.75)]
    af.paint <- unlist(Descendants(tree, node = af.edges, type = "children"))
    
    zf.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,6]>0.75)]
    zf.paint <- unlist(Descendants(tree, node = zf.edges, type = "children"))
    
    plist <- list(as.paint, zs.paint, am.paint, zm.paint, af.paint, zf.paint)
    nvec <- c("AS", "ZS", "AM", "ZM", "AF", "ZF")
    
    for (i in 1:length(plist))
    {
      tryCatch(
        {
          tree <- paintBranches(tree, edge = plist[[i]], state = nvec[i])
        }, error = function(e)
        {
          message("Some regimes could not be printed on tree.")
        })
    }
  }
  return(tree)
}

azPaint <- function(tree, anc, rate.cat)
{
  ### This function takes a tree with branch lengths, and a set of ancestral
  ### state reconstructions made either with ancRECON or from summarizing
  ### across all trees in subsample. ASR should be a datatable with cols
  ### labeled with rate categories AS, ZS, etc. The function sums across 
  ### all zoox states and azoox states to plot an asr of just zoox vs azoox
  ### The function returns a tree with branches painted according to whether
  ### their ancestral node is zoox or azoox
  
  anc <- copy(anc)
  
  if (rate.cat == 1)
  {
    z.edges <- unique(tree$edge[,1])[which(anc[,Z]>0.75)]
    z.paint <- unlist(Descendants(tree, node = z.edges, type = "children"))
    
    a.edges <- unique(tree$edge[,1])[which(anc[,A]>0.75)]
    a.paint <- unlist(Descendants(tree, node = a.edges, type = "children"))
    
    tree <- paintBranches(tree, edge = z.paint, state = "Z")
    tree <- paintBranches(tree, edge = a.paint, state = "A")
    
    outlist <- list(tree = tree, nodeprobs = anc[, .(A = A, Z = Z)])
    
  } else
  {
    anc[, zsum := rowSums(.SD), .SDcols = grep("Z", names(anc))]
    anc[, asum := rowSums(.SD), .SDcols = grep("A", names(anc))]
    
    z.edges <- unique(tree$edge[,1])[which(anc[,zsum]>0.75)]
    z.paint <- unlist(Descendants(tree, node = z.edges, type = "children"))
    
    a.edges <- unique(tree$edge[,1])[which(anc[,asum]>0.75)]
    a.paint <- unlist(Descendants(tree, node = a.edges, type = "children"))
    
    tree <- paintBranches(tree, edge = z.paint, state = "Z")
    tree <- paintBranches(tree, edge = a.paint, state = "A")
    
    outlist <- list(tree = tree, nodeprobs = anc[, .(A = asum, Z = zsum)])
    
  }
  return(outlist)
}

fastPaint <- function(tree, anc, rate.cat)
{
  
  anc <- copy(anc)
  
  if (rate.cat == 2)
  {
    anc[, fsum := rowSums(.SD), .SDcols = c("AF", "ZF")]
    
    f.edges <- unique(tree$edge[,1])[which(anc[,fsum]>0.75)]
    f.paint <- unlist(Descendants(tree, node = f.edges, type = "children"))
    
    tree <- paintBranches(tree, edge = f.paint, state = "F")
  } else
  {
    anc[, fsum := rowSums(.SD), .SDcols = c("ZS", "AM", "AF", "ZF")]
    
    f.edges <- unique(tree$edge[,1])[which(anc[,fsum]>0.75)]
    f.paint <- unlist(Descendants(tree, node = f.edges, type = "children"))
    
    tree <- paintBranches(tree, edge = f.paint, state = "F")
  }
  
  return(tree)
  
}

maxPaint <- function(tree, anc, rate.cat)
{
  anc <- copy(anc)
  anc$lik.anc.states[,max :=  names(.SD)[max.col(.SD)]]
  if (rate.cat == 1)
  {
    
    as.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,max] == "A")]
    as.paint <- unlist(Descendants(tree, node = as.edges, type = "children"))
    
    zs.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,max] == "Z")]
    zs.paint <- unlist(Descendants(tree, node = zs.edges, type = "children"))
    
    plist <- list(as.paint, zs.paint)
    nvec <- c("A", "Z")
    
    for (i in 1:length(plist))
    {
      tryCatch(
        {
          tree <- paintBranches(tree, edge = plist[[i]], state = nvec[i])
        }, error = function(e)
        {
          message("Some regimes could not be printed on tree.")
        })
    }
    
  }
  else if (rate.cat == 2)
  {
    as.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,max] == "AS")]
    as.paint <- unlist(Descendants(tree, node = as.edges, type = "children"))
    
    zs.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,max] == "ZS")]
    zs.paint <- unlist(Descendants(tree, node = zs.edges, type = "children"))
    
    af.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,max] == "AF")]
    af.paint <- unlist(Descendants(tree, node = af.edges, type = "children"))
    
    zf.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,max] == "ZF")]
    zf.paint <- unlist(Descendants(tree, node = zf.edges, type = "children"))
    
    plist <- list(as.paint, zs.paint, af.paint, zf.paint)
    nvec <- c("AS", "ZS", "AF", "ZF")
    
    for (i in 1:length(plist))
    {
      tryCatch(
        {
          tree <- paintBranches(tree, edge = plist[[i]], state = nvec[i])
        }, error = function(e)
        {
          message("Some regimes could not be printed on tree.")
        })
    }
  } else
  {
    # Designate edges that should be painted in each state
    as.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,max] == "AS")]
    as.paint <- unlist(Descendants(tree, node = as.edges, type = "children"))
    
    zs.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,max] == "ZS")]
    zs.paint <- unlist(Descendants(tree, node = zs.edges, type = "children"))
    
    am.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,max] == "AM")]
    am.paint <- unlist(Descendants(tree, node = am.edges, type = "children"))
    
    zm.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,max] == "ZM")]
    zm.paint <- unlist(Descendants(tree, node = zm.edges, type = "children"))
    
    af.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,max] == "AF")]
    af.paint <- unlist(Descendants(tree, node = af.edges, type = "children"))
    
    zf.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,max] == "ZF")]
    zf.paint <- unlist(Descendants(tree, node = zf.edges, type = "children"))
    
    plist <- list(as.paint, zs.paint, am.paint, zm.paint, af.paint, zf.paint)
    nvec <- c("AS", "ZS", "AM", "ZM", "AF", "ZF")
    
    for (i in 1:length(plist))
    {
      tryCatch(
        {
          tree <- paintBranches(tree, edge = plist[[i]], state = nvec[i])
        }, error = function(e)
        {
          message("Some regimes could not be printed on tree.")
        })
    }
  }
  return(tree)
}

ptreePaint <- function(tree, anc, rate.cat, p)
{
  if (rate.cat == 1)
  {
    as.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,1]>p)]
    as.paint <- unlist(Descendants(tree, node = as.edges, type = "children"))
    
    zs.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,2]>p)]
    zs.paint <- unlist(Descendants(tree, node = zs.edges, type = "children"))
    
    plist <- list(as.paint, zs.paint)
    nvec <- c("A", "Z")
    
    for (i in 1:length(plist))
    {
      tryCatch(
        {
          tree <- paintBranches(tree, edge = plist[[i]], state = nvec[i])
        }, error = function(e)
        {
          message("Some regimes could not be printed on tree.")
        })
    }
    
  }
  else if (rate.cat == 2)
  {
    as.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,1]>p)]
    as.paint <- unlist(Descendants(tree, node = as.edges, type = "children"))
    
    zs.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,2]>p)]
    zs.paint <- unlist(Descendants(tree, node = zs.edges, type = "children"))
    
    af.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,3]>p)]
    af.paint <- unlist(Descendants(tree, node = af.edges, type = "children"))
    
    zf.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,4]>p)]
    zf.paint <- unlist(Descendants(tree, node = zf.edges, type = "children"))
    
    plist <- list(as.paint, zs.paint, af.paint, zf.paint)
    nvec <- c("AS", "ZS", "AF", "ZF")
    
    for (i in 1:length(plist))
    {
      tryCatch(
        {
          tree <- paintBranches(tree, edge = plist[[i]], state = nvec[i])
        }, error = function(e)
        {
          message("Some regimes could not be printed on tree.")
        })
    }
  } else
  {
    # Designate edges that should be painted in each state
    as.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,1]>p)]
    as.paint <- unlist(Descendants(tree, node = as.edges, type = "children"))
    
    zs.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,2]>p)]
    zs.paint <- unlist(Descendants(tree, node = zs.edges, type = "children"))
    
    am.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,3]>p)]
    am.paint <- unlist(Descendants(tree, node = am.edges, type = "children"))
    
    zm.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,4]>p)]
    zm.paint <- unlist(Descendants(tree, node = zm.edges, type = "children"))
    
    af.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,5]>p)]
    af.paint <- unlist(Descendants(tree, node = af.edges, type = "children"))
    
    zf.edges <- unique(tree$edge[,1])[which(anc$lik.anc.states[,6]>p)]
    zf.paint <- unlist(Descendants(tree, node = zf.edges, type = "children"))
    
    plist <- list(as.paint, zs.paint, am.paint, zm.paint, af.paint, zf.paint)
    nvec <- c("AS", "ZS", "AM", "ZM", "AF", "ZF")
    
    for (i in 1:length(plist))
    {
      tryCatch(
        {
          tree <- paintBranches(tree, edge = plist[[i]], state = nvec[i])
        }, error = function(e)
        {
          message("Some regimes could not be printed on tree.")
        })
    }
  }
  return(tree)
}


maxFast <- function(tree, anc, rate.cat)
{
  
  anc <- copy(anc)
  
  if (rate.cat == 2)
  {
    
    anc[, S := rowSums(.SD), .SDcols = c("AS", "ZS")]
    anc[, L := rowSums(.SD), .SDcols = c("AF", "ZF")]
    anc <- anc[, .(S, L)]
    anc[,max :=  names(.SD)[max.col(.SD)]]
    
    s.edges <- unique(tree$edge[,1])[which(anc[, max] == "S")]
    s.paint <- unlist(Descendants(tree, node = s.edges, type = "children"))
    
    l.edges <- unique(tree$edge[,1])[which(anc[, max] == "L")]
    l.paint <- unlist(Descendants(tree, node = l.edges, type = "children"))
    
    tree <- paintBranches(tree, edge = s.paint, state = "Stable")
    tree <- paintBranches(tree, edge = l.paint, state = "Labile")
    
  } else
  {
    
    anc[, S := rowSums(.SD), .SDcols = c("AS", "ZM")]
    anc[, L := rowSums(.SD), .SDcols = c("ZS", "AM")]
    anc[, V := rowSums(.SD), .SDcols = c("AF", "ZF")]
    anc <- anc[, .(S, L, V)]
    anc[,max :=  names(.SD)[max.col(.SD)]]
    
    s.edges <- unique(tree$edge[,1])[which(anc[, max] == "S")]
    s.paint <- unlist(Descendants(tree, node = s.edges, type = "children"))
    
    l.edges <- unique(tree$edge[,1])[which(anc[, max] == "L")]
    l.paint <- unlist(Descendants(tree, node = l.edges, type = "children"))
    
    v.edges <- unique(tree$edge[,1])[which(anc[, max] == "V")]
    v.paint <- unlist(Descendants(tree, node = v.edges, type = "children"))
    
    tree <- paintBranches(tree, edge = s.paint, state = "Stable")
    tree <- paintBranches(tree, edge = l.paint, state = "Labile")
    tree <- paintBranches(tree, edge = v.paint, state = "Volatile")
  }
  
  return(tree)
  
}

# Read in tree and traits------------------------------------------------------
stree <- mtree <- read.nexus(here("data/updated_trees_traits/stree_traits", 
                                  "stree.trees"))

traits <- fread(here("data/updated_trees_traits/stree_traits",
                     "stree_traits_B_as_Z.csv"),
                header = FALSE, col.names = c("taxa", "state"))


# Format traits and tip labels-------------------------------------------------
# Drop taxa without data from tree
stree <- lapply(stree, drop.tip, tip = traits[state == "-", taxa])
class(stree) <- "multiPhylo"

# Remove taxa with missing data from traits
traits <- traits[state != "-"]


# Sample 100 trees used to estimate rates--------------------------------------
set.seed(1)
stree <- sample(stree, 100, replace = FALSE)


# Read in corHMM runs----------------------------------------------------------
# Read in asr summarized across all 100 trees
f.anc <- readRDS(here("analysis/stree_asr", "stree_subsample_nodeframes.rds"))

# Read in 100 individual asrs
anc3 <- readRDS(here("analysis/stree_asr", "stree_3rate_anc.rds"))

# Plot ancestral state reconstructions-----------------------------------------
# Make consensus tree
ctree <- consensus(stree, p = 0.95)
# Compute banch lengths of ctree for mapping states. The actual branch lengths
# don't matter at this point. They just need values or treePaint won't work. 
ctree <- compute.brlen(ctree)

### Lengthen some short branches for presentation
# Plot with edge labels to see which to change
plotTree(ctree, ftype = "off", type = "fan")
edgelabels(frame = "n")

# Root
ctree$edge.length[1] <- ctree$edge.length[1] + 0.3
ctree$edge.length[2122] <- ctree$edge.length[2122] + 0.3

### Robusta
ctree$edge.length[1041] <- ctree$edge.length[1041] + 0.01
ctree$edge.length[1040] <- ctree$edge.length[1040] - 0.01
ctree$edge.length[2121] <- ctree$edge.length[2121] + 0.01

ctree$edge.length[1043] <- ctree$edge.length[1043] + 0.01
ctree$edge.length[1042] <- ctree$edge.length[1042] - 0.01
ctree$edge.length[2091] <- ctree$edge.length[2091] + 0.01

ctree$edge.length[1050] <- ctree$edge.length[1050] + 0.01
ctree$edge.length[1049] <- ctree$edge.length[1049] - 0.01
ctree$edge.length[1480] <- ctree$edge.length[1480] + 0.01

ctree$edge.length[1052] <- ctree$edge.length[1052] + 0.01
ctree$edge.length[1051] <- ctree$edge.length[1051] - 0.01
ctree$edge.length[1425] <- ctree$edge.length[1425] + 0.01

ctree$edge.length[1053] <- ctree$edge.length[1053] + 0.005
ctree$edge.length[1052] <- ctree$edge.length[1052] - 0.005
ctree$edge.length[1424] <- ctree$edge.length[1424] + 0.005

ctree$edge.length[1046] <- ctree$edge.length[1046] + 0.005
ctree$edge.length[1045] <- ctree$edge.length[1045] - 0.005
ctree$edge.length[1865] <- ctree$edge.length[1865] + 0.005

ctree$edge.length[1788] <- ctree$edge.length[1788] + 0.01
ctree$edge.length[1787] <- ctree$edge.length[1787] - 0.01
ctree$edge.length[1858] <- ctree$edge.length[1858] + 0.01

### Complexa
ctree$edge.length[3] <- ctree$edge.length[3] + 0.01
ctree$edge.length[2] <- ctree$edge.length[2] - 0.01
ctree$edge.length[1039] <- ctree$edge.length[1039] + 0.01


ctree$edge.length[561] <- ctree$edge.length[561] + 0.01
ctree$edge.length[560] <- ctree$edge.length[560] - 0.01
ctree$edge.length[1036] <- ctree$edge.length[1036] + 0.01


ctree$edge.length[562] <- ctree$edge.length[562] + 0.01
ctree$edge.length[560] <- ctree$edge.length[560] - 0.01
ctree$edge.length[1036] <- ctree$edge.length[1036] + 0.01
ctree$edge.length[1023] <- ctree$edge.length[1023] + 0.01

ctree$edge.length[957] <- ctree$edge.length[957] + 0.01
ctree$edge.length[956] <- ctree$edge.length[956] - 0.01
ctree$edge.length[1016] <- ctree$edge.length[1016] + 0.01

ctree$edge.length[958] <- ctree$edge.length[958] + 0.01
ctree$edge.length[956] <- ctree$edge.length[956] - 0.01
ctree$edge.length[1015] <- ctree$edge.length[1015] + 0.01
ctree$edge.length[1016] <- ctree$edge.length[1016] + 0.01


# Make color vectors for plotting combined rate cats
# 1 rate
r1.cols <- c("grey", "#084594", "#E6AB02")
names(r1.cols) <- c(1, "A", "Z")

# f3.cols <- c("grey", "grey", "red", "red", "grey", "red", "red")
# names(f3.cols) <- c(1, "AS", "ZS", "AM", "ZM", "AF", "ZF")

# Fast categories
# rc.cols <- c("grey", "red")
# names(rc.cols) <- c(1, "F")
fast.cols <- c("#bdbdbd", "#fc9272", "#a50f15")
names(fast.cols) <- c("Stable", "Labile", "Volatile")


# Regular colors for plotting all rate cats
r3.cols <- c("grey", "#084594", "#F781BF", "#1B9E77", "#E6AB02", "#7570B3", "#E41A1C")
names(r3.cols) <- c(1, "AS", "ZS", "AM", "ZM", "AF", "ZF")

# Make vector of trait colors
trait.cols <- traits[, state]
names(trait.cols) <- traits[, taxa]
trait.cols[trait.cols == "A"] <- "#084594"
trait.cols[trait.cols == "Z"] <- "#E6AB02"



### Format f.anc to look like asrs generated from ancRECON so I can use the 
### treePaint function I wrote
fanc3 <- f.anc[[3]]
setnames(fanc3, c("AS", "ZS", "AM", "ZM", "AF", "ZF", "ctree.nodes"))
fanc3 <- fanc3[, .(AS = mean(AS), ZS = mean(ZS), 
                   AM = mean(AM), ZM = mean(ZM), 
                   AF = mean(AF), ZF = mean(ZF)), 
               by = ctree.nodes]
fanc3[, ctree.nodes := NULL]
fanc3 <- list(fanc3)
names(fanc3) <- "lik.anc.states"


### Paint trees with azoox/zoox recons
og3.asr <- azPaint(ctree, anc = fanc3$lik.anc.states, rate.cat = 3)

### Paint summed fast categories on trees
rc3.asr <- maxFast(ctree, fanc3$lik.anc.states, rate.cat = 3)

### Paint trees with rate categories
# 2 rate
f3.asr <- treePaint(tree = ctree, anc = fanc3, rate.cat = 3)

### Paint tree with rate cats based on max p at each node
m3.asr <- maxPaint(tree = ctree, anc = fanc3, rate.cat = 3)

### Paint tree with different cutoff values
f5.asr <- ptreePaint(tree = ctree, anc = fanc3, rate.cat = 3, p = 0.5)

### Plot-----------------------------------------------------------------------

# Plot binary asr
pdf(file = "binary_asr.pdf", width = 5, height = 5)
plotSimmap(og3.asr$tree, type = "fan", ftype = "off", r1.cols, lwd = 0.5)
dev.off()

# Plot binary with names for reference 
pdf(file = "binary_reference.pdf", width = 60.84, height = 40.63)
plotSimmap(og3.asr$tree, type = "fan", fsize = 0.3, r1.cols, lwd = 0.5)
nodelabels(pie = as.matrix(og3.asr$nodeprobs), cex = 0.05,
           piecol = r1.cols[-1])
tiplabels(pch = 19, cex = 0.5, col = trait.cols[ctree$tip.label])
dev.off()

# Plot summed fast cats
pdf(file = "fast_cats.pdf", width = 5, height = 5)
plotSimmap(rc3.asr, type = "fan", ftype = "off", fast.cols, lwd = 0.5)
dev.off()

# Plot all 4 categories
pdf(file = "all4_asr.pdf", width = 5, height = 5)
plotSimmap(f3.asr, type = "fan", ftype = "off", r3.cols, lwd = 0.5)
dev.off()

pdf(file = "all4_maxp_asr.pdf", width = 5, height = 5)
plotSimmap(m3.asr, type = "fan", ftype = "off", r3.cols, lwd = 0.5)
dev.off()

pdf(file = "all4_maxp_reference.pdf", width = 60.84, height = 40.63)
plotSimmap(m3.asr, type = "fan", fsize = 0.3, r3.cols, lwd = 1)
nodelabels(pie = as.matrix(fanc3$lik.anc.states), cex = 0.05,
           piecol = r3.cols[-1])
tiplabels(pch = 19, cex = 0.5, col = trait.cols[ctree$tip.label])
dev.off()


# Plot individual asrs---------------------------------------------------------

anc3 <- lapply(anc3, "[[", "lik.anc.states")
anc3 <- lapply(anc3, data.table)
anc3 <- lapply(anc3, setnames, c("AS", "ZS", "AM", "ZM", "AF", "ZF"))
#anc3 <- lapply(anc3, list)
names(anc3) <- rep("lik.anc.states", length(anc3))


all.paint <- list()
for (i in 1:length(stree)){
  all.paint[[i]] <- maxPaint(stree[[i]], anc3[i], rate.cat = 3)
}



# pdf(file = "3rate_asr1.pdf", width = 50, height = 50)
# layout(matrix(1:100, 5, 20, byrow = TRUE))
# lapply(all.paint, plotSimmap, ftype = "off", r3.cols)
# dev.off()

vec <- c(1:100)
pdf(file = "3rate_asr1.pdf", width = 50, height = 50)
layout(matrix(1:100, 5, 20, byrow = TRUE))
for(i in 1:length(all.paint)){
  plotSimmap(all.paint[[i]], ftype = "off", r3.cols, mar = c(0.1, 0.1, 1.5, 0.1))
  mtext(vec[i])
}
dev.off()




