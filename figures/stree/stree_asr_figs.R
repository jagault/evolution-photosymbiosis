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

treePaint <- function(tree, anc, rate.cat)
{
  ### This function takes a tree with branch lengths, and a set of ancestral
  ### state reconstructions made either with ancRECON or from summarizing
  ### across all trees in subsample. ASR should be a datatable with cols
  ### labeled with rate categories AS, ZS, etc. 
  ### The function returns a tree with branches painted according to the most 
  ### likely state at their ancestor node. Unlike the maxPaint function, this 
  ### uses a cutoff. If the probability of any given state is less than 0.75,
  ### the branch is painted grey to denote uncertainty. 
  
  ### tree is a tree of class phylo with branch lengths
  ### anc is a datatable that lists the probability at each internal node of 
  ### each state. 
  ### rate.cat specifies the number of rate categories that were fit to the tree
  ### to estimate the rates used for the ancestral state reconstruction
  
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

maxPaint <- function(tree, anc, rate.cat)
{
  
  ### This function takes a tree with branch lengths, and a set of ancestral
  ### state reconstructions made either with ancRECON or from summarizing
  ### across all trees in subsample. ASR should be a datatable with cols
  ### labeled with rate categories AS, ZS, etc. 
  ### The function returns a tree with branches painted according to the most 
  ### likely state at their ancestor node. 
  
  ### tree is a tree of class phylo with branch lengths
  ### anc is a datatable that lists the probability at each internal node of 
  ### each state. 
  ### rate.cat specifies the number of rate categories that were fit to the tree
  ### to estimate the rates used for the ancestral state reconstruction
  
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


### Paint trees with rate categories
# 2 rate
f3.asr <- treePaint(tree = ctree, anc = fanc3, rate.cat = 3)

### Paint tree with rate cats based on max p at each node
m3.asr <- maxPaint(tree = ctree, anc = fanc3, rate.cat = 3)


### Plot-----------------------------------------------------------------------

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

# Paint all 100 trees
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




