#Clear workspace
rm(list = ls())

# Set working directory
setwd("~/zoox/results/2019-04-05/")

# Load packages
library(phytools)
library(data.table)
library(parallel)


# Define functions-------------------------------------------------------------
# ctreeAnc takes a multiphylo object (trees), a corresponding set of 
# marginal asr's from applying ancRECON across trees (clist), using rates 
# estimated from corHMM, and a consensus
# tree (ctree). It makes a datatable of probability of being in each state
# for nodes that are bifurcating in the ctree. These are indexed by their 
# corresponding ctree node. 
ctreeAnc <- function (trees, ctree, clist){
  # Return nodes in trees that match with ctree nodes
  matched.nodes <- lapply(trees, function (x) matchNodes(ctree, x)[,2])
  # Make vector of nodes in ctree 
  ctree.nodes <- unique(ctree$edge[,1])
  # Make vectors of nodes in trees. This will be used along with matched.nodes
  # to create an index of state probs to select. State probs are listed in same
  # order as nodes on the tree. 
  trees.nodes <- lapply(trees, function(x) unique(x$edge[,1]))
  # Get index of nodes that will be used to select state probs
  state.index <- mapply(match, matched.nodes, trees.nodes, SIMPLIFY = FALSE)
  # Select state probs using state.index. These are in order of ctree nodes
  state.probs <- mapply(function(x,y) x$lik.anc.states[y,], clist, state.index, 
                        SIMPLIFY = FALSE)
  # Add column that gives corresponding ctree node for each state prob
  state.probs <- lapply(state.probs, cbind, ctree.nodes)
  # Convert list of node probs to datatable 
  nodeframe <- rbindlist(lapply(state.probs, data.table))
  return(nodeframe)
}

# Read in supertree tree and traits--------------------------------------------
mtree <- read.nexus("~/zoox/data/2017-12-29/mtree_traits/mtree.trees")

traits <- fread("~/zoox/data/2017-12-29/mtree_traits/mtree_traits_B_as_Z.csv",
                header = FALSE, col.names = c("taxa", "state"))


# Format traits and tip labels-------------------------------------------------
# Drop taxa without data from tree
mtree <- lapply(mtree, drop.tip, tip = traits[state == "-", taxa])
class(mtree) <- "multiPhylo"

# Remove taxa with missing data from traits
traits <- traits[state != "-"]


# Sample 100 trees used to estimate rates--------------------------------------
set.seed(12)
mtree <- sample(mtree, 100, replace = FALSE)

# Read in ancestral state reconstructions--------------------------------------
# Facultative coded as zoox
anc1 <- readRDS("~/zoox/results/2019-04-04/mtree_1rate_BasZ_anc.rds")
anc2 <- readRDS("~/zoox/results/2019-04-04/mtree_2rate_BasZ_anc.rds")
anc3 <- readRDS("~/zoox/results/2019-04-04/mtree_3rate_BasZ_anc.rds")

# Facultative coded as azoox
az.anc1 <- readRDS("~/zoox/results/2019-04-04/mtree_1rate_BasA_anc.rds")
az.anc2 <- readRDS("~/zoox/results/2019-04-04/mtree_2rate_BasA_anc.rds")
az.anc3 <- readRDS("~/zoox/results/2019-04-04/mtree_3rate_BasA_anc.rds")

# Prep data for summary and run------------------------------------------------

ctree <- consensus(mtree, p = 0.95)

anclist <- list(anc1, anc2, anc3, az.anc1, az.anc2, az.anc3)

nodeframes <- mclapply(anclist, ctreeAnc, trees = mtree, ctree = ctree,
                       mc.cores = 40)

saveRDS(nodeframes, "mtree_subsample_nodeframes.rds")




