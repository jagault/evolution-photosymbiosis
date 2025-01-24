# Load packages
library(phytools)
library(data.table)
library(parallel)
library(here)

# Set working directory
setwd(here("analysis/stree_asr"))

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

# Select trees used for analysis-----------------------------------------------
set.seed(1)
stree <- sample(stree, 100, replace = FALSE)


# Read in ancestral state reconstructions--------------------------------------
anc1 <- readRDS(here("analysis/stree_asr", "stree_1rate_anc.rds"))
anc2 <- readRDS(here("analysis/stree_asr", "stree_2rate_anc.rds"))
anc3 <- readRDS(here("analysis/stree_asr", "stree_3rate_anc.rds"))

# Prep data for summary and run------------------------------------------------
# Compute 95% consensus tree
ctree <- consensus(stree, p = 0.95)

# Make list of ancestral state reconstructions
anclist <- list(anc1, anc2, anc3)

# Summarize probability at selected nodes
nodeframes <- mclapply(anclist, ctreeAnc, trees = stree, ctree = ctree,
                       mc.cores = 40)

saveRDS(nodeframes, "stree_subsample_nodeframes.rds")