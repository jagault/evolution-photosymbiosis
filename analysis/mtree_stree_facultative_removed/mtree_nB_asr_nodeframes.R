# Load packages
library(phytools)
library(data.table)
library(corHMM)
library(parallel)

# Set working directory
setwd(here("analysis/mtree_stree_facultative_removed"))

# Define functions-------------------------------------------------------------
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

# Read in trees and traits-----------------------------------------------------
mtree <- read.nexus(here("data/updated_trees_traits/mtree_traits", 
                         "mtree.trees"))

traits <- fread(here("data/updated_trees_traits/mtree_traits",
                     "mtree_traits_B_as_Z.csv"),
                header = FALSE, col.names = c("taxa", "state"))

# Format traits and tip labels-------------------------------------------------
# Molecular tree
mtree <- lapply(mtree, drop.tip, tip = traits[state == "-", taxa])
class(mtree) <- "multiPhylo"

mtree <- lapply(mtree, drop.tip, tip = traits[state == "B", taxa])
class(mtree) <- "multiPhylo"

traits <- traits[state != "-"]
traits <- traits[state != "B"]

traits[state == "Z", Z := 1]
traits[state == "A", Z := 0]


# Randomly sample 100 trees----------------------------------------------------
set.seed(12)
mtree <- sample(mtree, 100, replace = FALSE)

# Read in corHMM runs----------------------------------------------------------
rate1 <- readRDS(here("analysis/mtree_stree_facultative_removed", 
                      "mtree_1rate_nb.rds"))
rate2 <- readRDS(here("analysis/mtree_stree_facultative_removed", 
                      "mtree_2rate_nb.rds"))
rate3 <- readRDS(here("analysis/mtree_stree_facultative_removed", 
                      "mtree_3rate_nb.rds"))

# Make rate vectors------------------------------------------------------------
# 1 rate
par1 <- sapply(rate1, "[", "solution")
par1 <- lapply(par1, as.vector)
par1 <- lapply(par1, na.omit)

# 2 rate
par2 <- sapply(rate2, "[", "solution")
par2 <- lapply(par2, as.vector)
par2 <- lapply(par2, na.omit)

# 3 rate
par3 <- sapply(rate3, "[", "solution")
par3 <- lapply(par3, as.vector)
par3 <- lapply(par3, na.omit)


# Estimate ancestral states----------------------------------------------------
# Use the mc.cores argument to specify how many cores to estimate ancestral 
# states in parallel

# 1 rate
anc1 <- mcmapply(ancRECON, mtree, p = par1, 
                 MoreArgs = list(data = traits[, .(taxa, Z)], 
                                 method = "marginal", ntraits = 1, 
                                 charnum = 1, rate.cat = 1,
                                 model = "ARD"), 
                 mc.cores = 40, SIMPLIFY = FALSE)

saveRDS(anc1, "mtree_nB_1rate_anc.rds")

# 2 rate
anc2 <- mcmapply(ancRECON, mtree, p = par2,
                 MoreArgs = list(data = traits[, .(taxa, Z)],
                                 method = "marginal", hrm = TRUE, 
                                 rate.cat = 2),
                 mc.cores = 40, SIMPLIFY = FALSE)

saveRDS(anc2, "mtree_nB_2rate_anc.rds")

# 3 rate
anc3 <- mcmapply(ancRECON, mtree, p = par3,
                 MoreArgs = list(data = traits[, .(taxa, Z)],
                                 method = "marginal", hrm = TRUE, 
                                 rate.cat = 3),
                 mc.cores = 40, SIMPLIFY = FALSE)

saveRDS(anc3, "mtree_nB_3rate_anc.rds")


# Summarize across all 100 trees-----------------------------------------------
# Compute 95% consensus tree
ctree <- consensus(mtree, p = 0.95)

# Make list of ancestral state reconstructions
anclist <- list(anc1, anc2, anc3)

# Summarize probability at selected nodes
nodeframes <- mclapply(anclist, ctreeAnc, trees = mtree, ctree = ctree,
                       mc.cores = 40)

saveRDS(nodeframes, "mtree_nB_nodeframes.rds")
