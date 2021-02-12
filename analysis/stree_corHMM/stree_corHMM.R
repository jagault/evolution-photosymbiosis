# Load packages
library(phytools)
library(data.table)
library(phangorn)
library(corHMM)
library(here)

# Set working directory
setwd(here("analysis/stree_corHMM"))

# Read in tree and traits------------------------------------------------------
stree <- mtree <- read.nexus(here("data/updated_trees_traits/stree_traits", 
                                  "stree.trees"))

traits <- fread(here("data/updated_trees_traits/stree_traits",
                     "stree_traits_B_as_Z.csv"),
                header = FALSE, col.names = c("taxa", "state"))

# Format traits and tip labels-------------------------------------------------
# Make vector of traits
tvec <- traits[, state]
names(tvec) <- traits[, taxa]

# Drop taxa without data from tree and traits
stree <- lapply(stree, drop.tip, tip = names(tvec[tvec == "-"]))
class(stree) <- "multiPhylo"

# Remove taxa with missing data
tvec <- tvec[!tvec == "-"]

# Make dataframe of traits for corHMM
tm <- to.matrix(tvec, c("A","Z"))
tframe <- data.frame(tm)
tframe <- data.table(tframe, keep.rownames = TRUE)

# Randomly sample 100 trees----------------------------------------------------
set.seed(1)

stree <- sample(stree, 100, replace = FALSE)

# Run corHMM-------------------------------------------------------------------
# corHMM has built-in capability for running analyses in parallel. Adjust the 
# n.cores argument to specify the number of cores to use. 

stree_1rate <- lapply(stree, corHMM, data = tframe[, .(rn, Z)], rate.cat = 1,
                      node.states = "none", nstarts = 120,
                      n.cores = 40)

saveRDS(stree_1rate, "stree_1rate.rds")


stree_2rate <- lapply(stree, corHMM, data = tframe[, .(rn, Z)], rate.cat = 2,
                      node.states = "none", nstarts = 120,
                      n.cores = 40)

saveRDS(stree_2rate, "stree_2rate.rds")


stree_3rate <- lapply(stree, corHMM, data = tframe[, .(rn, Z)], rate.cat = 3,
                      node.states = "none", nstarts = 120,
                      n.cores = 40)

saveRDS(stree_3rate, "stree_3rate.rds")


stree_4rate <- lapply(stree, corHMM, data = tframe[, .(rn, Z)], rate.cat = 4,
                      node.states = "none", nstarts = 120,
                      n.cores = 40)

saveRDS(stree_4rate, "stree_4rate.rds")