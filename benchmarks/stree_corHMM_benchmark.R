# Load packages
library(phytools)
library(data.table)
library(phangorn)
library(corHMM)
library(here)

# Set working directory
setwd(here("analysis/stree_corHMM"))

# Read in tree and traits------------------------------------------------------
stree <- read.nexus(here("data/updated_trees_traits/stree_traits", 
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

# Run corHMM-------------------------------------------------------------------
# corHMM has built-in capability for running analyses in parallel. Adjust the 
# n.cores argument to specify the number of cores to use. 

b1 <- system.time(
  stree_1rate <- corHMM(stree[[1]], data = tframe[, .(rn, Z)], rate.cat = 1,
                        node.states = "none", nstarts = 0,
                        n.cores = 1)
)

b2 <- system.time(
  stree_2rate <- corHMM(stree[[1]], data = tframe[, .(rn, Z)], rate.cat = 2,
                        node.states = "none", nstarts = 0,
                        n.cores = 1)
)

b3 <- system.time(
  stree_3rate <- corHMM(stree[[1]], data = tframe[, .(rn, Z)], rate.cat = 3,
                        node.states = "none", nstarts = 0,
                        n.cores = 1)
)