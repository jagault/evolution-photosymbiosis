# Load packages
library(ape)
library(data.table)
library(corHMM)
library(parallel)
library(here)

# Set working directory
setwd(here("analysis/mtree_stree_facultative_removed"))

# Read in trees and traits-----------------------------------------------------

# supertree
stree <- mtree <- read.nexus(here("data/updated_trees_traits/stree_traits", 
                                  "stree.trees"))

traits <- fread(here("data/updated_trees_traits/stree_traits",
                     "stree_traits_B_as_Z.csv"),
                header = FALSE, col.names = c("taxa", "state"))

# molecular tree
mtree <- read.nexus(here("data/updated_trees_traits/mtree_traits", 
                         "mtree.trees"))

traits <- fread(here("data/updated_trees_traits/mtree_traits",
                     "mtree_traits_B_as_Z.csv"),
                header = FALSE, col.names = c("taxa", "state"))

# Format traits and tip labels-------------------------------------------------

# Supertree
stree <- lapply(stree, drop.tip, tip = straits[state == "-", taxa])
class(stree) <- "multiPhylo"

stree <- lapply(stree, drop.tip, tip = straits[state == "B", taxa])
class(stree) <- "multiPhylo"

straits <- straits[state != "-"]
straits <- straits[state != "B"]

straits[state == "Z", Z := 1]
straits[state == "A", Z := 0]


# Molecular tree
mtree <- lapply(mtree, drop.tip, tip = mtraits[state == "-", taxa])
class(mtree) <- "multiPhylo"

mtree <- lapply(mtree, drop.tip, tip = mtraits[state == "B", taxa])
class(mtree) <- "multiPhylo"

mtraits <- mtraits[state != "-"]
mtraits <- mtraits[state != "B"]

mtraits[state == "Z", Z := 1]
mtraits[state == "A", Z := 0]


# Randomly sample 100 trees----------------------------------------------------
set.seed(1)
stree <- sample(stree, 100, replace = FALSE)

set.seed(12)
mtree <- sample(mtree, 100, replace = FALSE)

# Run analysis-----------------------------------------------------------------
# corHMM has built-in capability for running analyses in parallel. Adjust the 
# n.cores argument to specify the number of cores to use. 

### Molecular tree
mtree_1rate <- lapply(mtree, corHMM, data = mtraits[, .(taxa, Z)], 
                      rate.cat = 1, node.states = "none", nstarts = 100,
                      n.cores = 40)

saveRDS(mtree_1rate, "mtree_1rate_nb.rds")


mtree_2rate <- lapply(mtree, corHMM, data = mtraits[, .(taxa, Z)], 
                      rate.cat = 2, node.states = "none", nstarts = 100,
                      n.cores = 40)

saveRDS(mtree_2rate, "mtree_2rate_nb.rds")


mtree_3rate <- lapply(mtree, corHMM, data = mtraits[, .(taxa, Z)], 
                      rate.cat = 3, node.states = "none", nstarts = 100,
                      n.cores = 40)

saveRDS(mtree_3rate, "mtree_3rate_nb.rds")


### Supertree
stree_1rate <- lapply(stree, corHMM, data = straits[, .(taxa, Z)], 
                      rate.cat = 1, node.states = "none", nstarts = 100,
                      n.cores = 40)

saveRDS(stree_1rate, "stree_1rate_nb.rds")


stree_2rate <- lapply(stree, corHMM, data = straits[, .(taxa, Z)], 
                      rate.cat = 2, node.states = "none", nstarts = 100,
                      n.cores = 40)

saveRDS(stree_2rate, "stree_2rate_nb.rds")


stree_3rate <- lapply(stree, corHMM, data = straits[, .(taxa, Z)], 
                      rate.cat = 3, node.states = "none", nstarts = 100,
                      n.cores = 40)

saveRDS(stree_3rate, "stree_3rate_nb.rds")