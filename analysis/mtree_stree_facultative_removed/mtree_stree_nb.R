#Clear workspace
rm(list = ls())

# Set working directory
setwd("~/zoox/results/2019-05-13/")

# Load packages
library(ape)
library(data.table)
library(corHMM)
library(parallel)


# Read in trees and traits-----------------------------------------------------
stree <- read.nexus("~/zoox/data/2017-12-29/stree_traits/stree.trees")

straits <- fread("~/zoox/data/2017-12-29/stree_traits/stree_traits.csv",
                 header = FALSE, col.names = c("taxa", "state"))


mtree <- read.nexus("~/zoox/data/2017-12-29/mtree_traits/mtree.trees")

mtraits <- fread("~/zoox/data/2017-12-29/mtree_traits/mtree_traits.csv",
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



