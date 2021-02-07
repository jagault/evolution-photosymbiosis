#Clear workspace
rm(list = ls())

# Set working directory
setwd("~/zoox/results/2019-02-04/")

# Load packages
library(phytools)
library(data.table)
library(phangorn)
library(corHMM)
library(parallel)

# Read in molecular tree and traits--------------------------------------------
stree <- read.nexus("~/zoox/data/2017-12-29/stree_traits/stree.trees")

traits <- fread("~/zoox/data/2017-12-29/stree_traits/stree_traits_B_as_Z.csv",
                header = FALSE, col.names = c("taxa", "state"))


# Format traits and tip labels-------------------------------------------------
# Make vector of traits for make.simmap
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

# Run two rate regime across full posterior------------------------------------
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


stree_5rate <- lapply(stree, corHMM, data = tframe[, .(rn, Z)], rate.cat = 5,
                      node.states = "none", nstarts = 120,
                      n.cores = 40)

saveRDS(stree_5rate, "stree_5rate.rds")
