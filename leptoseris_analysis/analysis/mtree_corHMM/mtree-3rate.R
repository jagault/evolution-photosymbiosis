# Load packages
library(phytools)
library(data.table)
library(phangorn)
library(corHMM)

# Set working directory
# setwd("/home/analysis")

# Capture command line arguments
args <- commandArgs(trailingOnly = T)
# Interprets options as characters. Change to integer
args <- as.integer(args)

# Read in tree and traits------------------------------------------------------
mtree <- read.nexus("mtree_lt.trees")

traits <- fread(file = "mtree_traits_B_as_Z_lt.csv", 
                header = FALSE, col.names = c("taxa", "state"))

# Format traits and tip labels-------------------------------------------------
# Make vector of traits
tvec <- traits[, state]
names(tvec) <- traits[, taxa]

# Drop taxa without data from tree and traits
mtree <- lapply(mtree, drop.tip, tip = names(tvec[tvec == "-"]))
class(mtree) <- "multiPhylo"

# Remove taxa with missing data
tvec <- tvec[!tvec == "-"]

# Make dataframe of traits for corHMM
tm <- to.matrix(tvec, c("A","Z"))
tframe <- data.frame(tm)
tframe <- data.table(tframe, keep.rownames = TRUE)

# Run corHMM-------------------------------------------------------------------
mtree_3rate <- corHMM(mtree[[args[1]]], data = tframe[, .(rn, Z)], rate.cat = 3,
                      node.states = "none", nstarts = 1,
                      n.cores = 1)

# Write to file----------------------------------------------------------------
saveRDS(mtree_3rate,
        file = paste("mtree-3rate", "-", "t", args[1], "-", "r", args[2], 
                     ".rds", sep = ""))