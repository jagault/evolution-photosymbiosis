# Load packages
library(phytools)
library(data.table)
library(phangorn)
library(corHMM)

# Set working directory
setwd("/home/analysis")

# Capture command line arguments
args <- commandArgs(trailingOnly = T)
# Interprets options as characters. Change to integer
args <- as.integer(args)

# Read in tree and traits------------------------------------------------------
stree <- read.nexus("stree.trees")

traits <- fread(file = "stree_traits_B_as_Z.csv", 
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
stree_4rate <- corHMM(stree[[args[1]]], data = tframe[, .(rn, Z)], rate.cat = 4,
                      node.states = "none", nstarts = 1,
                      n.cores = 1)

# Write to file----------------------------------------------------------------
saveRDS(stree_4rate,
        file = paste("stree-4rate", "-", "t", args[1], "-", "r", args[2], 
                     ".rds", sep = ""))