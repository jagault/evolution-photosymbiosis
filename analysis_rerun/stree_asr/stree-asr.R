# Load packages
library(phytools)
library(data.table)
library(phangorn)
library(corHMM)
# library(here)

# Set working directory
# setwd("/home/analysis")

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

# Read in corHMM runs----------------------------------------------------------
rate3 <- readRDS("stree_3rate.rds")

# Make rate vectors------------------------------------------------------------# 3 rate
par3 <- sapply(rate3, "[", "solution")
par3 <- lapply(par3, as.vector)
par3 <- lapply(par3, na.omit)

# Estimate ancestral states----------------------------------------------------
anc3 <- ancRECON(stree[[args[1]]], p = par3[args[1]], 
                 data = traits[, .(taxa, Z)], method = "marginal", hrm = TRUE, 
                 rate.cat = 3)

# Write to file----------------------------------------------------------------
saveRDS(anc3, file = paste("stree-3rate-anc-t", args[1], ".rds", sep = ""))
