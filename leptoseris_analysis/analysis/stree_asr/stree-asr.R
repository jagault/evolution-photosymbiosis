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
stree <- read.nexus("stree_lt.trees")

traits <- fread(file = "stree_traits_B_as_Z_lt.csv", 
                header = FALSE, col.names = c("taxa", "state"))

# Format traits and tip labels-------------------------------------------------
# Drop taxa without data from tree
stree <- lapply(stree, drop.tip, tip = traits[state == "-", taxa])
class(stree) <- "multiPhylo"

# Remove taxa with missing data from traits
traits <- traits[state != "-"]

# Format trait data for corHMM
traits[state == "Z", Z := 1]
traits[state == "A", Z := 0]

# Read in corHMM runs----------------------------------------------------------
rate3 <- readRDS("stree-3rate.rds")

# Make rate vectors-------------------------------------------------------------
par3 <- sapply(rate3, "[", "solution")
par3 <- lapply(par3, as.vector)
par3 <- lapply(par3, na.omit)

# Estimate ancestral states----------------------------------------------------
anc3 <- ancRECON(stree[[args[1]]], p = par3[[args[1]]], 
                 data = traits[, .(taxa, Z)], method = "marginal", hrm = TRUE, 
                 rate.cat = 3)

# Write to file----------------------------------------------------------------
saveRDS(anc3, file = paste("stree-3rate-lt-anc-t", args[1], ".rds", sep = ""))