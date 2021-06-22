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
mtree <- read.nexus("mtree.trees")

traits <- fread(file = "mtree_traits_B_as_Z.csv", 
                header = FALSE, col.names = c("taxa", "state"))

# Format traits and tip labels-------------------------------------------------
# Drop taxa without data from tree
mtree <- lapply(mtree, drop.tip, tip = traits[state == "-", taxa])
class(mtree) <- "multiPhylo"

# Remove taxa with missing data from traits
traits <- traits[state != "-"]

# Format trait data for corHMM
traits[state == "Z", Z := 1]
traits[state == "A", Z := 0]

# Read in corHMM runs----------------------------------------------------------
rate2 <- readRDS("mtree-2rate.rds")

# Make rate vectors-------------------------------------------------------------
par2 <- sapply(rate2, "[", "solution")
par2 <- lapply(par2, as.vector)
par2 <- lapply(par2, na.omit)

# Estimate ancestral states----------------------------------------------------
anc2 <- ancRECON(mtree[[args[1]]], p = par2[[args[1]]], 
                 data = traits[, .(taxa, Z)], method = "marginal", hrm = TRUE, 
                 rate.cat = 2)

# Write to file----------------------------------------------------------------
saveRDS(anc2, file = paste("mtree-2rate-anc-t", args[1], ".rds", sep = ""))