# Load packages
library(phytools)
library(data.table)
library(parallel)
library(here)

# Set working directory
setwd(here("leptoseris_analysis/analysis/stree_asr"))

# Source helper functions
source(here("R/asr-summary.R"))


# Read in tree and traits-------------------------------------------------------
stree <- read.nexus(here("leptoseris_analysis/data/stree_leptoseris",
                         "stree_lt.trees"))

traits <- fread(file = here("leptoseris_analysis/data/stree_leptoseris",
                            "stree_traits_B_as_Z_lt.csv"),
                header = FALSE, col.names = c("taxa", "state"))


# Format traits and tip labels--------------------------------------------------
# Drop taxa without data from tree
stree <- lapply(stree, drop.tip, tip = traits[state == "-", taxa])
class(stree) <- "multiPhylo"


# Read in ancestral state reconstructions---------------------------------------
anc <- readRDS(here("leptoseris_analysis/analysis/stree_asr", 
                    "stree-asr-lt.rds"))


# Compute 95% consensus tree----------------------------------------------------
ctree <- consensus(stree, p = 0.95)


# Summarize asr's on single ctree-----------------------------------------------
nodeframes <- mcCtreeAnc(trees = stree, ctree = ctree, clist = anc)

# Write to file
saveRDS(nodeframes, "stree-lt-nodeframes.rds")