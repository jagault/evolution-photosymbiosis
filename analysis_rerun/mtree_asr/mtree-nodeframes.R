# Load packages
library(phytools)
library(data.table)
library(parallel)
library(here)

# Set working directory
setwd(here("analysis_rerun/mtree_asr"))

# Source helper functions
source(here("R/asr-summary.R"))


# Read in tree and traits-------------------------------------------------------
mtree <- read.nexus(here("data/updated_trees_traits/mtree_traits",
                         "mtree.trees"))

traits <- fread(file = here("data/updated_trees_traits/mtree_traits",
                            "mtree_traits_B_as_Z.csv"),
                header = FALSE, col.names = c("taxa", "state"))


# Format traits and tip labels--------------------------------------------------
# Drop taxa without data from tree
mtree <- lapply(mtree, drop.tip, tip = traits[state == "-", taxa])
class(mtree) <- "multiPhylo"


# Read in ancestral state reconstructions---------------------------------------
anc <- readRDS(here("analysis_rerun/mtree_asr", 
                    "mtree-asr.rds"))

# Compute 95% consensus tree----------------------------------------------------
ctree <- consensus(mtree, p = 0.95)


# Summarize asr's on single ctree-----------------------------------------------
nodeframes <- mcCtreeAnc(trees = mtree, ctree = ctree, clist = anc)

# Write to file
saveRDS(nodeframes, "mtree-nodeframes.rds")