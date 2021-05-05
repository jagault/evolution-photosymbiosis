# Load packages
library(phytools)
library(data.table)
library(phangorn)
library(corHMM)
library(here)

# Set working directory
setwd(here("analysis/stree_corHMM"))

# Read in tree and traits------------------------------------------------------
stree <- read.nexus(here("data/updated_trees_traits/stree_traits", 
                                  "stree.trees"))

traits <- fread(here("data/updated_trees_traits/stree_traits",
                     "stree_traits_B_as_Z.csv"),
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
# corHMM has built-in capability for running analyses in parallel. Adjust the 
# n.cores argument to specify the number of cores to use. 

t1 <- system.time(
  stree_1rate <- corHMM(stree[[1]], data = tframe[, .(rn, Z)], rate.cat = 1,
                        node.states = "none", nstarts = 1,
                        n.cores = 1)
)

t2 <- system.time(
  stree_2rate <- corHMM(stree[[1]], data = tframe[, .(rn, Z)], rate.cat = 2,
                        node.states = "none", nstarts = 1,
                        n.cores = 1)
)

t3 <- system.time(
  stree_3rate <- corHMM(stree[[1]], data = tframe[, .(rn, Z)], rate.cat = 3,
                        node.states = "none", nstarts = 1,
                        n.cores = 1)
)

# Test different approaches to random restarts---------------------------------

t1r <- system.time(
  stree_1rate_r <- corHMM(stree[[1]], data = tframe[, .(rn, Z)], rate.cat = 1,
                          node.states = "none", nstarts = 20,
                          n.cores = 1)
)

st1 <- list()
tloop1 <- system.time(
for (i in 1:19) {
  st1[[i]] <- corHMM(stree[[1]], data = tframe[, .(rn, Z)], rate.cat = 1,
                     node.states = "none", nstarts = 1,
                     n.cores = 1)
}
)




t2r <- system.time(
  stree_2rate_r <- corHMM(stree[[1]], data = tframe[, .(rn, Z)], rate.cat = 2,
                          node.states = "none", nstarts = 5,
                          n.cores = 1)
)

st2 <- list()
tloop2 <- system.time(
  for (i in 1:4) {
    st2[[i]] <- corHMM(stree[[1]], data = tframe[, .(rn, Z)], rate.cat = 2,
                       node.states = "none", nstarts = 1,
                       n.cores = 1)
  }
)