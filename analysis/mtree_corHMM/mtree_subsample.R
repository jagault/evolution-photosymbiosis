# Set working directory
setwd("~/zoox/results/2019-03-26/")

# Load packages
library(data.table)
library(corHMM)
library(phytools)


# Read in supertree tree and traits--------------------------------------------
mtree <- read.nexus("~/zoox/data/2017-12-29/mtree_traits/mtree.trees")

traits <- fread("~/zoox/data/2017-12-29/mtree_traits/mtree_traits_B_as_Z.csv",
                header = FALSE, col.names = c("taxa", "state"))

az.traits <- fread("~/zoox/data/2017-12-29/mtree_traits/mtree_traits_B_as_A.csv",
                   header = FALSE, col.names = c("taxa", "state"))

# Format traits and tip labels-------------------------------------------------
# Drop taxa without data from tree
mtree <- lapply(mtree, drop.tip, tip = traits[state == "-", taxa])
class(mtree) <- "multiPhylo"

# Remove taxa with missing data from traits
traits <- traits[state != "-"]
az.traits <- az.traits[state != "-"]


# Format trait data for corHMM
traits[state == "Z", Z := 1]
traits[state == "A", Z := 0]

az.traits[state == "Z", Z := 1]
az.traits[state == "A", Z := 0]

# Select 100 mtrees
set.seed(12)
mtree <- sample(mtree, 100, replace = FALSE)

# Run corHMM-------------------------------------------------------------------
mtree_1rate <- lapply(mtree, corHMM, data = traits[, .(taxa, Z)], rate.cat = 1,
                      node.states = "none", nstarts = 100, n.cores = 40,
                      diagn = TRUE)

saveRDS(mtree_1rate, "mtree_1rate.rds")


mtree_2rate <- lapply(mtree, corHMM, data = traits[, .(taxa, Z)], rate.cat = 2,
                      node.states = "none", nstarts = 100, n.cores = 40,
                      diagn = TRUE)

saveRDS(mtree_2rate, "mtree_2rate.rds")


mtree_3rate <- lapply(mtree, corHMM, data = traits[, .(taxa, Z)], 
                      rate.cat = 3, node.states = "none", 
                      nstarts = 100, n.cores = 40, diagn = TRUE)

saveRDS(mtree_3rate, "mtree_3rate.rds")


mtree_4rate <- lapply(mtree, corHMM, data = traits[, .(taxa, Z)], 
                      rate.cat = 4, node.states = "none", 
                      nstarts = 100, n.cores = 40, diagn = TRUE)

saveRDS(mtree_4rate, "mtree_4rate.rds")


### Facultative coded as azoox
az.mtree_1rate <- lapply(mtree, corHMM, data = az.traits[, .(taxa, Z)], 
                         rate.cat = 1, node.states = "none", nstarts = 100, 
                         n.cores = 40, diagn = TRUE)

saveRDS(az.mtree_1rate, "mtree_1rate_BasA.rds")


az.mtree_2rate <- lapply(mtree, corHMM, data = az.traits[, .(taxa, Z)], 
                         rate.cat = 2, node.states = "none", nstarts = 100, 
                         n.cores = 40, diagn = TRUE)

saveRDS(az.mtree_2rate, "mtree_2rate_BasA.rds")

az.mtree_3rate <- lapply(mtree, corHMM, data = az.traits[, .(taxa, Z)], 
                         rate.cat = 3, node.states = "none", nstarts = 100, 
                         n.cores = 40, diagn = TRUE)

saveRDS(az.mtree_3rate, "mtree_3rate_BasA.rds")



