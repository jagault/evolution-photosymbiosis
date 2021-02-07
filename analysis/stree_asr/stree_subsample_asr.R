#Clear workspace
rm(list = ls())

# Set working directory
setwd("~/zoox/results/2019-03-28/")

# Load packages
library(ape)
library(data.table)
library(corHMM)
library(parallel)


# Read in supertree tree and traits--------------------------------------------
stree <- read.nexus("~/zoox/data/2017-12-29/stree_traits/stree.trees")

traits <- fread("~/zoox/data/2017-12-29/stree_traits/stree_traits_B_as_Z.csv",
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

# Select tree used for analysis------------------------------------------------
set.seed(1)

stree <- sample(stree, 100, replace = FALSE)

# Read in corHMM runs----------------------------------------------------------
rate1 <- readRDS("~/zoox/results/2019-02-04/stree_1rate.rds")
rate2 <- readRDS("~/zoox/results/2019-02-04/stree_2rate.rds")
rate3 <- readRDS("~/zoox/results/2019-02-04/stree_3rate.rds")


# Make rate vectors------------------------------------------------------------
# 1 rate
par1 <- sapply(rate1, "[", "solution")
par1 <- lapply(par1, as.vector)
par1 <- lapply(par1, na.omit)

# 2 rate
par2 <- sapply(rate2, "[", "solution")
par2 <- lapply(par2, as.vector)
par2 <- lapply(par2, na.omit)

# 3 rate
par3 <- sapply(rate3, "[", "solution")
par3 <- lapply(par3, as.vector)
par3 <- lapply(par3, na.omit)

# Estimate ancestral states----------------------------------------------------
# 1 rate
anc1 <- mcmapply(ancRECON, stree, p = par1, 
                 MoreArgs = list(data = traits[, .(taxa, Z)], 
                                 method = "marginal", ntraits = 1, 
                                 charnum = 1, rate.cat = 1,
                                 model = "ARD"), 
                 mc.cores = 40, SIMPLIFY = FALSE)

saveRDS(anc1, "stree_1rate_anc.rds")

# 2 rate
anc2 <- mcmapply(ancRECON, stree, p = par2,
                    MoreArgs = list(data = traits[, .(taxa, Z)],
                                    method = "marginal", hrm = TRUE, 
                                    rate.cat = 2),
                 mc.cores = 40, SIMPLIFY = FALSE)

saveRDS(anc2, "stree_2rate_anc.rds")

# 3 rate
anc3 <- mcmapply(ancRECON, stree, p = par3,
                 MoreArgs = list(data = traits[, .(taxa, Z)],
                                 method = "marginal", hrm = TRUE, 
                                 rate.cat = 3),
                 mc.cores = 40, SIMPLIFY = FALSE)

saveRDS(anc3, "stree_3rate_anc.rds")

