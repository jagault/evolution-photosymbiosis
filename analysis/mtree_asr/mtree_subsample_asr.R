#Clear workspace
rm(list = ls())

# Set working directory
setwd("~/zoox/results/2019-04-04/")

# Load packages
library(ape)
library(data.table)
library(corHMM)
library(parallel)


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

# Read in corHMM runs----------------------------------------------------------
rate1 <- readRDS("~/zoox/results/2019-03-26/mtree_1rate.rds")
rate2 <- readRDS("~/zoox/results/2019-03-26/mtree_2rate.rds")
rate3 <- readRDS("~/zoox/results/2019-03-26/mtree_3rate.rds")

az.rate1 <- readRDS("~/zoox/results/2019-03-26/mtree_1rate_BasA.rds")
az.rate2 <- readRDS("~/zoox/results/2019-03-26/mtree_2rate_BasA.rds")
az.rate3 <- readRDS("~/zoox/results/2019-03-26/mtree_3rate_BasA.rds")


# Make rate vectors------------------------------------------------------------
### Facultative coded as zoox
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


### Facultative coded as azoox
# 1 rate
az.par1 <- sapply(az.rate1, "[", "solution")
az.par1 <- lapply(az.par1, as.vector)
az.par1 <- lapply(az.par1, na.omit)

# 2 rate
az.par2 <- sapply(az.rate2, "[", "solution")
az.par2 <- lapply(az.par2, as.vector)
az.par2 <- lapply(az.par2, na.omit)

# 3 rate
az.par3 <- sapply(az.rate3, "[", "solution")
az.par3 <- lapply(az.par3, as.vector)
az.par3 <- lapply(az.par3, na.omit)


# Estimate ancestral states----------------------------------------------------

### Facultative coded as zoox
# 1 rate
anc1 <- mcmapply(ancRECON, mtree, p = par1, 
                 MoreArgs = list(data = traits[, .(taxa, Z)], 
                                 method = "marginal", ntraits = 1, 
                                 charnum = 1, rate.cat = 1,
                                 model = "ARD"), 
                 mc.cores = 40, SIMPLIFY = FALSE)

saveRDS(anc1, "mtree_1rate_BasZ_anc.rds")

# 2 rate
anc2 <- mcmapply(ancRECON, mtree, p = par2,
                 MoreArgs = list(data = traits[, .(taxa, Z)],
                                 method = "marginal", hrm = TRUE, 
                                 rate.cat = 2),
                 mc.cores = 40, SIMPLIFY = FALSE)

saveRDS(anc2, "mtree_2rate_BasZ_anc.rds")

# 3 rate
anc3 <- mcmapply(ancRECON, mtree, p = par3,
                 MoreArgs = list(data = traits[, .(taxa, Z)],
                                 method = "marginal", hrm = TRUE, 
                                 rate.cat = 3),
                 mc.cores = 40, SIMPLIFY = FALSE)

saveRDS(anc3, "mtree_3rate_BasZ_anc.rds")

### Facultative coded as azoox
# 1 rate
az.anc1 <- mcmapply(ancRECON, mtree, p = az.par1, 
                    MoreArgs = list(data = az.traits[, .(taxa, Z)], 
                                    method = "marginal", ntraits = 1, 
                                    charnum = 1, rate.cat = 1,
                                    model = "ARD"), 
                    mc.cores = 40, SIMPLIFY = FALSE)

saveRDS(az.anc1, "mtree_1rate_BasA_anc.rds")


# 2 rate
az.anc2 <- mcmapply(ancRECON, mtree, p = az.par2,
                    MoreArgs = list(data = az.traits[, .(taxa, Z)],
                                    method = "marginal", hrm = TRUE, 
                                    rate.cat = 2),
                    mc.cores = 40, SIMPLIFY = FALSE)

saveRDS(az.anc2, "mtree_2rate_BasA_anc.rds")


# 3 rate
az.anc3 <- mcmapply(ancRECON, mtree, p = az.par3,
                    MoreArgs = list(data = az.traits[, .(taxa, Z)],
                                    method = "marginal", hrm = TRUE, 
                                    rate.cat = 3),
                    mc.cores = 40, SIMPLIFY = FALSE)

saveRDS(az.anc3, "mtree_3rate_BasA_anc.rds")
