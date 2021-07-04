# Load packages
library(data.table)
library(phytools)
library(phangorn)
library(here)

# Set working directory
setwd(here("analysis_rerun/results/stree/gains_losses"))

# Source helper functions
source(here("R/aic-rate-summary.R"))
source(here("R/gains-losses.R"))


# Read in tree and traits-------------------------------------------------------
stree <- read.nexus(here("data/updated_trees_traits/stree_traits",
                         "stree.trees"))

traits <- fread(file = here("data/updated_trees_traits/stree_traits",
                            "stree_traits_B_as_Z.csv"),
                header = FALSE, col.names = c("taxa", "state"))


# Format traits and tip labels--------------------------------------------------
# Drop taxa without data from tree
stree <- lapply(stree, drop.tip, tip = traits[state == "-", taxa])
class(stree) <- "multiPhylo"


# Read in corHMM runs----------------------------------------------------------
# Read in asr summarized across all 1000 trees
fanc3 <- readRDS(here("analysis_rerun/stree_asr", 
                      "stree-nodeframes.rds"))

# Read in 1000 individual asrs
anc3 <- readRDS(here("analysis_rerun/stree_asr", 
                     "stree-asr.rds"))


# Calculate consensus state probs at internal nodes and tips-------------------
# Make consensus tree
ctree <- consensus(stree, p = 0.95)
ctree <- compute.brlen(ctree)

# Get state probs at internal nodes
setnames(fanc3, c("AS", "ZS", "AM", "ZM", "AF", "ZF", "ctree.nodes"))
fanc3 <- fanc3[, .(AS = mean(AS), ZS = mean(ZS),
                   AM = mean(AM), ZM = mean(ZM),
                   AF = mean(AF), ZF = mean(ZF)), 
               by = ctree.nodes]
fanc3[, ctree.nodes := NULL]

### Get tip probs
tip.probs <- lapply(anc3, "[[", "lik.tip.states")
tip.probs <- lapply(tip.probs, data.table)

tip.index <- lapply(stree, "[[", "tip.label")
tip.index <- lapply(tip.index, data.table)
lapply(tip.index, setnames, "species")

tip.probs <- rbindlist(mapply(cbind, tip.probs, tip.index, SIMPLIFY = FALSE))
setnames(tip.probs, c("AS", "ZS", "AM", "ZM", "AF", "ZF", "species"))

tip.probs <- tip.probs[, .(AS = mean(AS), ZS = mean(ZS),
                           AM = mean(AM), ZM = mean(ZM),
                           AF = mean(AF), ZF = mean(ZF)), 
                       by = species]

# Count number of changes on consensus asr-------------------------------------

### All 4 states
tip.vec <- apply(tip.probs[,2:7], 1, which.max)
node.vec <- apply(fanc3, 1, which.max)
cat <- c(tip.vec, node.vec)

cc <- getChanges(ctree, cat, rate.cat = 3)
cc <- melt(cc, variable = "type", value = "ctree")


### Binary gains and losses
az.fanc3 <- copy(fanc3)
az.tip <- copy(tip.probs)

az.fanc3[, asum := rowSums(.SD), .SDcols = grep("A", names(az.fanc3))]
az.fanc3[, zsum := rowSums(.SD), .SDcols = grep("Z", names(az.fanc3))]

az.tip[, asum := rowSums(.SD), .SDcols = grep("A", names(az.tip))]
az.tip[, zsum := rowSums(.SD), .SDcols = grep("Z", names(az.tip))]


az.node.vec <- apply(az.fanc3[, .(asum, zsum)], 1, which.max)
az.tip.vec <- apply(az.tip[, .(asum, zsum)], 1, which.max)

az.cat <- c(az.tip.vec, az.node.vec)

az.cc <- getChanges(ctree, az.cat, rate.cat = 1)
az.cc <- melt(az.cc, variable = "type", value = "ctree")


# Count number of changes across posterior-------------------------------------
### All 4 states
ftips <- lapply(anc3, "[[", "lik.tip.states")
fnodes <- lapply(anc3, "[[", "lik.anc.states")

ftips <- lapply(ftips, function(x) apply(x, 1, which.max))
fnodes <- lapply(fnodes, function(x) apply(x, 1, which.max))
fcat <- mapply(function(x,y) c(x,y), ftips, fnodes, SIMPLIFY = FALSE)

fchanges <- rbindlist(mapply(getChanges, stree, fcat, rate.cat = 3, 
                             SIMPLIFY = FALSE))


cmelt <- melt(fchanges, variable.name = "type", value = "n")
csumm <- cmelt[,  .(mean = mean(n), median = median(n)), by = type]


# Get confidence intervals around median number of changes
cq <- apply(fchanges[,1:ncol(fchanges)], 2, bootMed, n = 100000)
cq <- data.table(t(cq), keep.rownames = T)
setnames(cq, c("type", "0.025", "median", "0.975"))

setkey(cq, type)
setkey(cc, type)
cq <- cq[cc]
setcolorder(cq, c("type", "ctree", "0.025", "median", "0.975"))

write.csv(cq, file = "all_changes.csv")


### Binary gains and losses
az.ftips <- lapply(anc3, "[[", "lik.tip.states")
az.fnodes <- lapply(anc3, "[[", "lik.anc.states")

az.ftips <- lapply(az.ftips, data.table)
az.fnodes <- lapply(az.fnodes, data.table)

lapply(az.ftips, setnames, c("AS", "ZS", "AM", "ZM", "AF", "ZF"))
lapply(az.fnodes, setnames, c("AS", "ZS", "AM", "ZM", "AF", "ZF"))

# lapply(az.ftips, function(x) 
#   x[, asum := rowSums(.SD), .SDcols = grep("A", names(x))])
# lapply(az.ftips, function(x) 
#   x[, zsum := rowSums(.SD), .SDcols = grep("Z", names(x))])
# 
# lapply(az.fnodes, function(x) 
#   x[, asum := rowSums(.SD), .SDcols = grep("A", names(x))])
# lapply(az.fnodes, function(x) 
#   x[, zsum := rowSums(.SD), .SDcols = grep("Z", names(x))])

az.tvec <- lapply(az.ftips, function(x) apply(x[,.(asum, zsum)], 1, which.max))
az.nvec <- lapply(az.fnodes, function(x) apply(x[,.(asum, zsum)],1, which.max))

az.fcat <- mapply(function(x,y) c(x,y), az.tvec, az.nvec, SIMPLIFY = FALSE)


az.fchanges <- rbindlist(mapply(getChanges, stree, az.fcat, rate.cat = 1,
                                SIMPLIFY = FALSE))


az.cmelt <- melt(az.fchanges, variable.name = "type", value = "n")
az.csumm <- az.cmelt[,  .(mean = mean(n), median = median(n)), by = type]
az.csumm

# Get confidence intervals around median number of changes
az.cq <- apply(az.fchanges[,1:ncol(az.fchanges)], 2, bootMed, n = 100000)
az.cq <- data.table(t(az.cq), keep.rownames = T)

setnames(az.cq, c("type", "0.025", "median", "0.975"))

setkey(az.cq, type)
setkey(az.cc, type)
az.cq <- az.cq[az.cc]
setcolorder(az.cq, c("type", "ctree", "0.025", "median", "0.975"))
az.cq
az.cq <- az.cq[-(1:2),]
write.csv(az.cq, file = "binary_changes.csv")