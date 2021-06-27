# Load packages
library(data.table)
library(phytools)
library(phangorn)

# Set working directory
setwd(here("analysis_rerun/results/mtree/gains_losses"))

# Source helper functions
source(here("R/aic-rate-summary.R"))
source(here("R/gains-losses.R"))


# Read in tree and traits------------------------------------------------------
mtree <- read.nexus(here("data/updated_trees_traits/mtree_traits", 
                         "mtree.trees"))

traits <- fread(here("data/updated_trees_traits/mtree_traits",
                     "mtree_traits_B_as_Z.csv"),
                header = FALSE, col.names = c("taxa", "state"))


# Format traits and tip labels-------------------------------------------------
# Drop taxa without data from tree
mtree <- lapply(mtree, drop.tip, tip = traits[state == "-", taxa])
class(mtree) <- "multiPhylo"

# Remove taxa with missing data from traits
traits <- traits[state != "-"]


# Read in corHMM runs----------------------------------------------------------
# Read in asr summarized across all 1000 trees
fanc2 <- readRDS(here("analysis_rerun/mtree_asr", 
                      "mtree-nodeframes.rds"))

# Read in 1000 individual asrs
anc2 <- readRDS(here("analysis_rerun/mtree_asr", 
                     "mtree-asr.rds"))


# Calculate consensus state probs at internal nodes and tips-------------------
# Make consensus tree
ctree <- consensus(mtree, p = 0.95)

# Drop outgroup
ctree <- drop.tip(ctree, c("Discosoma", "Rhodactis", "Ricordea_florida"))

# Get state probs at internal nodes
setnames(fanc2, c("AS", "ZS", "AF", "ZF", "ctree.nodes"))
fanc2 <- na.omit(fanc2)
fanc2 <- fanc2[, .(AS = mean(AS), ZS = mean(ZS), AF = mean(AF), ZF = mean(ZF)), 
               by = ctree.nodes]
# Reorder nodes after dropping nas
setkey(fanc2, ctree.nodes)
# Remove nodes that were dropped from tree
fanc2 <- fanc2[-which(ctree.nodes %in% c(579, 858, 859)),]
# Drop ctree nodes
fanc2[, ctree.nodes := NULL]

### Get tip probs
tip.probs <- lapply(anc2, "[[", "lik.tip.states")
tip.probs <- lapply(tip.probs, data.table)

tip.index <- lapply(mtree, "[[", "tip.label")
tip.index <- lapply(tip.index, data.table)
lapply(tip.index, setnames, "species")

tip.probs <- rbindlist(mapply(cbind, tip.probs, tip.index, SIMPLIFY = FALSE))
setnames(tip.probs, c("AS", "ZS", "AF", "ZF", "species"))

tip.probs <- tip.probs[, .(AS = mean(AS), ZS = mean(ZS), 
                           AF = mean(AF), ZF = mean(ZF)), 
                       by = species]
# Drop tips that were dropped from tree
tip.probs <- tip.probs[-which(species %in% 
                              c("Discosoma", "Rhodactis", "Ricordea_florida"))]
tip.probs

### 100 tree posterior~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Check that the same nodes need to be dropped for all mtrees
# # 579, 580, 581 for all
# unique(unlist(lapply(mtree, function(x) Ancestors(x, which(x$tip.label %in% 
#                     c("Discosoma", "Rhodactis", "Ricordea_florida")), "all")))) 
# # This gives same answer bc 1st 3 tips are always the outgroup
# unique(unlist(lapply(mtree, Ancestors, 1:3, "all")))
# # Again, see that 1st 3 tip labels are always the outgroup
# unique(unlist(lapply(mtree, function(x) x$tip.label[1:3])))
# # So based on the above I can just drop the 1st three entries from every 
# # matrix of tip probs because 1st 3 entries are always the outgroup. 
# # I can drop the 1st three entries for nodeprobs because the descendents
# # of the outgroup are always the 1st 3 nodes. 
tips <- lapply(anc2, "[[", "lik.tip.states")
nodes <- lapply(anc2, "[[", "lik.anc.states")

tips <- lapply(tips, function(x) x[-c(1:3),])
nodes <- lapply(nodes, function(x) x[-c(1:3),])


# Count number of changes on consensus asr-------------------------------------

### All 4 states
tip.vec <- apply(tip.probs[,2:5], 1, which.max)
node.vec <- apply(fanc2, 1, which.max)
cat <- c(tip.vec, node.vec)

cc <- getChanges(ctree, cat, rate.cat = 2)
cc <- melt(cc, variable = "type", value = "ctree")


### Binary gains and losses
az.fanc2 <- copy(fanc2)
az.tip <- copy(tip.probs)

az.fanc2[, asum := rowSums(.SD), .SDcols = grep("A", names(az.fanc2))]
az.fanc2[, zsum := rowSums(.SD), .SDcols = grep("Z", names(az.fanc2))]

az.tip[, asum := rowSums(.SD), .SDcols = grep("A", names(az.tip))]
az.tip[, zsum := rowSums(.SD), .SDcols = grep("Z", names(az.tip))]

az.node.vec <- apply(az.fanc2[, .(asum, zsum)], 1, which.max)
az.tip.vec <- apply(az.tip[, .(asum, zsum)], 1, which.max)

az.cat <- c(az.tip.vec, az.node.vec)

az.cc <- getChanges(ctree, az.cat, rate.cat = 1)
az.cc <- melt(az.cc, variable = "type", value = "ctree")




# Count number of changes across posterior-------------------------------------

### All 4 states
ftips <- lapply(anc2, "[[", "lik.tip.states")
fnodes <- lapply(anc2, "[[", "lik.anc.states")

ftips <- lapply(ftips, function(x) apply(x, 1, which.max))
fnodes <- lapply(fnodes, function(x) apply(x, 1, which.max))
fcat <- mapply(function(x,y) c(x,y), ftips, fnodes, SIMPLIFY = FALSE)

fchanges <- rbindlist(mapply(getChanges, mtree, fcat, rate.cat = 2, 
                             SIMPLIFY = FALSE))
fchanges 

cmelt <- melt(fchanges, variable.name = "type", value = "n")
csumm <- cmelt[,  .(mean = mean(n), median = median(n)), by = type]
csumm

# Get confidence intervals around median number of changes
cq <- apply(fchanges[,1:ncol(fchanges)], 2, bootMed, n = 100000)
cq <- data.table(t(cq), keep.rownames = T)

setnames(cq, c("type", "0.025", "median", "0.975"))

setkey(cq, type)
setkey(cc, type)
cq <- cq[cc]
setcolorder(cq, c("type", "ctree", "0.025", "median", "0.975"))
cq
write.csv(cq, file = "all_changes.csv")




### Binary gains and losses
az.ftips <- lapply(anc2, "[[", "lik.tip.states")
az.fnodes <- lapply(anc2, "[[", "lik.anc.states")

az.ftips <- lapply(az.ftips, data.table)
az.fnodes <- lapply(az.fnodes, data.table)

lapply(az.ftips, setnames, c("AS", "ZS", "AF", "ZF"))
lapply(az.fnodes, setnames, c("AS", "ZS", "AF", "ZF"))

lapply(az.ftips, function(x) 
  x[, asum := rowSums(.SD), .SDcols = grep("A", names(x))])
lapply(az.ftips, function(x) 
  x[, zsum := rowSums(.SD), .SDcols = grep("Z", names(x))])

lapply(az.fnodes, function(x) 
  x[, asum := rowSums(.SD), .SDcols = grep("A", names(x))])
lapply(az.fnodes, function(x) 
  x[, zsum := rowSums(.SD), .SDcols = grep("Z", names(x))])

az.tvec <- lapply(az.ftips, function(x) apply(x[,.(asum, zsum)], 1, which.max))
az.nvec <- lapply(az.fnodes, function(x) apply(x[,.(asum, zsum)],1, which.max))

az.fcat <- mapply(function(x,y) c(x,y), az.tvec, az.nvec, SIMPLIFY = FALSE)


az.fchanges <- rbindlist(mapply(getChanges, mtree, az.fcat, rate.cat = 1,
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
