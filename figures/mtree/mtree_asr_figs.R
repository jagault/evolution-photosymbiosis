# Load packages
library(data.table)
library(phytools)
library(phangorn)
library(here)

# Set working directory
setwd(here("figures/mtree"))

# Source helper functions for plotting
source(here("R/plotting_functions.R"))


# Read in tree and traits------------------------------------------------------
mtree <- read.nexus(here("data/updated_trees_traits/mtree_traits", 
                         "mtree.trees"))

traits <- fread(here("data/updated_trees_traits/mtree_traits",
                     "mtree_traits_B_as_Z.csv"),
                header = FALSE, col.names = c("taxa", "state"))

stree <- read.nexus(here("data/updated_trees_traits/stree_traits", 
                         "stree.trees"))


# Format traits and tip labels-------------------------------------------------
# Drop taxa without data from tree
mtree <- lapply(mtree, drop.tip, tip = traits[state == "-", taxa])
class(mtree) <- "multiPhylo"

# Remove taxa with missing data from traits
traits <- traits[state != "-"]


# Sample 100 trees used to estimate rates--------------------------------------
set.seed(12)
mtree <- sample(mtree, 100, replace = FALSE)

# Supertree
set.seed(1)
stree <- sample(stree, 100, replace = FALSE)

# Read in corHMM runs----------------------------------------------------------
# Read in asr summarized across all 100 trees
f.anc <- readRDS(here("analysis/mtree_asr", "mtree_subsample_nodeframes.rds"))


# Rotate clades in mtree to match stree----------------------------------------

# Make consensus tree of mtree and stree
mctree <- consensus(mtree, p = 0.95)
sctree <- consensus(stree, p = 0.95)

# Select index of tips from edge matrix
tips <- sctree$edge[which(sctree$edge[,2] <= length(sctree$tip.label)),2]
# Get tip labels using index of tips from edge matrix
tip.index <- sctree$tip.label[tips]
# Select only tips that are in the molecular tere
tip.index <- tip.index[which(tip.index %in% mctree$tip.label)]

# Make named vector of tips in desired order
tip.vec <- setNames(1:Ntip(mctree), tip.index)
# Drop the 10 NAs from vector that aren't present in stree
tip.vec <- tip.vec[1:568]

# Rotate mtree to match stree
rmctree <- minRotate(mctree, tip.vec)

# Reformat some tip labels that were changed when rotating tree
rmctree$tip.label <- sub("-", "(", rmctree$tip.label)
rmctree$tip.label <- sub("-", ")", rmctree$tip.label)

# Make matrix of old and new node levels
nodes <- matchNodes(mctree, rmctree, method = "descendants")

###  Make consensus tree to drop tips and lengthen branches
ctree <- rmctree
ctree <- compute.brlen(ctree)



# Drop selected tips from tree-------------------------------------------------

# # Plot tree to get nodes that need to be dropped from data
# plotTree(ctree, ftype = "off", type = "fan")
# nodelabels(frame = "n")
# # Drop 579, 858, 859

# Drop outgroup
ctree <- drop.tip(ctree, c("Discosoma", "Rhodactis", "Ricordea_florida"))

# Drop Anthemiphyllia 
ctree <- drop.tip(ctree, "Anthemiphyllia_patera")



# Lengthen branches for presentation-------------------------------------------

# plotTree(ctree, ftype = "off", type = "fan")
# edgelabels(frame = "n")

# Root
ctree$edge.length[1] <- ctree$edge.length[1] + 0.3
ctree$edge.length[844] <- ctree$edge.length[844] + 0.3

### Robusta
ctree$edge.length[744] <- ctree$edge.length[744] + 0.01
ctree$edge.length[743] <- ctree$edge.length[743] - 0.01
ctree$edge.length[774] <- ctree$edge.length[774] + 0.01

ctree$edge.length[617] <- ctree$edge.length[617] + 0.01
ctree$edge.length[616] <- ctree$edge.length[616] - 0.01
ctree$edge.length[696] <- ctree$edge.length[696] + 0.01

ctree$edge.length[388] <- ctree$edge.length[388] + 0.01
ctree$edge.length[386] <- ctree$edge.length[386] - 0.01
ctree$edge.length[607] <- ctree$edge.length[607] + 0.01
ctree$edge.length[610] <- ctree$edge.length[610] + 0.01

ctree$edge.length[387] <- ctree$edge.length[387] + 0.01
ctree$edge.length[386] <- ctree$edge.length[386] - 0.01
ctree$edge.length[610] <- ctree$edge.length[610] + 0.01


### Complexa
ctree$edge.length[171] <- ctree$edge.length[171] + 0.01
ctree$edge.length[170] <- ctree$edge.length[170] - 0.01
ctree$edge.length[377] <- ctree$edge.length[377] + 0.01

ctree$edge.length[172] <- ctree$edge.length[172] + 0.01
ctree$edge.length[170] <- ctree$edge.length[170] - 0.01
ctree$edge.length[365] <- ctree$edge.length[365] + 0.01
ctree$edge.length[377] <- ctree$edge.length[377] + 0.01

ctree$edge.length[337] <- ctree$edge.length[337] + 0.01
ctree$edge.length[335] <- ctree$edge.length[335] - 0.01
ctree$edge.length[360] <- ctree$edge.length[360] + 0.01
ctree$edge.length[361] <- ctree$edge.length[361] + 0.01

ctree$edge.length[336] <- ctree$edge.length[336] + 0.01
ctree$edge.length[335] <- ctree$edge.length[335] - 0.01
ctree$edge.length[361] <- ctree$edge.length[361] + 0.01

# Format asr for plotting------------------------------------------------------

### Replace ctree.nodes in f.anc[[2]] with updated nodes and drop data
# that was pruned from tree
fanc2 <- f.anc[[2]]
setnames(fanc2, c("AS", "ZS", "AF", "ZF", "ctree.nodes"))
fanc2 <- na.omit(fanc2)
fanc2 <- fanc2[, .(AS = mean(AS), ZS = mean(ZS), AF = mean(AF), ZF = mean(ZF)), 
               by = ctree.nodes]
# Removing NAs puts nodes out of order. Reorder them by setting key
setkey(fanc2, ctree.nodes)
# Update ctree.nodes
fanc2[, ctree.nodes := nodes[,2]]
# Reorder by ctree.nodes
setkey(fanc2, ctree.nodes)
# Remove nodes that were dropped from tree
fanc2 <- fanc2[-which(ctree.nodes %in% c(579, 858, 859)),]

### Format f.anc to look like asrs generated from ancRECON so I can use the 
### treePaint function I wrote
# 2 rate
fanc2[, ctree.nodes := NULL]
fanc2 <- list(fanc2)
names(fanc2) <- "lik.anc.states"



# Paint branches---------------------------------------------------------------

# Make color vectors for plotting combined rate cats
# 1 rate
r1.cols <- c("grey", "#084594", "#E6AB02")
names(r1.cols) <- c(1, "A", "Z")

# f2.cols <- c("grey", "grey", "grey", "red", "red")
# names(f2.cols) <- c(1, "AS", "ZS", "AF", "ZF")

# Fast categories
# rc.cols <- c("grey", "red")
# names(rc.cols) <- c(1, "F")
fast.cols <- c("#bdbdbd", "#fc9272")
names(fast.cols) <- c("Stable", "Labile")

# Regular colors for plotting all rate cats
# 2 rate
r2.cols <- c("grey", "#084594", "#E6AB02", "#1B9E77", "#F781BF")
names(r2.cols) <- c(1, "AS", "ZS", "AF", "ZF")

# Make vector of trait colors
trait.cols <- traits[, state]
names(trait.cols) <- traits[, taxa]
trait.cols[trait.cols == "A"] <- "#084594"
trait.cols[trait.cols == "Z"] <- "#E6AB02"


### Paint trees with azoox/zoox recons
# 2 rate
og2.asr <- azPaint(ctree, anc = fanc2$lik.anc.states, rate.cat = 2)

### Paint summed fast categories on trees
# 2 rate
rc2.asr <- maxFast(ctree, fanc2$lik.anc.states, rate.cat = 2)

### Paint trees with rate categories
# 2 rate
f2.asr <- treePaint(tree = ctree, anc = fanc2, rate.cat = 2)

### Paint tree with rate cats based on max p at each node
m2.asr <- maxPaint(tree = ctree, anc = fanc2, rate.cat = 2)

### Paint tree with different cutoff values
f5.asr <- ptreePaint(tree = ctree, anc = fanc2, rate.cat = 2, p = 0.5)

### Plot-----------------------------------------------------------------------

# Plot binary asr
pdf(file = "binary_asr.pdf", width = 5, height = 5)
plotSimmap(og2.asr$tree, type = "fan", ftype = "off", r1.cols, lwd = 0.5)
dev.off()

# Plot binary with names for reference 
pdf(file = "binary_reference.pdf", width = 60.84, height = 40.63)
plotSimmap(og2.asr$tree, type = "fan", fsize = 0.3, r1.cols, lwd = 0.5)
nodelabels(pie = as.matrix(og2.asr$nodeprobs), cex = 0.05,
           piecol = r1.cols[-1])
tiplabels(pch = 19, cex = 0.5, col = trait.cols[ctree$tip.label])
dev.off()

# Plot summed fast cats
pdf(file = "fast_cats.pdf", width = 5, height = 5)
plotSimmap(rc2.asr, type = "fan", ftype = "off", fast.cols, lwd = 0.5)
dev.off()

# Plot all 4 categories
pdf(file = "all4_asr.pdf", width = 5, height = 5)
plotSimmap(f2.asr, type = "fan", ftype = "off", r2.cols, lwd = 0.5)
dev.off()

pdf(file = "all4_maxp_asr.pdf", width = 5, height = 5)
plotSimmap(m2.asr, type = "fan", ftype = "off", r2.cols, lwd = 0.5)
dev.off()

pdf(file = "all4_maxp_reference.pdf", width = 60.84, height = 40.63)
plotSimmap(m2.asr, type = "fan", fsize = 0.3, r2.cols, lwd = 1)
nodelabels(pie = as.matrix(fanc2$lik.anc.states), cex = 0.05,
           piecol = r2.cols[-1])
tiplabels(pch = 19, cex = 0.5, col = trait.cols[ctree$tip.label])
dev.off()



# plotSimmap(f5.asr, type = "fan", fsize = 0.3, r2.cols, lwd = 1)
# nodelabels(pie = as.matrix(fanc2$lik.anc.states), cex = 0.05,
#            piecol = c("#0F4FA8", "#FF9F00", "#00A480", "red"))
# tiplabels(pch = 21, bg = trait.cols[ctree$tip.label], cex = 0.5, lwd = 0.00001)

