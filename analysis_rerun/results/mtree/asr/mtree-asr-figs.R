# Load packages
library(data.table)
library(phytools)
library(phangorn)
library(here)

# Set working directory
setwd(here("analysis_rerun/results/mtree/asr"))

# Source helper functions for plotting
source(here("R/asr-plot.R"))


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

# Read in corHMM runs----------------------------------------------------------
# Read in asr summarized across all 1000 trees
fanc2 <- readRDS(here("analysis_rerun/mtree_asr", 
                      "mtree-nodeframes.rds"))


# Rotate clades in mtree to match stree----------------------------------------

# Make consensus tree of mtree and stree
mctree <- consensus(mtree, p = 0.95)
sctree <- consensus(stree, p = 0.95)

# Select index of tips from edge matrix
tips <- sctree$edge[which(sctree$edge[,2] <= length(sctree$tip.label)),2]
# Get tip labels using index of tips from edge matrix
tip.index <- sctree$tip.label[tips]
# Select only tips that are in the molecular tree
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

# Drop outgroup
ctree <- drop.tip(ctree, c("Discosoma", "Rhodactis", "Ricordea_florida"))

# Drop Anthemiphyllia 
ctree <- drop.tip(ctree, "Anthemiphyllia_patera")

# Lengthen branches for presentation-------------------------------------------

# Plot with edge labels to see which to change
# plotTree(ctree, ftype = "off", type = "fan")
# edgelabels(frame = "n")

# Root
ctree$edge.length[1] <- ctree$edge.length[1] + 0.3
ctree$edge.length[835] <- ctree$edge.length[835] + 0.3

### Robusta
ctree$edge.length[735] <- ctree$edge.length[735] + 0.01
ctree$edge.length[734] <- ctree$edge.length[734] - 0.01
ctree$edge.length[764] <- ctree$edge.length[764] + 0.01

ctree$edge.length[611] <- ctree$edge.length[611] + 0.01
ctree$edge.length[610] <- ctree$edge.length[610] - 0.01
ctree$edge.length[688] <- ctree$edge.length[688] + 0.01

ctree$edge.length[385] <- ctree$edge.length[385] + 0.01
ctree$edge.length[383] <- ctree$edge.length[383] - 0.01
ctree$edge.length[601] <- ctree$edge.length[601] + 0.01
ctree$edge.length[604] <- ctree$edge.length[604] + 0.01

ctree$edge.length[384] <- ctree$edge.length[384] + 0.01
ctree$edge.length[383] <- ctree$edge.length[383] - 0.01
ctree$edge.length[604] <- ctree$edge.length[604] + 0.01


### Complexa
ctree$edge.length[169] <- ctree$edge.length[169] + 0.01
ctree$edge.length[168] <- ctree$edge.length[168] - 0.01
ctree$edge.length[374] <- ctree$edge.length[374] + 0.01

ctree$edge.length[170] <- ctree$edge.length[170] + 0.01
ctree$edge.length[168] <- ctree$edge.length[168] - 0.01
ctree$edge.length[362] <- ctree$edge.length[362] + 0.01
ctree$edge.length[374] <- ctree$edge.length[374] + 0.01

# ctree$edge.length[337] <- ctree$edge.length[337] + 0.01
# ctree$edge.length[335] <- ctree$edge.length[335] - 0.01
# ctree$edge.length[360] <- ctree$edge.length[360] + 0.01
# ctree$edge.length[361] <- ctree$edge.length[361] + 0.01

ctree$edge.length[334] <- ctree$edge.length[334] + 0.01
ctree$edge.length[333] <- ctree$edge.length[333] - 0.01
ctree$edge.length[358] <- ctree$edge.length[358] + 0.01


# Format asr for plotting------------------------------------------------------

### Replace ctree.nodes in f.anc[[2]] with updated nodes and drop data
# that was pruned from tree
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

### Format f.anc to look like asr's generated from ancRECON so I can use the 
### function I wrote
fanc2[, ctree.nodes := NULL]
fanc2 <- list(fanc2)
names(fanc2) <- "lik.anc.states"


# Paint branches---------------------------------------------------------------

### Paint tree with rate cats based on max p at each node
m2.asr <- maxPaint(tree = ctree, anc = fanc2, rate.cat = 2)


### Plot-----------------------------------------------------------------------

# Define branch colors for plotting
r2.cols <- c("grey", "#084594", "#E6AB02", "#1B9E77", "#F781BF")
names(r2.cols) <- c(1, "AS", "ZS", "AF", "ZF")

# Make vector of trait colors for tips
trait.cols <- traits[, state]
names(trait.cols) <- traits[, taxa]
trait.cols[trait.cols == "A"] <- "#084594"
trait.cols[trait.cols == "Z"] <- "#E6AB02"


# Plot based on maximum probability at nodes
pdf(file = "mtree_consensus_asr.pdf", width = 5, height = 5)
plotSimmap(m2.asr, type = "fan", ftype = "off", r2.cols, lwd = 0.5)
dev.off()

# Plot large figure with tip labels and pie charts at nodes
pdf(file = "mtree_consensus_asr_tiplabels.pdf", width = 60.84, height = 40.63)
plotSimmap(m2.asr, type = "fan", fsize = 0.3, r2.cols, lwd = 1)
nodelabels(pie = as.matrix(fanc2$lik.anc.states), cex = 0.05,
           piecol = r2.cols[-1])
tiplabels(pch = 19, cex = 0.5, col = trait.cols[ctree$tip.label])
dev.off()
