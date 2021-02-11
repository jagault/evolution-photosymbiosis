# Load packages
library(data.table)
library(phytools)
library(phangorn)
library(here)

# Set working directory
setwd(here("figures/stree"))

# Source helper functions for plotting
source(here("R/plotting_functions.R"))


# Read in tree and traits------------------------------------------------------
stree <- read.nexus(here("data/updated_trees_traits/stree_traits", 
                                  "stree.trees"))

traits <- fread(here("data/updated_trees_traits/stree_traits",
                     "stree_traits_B_as_Z.csv"),
                header = FALSE, col.names = c("taxa", "state"))


# Format traits and tip labels-------------------------------------------------
# Drop taxa without data from tree
stree <- lapply(stree, drop.tip, tip = traits[state == "-", taxa])
class(stree) <- "multiPhylo"

# Remove taxa with missing data from traits
traits <- traits[state != "-"]


# Sample 100 trees used to estimate rates--------------------------------------
set.seed(1)
stree <- sample(stree, 100, replace = FALSE)


# Read in corHMM runs----------------------------------------------------------
# Read in asr summarized across all 100 trees
f.anc <- readRDS(here("analysis/stree_asr", "stree_subsample_nodeframes.rds"))

# Read in 100 individual asrs
anc3 <- readRDS(here("analysis/stree_asr", "stree_3rate_anc.rds"))


# Prep consensus tree for plotting---------------------------------------------

# Make consensus tree
ctree <- consensus(stree, p = 0.95)

# Compute branch lengths of ctree for mapping states. The actual branch lengths
# don't matter. They just need values for helper functions to work.
ctree <- compute.brlen(ctree)

### Lengthen some short branches for presentation
# Plot with edge labels to see which to change
# plotTree(ctree, ftype = "off", type = "fan")
# edgelabels(frame = "n")

# Root
ctree$edge.length[1] <- ctree$edge.length[1] + 0.3
ctree$edge.length[2122] <- ctree$edge.length[2122] + 0.3

### Robusta
ctree$edge.length[1041] <- ctree$edge.length[1041] + 0.01
ctree$edge.length[1040] <- ctree$edge.length[1040] - 0.01
ctree$edge.length[2121] <- ctree$edge.length[2121] + 0.01

ctree$edge.length[1043] <- ctree$edge.length[1043] + 0.01
ctree$edge.length[1042] <- ctree$edge.length[1042] - 0.01
ctree$edge.length[2091] <- ctree$edge.length[2091] + 0.01

ctree$edge.length[1050] <- ctree$edge.length[1050] + 0.01
ctree$edge.length[1049] <- ctree$edge.length[1049] - 0.01
ctree$edge.length[1480] <- ctree$edge.length[1480] + 0.01

ctree$edge.length[1052] <- ctree$edge.length[1052] + 0.01
ctree$edge.length[1051] <- ctree$edge.length[1051] - 0.01
ctree$edge.length[1425] <- ctree$edge.length[1425] + 0.01

ctree$edge.length[1053] <- ctree$edge.length[1053] + 0.005
ctree$edge.length[1052] <- ctree$edge.length[1052] - 0.005
ctree$edge.length[1424] <- ctree$edge.length[1424] + 0.005

ctree$edge.length[1046] <- ctree$edge.length[1046] + 0.005
ctree$edge.length[1045] <- ctree$edge.length[1045] - 0.005
ctree$edge.length[1865] <- ctree$edge.length[1865] + 0.005

ctree$edge.length[1788] <- ctree$edge.length[1788] + 0.01
ctree$edge.length[1787] <- ctree$edge.length[1787] - 0.01
ctree$edge.length[1858] <- ctree$edge.length[1858] + 0.01


### Complexa
ctree$edge.length[3] <- ctree$edge.length[3] + 0.01
ctree$edge.length[2] <- ctree$edge.length[2] - 0.01
ctree$edge.length[1039] <- ctree$edge.length[1039] + 0.01

ctree$edge.length[561] <- ctree$edge.length[561] + 0.01
ctree$edge.length[560] <- ctree$edge.length[560] - 0.01
ctree$edge.length[1036] <- ctree$edge.length[1036] + 0.01

ctree$edge.length[562] <- ctree$edge.length[562] + 0.01
ctree$edge.length[560] <- ctree$edge.length[560] - 0.01
ctree$edge.length[1036] <- ctree$edge.length[1036] + 0.01
ctree$edge.length[1023] <- ctree$edge.length[1023] + 0.01

ctree$edge.length[957] <- ctree$edge.length[957] + 0.01
ctree$edge.length[956] <- ctree$edge.length[956] - 0.01
ctree$edge.length[1016] <- ctree$edge.length[1016] + 0.01

ctree$edge.length[958] <- ctree$edge.length[958] + 0.01
ctree$edge.length[956] <- ctree$edge.length[956] - 0.01
ctree$edge.length[1015] <- ctree$edge.length[1015] + 0.01
ctree$edge.length[1016] <- ctree$edge.length[1016] + 0.01


### Format ancestral state recon-----------------------------------------------

# Format f.anc to look like asrs generated from ancRECON for use with 
# ploting functions 
fanc3 <- f.anc[[3]]
setnames(fanc3, c("AS", "ZS", "AM", "ZM", "AF", "ZF", "ctree.nodes"))
fanc3 <- fanc3[, .(AS = mean(AS), ZS = mean(ZS), 
                   AM = mean(AM), ZM = mean(ZM), 
                   AF = mean(AF), ZF = mean(ZF)), 
               by = ctree.nodes]
fanc3[, ctree.nodes := NULL]
fanc3 <- list(fanc3)
names(fanc3) <- "lik.anc.states"


### Paint trees with rate categories-------------------------------------------

# Based on probability cutoff at nodes
f3.asr <- treePaint(tree = ctree, anc = fanc3, rate.cat = 3)

# Based on maximum probability at nodes
m3.asr <- maxPaint(tree = ctree, anc = fanc3, rate.cat = 3)


### Plot-----------------------------------------------------------------------

# Define branch colors for plotting
r3.cols <- c("grey", "#084594", "#F781BF", "#1B9E77", 
             "#E6AB02", "#7570B3", "#E41A1C")
names(r3.cols) <- c(1, "AS", "ZS", "AM", "ZM", "AF", "ZF")

# Make vector of trait colors for tips
trait.cols <- traits[, state]
names(trait.cols) <- traits[, taxa]
trait.cols[trait.cols == "A"] <- "#084594"
trait.cols[trait.cols == "Z"] <- "#E6AB02"

# Plot based on maximum probability at nodes
pdf(file = "stree_consensus_asr.pdf", width = 5, height = 5)
plotSimmap(m3.asr, type = "fan", ftype = "off", r3.cols, lwd = 0.5)
dev.off()

# Plot large figure with tip labels and pie charts at nodes 
pdf(file = "stree_consensus_asr_tiplabels.pdf", width = 60.84, height = 40.63)
plotSimmap(m3.asr, type = "fan", fsize = 0.3, r3.cols, lwd = 1)
nodelabels(pie = as.matrix(fanc3$lik.anc.states), cex = 0.05,
         piecol = r3.cols[-1])
tiplabels(pch = 19, cex = 0.5, col = trait.cols[ctree$tip.label])
dev.off()

# Plot individual asrs---------------------------------------------------------

# Format asr's for plotting function 
anc3 <- lapply(anc3, "[[", "lik.anc.states")
anc3 <- lapply(anc3, data.table)
anc3 <- lapply(anc3, setnames, c("AS", "ZS", "AM", "ZM", "AF", "ZF"))
names(anc3) <- rep("lik.anc.states", length(anc3))

# Paint all 100 trees
all.paint <- list()
for (i in 1:length(stree)){
  all.paint[[i]] <- maxPaint(stree[[i]], anc3[i], rate.cat = 3)
}

# Save as a pdf
vec <- c(1:100)
pdf(file = "stree_individual_asr.pdf", width = 50, height = 50)
layout(matrix(1:100, 5, 20, byrow = TRUE))
for(i in 1:length(all.paint)){
  plotSimmap(all.paint[[i]], ftype = "off", r3.cols, mar = c(0.1, 0.1, 1.5, 0.1))
  mtext(vec[i])
}
dev.off()