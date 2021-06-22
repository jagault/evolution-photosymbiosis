# Load packages
library(data.table)
library(phytools)
library(phangorn)
library(here)

# Set working directory
setwd(here("leptoseris_analysis/results/stree/asr"))

# Source helper functions for plotting
source(here("R/asr-plot.R"))


# Read in tree and traits-------------------------------------------------------
stree <- read.nexus(here("leptoseris_analysis/data/stree_leptoseris",
                         "stree_lt.trees"))

traits <- fread(file = here("leptoseris_analysis/data/stree_leptoseris",
                            "stree_traits_B_as_Z_lt.csv"),
                header = FALSE, col.names = c("taxa", "state"))


# Format traits and tip labels--------------------------------------------------
# Drop taxa without data from tree
stree <- lapply(stree, drop.tip, tip = traits[state == "-", taxa])
class(stree) <- "multiPhylo"


# Read in corHMM runs----------------------------------------------------------
# Read in asr summarized across all 1000 trees
fanc3 <- readRDS(here("leptoseris_analysis/analysis/stree_asr", 
                      "stree-lt-nodeframes.rds"))

# Read in 1000 individual asrs
anc3 <- readRDS(here("leptoseris_analysis/analysis/stree_asr", 
                     "stree-asr-lt.rds"))

# Prep consensus tree for plotting---------------------------------------------
# Make consensus tree
ctree <- consensus(stree, p = 0.95)

# Compute branch lengths of ctree for mapping states. The actual branch lengths
# don't matter. They just need values for helper functions to work.
ctree <- compute.brlen(ctree)


### Format ancestral state recon-----------------------------------------------

# Format f.anc to look like asrs generated from ancRECON for use with 
# ploting functions 
setnames(fanc3, c("AS", "ZS", "AM", "ZM", "AF", "ZF", "ctree.nodes"))
fanc3 <- fanc3[, .(AS = mean(AS), ZS = mean(ZS), 
                   AM = mean(AM), ZM = mean(ZM), 
                   AF = mean(AF), ZF = mean(ZF)), 
               by = ctree.nodes]
fanc3[, ctree.nodes := NULL]
fanc3 <- list(fanc3)
names(fanc3) <- "lik.anc.states"


### Paint trees with rate categories-------------------------------------------
# Based on maximum probability at nodes
m3.asr <- maxPaint(tree = ctree, anc = fanc3, rate.cat = 3)


### Plot-----------------------------------------------------------------------

# # Define branch colors for plotting
# r3.cols <- c("grey", "#084594", "#F781BF", "#1B9E77", 
#              "#E6AB02", "#7570B3", "#E41A1C")
# names(r3.cols) <- c(1, "AS", "ZS", "AM", "ZM", "AF", "ZF")

# Switch AS and AM to reflect that they have switched roles relative to the 
# original supertree analysis
# Define branch colors for plotting
r3.cols <- c("grey", "#084594", "#E6AB02", "#1B9E77", 
             "#F781BF", "#7570B3", "#E41A1C")
names(r3.cols) <- c(1, "AS", "ZS", "AM", "ZM", "AF", "ZF")


# Make vector of trait colors for tips
trait.cols <- traits[, state]
names(trait.cols) <- traits[, taxa]
trait.cols[trait.cols == "A"] <- "#084594"
trait.cols[trait.cols == "Z"] <- "#E6AB02"

# Plot based on maximum probability at nodes
pdf(file = "stree-lt-consensus-asr.pdf", width = 5, height = 5)
plotSimmap(m3.asr, type = "fan", ftype = "off", r3.cols, lwd = 0.5)
dev.off()

# Plot large figure with tip labels and pie charts at nodes 
pdf(file = "stree-lt-consensus-asr-tiplabels.pdf", width = 60.84, height = 40.63)
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
vec <- c(1:1000)
pdf(file = "stree_individual_asr.pdf", width = 50, height = 500)
layout(matrix(1:1000, 50, 20, byrow = TRUE))
for(i in 1:length(all.paint)){
  plotSimmap(all.paint[[i]], ftype = "off", r3.cols, 
             mar = c(0.1, 0.1, 1.5, 0.1))
  mtext(vec[i])
}
dev.off()