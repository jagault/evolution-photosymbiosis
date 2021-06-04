# Load packages
library(phytools)
library(data.table)
library(phangorn)
library(here)

# Set working directory
setwd(here("leptoseris_analysis/data/mtree_leptoseris"))

# Read in tree and traits------------------------------------------------------
mtree <- read.nexus(here("data/updated_trees_traits/mtree_traits", 
                         "mtree.trees"))

traits <- fread(here("data/updated_trees_traits/mtree_traits",
                     "mtree_traits_B_as_Z.csv"),
                header = FALSE, col.names = c("taxa", "state"))

# Add leptoseris troglodyta to tree--------------------------------------------
mtree <- .uncompressTipLabel(mtree)
set.seed(22)
for (i in 1:length(mtree)) {
  mtree[[i]] <- add.species.to.genus(mtree[[i]], 
                                     species = "Leptoseris_troglodyta",
                                     where = "random")
}

# Add Leptoseris troglodyta to traits-------------------------------------------
lt <- data.table("taxa" = c("Leptoseris_troglodyta"),
                   "state" = c("A"))
  
traits <- rbind(traits, lt)

# Test monophyly of Leptoseris--------------------------------------------------
# ll <- copy(traits)
# ll[, c("genus", "s1", "s2") := tstrsplit(taxa, "_")]
# 
# is.monophyletic(mtree[[1]], ll[genus == "Leptoseris", taxa], plot = T)
# 
# m <- c()
# for (i in 1:length(mtree)){
#   m[i] <- is.monophyletic(mtree[[i]], ll[genus == "Leptoseris", taxa])
# }
# m
# unique(m)
# 
# plot(mtree[[1]], cex = 0.2)
# edgelabels(frame = "n")

# Write to file-----------------------------------------------------------------
write.nexus(mtree, file = "mtree_lt.trees")

fwrite(traits, "mtree_traits_B_as_Z_lt.csv", col.names = FALSE, sep = "\t")
