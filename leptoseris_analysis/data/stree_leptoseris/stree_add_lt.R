# Load packages
library(phytools)
library(data.table)
library(phangorn)
library(here)

# Set working directory
setwd(here("leptoseris_analysis/data/stree_leptoseris"))


# Read in tree and traits------------------------------------------------------
stree <- read.nexus(here("data/updated_trees_traits/stree_traits", 
                         "stree.trees"))

traits <- fread(here("data/updated_trees_traits/stree_traits",
                     "stree_traits_B_as_Z.csv"),
                header = FALSE, col.names = c("taxa", "state"))


# Add leptoseris troglodyta to tree--------------------------------------------
stree <- .uncompressTipLabel(stree)
set.seed(22)
for (i in 1:length(stree)) {
  stree[[i]] <- add.species.to.genus(stree[[i]], 
                                     species = "Leptoseris_troglodyta",
                                     where = "random")
}


# Add Leptoseris troglodyta to traits-------------------------------------------
lt <- data.table("taxa" = c("Leptoseris_troglodyta"),
                 "state" = c("A"))

traits <- rbind(traits, lt)


# Write to file-----------------------------------------------------------------
write.nexus(stree, file = "stree_lt.trees")

fwrite(traits, "stree_traits_B_as_Z_lt.csv", col.names = FALSE, sep = "\t")