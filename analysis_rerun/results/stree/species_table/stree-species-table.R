# Load packages
library(data.table)
library(phytools)
library(phangorn)
library(here)

# Set working directory
setwd(here("analysis_rerun/results/stree/species_table"))

# Suppress scientific notation
options(scipen = 999)

# Read in tree and traits-------------------------------------------------------
stree <- read.nexus(here("data/updated_trees_traits/stree_traits",
                         "stree.trees"))

traits <- fread(file = here("data/updated_trees_traits/stree_traits",
                            "stree_traits.csv"),
                header = FALSE, col.names = c("taxa", "state"))


# Format traits and tip labels--------------------------------------------------
# Drop taxa without data from tree
stree <- lapply(stree, drop.tip, tip = traits[state == "-", taxa])
class(stree) <- "multiPhylo"


# Read in asr------------------------------------------------------------------
# Read in 1000 individual asrs
anc3 <- readRDS(here("analysis_rerun/stree_asr", 
                     "stree-asr.rds"))


# Read in updated list of taxa for family names
fkey <- fread(here("data/family_key", "family_key.csv"))


### Get tip probs--------------------------------------------------------------
tip.probs <- lapply(anc3, "[[", "lik.tip.states")
tip.probs <- lapply(tip.probs, data.table)

tip.index <- lapply(stree, "[[", "tip.label")
tip.index <- lapply(tip.index, data.table)
lapply(tip.index, setnames, "species")

tip.probs <- rbindlist(mapply(cbind, tip.probs, tip.index, SIMPLIFY = FALSE))
setnames(tip.probs, c("AS", "ZS", "AM", "ZM", "AF", "ZF", "species"))

tip.probs <- tip.probs[, .(AS = round(mean(AS), 4), ZS = round(mean(ZS), 4),
                           AM = round(mean(AM), 4), ZM = round(mean(ZM), 4),
                           AF = round(mean(AF), 4), ZF = round(mean(ZF), 4)), 
                       by = species]


# Add family names-------------------------------------------------------------
# Format fkey and drop outgroup from mtree
fkey <- fkey[, .(Family, Genus)]
setnames(fkey, c("family", "genus"))
fkey <- fkey[!(genus %in% c("Discosoma", "Ricordea", "Rhodactis"))]

# Add genus column to tip.probs
tip.probs[, genus := unlist(lapply(strsplit(tip.probs[, species], "_"), "[[", 1))]


# Merge tip.probs and taxa.list
setkey(fkey, genus)
setkey(tip.probs, genus)

tip.probs <- merge(tip.probs, fkey)


# Add trait column 
setnames(traits, c("species", "state"))
setkey(traits, species)
setkey(tip.probs, species)

tip.probs <- merge(tip.probs, traits)

tip.probs[state == "Z", state := "Z"]
tip.probs[state == "A", state := "AZ"]
tip.probs[state == "B", state := "F"]


# Format tip.probs
tip.probs[, genus := NULL]

setnames(tip.probs, c("Species", "Azoox Stable", "Zoox Labile", 
                      "Azoox Labile", "Zoox Stable", "Azoox Volatile", 
                      "Zoox Volatile", "Family", "Observed State"))

setcolorder(tip.probs, c("Family", "Species", "Observed State","Azoox Stable", 
                         "Zoox Stable", "Azoox Labile", "Zoox Labile", 
                         "Azoox Volatile", "Zoox Volatile"))

# Drop underscore from species names 
tip.probs[, Species := gsub("_", " ", tip.probs[, Species])] 


fwrite(tip.probs, file = "stree-species-table.csv")
