# Load packages
library(data.table)
library(phytools)
library(phangorn)
library(here)

# Set working directory
setwd(here("analysis_rerun/results/mtree/species_table"))

# Suppress scientific notation
options(scipen = 999)


# Read in tree and traits-------------------------------------------------------
mtree <- read.nexus(here("data/updated_trees_traits/mtree_traits",
                         "mtree.trees"))

traits <- fread(file = here("data/updated_trees_traits/mtree_traits",
                            "mtree_traits_B_as_Z.csv"),
                header = FALSE, col.names = c("taxa", "state"))


# Format traits and tip labels--------------------------------------------------
# Drop taxa without data from tree
mtree <- lapply(mtree, drop.tip, tip = traits[state == "-", taxa])
class(mtree) <- "multiPhylo"

# Remove taxa with missing data from traits
traits <- traits[state != "-"]


# Read in asr for each tree-----------------------------------------------------
anc2 <- readRDS(here("analysis_rerun/mtree_asr", 
                     "mtree-asr.rds"))

# Read in family key
fkey <- fread(here("data/family_key", "family_key.csv"))


### Get tip probs--------------------------------------------------------------
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

tip.probs[state == "Z", state := "Zooxanthellate"]
tip.probs[state == "A", state := "Azooxanthellate"]
tip.probs[state == "B", state := "Facultative"]


# Format tip.probs
tip.probs[, genus := NULL]

setnames(tip.probs, c("Species", "Azoox Stable", "Zoox Stable", 
                      "Azoox Labile", "Zoox Labile", 
                      "Family", "Observed State"))

setcolorder(tip.probs, c("Family", "Species", "Observed State","Azoox Stable", 
                         "Zoox Stable", "Azoox Labile", "Zoox Labile"))

# Drop underscore from species names 
tip.probs[, Species := gsub("_", " ", tip.probs[, Species])] 

fwrite(tip.probs, file = "mtree_2rate_tip_probs.csv")