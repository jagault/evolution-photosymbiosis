# Load packages
library(phytools)
library(data.table)
library(phangorn)
library(here)

# Set working directory
setwd(here("leptoseris_analysis/analysis/stree_corHMM/"))

# Define path to results, will be used to read in files in for-loop below
path <- c("~/Dropbox/osg_results/leptoseris/analysis/stree_corHMM/stree_2rate")

# Initialize vector for number of trees
t <- c(1:1000)

# Create empty vectors and lists to hold results
files <- c()
results <- list()
maxl <- list()
liks <- c()
reps <- c()

# Loop through results and save max-likelihood results for each tree
for (i in 1:length(t)){
  files <- list.files(path = paste(path, "/2-rate_", i, sep = ""), 
                      pattern = ".rds", full.names = T)
  
  results <- lapply(files, readRDS)
  
  liks <- sapply(results, "[", "loglik", simplify = T)
  
  maxl[[i]] <- results[[which.max(liks)]]
  reps[i] <- length(liks)
}

# Save list of results 
saveRDS(maxl, "stree-2rate.rds")