# Load packages
library(phytools)
library(data.table)
library(phangorn)
library(here)

# Set working directory
setwd(here("analysis_rerun/stree_corHMM/stree_2rate/"))

# Define path to results, will be used to read in files in for-loop below
path <- c("~/Dropbox/osg_results/stree/2-rate/")

# Initialize vector for number of trees
t <- c(1:1000)
# Initialize vector for number of replicates per tree
r <- c(1:99)

# Create empty vectors and lists to hold results
files <- c()
results <- list()
maxl <- list()
liks <- c()
reps <- c()

# Loop through results and save max-likelihood results for each tree
for (i in 1:length(t)){
  for (j in 1:length(r)){
    
    files[j] <- paste(path, "stree-2rate-t", i, "-r", j, ".rds", sep = "")
    results[[j]] <- readRDS(files[j])
    liks[j] <- results[[j]]$loglik
  }
  
  maxl[[i]] <- results[[which.max(liks)]]
  reps[i] <- length(liks)
}

# Save list of results 
saveRDS(maxl, "stree-2rate.rds")