# Load packages
library(phytools)
library(data.table)
library(phangorn)
library(here)

# Set working directory
setwd(here("analysis_rerun/mtree_corHMM/mtree_1rate/"))

# Define path to results, will be used to read in files in for-loop below
path <- c("~/Dropbox/osg_results/stree/mtree/mtree_1rate")

# Initialize vector for number of trees
t <- c(1:3361)

# Create empty vectors and lists to hold results
files <- c()
results <- list()
maxl <- list()
liks <- c()
reps <- c()

# Loop through results and save max-likelihood results for each tree
for (i in 1:length(t)){
  files <- list.files(path = paste(path, "/1-rate_", i, sep = ""), 
                      pattern = ".rds", full.names = T)
  
  results <- lapply(files, readRDS)
  
  liks <- sapply(results, "[", "loglik", simplify = T)
  
  maxl[[i]] <- results[[which.max(liks)]]
  reps[i] <- length(liks)
}

# Save list of results 
saveRDS(maxl, "mtree-1rate.rds")