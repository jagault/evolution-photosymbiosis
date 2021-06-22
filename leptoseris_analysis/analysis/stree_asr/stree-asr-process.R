# Load packages
library(here)

# Set working directory
setwd(here("leptoseris_analysis/analysis/stree_asr/"))

# Define path to results, will be used to read in files in for-loop below
path <- c("~/Dropbox/osg_results/ASR/leptoseris_asr/")

# Initialize vector for number of trees
t <- c(1:1000)

# Create empty vectors and lists to hold results
files <- c()
results <- list()

# Loop through asr results and add them to a single list
for (i in 1:length(t)){
  files <- list.files(path = paste(path, "lepto_stree_asr_", i, sep = ""),
                      pattern = ".rds", full.names = T)
  
  results[[i]] <- readRDS(files)
}

# Save list of results 
saveRDS(results, "stree-asr-lt.rds")