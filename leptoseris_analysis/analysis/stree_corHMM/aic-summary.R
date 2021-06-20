# Load packages
library(phytools)
library(data.table)
library(phangorn)
library(here)

# Set working directory
setwd(here("leptoseris_analysis/analysis/stree_corHMM/"))

# Read in files
stree.1rate <- readRDS("stree-1rate.rds")
stree.2rate <- readRDS("stree-2rate.rds")
stree.3rate <- readRDS("stree-3rate.rds")
stree.4rate <- readRDS("stree-4rate.rds")

# Select AICc scores
aic.1rate <- unlist(sapply(stree.1rate, "[", "AICc"))
aic.2rate <- unlist(sapply(stree.2rate, "[", "AICc"))
aic.3rate <- unlist(sapply(stree.3rate, "[", "AICc"))
aic.4rate <- unlist(sapply(stree.4rate, "[", "AICc"))

# Summarize
summary(aic.1rate)
summary(aic.2rate)
summary(aic.3rate)
summary(aic.4rate)

# Calculate standard error
std <- function(x) sd(x)/sqrt(length(x))

std(aic.1rate)
std(aic.2rate)
std(aic.3rate)
std(aic.4rate)

# Boxplot
aic <- cbind(aic.1rate, aic.2rate, aic.3rate, aic.4rate)
boxplot(aic)
