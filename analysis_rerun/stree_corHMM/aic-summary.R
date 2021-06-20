# Load packages
library(phytools)
library(data.table)
library(phangorn)
library(here)

# Set working directory
setwd(here("analysis_rerun/stree_corHMM/"))

# Read in files
stree.1rate <- readRDS("stree_1rate/stree-1rate.rds")
stree.2rate <- readRDS("stree_2rate/stree-2rate.rds")
stree.3rate <- readRDS("stree_3rate/stree-3rate.rds")
stree.4rate <- readRDS("stree_4rate/stree-4rate.rds")

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

# AICc weights
daic1 <- 335.6 - 277.5
daic2 <- 286.6 - 277.5
daic3 <- 277.5 - 277.5
daic4 <- 287.5 - 277.5

rl1 <- exp(-.5*daic1)
rl2 <- exp(-.5*daic2)
rl3 <- exp(-.5*daic3)
rl4 <- exp(-.5*daic4)

aicw1 <- rl1/sum(rl1, rl2, rl3, rl4)
aicw2 <- rl2/sum(rl1, rl2, rl3, rl4)
aicw3 <- rl3/sum(rl1, rl2, rl3, rl4)
aicw4 <- rl4/sum(rl1, rl2, rl3, rl4)

# Summary table
aicw1
aicw2
aicw3
aicw4