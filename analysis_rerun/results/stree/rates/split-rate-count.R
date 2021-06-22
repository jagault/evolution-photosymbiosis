# This script counts the number of trees for which the rate estimate of 
# ZS to AS is 0. This should give the number of trees for which ZS is the 
# absorbing state vs ZM. 

# Load packages
library(here)

# Set working directory
setwd(here("analysis_rerun/results/stree/rates/"))

# Suppress scientific notation
options(scipen = 999)

# Source helper functions
source(here("R/aic-rate-summary.R"))


# Read in corHMM runs----------------------------------------------------------
r3 <- readRDS(here("analysis_rerun/stree_corHMM/stree_3rate/stree-3rate.rds"))

# Make rate vectors-------------------------------------------------------------
par3 <- sapply(r3, "[", "solution")

zs.to.as <- sapply(par3, "[", c(2, 1))
zs.to.as <- na.omit(zs.to.as)

sum(zs.to.as < 0.0000001)
# About 16% of the trees have ZS as absorbing state. This is close to the 
# 19% from the original run. 