# Load packages
library(data.table)
library(ggplot2)
library(here)

# Set working directory
setwd(here("analysis_rerun/results/stree/rates/"))

# Suppress scientific notation
options(scipen = 999)

# Source helper functions
source(here("R/aic-rate-summary.R"))


# Read in corHMM runs----------------------------------------------------------
r1 <- readRDS(here("analysis_rerun/stree_corHMM/stree_1rate/stree-1rate.rds"))
r2 <- readRDS(here("analysis_rerun/stree_corHMM/stree_2rate/stree-2rate.rds"))
r3 <- readRDS(here("analysis_rerun/stree_corHMM/stree_3rate/stree-3rate.rds"))
r4 <- readRDS(here("analysis_rerun/stree_corHMM/stree_4rate/stree-4rate.rds"))


# Summarize rates and model selection------------------------------------------
sr1 <- corSumm(r1, rate.cat = 1)
sr2 <- corSumm(r2, rate.cat = 2)
sr3 <- corSumm(r3, rate.cat = 3)
sr4 <- corSumm(r4, rate.cat = 4)


### Model selection~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
aic.table <- data.table(model = c("TH", "HRM+2", "HRM+3", "HRM+4"),
                        num.rate.class = c(1, 2, 3, 4),
                        num.param = c(2, 8, 14, 20),
                        mean.AICc = c(round(sr1$aic.summ[, mean], 3),
                                      round(sr2$aic.summ[, mean], 3),
                                      round(sr3$aic.summ[, mean], 3),
                                      round(sr4$aic.summ[, mean], 3)),
                        AICc.weight.percent = c(round(aicw1*100, 3), 
                                                round(aicw2*100, 3), 
                                                round(aicw3*100, 3),
                                                round(aicw4*100, 3)))
aic.table
fwrite(aic.table, file = "aic_weights.csv")


### Model selection plots-------------------------------------------------------
aics <- cbind(sr1$aic, sr2$aic, sr3$aic, sr4$aic)
setnames(aics, c("TH", "HRM+2", "HRM+3", "HRM+4"))
aics <- melt(aics, variable = "model", value = "AICc")


# pdf(file = "aic_histogram.pdf", width = 6.83, height = 4.31)
ggplot(aics) + geom_histogram(aes(x = AICc, fill = model), 
                              alpha = .5, binwidth = 2)
# dev.off()


### Line plot of mean aics
maic <- data.table(aic = c(sr1$aic.summ[, mean], sr2$aic.summ[, mean], 
                           sr3$aic.summ[, mean], sr4$aic.summ[, mean]),
                   stderr = c(sr1$aic.summ[, std.err], sr2$aic.summ[, std.err],
                              sr3$aic.summ[, std.err], sr4$aic.summ[, std.err]),
                   model = c("TH", "HRM+2", "HRM+3", "HRM+4"))



# pdf(file = "aic_curve.pdf", width = 6.83, height = 4.31)
ggplot(maic, aes(x = model, y = -aic, group = 1)) + geom_point() + geom_path() +
  scale_x_discrete(limits= c("TH", "HRM+2", "HRM+3", "HRM+4")) + 
  geom_errorbar(aes(ymin = -aic - stderr, ymax = -aic + stderr), width = 0.1) +
  xlab("Model") + 
  ylab("-AIC") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# dev.off()



### Line plot for log-likelihoods
lik1 <- unlist(sapply(rate1, "[", "loglik"))
lik2 <- unlist(sapply(rate2, "[", "loglik"))
lik3 <- unlist(sapply(rate3, "[", "loglik"))
lik4 <- unlist(sapply(rate4, "[", "loglik"))


loglik <- data.table("TH" = lik1, "HRM+2" = lik2, "HRM+3" = lik3, 
                     "HRM+4" = lik4)

loglik <- melt(loglik, value = "loglik", variable.name = "model")

loglik <- loglik[, .(loglik = mean(loglik), std.err = sd(loglik)/.N), 
                 by = "model"]

# pdf(file = "likelihood_curve.pdf", width = 6.83, height = 4.31)
ggplot(loglik, aes(x = model, y = loglik, group = 1)) + geom_point() + 
  geom_path() +
  scale_x_discrete(limits= c("TH", "HRM+2", "HRM+3", "HRM+4")) + 
  geom_errorbar(aes(ymin = loglik - std.err, ymax = loglik + std.err), 
                width = 0.1) + 
  ylab("-LnLik") + 
  xlab("Model") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# dev.off()


### Rates-----------------------------------------------------------------------
sr1q <- apply(sr1$par[,2:ncol(sr1$par)], 2, bootMed, n = 100000)
sr1q <- data.table(t(sr1q), keep.rownames = T)


sr2q <- apply(sr2$par[,2:ncol(sr2$par)], 2, bootMed, n = 100000)
sr2q <- data.table(t(sr2q), keep.rownames = T)


sr3q <- apply(sr3$par[,2:ncol(sr3$par)], 2, bootMed, n = 100000)
sr3q <- data.table(t(sr3q), keep.rownames = T)


sr3q <- apply(sr3$par[,2:ncol(sr3$par)], 2, bootMed, n = 100000)
sr3q <- data.table(t(sr3q), keep.rownames = T)


sr4q <- apply(sr4$par[,2:ncol(sr4$par)], 2, bootMed, n = 100000)
sr4q <- data.table(t(sr4q), keep.rownames = T)


setnames(sr1q, c("transition", "lci", "median", "uci"))
setnames(sr2q, c("transition", "lci", "median", "uci"))
setnames(sr3q, c("transition", "lci", "median", "uci"))
setnames(sr4q, c("transition", "lci", "median", "uci"))


sr1q[, c("lci", "median", "uci") := list(round(lci*1000, 3), 
                                         round(median*1000, 3), 
                                         round(uci*1000, 3))]

sr2q[, c("lci", "median", "uci") := list(round(lci*1000, 3), 
                                         round(median*1000, 3), 
                                         round(uci*1000, 3))]

sr3q[, c("lci", "median", "uci") := list(round(lci*1000, 3), 
                                         round(median*1000, 3), 
                                         round(uci*1000, 3))]

sr4q[, c("lci", "median", "uci") := list(round(lci*1000, 3), 
                                         round(median*1000, 3), 
                                         round(uci*1000, 3))]

sr1q[, c("minus", "plus") := list(round(median-lci, 3), round(uci-median, 3))]
sr2q[, c("minus", "plus") := list(round(median-lci, 3), round(uci-median, 3))]
sr3q[, c("minus", "plus") := list(round(median-lci, 3), round(uci-median, 3))]
sr4q[, c("minus", "plus") := list(round(median-lci, 3), round(uci-median, 3))]

fwrite(sr1q, file = "1rate_rates.csv")
fwrite(sr2q, file = "2rate_rates.csv")
fwrite(sr3q, file = "3rate_rates.csv")
fwrite(sr4q, file = "4rate_rates.csv")