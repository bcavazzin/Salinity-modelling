
# run BUGS model script for two-step model
# and joint prior on alpha, beta
# and forest plots

library(R2jags)
library(dplyr)
library(reshape2)
library(mcmc)
library(coda)
library(lattice)
library(R2WinBUGS)
library(mcmcplots)
library(bayesplot)
library(tidyverse)


# load data set
load("raw data/dat_long.RData")

dat_intens <- read.csv(here::here("raw data", "data_signal_intensity.csv"),
                       check.names = FALSE, header = TRUE)

names(dat_long)[names(dat_long) == "LakeID"] <- "lake_label"

bac_keep <- c('Ib', 'IIa', 'IIb', 'IIb.', 'IIIa', 'IIIc')

dat_long <- dat_long[dat_long$bacteria %in% bac_keep, ]
dat_long$bacteria <- as.factor(as.character(dat_long$bacteria))

names(dat_intens)[names(dat_intens) == "LakeID"] <- "lake_label"
dat_intens <- dat_intens[, c("lake_label", bac_keep, "MAT", "Wsalinity")]
dat_intens <- dat_intens[dat_intens$Wsalinity != 0, ]


# jags set-up

filein <- "BUGS/model_predict_lm_joint.txt"

params <- c("alpha", "beta", "mu_alpha", "sd_alpha", "mu_beta",
            "sd_beta","pred_mean_lsalinity", "pred_lsalinity",
            "mu", "intens_pred")

n.iter <- 20000
n.burnin <- 1000
n.thin <- floor((n.iter - n.burnin)/500)
n_dat <- nrow(dat_long)

dat_total <- dat_long |> arrange(lake_label)
  
lake_labels <- sort(unique(dat_total$lake_label))

n_lakes <- length(lake_labels)

intens_mat <- dat_intens %>%
  select("lake_label", bac_keep) %>%
  mutate(across(everything(), ~replace(., . == 0, NA)),
         across(.cols = bac_keep, .fn = log)) |> 
  arrange(lake_label) |> 
  select(all_of(bac_keep))

dataJags <-
  list(bac_id = as.numeric(as.factor(dat_total$bacteria)),
       log_salinity = log(dat_intens$Wsalinity), 
       lake_id = as.numeric(as.factor(dat_total$lake_label)),
       intensity = intens_mat,
       n_bac = length(bac_keep),
       n_dat = n_dat,
       n_lake = n_lakes)

res_bugs <-
  jags(data = dataJags,
       inits = NULL,
       parameters.to.save = params,
       model.file = filein,
       n.chains = 1,
       n.iter,
       n.burnin,
       n.thin,
       DIC = TRUE)

print(res_bugs)
mcmcplot(res_bugs)
# plots <- traceplot(res_bugs)

R2WinBUGS::attach.bugs(res_bugs$BUGSoutput)

output <- res_bugs$BUGSoutput

x <- output$sims.matrix

library(ggplot2)

mcmc_areas(x, regex_pars = "^alpha") #pars = c("alpha[1]","alpha[2]","alpha[3]"))
mcmc_areas(x, regex_pars = "^beta") #pars = c("beta[1]","beta[2]","beta[3]"))
mcmc_areas(x, pars = c("missing")) #+ xlim(0, 50)
mcmc_areas(x, regex_pars = "intens_pred\\[1]") #+ xlim(0, 50)
mcmc_areas(x, pars = c("intens_pred[1071]","intens_pred[1072]","intens_pred[1073]","intens_pred[1076]"))

save(output, file = "output_data/BUGS_output_missing.RData")

