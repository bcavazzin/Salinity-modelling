
# run BUGS model script
# and forest plots

# library ####
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

# load data set ####
load("raw data/dat_long.RData")
dat <- read.csv(here::here("raw data", "data_signal_intensity.csv"),
                check.names = FALSE, header = TRUE)
dat <- dat[dat$Wsalinity != 0, ]


# jags set-up ####

filein <- "BUGS/model_missing_lm_long.txt"
params <- c("alpha", "beta", "mu_alpha", "sd_alpha", "mu_beta",
            "sd_beta", "missing", "log_missing", "mu.x", "p.x", "log_salinity")#,
            #  "pred_mean_lsalinity", #"pred_lsalinity")
            # "mu", "intens_pred")

n.iter <- 20000
n.burnin <- 1000
n.thin <- floor((n.iter - n.burnin)/500)
n_dat <- nrow(dat_long)

missing_lake_id <- 935
missing_lake_dat <-
  dat_long |> 
  filter(LakeID == missing_lake_id)

# Lake ID and bacteria gorup matrix
intens_mat <- dat %>% select("LakeID", "IIIa":"Ic")
intens_mat[, 2:16] <- log(intens_mat[2:16],)
intens_mat[intens_mat == -Inf] <- NA

# order in same way as dat_long
intens_matrix <- intens_mat[,-1]
rownames(intens_matrix) <- intens_mat[,1]

intens_matrix <- intens_matrix[levels(as.factor(dat_long$LakeID)), ]
intens_matrix <- intens_matrix[, levels(dat_long$bacteria)] 

intens_matrix[nrow(intens_matrix) + 1, ] <- intens_matrix[as.character(missing_lake_id), ]

##TODO: could just append to bottom of dat_long?

dataJags <-
  list(bac_id = as.numeric(as.factor(c(dat_long$bacteria, missing_lake_dat$bacteria))),
       n_bac = length(unique(dat_long$bacteria)),
       log_salinity = c(unique(dat_long$log_salinity), rep(NA, nrow(missing_lake_dat))),
       lake = c(as.numeric(as.factor(dat_long$LakeID)), 
                rep(max(as.numeric(as.factor(dat_long$LakeID))) + 1, nrow(missing_lake_dat))),
       intensity = intens_matrix,  ### c(dat_long$log_intensity, missing_lake_dat$log_intensity)
       n_dat = n_dat + nrow(missing_lake_dat),
       n_obs = n_dat,
       n_miss = nrow(missing_lake_dat),
       n_lake_miss = length(unique(missing_lake_dat$LakeID)),
       n_lake_obs = length(unique(dat_long$LakeID)))

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

