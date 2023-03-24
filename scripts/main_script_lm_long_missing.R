
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
dat_intens <- read.csv(here::here("raw data", "data_signal_intensity.csv"),
                check.names = FALSE, header = TRUE)
dat_intens <- dat_intens[dat_intens$Wsalinity != 0, ]

# jags set-up ####

filein <- "BUGS/model_missing_lm_long2.txt"
params <- c("alpha", "beta", "mu_alpha", "sd_alpha", "mu_beta",
            "sd_beta", "missing", "log_missing", "mu.x", "p.x", "log_salinity")#,
            #  "pred_mean_lsalinity", #"pred_lsalinity")
            # "mu", "intens_pred")

n.iter <- 20000
n.burnin <- 1000
n.thin <- floor((n.iter - n.burnin)/500)
n_dat <- nrow(dat_long)

missing_lake_id <- 861

missing_lake_dat <-
  filter(dat_long, LakeID == missing_lake_id) |> 
  mutate(log_salinity = NA,
         LakeID = 10000)

dat_total <-
  bind_rows(dat_long, missing_lake_dat) |> 
  arrange(LakeID)

n_missing_dat <- nrow(missing_lake_dat)
n_missing_lake <- length(missing_lake_id)
  
bac_names <- levels(dat_total$bacteria)
lakeIDs <- sort(unique(dat_total$LakeID))
n_lakes <- length(lakeIDs)

sal_dat <- dat_intens[order(dat_intens$LakeID),]
salinity_dat <- log(sal_dat$Wsalinity)
salinity_dat <- append(salinity_dat,'NA')

intens_mat <- dat_intens %>%
  select("LakeID", "IIIa":"Ic") %>% 
  mutate(across(everything(), ~replace(., . == 0, NA)),
         across(
           .cols = "IIIa":"Ic",
           .fn = log)) |> 
  arrange(LakeID) |> 
  select(all_of(bac_names))

intens_mat[n_lakes, ] <- intens_mat[which(lakeIDs == missing_lake_id), ]

dataJags <-
  list(bac_id = as.numeric(as.factor(dat_total$bacteria)),
       n_bac = length(bac_names),
       log_salinity = salinity_dat, #dat_total$log_salinity
       lake = as.numeric(as.factor(dat_total$LakeID)),
       intensity = intens_mat,
       n_dat = n_dat + n_missing_dat,
       n_obs = n_dat,
       n_miss = n_missing_dat,
       n_lake_miss = n_missing_lake,
       n_lake = n_lakes,
       n_lake_obs = n_lakes-n_missing_lake)

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

