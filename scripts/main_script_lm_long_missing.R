
# run BUGS model script
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

load("raw data/dat_long.RData")

# jags set-up

filein <- "BUGS/model_missing_lm_long.txt"
params <- c("alpha", "beta", "mu_alpha", "sd_alpha", "mu_beta",
            "sd_beta", "missing", "log_missing", "mu.x", "p.x")#,
            #  "pred_mean_lsalinity", #"pred_lsalinity")
            # "mu", "intens_pred")

#dat_long <- subset(dat_long, bacteria == c("Ib", "IIa", "IIb" ,"IIIa"))

n.iter <- 20000
n.burnin <- 1000
n.thin <- floor((n.iter - n.burnin)/500)
n_dat <- nrow(dat_long)

missing_lake_dat <-
  dat_long |> 
  filter(LakeID == 926)

##TODO: could just append to bottom of dat_long?

dataJags <-
  list(bac_id = as.numeric(as.factor(c(dat_long$bacteria, missing_lake_dat$bacteria))),
       n_bac = length(unique(dat_long$bacteria)),
       log_salinity = c(unique(dat_long$log_salinity), rep(NA, nrow(missing_lake_dat))),
       lake = c(as.numeric(as.factor(dat_long$LakeID)), 
                rep(max(as.numeric(as.factor(dat_long$LakeID))) + 1, nrow(missing_lake_dat))),
       intensity = c(dat_long$log_intensity, missing_lake_dat$log_intensity),
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
mcmc_areas(x, regex_pars = "intens_pred\\[1") #+ xlim(0, 50)
mcmc_areas(x, pars = c("intens_pred[1071]","intens_pred[1072]","intens_pred[1073]","intens_pred[1076]"))

save(output, file = "output_data/BUGS_output_missing.RData")

