
# run BUGS model script
# for linear regression missing data model
# and forest plots


library(R2jags)
library(dplyr)
library(reshape2)
library(R2WinBUGS)


load("raw data/dat_long.RData")
load("raw data/intens_mat.RData")


#################
# jags set-up

filein <- "BUGS/model_missing_lm_long.txt"

params <- c("alpha", "beta", "mu_alpha", "sd_alpha", "mu_beta", "sd_beta",
            "missing", "log_missing", "mu.x", "p.x", "log_salinity", 
            "beta_s",
            "pred_mean_lsalinity", "pred_mean_salinity",
            "pred_intens_lsalinity", "pred_intens_salinity",
            "pred_lsalinity", "pred_salinity")

n.iter <- 20000
n.burnin <- 1000
n.thin <- floor((n.iter - n.burnin)/500)
n_dat <- nrow(dat_long)

## define inits?

####################
# prep missing data

# select a lake
missing_lake_name <- 908

missing_lake_dat <-
  filter(dat_long, Lake_name == missing_lake_name) |> 
  mutate(log_salinity = NA,
         Lake_name = 10000) #1000 lake name same as for ex. 854

n_missing_dat <- nrow(missing_lake_dat)
n_missing_lake <- length(missing_lake_name)

dat_total <-
  bind_rows(dat_long, missing_lake_dat) |> 
  arrange(Lake_name)

lakeNames <- sort(unique(dat_total$Lake_name)) # previously called LakeIDs
n_lakes <- length(lakeNames)

log_salinity <- dat_total |>
  group_by(Lake_name) |>
  arrange(Lake_name) |> 
  summarise(log_salinity = first(log_salinity)) |> 
  pull()

intens_mat[n_lakes, ] <- intens_mat[which(lakeNames == missing_lake_name), ]


#############
# run model

dataJags <-
  list(bac_id = as.numeric(as.factor(dat_total$bacteria)),
       n_bac = nlevels(dat_total$bacteria),
       log_salinity = log_salinity,
       lake_id = as.numeric(as.factor(dat_total$Lake_name)),
       intensity = intens_mat, 
       n_obs = n_dat,
       n_miss = n_missing_dat,
       n_dat = n_dat + n_missing_dat,
       n_lake_miss = n_missing_lake,
       n_lake = n_lakes,
       n_lake_obs = n_lakes - n_missing_lake)

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
R2WinBUGS::attach.bugs(res_bugs$BUGSoutput)

output <- res_bugs$BUGSoutput
x <- output$sims.matrix

save(output, file = "output_data/BUGS_output_missing.RData")


##########
# plots

library(ggplot2)
library(bayesplot)

mcmcplot(res_bugs)
# plots <- traceplot(res_bugs)

mcmc_areas(x, regex_pars = "^alpha") #pars = c("alpha[1]","alpha[2]","alpha[3]"))
mcmc_areas(x, regex_pars = "^beta") #pars = c("beta[1]","beta[2]","beta[3]"))
mcmc_areas(x, pars = c("missing")) #+ xlim(0, 50)
mcmc_areas(x, regex_pars = "intens_pred\\[1]") #+ xlim(0, 50)
mcmc_areas(x, pars = c("intens_pred[1071]",
                       "intens_pred[1072]",
                       "intens_pred[1073]",
                       "intens_pred[1076]"))

mcmc_areas(x, regex_pars = "pred_mean_lsalinity\\[107")
mcmc_areas(x, regex_pars = "pred_mean_salinity\\[107") + xlim(0, 10)

mcmc_areas(x, regex_pars = "pred_lsalinity\\[107")
mcmc_areas(x, regex_pars = "pred_salinity\\[107")

mcmc_areas(x, regex_pars = "pred_intens_lsalinity\\[107")
mcmc_areas(x, regex_pars = "pred_intens_salinity\\[107") + xlim(0,10)
