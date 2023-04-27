
# run BUGS model script
# for linear regression missing data model
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


load("raw data/dat_long.RData")
load("raw data/intens_mat.RData")


#################
# jags set-up

filein <- "BUGS/model_missing_lm_long.txt"

params <- c("alpha", "beta", "mu_alpha", "sd_alpha", "mu_beta", "sd_beta",
            "missing", "log_missing", "mu.x", "p.x", "log_salinity", 
            "beta_s", "gamma")#,
#  "pred_mean_lsalinity", #"pred_lsalinity")
# "mu", "intens_pred")

n.iter <- 20000
n.burnin <- 1000
n.thin <- floor((n.iter - n.burnin)/500)
n_dat <- nrow(dat_long)


####################
# prep missing data

# select a lake
missing_lake_name <- 908

missing_lake_dat <-
  filter(dat_long, Lake_name == missing_lake_name) |> 
  mutate(log_salinity = NA,
         Lake_name = 10000) #1000 lake name same as for ex. 854

dat_total <-
  bind_rows(dat_long, missing_lake_dat) |> 
  arrange(Lake_name)

n_missing_dat <- nrow(missing_lake_dat)
n_missing_lake <- length(missing_lake_name)

bac_names <- levels(dat_total$bacteria)
lakeNames <- sort(unique(dat_total$Lake_name)) # previously called LakeIDs
n_lakes <- length(lakeNames)

sal_dat <- dat_intens[order(dat_intens$Lake_name), ]
salinity_dat <- append(log(sal_dat$Wsalinity), NA)

MAT_miss <- dat_intens$MAT[dat_intens$Lake_name == missing_lake_name]
MAT <- append(dat_intens$MAT, MAT_miss)

#
intens_mat[n_lakes, ] <- intens_mat[which(lakeNames == missing_lake_name), ]


#############
# run model

dataJags <-
  list(bac_id = as.numeric(as.factor(dat_total$bacteria)),
       n_bac = length(bac_names),
       log_salinity = as.numeric(salinity_dat),
       lake_id = as.numeric(as.factor(dat_total$Lake_name)),
       intensity = intens_mat, 
       #MAT = as.numeric(MAT),
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
R2WinBUGS::attach.bugs(res_bugs$BUGSoutput)

output <- res_bugs$BUGSoutput
x <- output$sims.matrix

save(output, file = "output_data/BUGS_output_missing.RData")


##########
# plots

library(ggplot2)

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

##

for (i in 1:dataJags$n_dat) {
  x <- dataJags$log_salinity[dataJags$lake_id[i]] 
  print(x)
}
