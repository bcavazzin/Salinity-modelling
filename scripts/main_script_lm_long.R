
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
params <- c("mu", "alpha", "beta", "missing")#, "pred")

n.iter <- 20000
n.burnin <- 1000
n.thin <- floor((n.iter - n.burnin)/500)

dataJags <-
  list(bac_id = c(as.numeric(as.factor(dat_long$bacteria)), 1),
       n_bac = length(unique(dat_long$bacteria)),
       log_salinity = c(dat_long$log_salinity, NA),
       intensity = c(dat_long$log_intensity, 13),
       n_dat = nrow(dat_long) + 1)

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
mcmc_areas(x, pars = c("alpha[1]","alpha[2]","alpha[3]"))
mcmc_areas(x, pars = c("beta[1]","beta[2]","beta[3]"))
mcmc_areas(x, pars = c("missing")) #+ xlim(0, 50)


save(output, file = "output_data/BUGS_output.RData")

