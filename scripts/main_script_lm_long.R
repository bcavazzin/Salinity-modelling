
# run BUGS model script
# version _without_ missingness model
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
dat <- read.csv(here::here("raw data", "data_signal_intensity.csv"),
                check.names = FALSE, header = TRUE)
dat <- dat[dat$Wsalinity != 0, ]

# jags set-up

filein <- "BUGS/model_lm_long.txt"
params <- c("alpha", "beta", "mu_alpha", "sd_alpha","mu_beta", "sd_beta")

n.iter <- 50000
n.burnin <- 1000
n.thin <- floor((n.iter - n.burnin)/500)
n_dat <- nrow(dat_long)

# Lake ID and bacteria gorup matrix
intens_mat <- dat %>% select("LakeID", "IIIa":"Ic")
intens_mat[, 2:16] <- log(intens_mat[2:16],)
intens_mat[intens_mat == -Inf] <- NA

# order in same way as dat_long
intens_matrix <- intens_mat[,-1]
rownames(intens_matrix) <- intens_mat[,1]
intens_matrix <- intens_matrix[levels(as.factor(dat_long$LakeID)), ]
intens_matrix <- intens_matrix[, levels(dat_long$bacteria)] 

# dataJags <-
#   list(bac_id = as.numeric(as.factor(dat_long$bacteria)),
#        n_bac = length(unique(dat_long$bacteria)),
#        log_salinity = dat_long$log_salinity,
#        intensity = dat_long$log_intensity,
#        n_dat = nrow(dat_long))

dataJags <-
  list(bac_id = as.numeric(as.factor(dat_long$bacteria)),
       n_bac = length(unique(dat_long$bacteria)),
       log_salinity = unique(dat_long$log_salinity),
       lake = as.numeric(as.factor(dat_long$LakeID)),
       intensity = intens_matrix, #dat_long$log_intensity,
       n_lake = length(unique(dat_long$log_salinity)),
       n_dat = n_dat)

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
mcmc_areas(x, regex_pars = "alpha") #pars = c("alpha[1]","alpha[2]","alpha[3]"))
mcmc_areas(x, regex_pars = "beta") #pars = c("beta[1]","beta[2]","beta[3]"))

save(output, file = "output_data/BUGS_output.RData")

