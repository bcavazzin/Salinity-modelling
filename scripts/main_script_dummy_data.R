
# run BUGS model script
# for linear regression missing data model
# with dummy data


library(R2jags)
library(dplyr)
library(reshape2)
library(R2WinBUGS)


#################
# jags set-up

filein <- "BUGS/model_missing_lm_long.txt"

params <- c("alpha", "beta", "mu_alpha", "sd_alpha", "mu_beta", "sd_beta",
            "missing", "log_missing", "mu.x", "p.x", "log_salinity", 
            "pred_mean_lsalinity", "pred_mean_salinity",
            "pred_intens_lsalinity", "pred_intens_salinity",
            "pred_lsalinity", "pred_salinity")

n.iter <- 20000
n.burnin <- 1000
n.thin <- floor((n.iter - n.burnin)/500)

## define inits?

##############
# dummy data

n_sample <- 5
alpha <- 10
beta <- -1

dat_total <- 
  data.frame(Lake_name = "1",
             bacteria = as.factor(LETTERS[1:n_sample])) |> 
  mutate(log_salinity = rnorm(1, 0, 1),
         log_intensity = rnorm(n_sample, alpha + beta*log_salinity, 3))

intens_mat <- t(dat_total$log_intensity)

################
# missing data

dat_total <- rbind(dat_total,
                   c("2", "A", NA_real_, 10))
intens_mat <- rbind(intens_mat,
                    c(10, NA, NA, NA, NA))

#############
# run model

dataJags <-
  list(bac_id = as.numeric(dat_total$bacteria),
       n_bac = nlevels(dat_total$bacteria),
       log_salinity = as.numeric(dat_total$log_salinity),
       lake_id = as.numeric(as.factor(dat_total$Lake_name)),
       intensity = intens_mat, 
       n_dat = nrow(dat_total),
       n_lake_miss = 1,
       n_lake_obs = 1,
       n_lake = 2)

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

save(output, file = "output_data/BUGS_output_missing_dummy.RData")


##########
# plots

library(ggplot2)
library(bayesplot)

mcmc_areas(x, regex_pars = "^alpha")
mcmc_areas(x, regex_pars = "^beta")
mcmc_areas(x, pars = c("missing"))
mcmc_areas(x, regex_pars = "intens_pred\\[1]")

mcmc_areas(x, regex_pars = "pred_mean_lsalinity")
mcmc_areas(x, regex_pars = "pred_mean_salinity")

mcmc_areas(x, regex_pars = "pred_lsalinity")
mcmc_areas(x, regex_pars = "pred_salinity")

mcmc_areas(x, regex_pars = "pred_intens_lsalinity")
mcmc_areas(x, regex_pars = "pred_intens_salinity")
