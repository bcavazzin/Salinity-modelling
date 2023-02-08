
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


##TODO
# lake random effect?

# mg/g sediment
dat <- readr::read_csv("raw data/Dataset_concentration.csv")

# clean data

# remove dots
# names(dat) <- gsub("\\.", "", names(dat))

bac_names <- c("IIIa","IIIa.","IIIb","IIIb.","IIIc","IIIc.","IIa",
               "IIa.","IIb","IIb.","IIc","IIc.",
               "Ia","Ib","Ic")

bac_names_log <- paste0(bac_names, "_log")

# remove 0 entries
dat[dat == 0] <- NA
dat <- dat[complete.cases(dat), ]

dat <- dat[dat$Wsalinity != 0, ]

dat$log_salinity <- log(dat$Wsalinity)

dat <- 
  dat %>% mutate(across(all_of(bac_names), list(log = log)))


# jags set-up

filein <- "BUGS/model_missing_lm.txt"
params <- c("mu", "alpha", "beta", "missing")#, "pred")

n.iter <- 20000
n.burnin <- 1000
n.thin <- floor((n.iter - n.burnin)/500)

dataJags <-
  list(bac = as.matrix(dat[, bac_names_log]) |>
         # rbind(rep(10, length(bac_names))),  # missing data bacteria
         rbind(dat[nrow(dat), bac_names_log]),
       n_bac = length(bac_names),
       envir = c(dat$log_salinity, NA),
       n_dat = nrow(dat) + 1)

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

