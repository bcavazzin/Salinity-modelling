
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

# dataset ####

dat <- read.csv("Dataset_Bianca.csv", check.names = FALSE, header = TRUE)

dat$log_salinity <- log(dat$Wsalinity)
dat$log_WpH <- log(dat$WpH)
dat$log_lakeArea <- log(dat$`Lake area`)
dat <- dat[dat$Wsalinity != 0, ]
dat <- dat[dat$`Lake area` != 0, ]
dat$Wsalinity <- ifelse(dat$Wsalinity <= 0, 0.0001, dat$Wsalinity)

bac_names <- c("IIIa","IIIa.","IIIb","IIIb.","IIIc","IIIc.","IIa","IIa.",
               "IIb","IIb.","IIc","IIc.","Ia","Ib","Ic")
n_bac <- length(bac_names)

# # standardize covariates
# dat <- mutate(dat, across(Lat:log_lakeArea, ~ (.x - mean(.x))/sd(.x)))


filein <- "model_missing.txt"
params <- c("mu", "alpha", "beta", "log_salinity")#, "pred")

n.iter <- 20000
n.burnin <- 1000
n.thin <- floor((n.iter - n.burnin)/500)


# # no missing covariates
# dataJags <-
#   list(bac = matrix(c(dat$Ic, dat$IIc, dat$IIIc),
#                     ncol = 3, byrow = FALSE),
#        n_bac = 3,
#        envir = c(dat$Wsalinity, NA),
#        n_dat = nrow(dat))

# # only one bacteria
# dataJags <-
#   list(bac = c(dat$Ic, 3),
#        envir = c(dat$Wsalinity, NA),
#        n_dat = nrow(dat) + 1)
# filein <- "model_one_bac.txt"
# params <- c("mu", "alpha", "beta")

# # all bacteria
##TODO: supply separate log_salinity for each bacteria?
dataJags <-
  list(bac = as.matrix(dat[, bac_names])|>
         rbind(c(16,19,1,2,0,15,17,0,6,10,0,8,0,6,1)),
       n_bac = n_bac,
       log_salinity = c(dat$log_salinity, NA),
       WT = c(dat$WT, NA),
       MAT = c(dat$MAT, NA),
       log_lakeArea = c(dat$log_lakeArea, NA),
       log_WpH = c(dat$log_WpH, NA),
       n_dat = nrow(dat) + 1)
filein <- "model_missing.txt"

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

# print(res_bugs)
# mcmcplot(res_bugs)
# plots <- traceplot(res_bugs)

R2WinBUGS::attach.bugs(res_bugs$BUGSoutput)

output <- res_bugs$BUGSoutput

x <- output$sims.matrix
mcmc_areas(x, pars = c("beta[1]","beta[2]","beta[3]"))
mcmc_areas(x, pars = c("envir[92]"))


save(output, file = "BUGS_output.RData")

