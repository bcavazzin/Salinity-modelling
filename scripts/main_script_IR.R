
# run BUGS model script
# and forest plots

library(dplyr)
library(reshape2)
library(lattice)
library(R2jags)
library(R2WinBUGS)
library(mcmc)
library(coda)
library(mcmcplots)
library(bayesplot)
library(tidyr)
library(ggplot2)

# dataset ####

# counts
# including 0s
dat <- read.csv("raw data/data_signal_intensity.csv",
                check.names = FALSE, header = TRUE)

dat$log_salinity <- log(dat$Wsalinity)
dat$log_WpH <- log(dat$WpH)
dat$log_lakeArea <- log(dat$`Lake area`)

### IR6me scatter
plot(dat$log_salinity, dat$IR6) # ratio
plot(dat$log_salinity, log(dat$IR6.Num)) # numer
plot(dat$log_salinity, log(dat$IR6.den)) # denom

## bacteria vs sal
dat.2 <- dat %>%
  pivot_longer(IIIa:Ic, names_to = "bacteria", values_to = "intensity")
dat.2$log_intensity <- log(dat.2$intensity)
  
p1 <- ggplot(dat.2, aes(x = log_salinity, y = log_intensity)) +
  geom_point() +
  geom_smooth(method = lm) +
  facet_wrap(vars(bacteria), ncol = 4) +
  labs(x = "Salinity (log)", y = "Signal intensity (log)") +
  theme_bw()
p1


dat <- dat[dat$Wsalinity != 0, ]
dat <- dat[dat$`Lake area` != 0, ]
# dat <- dat[dat$bac != 0, ]

dat$Wsalinity <- ifelse(dat$Wsalinity <= 0, 0.0001, dat$Wsalinity)

bac_names <- c("IIIa","IIIa.","IIIb","IIIb.","IIIc","IIIc.",
               "IIa","IIa.", "IIb","IIb.","IIc","IIc.",
               "Ia","Ib","Ic")
n_bac <- length(bac_names)

     # # standardize covariates
# dat <- mutate(dat, across(Lat:log_lakeArea, ~ (.x - mean(.x))/sd(.x)))


filein <- "model_missing.txt"
params <- c("r", "mu", "alpha", "beta", "log_salinity", "salinity_mis")#, "pred")

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

bac <- as.matrix(dat[, bac_names])
bac <- rbind(bac, bac[23, ])
mode(bac) <- "integer"


# # full model
# dataJags <-
#   list(bac = bac,
#        n_bac = n_bac,
#        log_salinity = c(dat$log_salinity, NA),
#        WT = c(dat$WT, NA),
#        MAT = c(dat$MAT, NA),
#        log_lakeArea = c(dat$log_lakeArea, NA),
#        log_WpH = c(dat$log_WpH, NA),
#        n_dat = nrow(bac))

# salinity only
dataJags <-
  list(bac = bac,
       n_bac = n_bac,
       log_salinity = c(dat$log_salinity, NA),
       n_dat = nrow(bac))

filein <- "BUGS/model_missing.txt"

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
mcmc_areas(x, pars = c("beta[1,1]","beta[1,2]","beta[1,3]"))
mcmc_areas(x, pars = c("log_salinity[88]"))
mcmc_areas(x, pars = c("log_salinity[88]"), transformations = "exp")
mcmc_areas(x, pars = c("salinity_mis")) +
  # xlim(0.5,100.5)
  coord_cartesian(xlim = c(0, 2000))

mcmc_areas(x, regex_pars = "alpha")
mcmc_areas(x, regex_pars = "beta")


save(output, file = "BUGS_output.RData")

