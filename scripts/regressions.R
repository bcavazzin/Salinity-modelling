
# regressions

library(dplyr)
library("rstan")
library("rstanarm")
library("bayesplot")
library("ggplot2")
library("broom")


dat_total <- readr::read_csv("Dataset_concentration.csv")

# select bacteria
var_names <- c("Lat","Long","Sample Depth","Lake area","MAT","Soi_pH","WT","Surface WT","Wsalinity","WpH")   
bac_name <- "Ia"
# bac_name <- "IIIa"

dat <- dat_total[, c(var_names, bac_name)]

dat <- rename(dat, "bac" = bac_name) |> 
  mutate(bac = as.integer(round(bac, 0)))

# dat$Ia <- floor(dat$Ia)  # for poisson model

dat$log_bac <- log(dat$bac)
dat$log_salinity <- log(dat$Wsalinity)
dat$log_WpH <- log(dat$WpH)
dat$log_lakeArea <- log(dat$`Lake area`)
dat <- dat[dat$Wsalinity != 0, ]
dat <- dat[dat$bac != 0, ]
dat <- dat[dat$`Lake area` != 0, ]

# remove outliers
# dat <- dat[!(dat$log_salinity > 2.5 & dat$log_bac > 3.8), ]
# dat <- dat[dat$bac < 400, ]

# change point variable
dat$log_salinity_star <- (dat$log_salinity > 1) * (dat$log_salinity - 1)  

# standardize covariates
dat <- 
  mutate(dat, across(Lat:WpH, ~ (.x - mean(.x))/sd(.x)))


############################
# basic plots

# pairs(dat)
# GGally::ggpairs(dat, columns = 20:34)


################
# normal fits

y <- dat$bac

lm_null <- stan_glm(log_bac ~ 1,
                    family = "gaussian",
                    data = dat)
lm11 <- stan_glm(log_bac ~ log_salinity,
                 family = "gaussian",
                 data = dat)
lm22 <- stan_glm(log_bac ~ log_salinity + WT + MAT + log_lakeArea + log_WpH,
                 family = "gaussian",
                 data = dat)

# change point model
lm_cp <- stan_glm(log_bac ~ log_salinity + log_salinity_star,
                  family = "gaussian",
                  data = dat)

# lm1 <- stan_lm(log_bac ~ log_salinity,
#                prior = NULL,
#                data = dat)
# lm2 <- stan_lm(log_bac ~ log_salinity + WT + MAT + log_lakeArea,
#                prior = NULL,
#                data = dat)

yrep_lm_null <- posterior_predict(lm_null)
yrep_lm11 <- posterior_predict(lm11)
yrep_lm22 <- posterior_predict(lm22)
yrep_lm_cp <- posterior_predict(lm_cp)
# yrep_lm <- posterior_predict(lm1)
# yrep_lm2 <- posterior_predict(lm2)

epred_lm22 <- posterior_epred(lm22)


##########
# plots

color_scheme_set("brightblue")
# ppc_dens_overlay(y, yrep[1:5, ])
# ppc_dens_overlay(y, yrep_null[1:5, ])
ppc_dens_overlay(dat$log_bac, yrep_lm11[1:5, ])
ppc_dens_overlay(dat$log_bac, yrep_lm_null[1:5, ])

ppc_dens_overlay(dat$bac, exp(yrep_lm_cp[1:100, ])) + xlim(0, 500)

# ppc_dens_overlay(dat$Ia, exp(yrep_lm2[1:100, ])) + xlim(0, 500)
ppc_dens_overlay(dat$bac, exp(yrep_lm11[1:100, ])) + xlim(0, 500) + ylim(0,0.04)
ppc_dens_overlay(dat$bac, exp(yrep_lm22[1:100, ])) + xlim(0, 500) + ylim(0,0.04)
ppc_dens_overlay(dat$bac, exp(yrep_lm_null[1:100, ])) + xlim(0, 500) + ylim(0,0.04)

ppc_dens_overlay(dat$bac, exp(epred_lm22[1:100, ])) + xlim(0, 500) + ylim(0,0.04)

ppc_dens_overlay(dat$log_bac, yrep_lm11[1:100, ]) + xlim(-5, 10) + ylim(0,0.3)
ppc_dens_overlay(dat$log_bac, yrep_lm22[1:100, ]) + xlim(-5, 10) + ylim(0,0.3)
ppc_dens_overlay(dat$log_bac, yrep_lm_null[1:100, ]) + xlim(-5, 10) + ylim(0,0.3)
ppc_dens_overlay(dat$log_bac, yrep_lm_cp[1:100, ]) + xlim(-5, 10) + ylim(0,0.3)


plot(dat$log_salinity, dat$log_bac)
abline(coef(lm11)[1], coef(lm11)[2])
segments(x0 = -2, y0 = coef(lm_cp)[1] - 2*coef(lm_cp)[2],
         x1 =  1, y1 = coef(lm_cp)[1] + coef(lm_cp)[2], col = "blue")
segments(x0 =  1, y0 = coef(lm_cp)[1] + coef(lm_cp)[2],
         x1 =  4, y1 = coef(lm_cp)[1] + 4*(coef(lm_cp)[2] + coef(lm_cp)[3]), col = "blue")

# frequentist

freq_lm <- lm(log_bac ~ log_salinity + log_salinity_star,
              data = dat)

freq_lm_pred <- predict(freq_lm)


plot(dat$log_salinity, dat$log_bac)
segments(x0 = -2, y0 = coef(freq_lm)[1] - 2*coef(freq_lm)[2],
         x1 =  1, y1 = coef(freq_lm)[1] + coef(freq_lm)[2], col = "blue")
segments(x0 =  1, y0 = coef(freq_lm)[1] + coef(freq_lm)[2],
         x1 =  4, y1 = coef(freq_lm)[1] + 4*(coef(freq_lm)[2] + coef(freq_lm)[3]), col = "blue")


###################
# poisson fits

library(brms)

fit1 <- brm(bac ~ log_salinity,
            family = "poisson",
            data = dat)

fit1

yrep_fit1 <- posterior_predict(fit1)

ppc_dens_overlay(dat$bac, yrep_fit1[1:50, ])


######################
# negative binomial

fit_nb0 <- brm(bac ~ 1,
               family = negbinomial(),
               data = dat)
# full model
fit2 <- brm(bac ~ log_salinity + WT + MAT + log_lakeArea + log_WpH,
            family = negbinomial(),
            data = dat)

fit3 <- brm(bac ~ log_salinity,
            family = zero_inflated_negbinomial(),
            data = dat)

fit2
fit3

yrep_nb0 <- posterior_predict(fit_nb0)
yrep_fit2 <- posterior_predict(fit2)
yrep_fit3 <- posterior_predict(fit3)

epred_nb0 <- posterior_epred(fit_nb0)
epred_fit2 <- posterior_epred(fit2)
epred_fit3 <- posterior_epred(fit3)

ppc_dens_overlay(dat$bac, yrep_fit2[1:100, ]) + xlim(0, 800)
ppc_dens_overlay(dat$bac, yrep_fit3[1:100, ])
ppc_dens_overlay(dat$bac, yrep_nb0[1:100, ])

ppc_dens_overlay(dat$log_bac, log(yrep_fit2[1:100, ]))
ppc_dens_overlay(dat$log_bac, log(yrep_fit3[1:100, ]))

ppc_dens_overlay(dat$bac, epred_fit2[1:100, ]) + xlim(0, 800)
ppc_dens_overlay(dat$bac, epred_fit3[1:100, ]) + xlim(0, 800)
ppc_dens_overlay(dat$bac, epred_nb0[1:100, ]) + xlim(0, 800)

ppc_dens_overlay(dat$log_bac, log(epred_fit2[1:100, ]))

