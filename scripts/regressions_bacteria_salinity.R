
# regressions
# bacteria measurement against salinity 

library(dplyr)
library(broom)
library(rstan)
library(rstanarm)
library(brms)
library(bayesplot)
library(ggplot2)
library(lme4) 
library(brms)
library(tidyr)
library(tidybayes)
library(modelr)
library(ggdist)
library(magrittr)
library(tibble)

# dataset ####

# counts
# including 0s
dat <- read.csv("raw data/data_signal_intensity.csv",
                check.names = FALSE, header = TRUE)

dat <- dat[dat$Wsalinity != 0, ]

dat$log_salinity <- log(dat$Wsalinity)
dat$log_WpH <- log(dat$WpH)
dat$log_lakeArea <- log(dat$`Lake area`)

dat <- dat[,-c(8,10,14)]

## bacteria vs sal
dat.2 <- dat %>%
  pivot_longer(IIIa:Ic, names_to = "bacteria", values_to = "intensity")
dat.2 <- dat.2[dat.2$intensity != 0, ]

dat.2$log_intensity <- log(dat.2$intensity)
dat.2$bacteria <- as.factor(dat.2$bacteria)

############################
# basic plots

# pairs(dat)
# GGally::ggpairs(dat, columns = 20:34)


################
# normal fits

lm_null <- stan_glm(log_intensity ~ 1,
                    family = "gaussian",
                    data = dat.2)

yrep_lm_null <- posterior_predict(lm_null)

lm11 <- stan_glmer(log_intensity ~ log_salinity + (log_salinity|bacteria),
                 family = "gaussian",
                 data = dat.2)

yrep_lm11 <- posterior_predict(lm11)

lm12 <- stan_glmer(log_intensity ~ 1 + (1 |bacteria),
                   family = "gaussian",
                   data = dat.2)

yrep_lm12 <- posterior_predict(lm12)

##########
# plots

color_scheme_set("brightblue")
ppc_dens_overlay(dat.2$log_intensity, yrep_lm11[1:50, ])
ppc_dens_overlay(dat.2$log_intensity, yrep_lm_null[1:50, ])
ppc_dens_overlay(dat.2$log_intensity, yrep_lm12[1:50, ])

ppc_dens_overlay(dat.2$log_intensity, exp(yrep_lm11[1:100, ]))
ppc_dens_overlay(dat.2$log_intensity, exp(yrep_lm_null[1:100, ])) 

plot(dat.2$log_salinity, dat.2$log_intensity)
abline(coef(lm11)[1], coef(lm11)[2])

# frequentist

freq_lm <- lmer(log_intensity ~ log_salinity + (log_salinity|bacteria),
              data = dat.2)

freq_lm_pred <- predict(freq_lm)

###################
# 
m_no_pooling <- lmList(log_intensity ~ log_salinity | bacteria, dat.2) 

df_no_pooling <- tibble(
  Model = "No pooling",
  Subject = unique(dat.2$bacteria),
  Intercept = coef(m_no_pooling)[1], 
  Slope_Salinity = coef(m_no_pooling)[2]
)

m_pooled <- lm(log_intensity ~ log_salinity, dat.2) 

# Repeat the intercept and slope terms for each participant
df_pooled <- tibble(
  Model = "Complete pooling",
  Subject = unique(dat.2$bacteria),
  Intercept = coef(m_pooled)[1], 
  Slope_Salinity = coef(m_pooled)[2]
)

head(df_pooled)

# join raw data
df_models <- bind_rows(df_pooled, df_no_pooling) %>% 
  left_join(dat.2, by = "bacteria")

p_model_comparison <- ggplot(df_models) + 
  aes(x = Days, y = Reaction) + 
  # Set the color mapping in this layer so the points don't get a color
  geom_abline(
    aes(intercept = Intercept, slope = Slope_Days, color = Model),
    size = .75
  ) + 
  geom_point() +
  facet_wrap("Subject") +
  labs(x = xlab, y = ylab) + 
  scale_x_continuous(breaks = 0:4 * 2) + 
  # Fix the color palette 
  scale_color_brewer(palette = "Dark2") + 
  theme(legend.position = "top", legend.justification = "left")

p_model_comparison

mod.1 <- lmer(log_intensity ~ 1 + log_salinity + (1 + log_salinity | bacteria), dat.2)
mod.1

df_partial_pooling <- coef(mod.1)[["bacteria"]] %>% 
  rownames_to_column("bacteria") %>% 
  as_tibble() %>% 
  rename(Intercept = `(Intercept)`, Slope_Salinity = log_salinity) %>% 
  add_column(Model = "Partial pooling")

head(df_partial_pooling)

p_mod.1 <- ggplot(df_partial_pooling) + 
  aes(x = log_salinity, y = log_intensity) + 
  # Set the color mapping in this layer so the points don't get a color
  geom_abline(
    aes(intercept = Intercept, slope = Slope_Salinity),
    size = .75
  ) + 
  geom_point() +
  facet_wrap("bacteria") +
  #labs(x = Salinity (Log), y = Bacteria (log)) + 
  #scale_x_continuous() + 
  # Fix the color palette 
  scale_color_brewer(palette = "Dark2") + 
  theme(legend.position = "top", legend.justification = "left")

p_mod.1 

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

ppc_dens_overlay(dat$bac, yrep_fit2[1:100, ])
ppc_dens_overlay(dat$bac, yrep_fit3[1:100, ])
ppc_dens_overlay(dat$bac, yrep_nb0[1:100, ])

ppc_dens_overlay(dat$log_bac, log(yrep_fit2[1:100, ]))
ppc_dens_overlay(dat$log_bac, log(yrep_fit3[1:100, ]))

ppc_dens_overlay(dat$bac, epred_fit2[1:100, ]) + xlim(0, 800)
ppc_dens_overlay(dat$bac, epred_fit3[1:100, ]) + xlim(0, 800)
ppc_dens_overlay(dat$bac, epred_nb0[1:100, ]) + xlim(0, 800)

ppc_dens_overlay(dat$log_bac, log(epred_fit2[1:100, ]))

