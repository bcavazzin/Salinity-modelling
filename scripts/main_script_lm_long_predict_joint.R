
# run BUGS model script for two-step model
# and joint prior on alpha, beta
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


# load data set
load("raw data/dat_long.RData")

dat_intens <- read.csv(here::here("raw data", "data_signal_intensity.csv"),
                       check.names = FALSE, header = TRUE)

names(dat_long)[names(dat_long) == "LakeID"] <- "lake_label"

bac_keep <- c('Ib', 'Ic', 'IIa', 'IIb', 'IIb.', 'IIc','IIIa', 'IIIa.',
              'IIIb','IIIb.','IIIc')

dat_long <- dat_long[dat_long$bacteria %in% bac_keep, ]
dat_long$bacteria <- as.factor(as.character(dat_long$bacteria))

names(dat_intens)[names(dat_intens) == "LakeID"] <- "lake_label"
dat_intens <- dat_intens[, c("lake_label", bac_keep, "MAT", "Wsalinity")]
dat_intens <- dat_intens[dat_intens$Wsalinity != 0, ]

# jags set-up

filein <- "BUGS/model_missing_lm_long.txt"

params <- c("alpha", "beta", "mu_alpha", "sd_alpha", "mu_beta",
            "sd_beta","pred_mean_lsalinity", "pred_lsalinity",
            "mu", "intens_pred", "mu.x", "p.x")

n.iter <- 20000
n.burnin <- 1000
n.thin <- floor((n.iter - n.burnin)/500)
n_dat <- nrow(dat_long)

dat_total <- dat_long |> arrange(lake_label)
  
lake_labels <- sort(unique(dat_total$lake_label))

n_lakes <- length(lake_labels)

intens_mat <- dat_intens %>%
  select("lake_label", bac_keep) %>%
  mutate(across(everything(), ~replace(., . == 0, NA)),
         across(.cols = bac_keep, .fn = log)) |> 
  arrange(lake_label) |> 
  select(all_of(bac_keep))

dataJags <-
  list(bac_id = as.numeric(as.factor(dat_total$bacteria)),
       log_salinity = log(dat_intens$Wsalinity), 
       lake_id = as.numeric(as.factor(dat_total$lake_label)),
       intensity = intens_mat,
       MAT = as.numeric(dat_intens$MAT),
       n_bac = length(bac_keep),
       n_dat = n_dat,
       n_lake = n_lakes
       )

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

#options(max.print=100000)
print(res_bugs)

#mcmcplot(res_bugs)
# plots <- traceplot(res_bugs)

R2WinBUGS::attach.bugs(res_bugs$BUGSoutput)
output <- res_bugs$BUGSoutput

# BUGS output into data frame
x <- output$sims.matrix
x.dat <- as.data.frame(x)

# create long dataset of model output alphas and betas for each bacteria
alpha.poster <- x.dat %>%
  as.data.frame()%>%
  select(starts_with("alpha")) %>%
  pivot_longer(everything(),
               names_to = "bacteria.alpha",
               values_to = "Alpha")
beta.poster <- x.dat %>%
  as.data.frame()%>%
  select(starts_with("beta")) %>%
  pivot_longer(everything(),
               names_to = "bacteria.beta",
               values_to = "Beta")
# merge the two
dat.post <- cbind(alpha.poster, beta.poster)

dat.post <- dat.post %>% 
  separate(bacteria.alpha, c('param', 'bacteria')) %>%
  mutate(
    bacteria = case_match(bacteria,
                          "1" ~ "Ib",
                          "2" ~ "Ic",
                          "3" ~ "IIa",
                          "4" ~ "IIb",
                          "5"~ "IIb.",
                          "6"~ "IIc",
                          "7"~ "IIIa",
                          "8"~ "IIIa.",
                          "9"~ "IIIb",
                          "10" ~ "IIIb.",
                          "11" ~ "IIIc"))
  

# 2D density plot of the alphas and betas 
ggplot(dat.post) + 
  aes(x = Alpha, y = Beta) + 
  stat_density_2d(
    aes(fill = after_stat(level)), 
    geom = "polygon", 
    # normalized density so all colors appear in each plot
    contour_var = "ndensity"
  ) +
  facet_wrap("bacteria") + 
  xlab("Intercept estimate") + 
  ylab("Slope estimate") +
  guides(fill = "none") 

#select only alphas and beta per bacteria (remove mu)
x.dat2 <- x.dat %>%
  select("alpha[1]":"beta[11]")

xx_long <- x.dat2 %>%
  select_all() %>%
  rowid_to_column("draw") %>%
  tidyr::pivot_longer(
    cols = c(-draw),
    # when we make new columns with pivot_ functions, the
    # they get quotes
    names_to = "Parameter", 
    values_to = "Value"
  ) %>%
  separate(Parameter, c('param', 'bacteria')) %>%
  mutate(
    bacteria = case_match(bacteria,
      "1" ~ "Ib",
      "2" ~ "Ic",
      "3" ~ "IIa",
      "4" ~ "IIb",
      "5"~ "IIb.",
      "6"~ "IIc",
      "7"~ "IIIa",
      "8"~ "IIIa.",
      "9"~ "IIIb",
      "10" ~ "IIIb.",
      "11" ~ "IIIc"
    )
  )

set.seed(20220330)

xx_samples <- xx_long %>%
  filter(draw %in% sample(1:1000, size = 50)) %>%
  tidyr::pivot_wider(names_from = param, values_from = Value)

xx_samples

ggplot(xx_samples) +
  facet_wrap("bacteria") +
  geom_point(data = dat_long, aes(x = log_salinity, y = log_intensity), inherit.aes = FALSE) +
  geom_abline(aes(intercept = alpha, slope = -beta), size = 0.75, #-beta because the BUGS code is positive relationship 
              color = "#3366FF", 
              alpha = .1)

library(ggplot2)
mcmc_areas(x, regex_pars = "^alpha") #pars = c("alpha[1]","alpha[2]","alpha[3]"))
mcmc_areas(x, regex_pars = "^beta") #pars = c("beta[1]","beta[2]","beta[3]"))
mcmc_areas(x, pars = c("missing")) #+ xlim(0, 50)
mcmc_areas(x, regex_pars = "intens_pred\\[1]") #+ xlim(0, 50)
mcmc_areas(x, pars = c("intens_pred[1071]","intens_pred[1072]","intens_pred[1073]","intens_pred[1076]"))

save(output, file = "output_data/BUGS_output_missing.RData")


########################
# from GBs book
# BUGS:
#   model {
#   for (i in 1:N) {
#     y[i] ~ dnorm(mu[i],tau)
#     mu[i] <- alpha + beta*X[i]
#   }
#   # Blocking of coefficients
#   coef[1:2] ~ dmnorm(m[1:2],P[1:2,1:2])
#   alpha <- coef[1]
#   beta <- coef[2]
#   lsigma ~ dunif(-k,k)
#   sigma <- exp(lsigma)
#   tau <- pow(sigma,-2)
# }
# 
# R: m <- c(0,0)
# prec <- (1/h^2)*diag(2)