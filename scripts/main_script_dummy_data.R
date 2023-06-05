
# run BUGS model script
# for linear regression missing data model
# with dummy data


library(R2jags)
library(dplyr)
library(reshape2)
library(R2WinBUGS)
library(tidyr)
library(tibble)


#################
# jags set-up

filein <- "BUGS/model_missing_dummy.txt"

params <- c("alpha", "beta", "mu_alpha", "sd_alpha", "mu_beta", "sd_beta",
            "missing", "log_missing", "m.missing",
            "mu.x", "p.x", "log_salinity", 
            "pred_mean_lsalinity", "pred_mean_salinity",
            "pred_intens_lsalinity", "pred_intens_salinity",
            "pred_lsalinity", "pred_salinity")

n.iter <- 20000
n.burnin <- 1000
n.thin <- floor((n.iter - n.burnin)/500)

##TODO: what values?
jags.inits <- function() {
  list("mu.x" = rnorm(1),
       "sd.x" = runif(1),
       "sd_alpha" = rnorm(1),
       "sd_beta" = runif(1),
       "mu_alpha" = runif(1),
       "mu_beta" = runif(1),
       "sd.tau" = rnorm(1))
}

##############
# dummy data 20 lakes, 5 bacteria per lake

n_lakes <- 20
n_entires <- 100
n_bacteria <- 5
alpha0 <- 10
beta0 <- -1

LakeName <- rep(1:20, each = 5)
bacteria <- rep(as.factor(LETTERS[1:n_bacteria]), 20)

dummy <- tibble(LakeName, bacteria) %>%
  mutate(log_salinity = rep(rnorm(n_lakes, 0, 1), each = 5))|>
  group_by(bacteria) |>
  mutate(alpha_bac = rnorm(1, alpha0, 1), 
         beta_bac = rnorm(1, beta0, 1)) %>% 
  ungroup() %>% 
  mutate(log_intensity = rnorm(n_entires, alpha_bac + beta_bac*log_salinity, 1))

# Long dataset with last lake missing salinity
dat_total <- select(dummy, -alpha_bac, -beta_bac) |>
    mutate(log_salinity = ifelse(LakeName == 20, NA, log_salinity))

dat_total$log_salinity = as.numeric(dat_total$log_salinity)
dat_total$log_intensity = as.numeric(dat_total$log_intensity)

intens_mat <- dat_total %>%
  pivot_wider(names_from = bacteria, values_from = log_intensity) |>
  select(-log_salinity, -LakeName)
  # missing salinity for lake 20 
  # mutate(log_salinity = ifelse(LakeName == 20, NA, LakeName))

salinity_dat <- dat_total %>%
  pivot_wider(names_from = bacteria, values_from = log_intensity) |>
  select(log_salinity)

missing_lake <- dat_total[dat_total$LakeName == "20", ]
observed_dat <- dat_total[dat_total$LakeName < 20, ] 

missing_lake_name <- 20
n_missing_lake <- length(missing_lake_name)
lakeNames <- sort(unique(dat_total$LakeName))

#############
# run model

dataJags <-
  list(bac_id = as.numeric(dat_total$bacteria),
       n_bac = nlevels(dat_total$bacteria),
       log_salinity = salinity_dat$log_salinity,
       lake_id = as.numeric(as.factor(dat_total$LakeName)),
       intensity = intens_mat, 
       #n_obs = nrow(observed_dat),
       #n_miss = nrow(missing_lake),
       n_dat = nrow(dat_total),
       n_lake_miss = n_missing_lake,
       n_lake = length(lakeNames), 
       n_lake_obs = n_lakes - n_missing_lake
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

##########
# scatter plots

x.dat <- as.data.frame(x)

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

x.dat2 <- x.dat %>%
  select("alpha[1]":"beta[5]")

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
                          "1" ~ "A",
                          "2" ~ "B",
                          "3" ~ "C",
                          "4" ~ "D",
                          "5"~ "E")
  )

set.seed(20220330)

xx_samples <- xx_long %>%
  filter(draw %in% sample(1:40000, size = 5000)) %>%
  tidyr::pivot_wider(names_from = param, values_from = Value)

xx_samples

#plot_data <- dat_total[-16,]

ggplot(xx_samples) +
  facet_wrap("bacteria") +
  geom_point(data = dat_total, aes(x = log_salinity, y = log_intensity), inherit.aes = FALSE) +
  geom_abline(aes(intercept = 10, slope = -1), col="red") +
  geom_abline(aes(intercept = alpha, slope = -beta), size = 0.75, #-beta because the BUGS code is positive relationship 
              color = "#3366FF", 
              alpha = .1) +
  xlim(-10,10)


