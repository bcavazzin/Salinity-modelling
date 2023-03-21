
# regressions
# bacteria measurement against salinity 

# libraries ####
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
library(tidyverse)
library(ggpubr)


# dataset counts
# including 0s
dat <- read.csv(here::here("raw data", "data_signal_intensity.csv"),
                check.names = FALSE, header = TRUE)

dat <- dat[dat$Wsalinity != 0, ]

dat$log_salinity <- log(dat$Wsalinity)
dat$log_WpH <- scale(log(dat$WpH), scale = FALSE)
dat$log_lakeArea <- log(dat$`Lake area`)
dat$MAT <- scale(dat$MAT, scale = FALSE)

# dat_long <- dat[, -c(8,10,14)]  # remove IIIc', IIa' and IIc'

# melt data
dat_long <- dat %>%
  pivot_longer(IIIa:Ic, names_to = "bacteria", values_to = "intensity") |> 
  mutate(log_intensity = log(intensity),
         bacteria = as.factor(bacteria)) |> 
  filter(intensity != 0) |>              # remove with signal intensity = 0
  select("LakeID", "bacteria", "log_salinity", "log_intensity")

#dat.3 <- dat_long[, c("LakeName", "bacteria", "log_salinity", "log_intensity")] #for Section 3

save(dat_long, file = "raw data/dat_long.RData")


############################
# basic plots

p1 <- ggplot(dat_long, aes(x = log_salinity, y = log_intensity)) +
  geom_point() +
  geom_smooth(method = lm) +
  facet_wrap(vars(bacteria), ncol = 4) +
  labs(x = "Salinity (log)", y = "Signal intensity (log)") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    p.accuracy = 0.001, r.accuracy = 0.01) + 
  theme_bw()
p1


# pairs(dat)
# GGally::ggpairs(dat, columns = 20:34)


## Frequentist no pool, pool, partial pool #################

names_bacteria <- levels(dat_long$bacteria) 

m_no_pooling <- list()

for (i in names_bacteria) {
  m_no_pooling[[i]] <- lm(log_intensity ~ log_salinity, filter(dat_long, bacteria == i)) 
}

df_no_pooling <- tibble(
  Model = "No pooling",
  bacteria = levels(dat_long$bacteria),
  Intercept = sapply(m_no_pooling, coef)[1, ], 
  Slope_Salinity = sapply(m_no_pooling, coef)[2, ]
)

#write.csv(df_no_pooling, file = "df_no_pool.csv")
#read.csv("df_no_pool.csv", row.names = 1)

m_pooled <- lm(log_intensity ~ log_salinity, dat_long) 

# Repeat the intercept and slope terms for each participant
df_pooled <- tibble(
  Model = "Complete pooling",
  bacteria = levels(dat_long$bacteria),
  Intercept = coef(m_pooled)[1], 
  Slope_Salinity = coef(m_pooled)[2]
)

head(df_pooled)

# join raw data
df_models.1 <- bind_rows(df_pooled, df_no_pooling) 

ggplot(df_models.1) +
  facet_wrap("bacteria") +
  geom_point(data = dat_long, aes(x = log_salinity, y = log_intensity), inherit.aes = FALSE) +
  geom_abline(aes(intercept = Intercept, slope = Slope_Salinity, color = Model), size = 0.75) +
  theme_bw()


mod.1 <- lmer(log_intensity ~ 1 + log_salinity + (1 + log_salinity | bacteria), dat_long)
mod.1

df_partial_pooling <- coef(mod.1)[["bacteria"]] %>% 
  rownames_to_column("bacteria") %>% 
  as_tibble() %>% 
  rename(Intercept = `(Intercept)`, Slope_Salinity = log_salinity) %>% 
  add_column(Model = "Partial pooling")

head(df_partial_pooling)

df_models.2 <- bind_rows(df_models.1, df_partial_pooling) 

ggplot(df_models.2) +
  facet_wrap("bacteria") +
  geom_point(data = dat_long, aes(x = log_salinity, y = log_intensity), inherit.aes = FALSE) +
  geom_abline(aes(intercept = Intercept, slope = Slope_Salinity, color = Model), size = 0.75) +
  theme_bw()

## Bayesian no pool, pool, partial pool####################


df_pulled <- bind_rows(df_no_pooling, df_partial_pooling)

b <- stan_glmer(
  log_intensity ~ log_salinity + (1 + log_salinity | bacteria),
  family = gaussian(),
  data = dat_long,
  # prior = normal(0, 2, autoscale = TRUE),
  # prior_intercept = normal(0, 5, autoscale = TRUE),
  # prior_covariance = decov(regularization = 2),
  # prior_aux = cauchy(0, 1, autoscale = TRUE),
  seed = 20211116
)

b
prior_summary(b) #http://mc-stan.org/rstanarm/articles/priors.html
xx <- extract(b)$stan_summary

df_posterior <- b %>% 
  as.data.frame() %>% 
  as_tibble()

## Average intercept and slope
ggplot(df_posterior) + 
  aes(x = `(Intercept)`, y = `log_salinity`) + 
  # Calculate the density
  stat_density_2d(aes(fill = stat(level)), geom = "polygon") +
  ggtitle("Where's the average intercept and slope?") + 
  xlab("Estimate for average intercept") + 
  ylab("Estimate for average slope") +
  # Use the same coordinate limits as last plot
  coord_cartesian(
    xlim = range(df_pulled$Intercept), 
    ylim = range(df_pulled$Slope_Salinity),
    expand = TRUE
  ) + 
  guides(fill = "none")


df_effects <- df_posterior %>%
  mutate(
    # Find all the columns with the pattern "b[(Intercept". Add the column
    # `(Intercept)` to each of those columns.
    across(
      .cols = matches("b\\[\\(Intercept"), 
      .fns = ~ . + `(Intercept)`
    ),
    # Again for slope
    across(
      .cols = matches("b\\[log_salinity"), 
      .fns = ~ . + log_salinity
    )
  )

# Convert to a long format
df_long_effects <- df_effects %>%
  select(matches("b\\[")) %>%
  rowid_to_column("draw") %>%
  tidyr::pivot_longer(
    cols = c(-draw),
    # when we make new columns with pivot_ functions, the
    # they get quotes
    names_to = "Parameter", 
    values_to = "Value"
  )

# df_long_effects <-
#   df_long_effects %>% 
#   mutate(
#     Effect = Parameter %>% 
#       stringr::str_detect("Intercept") %>%
#       ifelse(., "Intercept", "Slope_Salinity"),
#     Subject = Parameter %>%
#       stringr::str_extract("\\d\\d\\d")
#   ) %>% 
#   select(draw, bacteria, Effect, Value)

df_long_effects2 <-
  df_long_effects %>%
tidyr::separate_wider_delim(Parameter,
                            delim = " ", names = c("param", "bacteria")) %>% 
  mutate(bacteria = gsub("bacteria:", "", bacteria),
         bacteria = gsub("]", "", bacteria),
         param = gsub("b\\[", "", param),
         param = gsub("\\(", "", param),
         param = gsub(")", "", param))

df_long_effects2

set.seed(20220330)

df_samples <- df_long_effects2 %>%
  filter(draw %in% sample(1:4000, size = 50)) %>%
  tidyr::pivot_wider(names_from = param, values_from = Value)

df_samples

ggplot(df_samples) +
  facet_wrap("bacteria") +
  geom_point(data = dat_long, aes(x = log_salinity, y = log_intensity), inherit.aes = FALSE) +
  geom_abline(aes(intercept = Intercept, slope = log_salinity), size = 0.75,
              color = "#3366FF", 
              alpha = .1)

b %>%
  spread_draws(`(Intercept)`, b[,group]) %>%
  median_qi(condition_mean = `(Intercept)` + b, .width = c(.95, .8, .5))


# forest plot intercept
p1 <- b %>%
  spread_draws(`(Intercept)`, b[,group]) %>%
  median_qi(condition_mean = `(Intercept)` + b, .width = c(.95, .66)) %>%
  add_row(group = "Pulled Interc", condition_mean = 11.2, .lower = 11.00, 
          .upper = 11.33, .width = 0.66, 
          .point = "median", .interval = "qi") %>%
  add_row(group = "Pulled Interc", condition_mean = 11.2, .lower = 10.67, 
          .upper = 11.70, .width = 0.95, 
          .point = "median", .interval = "qi") %>%
  ggplot(aes(y = group, x = condition_mean, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  labs(x ="Mean Intercept", y = "") +
  scale_y_discrete(labels = c("Ia", "Ib", "Ic", "IIa", "IIa'", "IIb",
                              "IIb'", "IIc","IIc'", "IIIa", "IIIa'",
                              "IIIb", "IIIb'", "IIIc","IIIc'","Model av.")) + # "IIa'", "IIc'", "IIIc'"
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

# forest plot salinity slope
p2 <- b %>%
  spread_draws(`log_salinity`, b[,group]) %>%
  median_qi(condition_mean = `log_salinity` + b, .width = c(.95, .66)) %>%
  add_row(group = "Pulled slope", condition_mean = -0.453, .lower = -0.511, 
          .upper = -0.396, .width = 0.66, 
          .point = "median", .interval = "qi") %>%
  add_row(group = "Pulled slope", condition_mean = -0.453, .lower = -0.634, 
          .upper = -0.283, .width = 0.95, 
          .point = "median", .interval = "qi") %>%
  ggplot(aes(y = group, x = condition_mean, xmin = .lower, xmax = .upper)) +
  geom_pointinterval()  +
  labs(x ="Mean Slope", y = "") +
  scale_y_discrete(labels = c("Ia", "Ib", "Ic", "IIa", "IIa'", "IIb",
                              "IIb'", "IIc","IIc'", "IIIa", "IIIa'",
                              "IIIb", "IIIb'", "IIIc","IIIc'","Model av.")) + # "IIa'", "IIc'", "IIIc'"
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

#library(ggpubr)
ggarrange(p1, p2,
          ncol = 2, nrow = 1)

## Section 3 Individual lakes regression salinity vs counts ####

p.bar <- ggplot(dat.3, aes(x = bacteria , y = log_intensity)) +
  geom_bar(stat="identity") + 
  facet_wrap(vars(LakeName), ncol = 8) +
  #labs(x = "Salinity (log)", y = "Signal intensity (log)") +
  theme_bw()
p.bar

p.all <- ggplot(dat.3, aes(x = log_salinity, y = log_intensity, color=bacteria)) +
  geom_point() +
  #geom_smooth(method = lm) +
  labs(x = "Salinity (log)", y = "Signal intensity (log)") +
  theme_bw()
p.all

p.box <- ggplot(dat.3, aes(x = log_salinity, y = log_intensity, color=LakeName)) +
  geom_boxplot() +
  #geom_smooth(method = lm) +
  labs(x = "Salinity (log)", y = "Signal intensity (log)") +
  theme(legend.position = "none")
p.box

## diagnostics #####

help('pareto-k-diagnostic')
loo::pareto_k_table(b)

loo::pareto_k_table(loo(b))
loo::plot(loo(b))
plot(loo(b))

## stan_glm OLD ##############

lm_null <- stan_glm(log_intensity ~ 1,
                    family = "gaussian",
                    data = dat_long)

yrep_lm_null <- posterior_predict(lm_null)

lm11 <- stan_glmer(log_intensity ~ log_salinity + (log_salinity|bacteria),
                   family = "gaussian",
                   data = dat_long)

yrep_lm11 <- posterior_predict(lm11)

# plots

color_scheme_set("brightblue")
ppc_dens_overlay(dat_long$log_intensity, yrep_lm11[1:50, ])
ppc_dens_overlay(dat_long$log_intensity, yrep_lm_null[1:50, ])

ppc_dens_overlay(dat_long$log_intensity, exp(yrep_lm11[1:100, ]))
ppc_dens_overlay(dat_long$log_intensity, exp(yrep_lm_null[1:100, ])) 

plot(dat_long$log_salinity, dat_long$log_intensity)
abline(coef(lm11)[1], coef(lm11)[2])

# frequentist

freq_lm <- lmer(log_intensity ~ log_salinity + (log_salinity|bacteria),
                data = dat_long)

freq_lm_pred <- predict(freq_lm)
