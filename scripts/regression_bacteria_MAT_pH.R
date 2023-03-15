
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

# dataset ####

# counts
# including 0s
dat <- read.csv("raw data/data_signal_intensity.csv",
                check.names = FALSE, header = TRUE)

dat <- dat[dat$Wsalinity != 0, ]

dat$log_salinity <- log(dat$Wsalinity)
dat$log_WpH <- scale(log(dat$WpH), scale = FALSE)
dat$log_lakeArea <- log(dat$`Lake area`)
dat$MAT <- scale(dat$MAT, scale = FALSE)
dat$Soi_pH <- scale(dat$Soi_pH, scale = FALSE) 

dat.2 <- dat

# dat.2 <- dat[,-c(8,10,14)] #remove IIIc', IIa' and IIc'

dat.2 <- dat.2 %>%
  pivot_longer(IIIa:Ic, names_to = "bacteria", values_to = "intensity")
dat.2 <- dat.2[dat.2$intensity != 0, ] #remove  with signal intensity = 0

dat.2$log_intensity <- log(dat.2$intensity)
dat.2$bacteria <- as.factor(dat.2$bacteria)
#dat.3 <- dat.2[, c("LakeName", "bacteria", "log_salinity", "log_intensity")] #for Section 3
dat.2 <- dat.2[, c("bacteria", "log_salinity", "log_intensity", "MAT", "log_WpH", "Soi_pH", "IR6")] #for section 1 and 2

# regression + plots ####

p.soilpH <- ggplot(dat.2, aes(x = Soi_pH, y = log_intensity)) +
  geom_point() +
  geom_smooth(method = lm) +
  facet_wrap(vars(bacteria), ncol = 4) +
  labs(x = "Soil pH", y = "Signal intensity (log)") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           p.accuracy = 0.001, r.accuracy = 0.01) + 
  theme_bw()
p.soilpH

p.waterpH <- ggplot(dat.2, aes(x = log_WpH, y = log_intensity)) +
  geom_point() +
  geom_smooth(method = lm) +
  facet_wrap(vars(bacteria), ncol = 4) +
  labs(x = "Water pH", y = "Signal intensity (log)") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           p.accuracy = 0.001, r.accuracy = 0.01) + 
  theme_bw()
p.waterpH

p.MAT <- ggplot(dat.2, aes(x = MAT, y = log_intensity)) +
  geom_point() +
  geom_smooth(method = lm) +
  facet_wrap(vars(bacteria), ncol = 4) +
  labs(x = "MAT", y = "Signal intensity (log)") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           p.accuracy = 0.001, r.accuracy = 0.01) + 
  theme_bw()
p.MAT

p.IR6ME <- ggplot(dat.2, aes(x = log_salinity, y = IR6)) +
  geom_point() +
  geom_smooth(method = lm) +
  labs(x = "Salinity (log)", y = "IR6ME") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           p.accuracy = 0.001, r.accuracy = 0.01) + 
  theme_bw()
p.IR6ME


p.MBT <- ggplot(dat.2, aes(x = MAT, y = log_intensity)) +
  geom_point() +
  geom_smooth(method = lm) +
  facet_wrap(vars(bacteria), ncol = 4) +
  labs(x = "MAT", y = "Signal intensity (log)") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           p.accuracy = 0.001, r.accuracy = 0.01) + 
  theme_bw()
p.MBT

