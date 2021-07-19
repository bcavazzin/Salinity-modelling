###########################################################
###### Non-linear censored model - Bayesian approach ######
###########################################################

# http://mc-stan.org/misc/warnings.html

setwd("C:/Users/BCavazzin/OneDrive - British Museum/Canada_autumn2017/Canada_autumn2017/Spatial analysis/GDGT_%")

## Load packages ####
library('mgcv')
library("brms")
library("ggplot2")
library("bayesplot")
theme_set(theme_bw())
library("tidyr")
library("readxl")
library('viridis')
library(ggmap)
library(maps)
library(mapdata)
library(rgdal) #for map
library(cowplot)
library(dplyr)

## load data ####
survey <- read_excel("ENTIRE.DATASET_2b.xlsx")

## Refit?
CORES <- 4
REFIT <- FALSE

## Load Saskatchewan map ####

if (!file.exists("./src/ref/ne_50m_admin_1_states_provinces_lakes/ne_50m_admin_1_states_provinces_lakes.dbf")){
  download.file(file.path('http://www.naturalearthdata.com/http/',
                          'www.naturalearthdata.com/download/50m/cultural',
                          'ne_50m_admin_1_states_provinces_lakes.zip'), 
                f <- tempfile())
  unzip(f, exdir = "./src/ref/ne_50m_admin_1_states_provinces_lakes")
  rm(f)
}
region <- readOGR("./src/ref/ne_50m_admin_1_states_provinces_lakes", 'ne_50m_admin_1_states_provinces_lakes', encoding='UTF-8')
# region is defined in the first part of the code (see above)
regions <- subset(region, name %in% "Saskatchewan") 
plot(regions)



##### iso 0 #####
survey <- transform(survey,
                    Censored = ifelse(iso0 < 0.1, -1, 0))

# transform 0 values in <0.1 det. limit to run the gamma model
survey$iso0[survey$iso0 < 0.1] <- 0.05

if (REFIT || !file.exists("spatial-t2-censored-gamma-iso0.rds")) {
  fit.iso0 <- brm(bf(iso0 | cens(Censored) ~  t2(SampleLat, SampleLong, k = 8) +
                       offset(log(sum.iso))), # bf = bayes factors for non-linear models
                  data = survey,
                  family = Gamma(link = "log"),
                  warmup = 1500, iter = 5500, chains = 4, cores = CORES,
                  control = list(adapt_delta = 0.999, max_treedepth = 50))
  # warmup iter (sampling iteratuins) = evaluating multiple chains to give you some 
  # confidence that you've converged to the right location.
  # e.g. in fit5 we used 4 chains, each with 5500 iterations of which 
  # the first 1500 are warmup to calibrate the sampler, 
  # leading to a total of 16000 posterior samples
  # adapt_delta = target average proposal acceptance probability during Stan’s adaptation period, 
  # and increasing it will force Stan to take smaller steps.
  # max_treedepth for NUTS (No-U-Turn Sampler) putting a cap on the depth of the 
  # trees that it evaluates during each iteration
  
  saveRDS(fit.iso0, "spatial-t2-censored-gamma-iso0.rds")
  mod.iso0.ms <- marginal_smooths(fit.iso0)
  saveRDS(mod.iso0.ms, "spatial-t2-censored-gamma-iso0-marginal-smooths.rds")
  
} else {
  mod.iso0 <- readRDS("spatial-t2-censored-gamma-iso0.rds")
  mod.iso0.ms <- readRDS("spatial-t2-censored-gamma-iso0-marginal-smooths.rds")
}
summary(mod.iso0, waic = FALSE) #widely applicable information criterion (WAIC) based on the posterior likelihood
plot(mod.iso0, stype = "raster", theme = theme_bw() + theme(legend.position = "top"))
plot(marginal_effects(fit.iso0), points = TRUE)

# make a grid with dimentions of 100x100 with points from survey dataset
newlocs <- with(survey, expand.grid(SampleLong = seq(min(SampleLong),
                                                     max(SampleLong), length = 100),
                                    SampleLat  = seq(min(SampleLat),
                                                     max(SampleLat), length = 100),
                                    sum.iso = 100)) 
# sum.iso = 100 cause we are calculating GDGT proporion over 100 ug and not as an actual % 
# of the total isoprenod. The result is the same

# extract fitted values of the model fit
pred.iso0 <- fitted(mod.iso0, newdata = newlocs, scale = "response")

# combine the grid and the fitted values
pdata <- cbind(newlocs, pred.iso0)
# are any grid points too far from the data..?
ind <- exclude.too.far(newlocs$SampleLong, newlocs$SampleLat,
                       survey$SampleLong, survey$SampleLat, dist = 0.09)
# ...set them to NA
pdata[ind, -(1:2)] <- NA

# heatmap of Estimate fitted values
iso0.est <- ggplot(pdata, aes(x = SampleLong, y = SampleLat)) +
  geom_raster(aes(fill = log(Estimate))) +
  geom_contour(aes(z = log(Estimate))) +
  geom_point(data = survey) +
  coord_equal() +
  theme(legend.position = "top") +
  scale_fill_viridis(option = 'plasma', na.value = 'transparent') +
  geom_polygon(data=regions, aes(x=long, y=lat), colour = "gray28", 
               fill = "transparent") +
  coord_fixed(ylim = c(48.6, 53.3), ratio = 1) + # ratio beween lat and long
  labs(x = "Longitude", y="Latitude")
# + theme(legend.text = element_text(size=10))

# heatmap of Error fitted values
iso0.err <- ggplot(pdata, aes(x = SampleLong, y = SampleLat)) +
  geom_raster(aes(fill = log(Est.Error))) +
  geom_contour(aes(z = log(Est.Error))) +
  geom_point(data = survey) +
  coord_equal() +
  theme(legend.position = "top") +
  scale_fill_viridis(option = 'plasma', na.value = 'transparent') +
  geom_polygon(data=regions, aes(x=long, y=lat), colour = "gray28", 
               fill = "transparent") +
  coord_fixed(ylim = c(48.6, 53.3), ratio = 1) + # ratio beween lat and long
  labs(x = "Longitude", y="Latitude") 

iso0.p <- plot_grid(iso0.est, iso0.err, labels = c("A", "B"), nrow = 2, align = "v") 
plot(iso0.p)
pdf("iso0.p.pdf")
dev.off()