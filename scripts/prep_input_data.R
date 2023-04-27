# prep_input_data


#######################
# dataset of counts and salinity
# long format
# including 0s

dat <- read.csv(here::here("raw data", "data_signal_intensity.csv"),
                check.names = FALSE, header = TRUE)

dat <- dat[dat$Wsalinity != 0, ]

dat$log_salinity <- log(dat$Wsalinity)
dat$log_WpH <- scale(log(dat$WpH), scale = FALSE)
dat$log_lakeArea <- log(dat$`Lake area`)
dat$MAT <- scale(dat$MAT, scale = FALSE)

# melt data
dat_long <- dat %>%
  pivot_longer(IIIa:Ic,
               names_to = "bacteria",
               values_to = "intensity") |> 
  mutate(log_intensity = log(intensity),
         bacteria = as.factor(bacteria)) |> 
  filter(intensity != 0) |>              # remove with signal intensity = 0
  select("LakeID", "bacteria", "log_salinity", "log_intensity")

names(dat_long)[names(dat_long) == "LakeID"] <- "Lake_name"

save(dat_long, file = "raw data/dat_long.RData")


#################
# intensity
# wide format

dat_intens <- read.csv(here::here("raw data", "data_signal_intensity.csv"),
                       check.names = FALSE, header = TRUE)

dat_intens <- dat_intens[, names(dat_intens)!="LakeName"]
dat_intens <- dat_intens[dat_intens$Wsalinity != 0, ]
names(dat_intens)[names(dat_intens) == "LakeID"] <- "Lake_name"

bac_names <- c("Ia",          "Ib",          "Ic",
               "IIa", "IIa.", "IIb", "IIb.", "IIc", "IIc.",
               "IIIa","IIIa.","IIIb","IIIb.","IIIc","IIIc.")

intens_mat <-
  dat_intens %>%
  select("Lake_name", "IIIa":"Ic") %>% 
  mutate(across(everything(), ~replace(., . == 0, NA)),
         across(
           .cols = "IIIa":"Ic",
           .fn = log)) |> 
  arrange(Lake_name) |> 
  select(all_of(bac_names))

save(intens_mat, file = "raw data/intens_mat.RData")

