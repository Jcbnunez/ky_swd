#### Haldane Analysis
library(tidyverse)
library(magrittr)
library(foreach)
library(doMC)
registerDoMC(cores = 4)
library(gmodels)
library(data.table)
library(stringr)
library(reshape2)
library(lubridate)
library(ggdist)
library(ggbeeswarm)

### load data
pheno <- fread("Mean_CT_SD_CT_haldane.txt")

# ---------------------------
# 1. Parse dates
# ---------------------------

pheno$Date <- as.Date(pheno$Date, format = "%m/%d/%y")
pheno2 <- pheno %>%
  mutate(Date = as.Date(Date))

# ---------------------------
# 2. Create breeding season clock
#    May 1 - Dec 31 only
#    16 generations per year
# ---------------------------

first_year <- min(year(pheno2$Date))

pheno2 <- pheno2 %>%
  mutate(
    yr = year(Date),
    
    season_start = make_date(yr, 5, 1),
    season_end   = make_date(yr, 12, 31),
    
    breeding_day_within_year = case_when(
      Date < season_start ~ 0,
      Date > season_end   ~ 245,
      TRUE ~ as.numeric(Date - season_start)
    ),
    
    breeding_day =
      (yr - first_year) * 245 +
      breeding_day_within_year
  )

# ---------------------------
# 3. Calculate Haldanes
# ---------------------------

haldanes <- pheno2 %>%
  arrange(Locality, Date) %>%
  group_by(Locality) %>%
  mutate(
    
    prev_mean = lag(Mean_CTmin),
    prev_sd   = lag(SD_CTmin),
    prev_n    = lag(N),
    prev_day  = lag(breeding_day),
    prev_year = lag(yr),
    
    season_transition = case_when(
      is.na(prev_year) ~ NA_character_,
      yr == prev_year ~ "Within growing season",
      TRUE ~ "Across growing seasons"
    ),
    
    generations =
      (breeding_day - prev_day) * (16/245),
    
    sd_pooled = sqrt(
      ((N - 1) * SD_CTmin^2 +
         (prev_n - 1) * prev_sd^2) /
        (N + prev_n - 2)
    ),
    
    haldane =
      (Mean_CTmin - prev_mean) /
      (sd_pooled * generations)
    
  ) %>%
  ungroup()

###
haldanes %>%
  filter(!is.na(season_transition)) %>%
  group_by(season_transition) %>%
  summarise(mh = mean(abs(haldane)))

haldanes %>%
  filter(!is.na(season_transition)) %>%
ggplot(
       aes(season_transition,
           abs(haldane),
           fill = season_transition)) +
  
  geom_boxplot(
    width = 0.3,
    alpha = 0.5,
    outlier.shape = NA
  ) +
  geom_quasirandom(
    width = 0.15,
    size = 3,
    alpha = 0.8
  ) +
  theme_classic() +
  theme(legend.position = "none")-> haldane_plot

ggsave(haldane_plot, file = "haldane_plot.pdf",
       h= 4, w = 4)

