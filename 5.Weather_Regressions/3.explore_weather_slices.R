library(tidyverse)
library(magrittr)
library(data.table)
library(reshape2)
library(foreach)
library(forcats)
library(lme4)

###load data
dat <- fread("/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/ky_swd/5.Weather_Regressions/data/SlicedWeatherData.mgarvin.Feb25.csv")

dat$win %>% table
dat$stat %>% table
dat$variable %>% table

dat %>%
  filter(variable == "T2M") %>%
  filter(stat == "mean") %>%
  filter(Sample == "KY8") %>%
  ggplot(aes(
    x=windowstart,
    y=value,
    color=as.character(win)
  )) + geom_point(size = 3) #+ scale_color_gradient2(midpoint= 15)


