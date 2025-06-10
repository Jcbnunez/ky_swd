#########################
### Loading Libraries ###
#########################

library(tidyverse)
library(magrittr)
library(data.table)
library(reshape2)
library(poolfstat)
library(FactoMineR)
library(reshape2)
library(foreach)
library(forcats)
library(lme4)
library(lmerTest)

load("/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/ky_swd/5.Weather_Regressions/data/SlicedWeatherData.JCBN_MOG.Jun5_2025.Rdata")
load("/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/ky_swd/5.Weather_Regressions/data/GLM_models_output_JCBN_Jun6.Rdata")
traits <- fread("/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/ky_swd/5.Weather_Regressions/data/Means_CT_min_FINAL.csv")
traits %<>%
  mutate(Sample = paste(Generation, 
                               Time.point, 
                               ifelse(site == "Berea", "Berea", "Lexington"), 
                               year, sep = "_"))


setwd("/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/ky_swd/5.Weather_Regressions/plot")
#####
### Raw phenotype analyses
lmer(ctmin ~ (1 | site) + (1 | `Fruit type`) + (1 | year) +
       `Juliean Date`, 
     data = filter(traits, Generation == "Parental")) -> Pmodel
summary(Pmodel)
anova(Pmodel)


lmer(ctmin ~ (1 | site) + (1 | `Fruit type`)  + (1 | year) + 
       `Juliean Date`, 
     data = filter(traits, Generation == "f4")) -> Fmodel
summary(Fmodel)
anova(Fmodel)


#####

o %>%
  filter(variable == "T2M") %>%
  filter(stat == "prop. max") %>%
  filter(win == "c(0, 7)") ->
  prop.max

ggplot() +
  geom_boxplot(data=filter(prop.max, p > 0),
              aes(x=factor(Generation, levels = c("Parental", "f4")),
                  y=-log10(p_lrt))) +
  geom_point(data=filter(prop.max, p == 0),
              aes(x=factor(Generation, levels = c("Parental", "f4")),,
                  y=-log10(p_lrt)), 
             shape = 23, 
             size = 6, fill = "grey") +
  theme_bw()  -> p1


ggsave(p1, file = "p1.pdf",
       h = 3, w= 4)


####
results %>%
  filter(win == "c(0, 7)") %>%
  filter(variable == "T2M") %>%
  filter(stat == "prop. max") %>%
  left_join(traits) %>%
  mutate(site_year = paste(site, year, sep ="_"))->
  traits_envvar

traits_envvar %>%
  ggplot(aes(
    x=value,
    y=ctmin,
  )) + geom_point(size = 3.5, aes(fill=month(as.Date(Date, format = "%m/%d/%Y")), 
                                  shape=as.character(site)))+
  scale_shape_manual(values  = 21:22) + scale_fill_gradient2(midpoint = 9,
                                                             low="firebrick",
                                                             high="powderblue") +
  facet_grid(factor(Generation, levels = c("Parental","F4"))~.) + geom_smooth(method = "lm", se = F, 
                                         color = "black",
                                         linetype = "dashed",
                                         ) +
  theme_bw() -> propplot
propplot

ggsave(propplot, file = "propplot.pdf",
       h = 6, w= 7.5)


lmer(ctmin ~ (1 | site_year) + (1 | `Fruit type`)  +
       value, 
     data = filter(traits_envvar, Generation == "Parental")) -> Pmodel2
summary(Pmodel2)
anova(Pmodel2)

lmer(ctmin ~ (1 | site_year) + (1 | `Fruit type`)  +
       value, 
     data = filter(traits_envvar, Generation == "f4")) -> fmodel2
summary(fmodel2)
anova(fmodel2)
