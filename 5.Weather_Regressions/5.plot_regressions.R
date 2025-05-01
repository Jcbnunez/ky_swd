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

load("/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/ky_swd/5.Weather_Regressions/data/CTmin_Weather_Comparison_Unfiltered_MOG.Rdata")
load("/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/ky_swd/5.Weather_Regressions/data/GLM_models_output_MOG.Rdata")

setwd("/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/ky_swd/5.Weather_Regressions/plot")
#####

o %>%
  filter(variable == "T2M") %>%
  filter(stat == "prop. max") ->
  prop.max

ggplot() +
  geom_violin(data=filter(prop.max, p > 0),
              aes(x=win,
                  y=-log10(p_lrt))) +
  geom_point(data=filter(prop.max, p == 0),
              aes(x=win,
                  y=-log10(p_lrt)), 
             shape = 23, 
             size = 6, fill = "grey") +
  theme_bw() -> p1

####
o %>%
  filter(variable == "PRECTOTCORR") %>%
  filter(stat == "mean") ->
  mean

ggplot() +
  geom_violin(data=filter(mean, p > 0),
              aes(x=win,
                  y=-log10(p_lrt))) +
  geom_point(data=filter(mean, p == 0),
             aes(x=win,
                 y=-log10(p_lrt)), 
             shape = 23, 
             size = 6, fill = "grey") +
  theme_bw() -> p2

####
o %>%
  filter(variable == "PRECTOTCORR") %>%
  filter(stat == "maximum") ->
  maxim

ggplot() +
  geom_violin(data=filter(maxim, p > 0),
              aes(x=win,
                  y=-log10(p_lrt))) +
  geom_point(data=filter(maxim, p == 0),
             aes(x=win,
                 y=-log10(p_lrt)), 
             shape = 23, 
             size = 6, fill = "grey") +
  theme_bw() -> p3

ggsave(p1, file = "p1.pdf",
       h = 3, w= 4)
ggsave(p2, file = "p2.pdf",
       h = 3, w= 4)

####
pmax_ct<-get(load("/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/ky_swd/5.Weather_Regressions/data/CTmin_Weather_Comparison_PropMaxFilter_MOG.Rdata"))
write.table(pmax_ct, 
            file = "pmax_ct.txt", append = FALSE, 
            quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,
            col.names = TRUE, 
            qmethod = c("escape", "double"),
            fileEncoding = "")

pmax_ct %>%
  filter(!is.na(site)) %>%
  ggplot(aes(
    x=`prop. max`,
    y=ctmin,
    color=as.character(Time.point),
  )) + geom_point(size = 3.1) +
  geom_smooth(method = "lm", se = F, color = "black") +
  theme_bw()+
  scale_color_manual(values = c("steelblue","navyblue")) +
  theme(legend.position = "none") -> propplot

premen_ct<-get(load("/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/ky_swd/5.Weather_Regressions/data/CTmin_Weather_Comparison_PRECTOTCORR_meanfilter_MOG.Rdata"))
write.table(premen_ct, 
            file = "premen_ct.txt", append = FALSE, 
            quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,
            col.names = TRUE, 
            qmethod = c("escape", "double"),
            fileEncoding = "")

premen_ct %>%
  filter(!is.na(site)) %>%
  ggplot(aes(
    x=mean,
    y=ctmin,
    color=as.character(Time.point),
  )) + geom_point(size = 3.1) +
  geom_smooth(method = "lm", se = F, color = "black") +
  theme_bw()+
  scale_color_manual(values = c("steelblue","navyblue"))+
  theme(legend.position = "none") -> preceplot

ggsave(propplot, file = "propplot.pdf",
       h = 3, w= 3)
ggsave(preceplot, file = "preceplot.pdf",
       h = 3, w= 3)

####
model0 <- lmer(ctmin ~ (1 | site) + (1 | `Fruit type`), 
               data = pmax_ct)  
model1 <- lmer(ctmin ~ (1 | site) + (1 | `Fruit type`) + `prop. max`, 
               data = pmax_ct)  

p_lrt=anova(model1, model0, test="Chisq")[2,8]

##
cor.test(pmax_ct$ctmin, pmax_ct$`prop. max`)

cor.test(premen_ct$ctmin, premen_ct$mean)
