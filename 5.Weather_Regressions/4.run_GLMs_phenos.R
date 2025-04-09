library(tidyverse)
library(magrittr)
library(data.table)
library(reshape2)
library(foreach)
library(forcats)
library(lme4)

###load data
weather_slice <- fread("/gpfs2/scratch/jcnunez/Dsu.prelim.data/Spotty_fly_1_regressions/SlicedWeatherData.mgarvin.Feb25.csv")
weather_slice %<>%
  mutate(win = factor(win, levels = unique(weather_slice$win)))

traits <- fread("/gpfs2/scratch/jcnunez/Dsu.prelim.data/Spotty_fly_1_regressions/traits_for_seq_data.txt") 
meta <- fread("/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023/2024_Joaquin_KY/swd.metadata.txt")

names(weather_slice)[1] = "sampleId_orig"

#### Run models
#i = weather_slice$win[1]
#j = weather_slice$variable[1]
#k = weather_slice$stat[1]

o =
  foreach( k = unique(weather_slice$stat)[4],
           #k = "prop. max",
           .combine = "rbind",
           .errorhandling = "remove")%do%{
  foreach( j = unique(weather_slice$variable)[1],
           #j =T2M,
           .combine = "rbind",
           .errorhandling = "remove")%do%{
foreach( i = unique(weather_slice$win),
         .combine = "rbind",
         .errorhandling = "remove")%do%{
           
  message(paste(p, i, j, k,sep = " ") ) 
           
  weather_slice %>%
    filter(win == i) %>%
    filter(variable == j) %>%
    filter(stat == k) %>%
    filter( sampleId_orig != "KY20") %>%
    left_join(traits[,-2]) %>%
    left_join(meta) -> tmp.obj
  
  foreach( p = 0:100,
           .combine = "rbind",
           .errorhandling = "remove")%do%{
             
  if(p == 0){
    #model0 <- lmer(Ctmin ~ city + (1 | fruit_type), data = tmp.obj)  
    #model1 <- lmer(Ctmin ~ city + (1 | fruit_type) + value, data = tmp.obj)  
    model0 <- lmer(Ctmin ~ (1 | city) + (1 | fruit_type), data = tmp.obj)  
    model1 <- lmer(Ctmin ~ (1 | city) + (1 | fruit_type) + value, data = tmp.obj)  
    #model0 <- lm(Ctmin ~ city +  fruit_type, data = tmp.obj)  
    #model1 <- lm(Ctmin ~ city +  fruit_type + value, data = tmp.obj)  
    
    p_lrt=anova(model1, model0, test="Chisq")[2,8]
    #p_lrt=anova(model1, model0, test="Chisq")[2,5]
  } else if(p != 0){
   
   ## Randomize traits with fruit and city
   tmp.obj[sample(dim(tmp.obj)[1]),] %>%
     dplyr::select(Ctmin, city, fruit_type) ->
     randomized_sample
   
   tmp.obj %>% 
     dplyr::select(value) ->
     true_value
   
   ran.tmp = cbind(randomized_sample, true_value)
   
   model0 <- lmer(Ctmin ~ (1 | city) + (1 | fruit_type), 
                  data = ran.tmp)  
   model1 <- lmer(Ctmin ~ (1 | city) + (1 | fruit_type) + value, 
                  data = ran.tmp)  
   p_lrt=anova(model1, model0, test="Chisq")[2,8]
   
 }
             
  data.frame(
    win= i,
    p=p,
    variable = j,
    stat= k,
    p_lrt=p_lrt
  ) 
}}}}
o

o %<>%
  mutate(win = factor(win, levels = unique(weather_slice$win)))

o %>%
  group_by(p == 0,
           win) %>%
  summarise(m.Plrt = quantile(p_lrt, 0.05)) %>%
  dcast(win~`p == 0`) %>%
  mutate(TEST = `FALSE` > `TRUE`) ->
  significant_test

o %<>% left_join(significant_test)
  
ggplot() +
  geom_violin(data = filter(o, p > 0),
              aes(x=win,
                  y=-log10(p_lrt))
              ) +
  geom_point(data = filter(o, p == 0),
              aes(x=win,
                  y=-log10(p_lrt),
                  shape = TEST ),
             size = 3, fill = "red"
  ) + scale_shape_manual(values = 23:24) + 
  ggtitle("Days above 35°C") -> test.plot
ggsave(test.plot, file = "test.plot.pdf")
