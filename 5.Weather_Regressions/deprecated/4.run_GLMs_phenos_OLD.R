library(tidyverse)
library(magrittr)
library(data.table)
library(reshape2)
library(foreach)
library(forcats)
library(lme4)

###Filter


###load data
weather_slice <- get(load("/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023/Ellie_MS/Weather_data/SlicedWeatherData.JCBN_MOG.Jun5_2025.Rdata"))
weather_slice %<>%
  mutate(win = factor(win, levels = unique(weather_slice$win)))

traits <- fread("/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023/Ellie_MS/Phenotype_data/Means_CT_min_FINAL.csv") 
traits %<>%
  mutate(sampleId_orig = paste(Generation, 
                               Time.point, 
                               ifelse(site == "Berea", "Berea", "Lexington"), 
                               year, sep = "_"))

meta <- fread("/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023/2024_Joaquin_KY/swd.metadata.txt")

names(weather_slice)[1] = "sampleId_orig"

#### Run models
#i = weather_slice$win[1]; j = weather_slice$variable[1]; k = weather_slice$stat[1]; g = "f4"

o =
  foreach( g = c("f4","Parental"),
           .combine = "rbind",
           .errorhandling = "remove")%do%{
  foreach( k = unique(weather_slice$stat),
           .combine = "rbind",
           .errorhandling = "remove")%do%{
  foreach( j = unique(weather_slice$variable),
           .combine = "rbind",
           .errorhandling = "remove")%do%{
foreach( i = unique(weather_slice$win),
         .combine = "rbind",
         .errorhandling = "remove")%do%{
  
  #i=unique(weather_slice$stat)[4][1]; j = unique(weather_slice$variable)[1][1]; i = unique(weather_slice$win)[1]          
  message(paste(g, i, j, k,sep = " ") ) 
           
  weather_slice %>%
    filter(win == i) %>%
    filter(variable == j) %>%
    filter(stat == k) %>%
    left_join(traits, by = "sampleId_orig") %>%
    filter(Generation == g) -> tmp.obj
  
  foreach( p = 0:100,
           .combine = "rbind",
           .errorhandling = "remove")%do%{
             
  if(p == 0){
    #model0 <- lmer(Ctmin ~ city + (1 | fruit_type), data = tmp.obj)  
    #model1 <- lmer(Ctmin ~ city + (1 | fruit_type) + value, data = tmp.obj)  
    model0 <- lmer(ctmin ~ (1 | site) + (1 | `Fruit type`), data = tmp.obj)  
    model1 <- lmer(ctmin ~ (1 | site) + (1 | `Fruit type`) + value, data = tmp.obj)  
    #model0 <- lm(Ctmin ~ city +  fruit_type, data = tmp.obj)  
    #model1 <- lm(Ctmin ~ city +  fruit_type + value, data = tmp.obj)  
    
    p_lrt=anova(model1, model0, test="Chisq")[2,8]
    #p_lrt=anova(model1, model0, test="Chisq")[2,5]
  } else if(p != 0){
   
   ## Randomize traits with fruit and city
   tmp.obj[sample(dim(tmp.obj)[1]),] %>%
     dplyr::select(ctmin, site, `Fruit type`) ->
     randomized_sample
   
   tmp.obj %>% 
     dplyr::select(value) ->
     true_value
   
   ran.tmp = cbind(randomized_sample, true_value)
   
   model0 <- lmer(ctmin ~ (1 | site) + (1 | `Fruit type`), 
                  data = ran.tmp)  
   model1 <- lmer(ctmin ~ (1 | site) + (1 | `Fruit type`) + value, 
                  data = ran.tmp)  
   p_lrt=anova(model1, model0, test="Chisq")[2,8]
   
 }
             
  data.frame(
    win= i,
    p=p,
    Generation=g,
    variable = j,
    stat= k,
    p_lrt=p_lrt
  ) 
}}}}}
o

o %<>%
  mutate(win = factor(win, levels = unique(weather_slice$win)))

o %<>%
  mutate(Perm_type = case_when(
    p == 0 ~ "real",
    p != 0 ~ "perm"
  ))
###here save o and sent to JCBN
save(o, file = "GLM_models_output_JCBN_Jun6.Rdata")

