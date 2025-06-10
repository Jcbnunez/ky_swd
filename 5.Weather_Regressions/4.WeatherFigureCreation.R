library(tidyverse)
library(magrittr)
library(data.table)
library(reshape2)
library(foreach)
library(forcats)
library(lme4)

#####load data
weather_slice <- get(load("/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023/Ellie_MS/Weather_data/SlicedWeatherData.JCBN_MOG.Jun5_2025.Rdata"))
weather_slice %<>%
  mutate(win = factor(win, levels = unique(weather_slice$win))) %>%
  mutate(sampleId_orig = Sample)


traits <- fread("/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023/Ellie_MS/Phenotype_data/Means_CT_min_FINAL.csv") 
traits %<>%
  mutate(sampleId_orig = paste(Generation, 
                               Time.point, 
                               ifelse(site == "Berea", "Berea", "Lexington"), 
                               year, sep = "_"))

##view(traits)
#### Run models
#i = weather_slice$win[1];j = weather_slice$variable[1];k = weather_slice$stat[1]
#i=1;j=1;k=1;g="f4"

#### STOP there is a check point!

o =
  foreach( g = c("f4","Parental"),
           .combine = "rbind",
           .errorhandling = "remove")%do%{
  foreach( k = unique(weather_slice$stat),
           .combine = "rbind",
           .errorhandling = "remove")%do%{
             foreach( j = unique(weather_slice$variable),
                      #j =T2M,
                      .combine = "rbind",
                      .errorhandling = "remove")%do%{
                        #print("testing checkpoint 1")
                        #if (!(k %in% c("prop. min", "prop. max") && j != "T2M")) {
                        foreach( i = unique(weather_slice$win),
                                 .combine = "rbind",
                                 .errorhandling = "remove")%do%{
                                   #g="f4";k="prop. max";j="T2M";i="c(0, 7)"
                                   message(paste(g, i, j, k,sep = " ")) 
                                   
                                   weather_slice %>%
                                     filter(win == i) %>%
                                     filter(variable == j) %>%
                                     filter(stat == k) %>%
                                     left_join(traits[,c("sampleId_orig","ctmin","Generation",
                                                         "Fruit type","site", "year")]) %>%
                                     filter(Generation == g) %>%
                                     mutate(site_year = paste(site, year, sep ="_"))-> tmp.obj
                                   
                                   foreach( p = 0:100,
                                            .combine = "rbind",
                                            .errorhandling = "remove")%do%{
                                              
                                              #message("testing checkpoint 3??")
                                              if(p == 0){
                                                #model0 <- lmer(Ctmin ~ city + (1 | fruit_type), data = tmp.obj)  
                                                #model1 <- lmer(Ctmin ~ city + (1 | fruit_type) + value, data = tmp.obj)  
                                                model0 <- lmer(ctmin ~ (1 | site_year) + (1 | `Fruit type`), data = tmp.obj)  
                                                model1 <- lmer(ctmin ~ (1 | site_year) + (1 | `Fruit type`) + value, data = tmp.obj)  
                                                #model0 <- lm(Ctmin ~ city +  fruit_type, data = tmp.obj)  
                                                #model1 <- lm(Ctmin ~ city +  fruit_type + value, data = tmp.obj)  
                                                
                                                #message(head(model1))
                                                
                                                p_lrt=anova(model1, model0, test="Chisq")[2,8]
                                                #p_lrt=anova(model1, model0, test="Chisq")[2,5]
                                              } else if(p != 0){
                                                
                                                ## Randomize traits with fruit and city
                                                tmp.obj[sample(dim(tmp.obj)[1]),] %>%
                                                  dplyr::select(ctmin, site_year, `Fruit type`) ->
                                                  randomized_sample
                                                
                                                tmp.obj %>% 
                                                  dplyr::select(value) ->
                                                  true_value
                                                
                                                ran.tmp = cbind(randomized_sample, true_value)
                                                
                                                model0 <- lmer(ctmin ~ (1 | site_year) + (1 | `Fruit type`), 
                                                               data = ran.tmp)  
                                                model1 <- lmer(ctmin ~ (1 | site_year) + (1 | `Fruit type`) + value, 
                                                               data = ran.tmp)  
                                                p_lrt=anova(model1, model0, test="Chisq")[2,8]
                                                
                                              }
                                              
                                              #message(head(p_lrt))
                                              #message("testing checkpoint 4")
                                              
                                              data.frame(
                                                win= i,
                                                p=p,
                                                Generation = g,
                                                variable = j,
                                                stat= k,
                                                p_lrt=p_lrt
                                              ) 
                                              
                                            }
                                   
                                 }
                          
                        }
                        
                      }
             
             
          # }
             
             
             }

save(o, file = "GLM_models_output_JCBN_Jun6.Rdata")

##head(p_lrt)
#check if exists
#o 


#o <- output.file #only run if an output file is already created!!
#
#o %<>%
#  mutate(win = factor(win, levels = unique(weather_slice$win)))
#
#o %<>%
#  mutate(Perm_type = case_when(
#    p == 0 ~ "real",
#    p != 0 ~ "perm"
#  ))
####here save o and sent to JCBN
#
#
#####
o <- get(load("GLM_models_output_JCBN_Jun6.Rdata"))

o %<>% mutate(Perm_type = p==0)

o %<>%
  filter(!is.na(p_lrt))%>%
  group_by(Generation, Perm_type, win ,variable, stat) %>%
  summarise(uci = quantile(p_lrt, 0.05)) %>%
  dcast(win + variable + stat + Generation ~ Perm_type) %>%
  mutate(sig =`TRUE` < `FALSE`) 

o %>% 
  filter(Generation == "f4") %>%
  filter(sig == T) %>%
  mutate(model_set = paste(win,variable,stat,sep ="_")) ->
  f4sigs
o %>% 
  filter(Generation == "Parental") %>%
  filter(sig == T) %>%
  mutate(model_set = paste(win,variable,stat,sep ="_"))->
  Parsigs

Parsigs$model_set[which(Parsigs$model_set %in% f4sigs$model_set)] ->
  models_significant

o %>% 
  mutate(model_set = paste(win,variable,stat,sep ="_")) %>%
  filter(model_set %in% models_significant[1])


weather_slice %>%
  filter(win == "c(0, 7)") %>%
  filter(variable == "T2M") %>%
  filter(stat == "prop. max") %>%
  left_join(traits) -> best.model



