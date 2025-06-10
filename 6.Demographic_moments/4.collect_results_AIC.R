#### Moments analysis for SWD

library(tidyverse)
library(magrittr)
library(data.table)
library(reshape2)
library(foreach)

nullf <- system("ls null/*_output.null.txt", intern = T)

split_migf <- system("ls split_mig/*_output.migsym.txt", intern = T)

split_migasymf <- system("ls split_migasym/*_output.migasymm.txt", intern = T)

o1 = foreach(i=nullf, .combine = "rbind")%do%{
  tmp <- fread(i)
  
  data.table(Pair_name = tmp$Pair_name,
             AIC = tmp$AIC) %>%
    mutate(model = "null")
}

o2 = foreach(i=split_migf, .combine = "rbind")%do%{
  tmp <- fread(i)
  
  data.table(Pair_name = tmp$Pair_name,
             AIC = tmp$AIC) %>%
    mutate(model = "split_migf")
}

o3 = foreach(i=split_migasymf, .combine = "rbind")%do%{
  tmp <- fread(i)
  
  data.table(Pair_name = tmp$Pair_name,
             AIC = tmp$AIC) %>%
    mutate(model = "split_migasymf")
}

rbind(o1, o2, o3) ->
  AIC_models

AIC_models %>%
  group_by(model) %>%
  summarise(m = mean(log10(AIC)),
            sd = sd(log10(AIC)))
  

AIC_models %>%
  group_by(Pair_name, model) %>%
  slice_min(AIC) %>%
  ggplot(aes(
    x=model,
    y=log10(AIC)
  )) + geom_boxplot() ->
  AICplots

ggsave(AICplots, file = "AICplots.pdf")

####
bets_model = foreach(i=split_migf, .combine = "rbind")%do%{
  tmp <- fread(i)
  tmp
}

bets_model2 = foreach(i=split_migasymf, .combine = "rbind")%do%{
  tmp <- fread(i)
  tmp
}
#divergence_time = 2*Nref*Ts*g
#Nref= theta/(4*mu*L)

bets_model  %>%
  group_by(Pair_name) %>%
  slice_min(AIC) %>%
  mutate(Nref2 = theta/(4*2.8e-9*146738399)) %>%
  mutate(divT2 = 2*Nref2*Ts) %>%
  mutate(migR2 = m12/(2*Nref2)) %>%
  mutate(mKY_VA = migR2*(nu1*Nref2),
         mVA_KY = migR2*(nu2*Nref2),
  ) %>%
  mutate(NeKY = nu1*Nref2,
         NeVA = nu2*Nref2) ->
  updated_Estimates

###
updated_Estimates %>%
  ungroup() %>%
  summarise(max = max(divT2),
            mean = mean(divT2),
            med = median(divT2),
            min = min(divT2))

updated_Estimates %>%
  ungroup() %>%
  slice_min(AIC) %>% .$divT2

####
updated_Estimates %>%
  ungroup() %>%
  summarise(max = max(mKY_VA),
            mean = mean(mKY_VA),
            med = median(mKY_VA),
            min = min(mKY_VA))

updated_Estimates %>%
  ungroup() %>%
  slice_min(AIC) %>% .$mKY_VA

####
updated_Estimates %>%
  ungroup() %>%
  summarise(max = max(mVA_KY),
            mean = mean(mVA_KY),
            med = median(mVA_KY),
            min = min(mVA_KY))

updated_Estimates %>%
  ungroup() %>%
  slice_min(AIC) %>% .$mVA_KY

####
updated_Estimates %>%
  ungroup() %>%
  summarise(max = max(NeKY),
            mean = mean(NeKY),
            med = median(NeKY),
            min = min(NeKY))

updated_Estimates %>%
  ungroup() %>%
  slice_min(AIC) %>% .$NeKY

####
updated_Estimates %>%
  ungroup() %>%
  summarise(max = max(NeVA),
            mean = mean(NeVA),
            med = median(NeVA),
            min = min(NeVA))

updated_Estimates %>%
  ungroup() %>%
  slice_min(AIC) %>% .$NeVA
