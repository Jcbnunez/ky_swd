library(tidyverse)
library(magrittr)
library(data.table)
library(reshape2)
library(foreach)
library(forcats)
library(lme4)

###load data
traits <- fread("/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023/Ellie_MS/Phenotype_data/Means.CT.min.csv", header = T) 
traits %<>%
  mutate(city = case_when(site == "Berea" ~ "Berea",
                          site != "Berea" ~ "Lexington"
  )) %>%
  mutate(Sample = paste(Generation, Time.point, city,  year, sep = "_" ))

weather_slice <- get(load("/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023/Ellie_MS/Weather_data/SlicedWeatherData.JCBN_MOG.Apr16.Rdata"))
weather_slice %<>%
  mutate(win = factor(win, levels = unique(weather_slice$win)))

meta <- fread("/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023/MASTER_SWD_METADATA-METADATA_v2_MG.csv")
#names(weather_slice)[1] = "sampleId_orig"
#output.file <- fread("/gpfs2/scratch/mgarvin1/kellerflycollaboration/output.csv")

##view(traits)
#### Run models
#i = weather_slice$win[1];j = weather_slice$variable[1];k = weather_slice$stat[1]
#i=1
#j=1
#k=1


o =
  foreach( k = unique(weather_slice$stat),
           .combine = "rbind",
           .errorhandling = "remove")%do%{
             foreach( j = unique(weather_slice$variable),
                      #j =T2M,
                      .combine = "rbind",
                      .errorhandling = "remove")%do%{
                        #print("testing checkpoint 1")
                        if (!(k %in% c("prop. min", "prop. max") && j != "T2M")) {
                        foreach( i = unique(weather_slice$win),
                                 .combine = "rbind",
                                 .errorhandling = "remove")%do%{
                                   #print("testing checkpoint 2")
                                   
                                   message(paste(i, j, k,sep = " ") ) 
                                  
                                   weather_slice %>%
                                     filter(win == i) %>%
                                     filter(variable == j) %>%
                                     filter(stat == k) %>%
                                     #filter( sampleId_orig != "KY20") %>%
                                     left_join(traits[,c("Sample","ctmin")]) -> tmp.obj
                                   
                                   #message(head(tmp.obj))
                                   
                                   foreach( p = 0:100,
                                            
                                            .combine = "rbind",
                                            .errorhandling = "remove")%do%{
                                              
                                              #message("testing checkpoint 3??")
                                              if(p == 0){
                                                #model0 <- lmer(Ctmin ~ city + (1 | fruit_type), data = tmp.obj)  
                                                #model1 <- lmer(Ctmin ~ city + (1 | fruit_type) + value, data = tmp.obj)  
                                                model0 <- lmer(ctmin ~ (1 | city) + (1 | fruit_type), data = tmp.obj)  
                                                model1 <- lmer(ctmin ~ (1 | city) + (1 | fruit_type) + value, data = tmp.obj)  
                                                #model0 <- lm(Ctmin ~ city +  fruit_type, data = tmp.obj)  
                                                #model1 <- lm(Ctmin ~ city +  fruit_type + value, data = tmp.obj)  
                                                
                                                #message(head(model1))
                                                
                                                p_lrt=anova(model1, model0, test="Chisq")[2,8]
                                                #p_lrt=anova(model1, model0, test="Chisq")[2,5]
                                              } else if(p != 0){
                                                
                                                ## Randomize traits with fruit and city
                                                tmp.obj[sample(dim(tmp.obj)[1]),] %>%
                                                  dplyr::select(CTmin, city, fruit_type) ->
                                                  randomized_sample
                                                
                                                tmp.obj %>% 
                                                  dplyr::select(value) ->
                                                  true_value
                                                
                                                ran.tmp = cbind(randomized_sample, true_value)
                                                
                                                model0 <- lmer(CTmin ~ (1 | city) + (1 | fruit_type), 
                                                               data = ran.tmp)  
                                                model1 <- lmer(CTmin ~ (1 | city) + (1 | fruit_type) + value, 
                                                               data = ran.tmp)  
                                                p_lrt=anova(model1, model0, test="Chisq")[2,8]
                                                
                                              }
                                              
                                              #message(head(p_lrt))
                                              #message("testing checkpoint 4")
                                              
                                              data.frame(
                                                win= i,
                                                p=p,
                                                variable = j,
                                                stat= k,
                                                p_lrt=p_lrt
                                              ) 
                                            }}}}}

##head(p_lrt)
#check if exists
#o 


o <- output.file #only run if an output file is already created!!

o %<>%
  mutate(win = factor(win, levels = unique(weather_slice$win)))

o %<>%
  mutate(Perm_type = case_when(
    p == 0 ~ "real",
    p != 0 ~ "perm"
  ))
###here save o and sent to JCBN
save(o, file = "GLM_models_output_MOG.Rdata")

#####
o %>%
  filter(!is.na(p_lrt))%>%
  group_by(Perm_type, win ,variable, stat) %>%
  summarise(uci = quantile(p_lrt, 0.05)) %>%
  dcast(win + variable + stat ~ Perm_type)

significant_test <- o %>%
  filter(!is.na(p_lrt))%>%
  mutate(p0 = (p == 0)) %>%         # Create a column with TRUE/FALSE for p==0
  group_by(p0, win) %>%             # Group by p0 and win
  summarise(m_Plrt = quantile(p_lrt, 0.05)) %>%
  dcast(win ~ p0, value.var = "m_Plrt") %>% 
  mutate(TEST = `FALSE` > `TRUE`)

head(o)

o %<>% left_join(significant_test)


##Violin plot!
violin.plot <- ggplot() +
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

violin.plot
ggsave(violin.plot, file = "violin.plot.pdf")

#create observed p value
observed <- o %>%
  filter(p == 0) %>%
  select(win, variable, stat, obs_p = p_lrt)

#create emperical p value
empirical <- o %>%
  filter(p != 0) %>%
  left_join(observed, by = c("win", "variable", "stat")) %>%
  group_by(win, variable, stat, obs_p) %>%
  summarise(emp_p = mean(p_lrt <= obs_p), .groups = "drop")


# Create the birds-eye view plot!!
birdseye_plot <- ggplot(empirical, aes(x = win, y = emp_p)) + #change to obs_p?? ask jcbn
  geom_point(size = 3, aes(color = emp_p <= 0.05)) + #<-- identifies/sorts color
  scale_color_manual(name = "Significance",
                     values = c("TRUE" = "red", "FALSE" = "black")) + #makes colors!!
  facet_grid(stat ~ variable, scales = "free_y") + #not working??? #now working
  labs(title = "Empirical p–values by Window",
       x = "Window",
       y = "Empirical p–value") +
  scale_y_continuous(limits = c(0, 1)) + #sets the ylim to be 0-1
  theme_minimal() +
  theme(strip.text = element_text(size = 10),
        axis.text = element_text(size = 8),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),)


birdseye_plot


#head(empirical)
##### Create weather slices for sequenced data

