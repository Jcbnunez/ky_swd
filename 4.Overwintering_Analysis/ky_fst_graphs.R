#########################
### Loading Libraries ###
#########################

library(tidyverse)
library(magrittr)
library(data.table)
library(reshape2)
library(poolfstat)
library(FactoMineR)
require(foreach)
library(rnaturalearth)
library(rnaturalearthdata)
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

setwd("/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/ky_swd/4.Overwintering_Analysis/plots")

root="/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/ky_swd/4.Overwintering_Analysis/plots"
#require(gtools)
#library(ggpmisc)
#library(gt)

load("/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/ky_swd/4.Overwintering_Analysis/fst.matrix.meta.ky.exc201711.Rdata")
load("/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/ky_swd/4.Overwintering_Analysis/fst.matrix.meta.kyva.exc201711.Rdata")
load("/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/ky_swd/4.Overwintering_Analysis/PCA.coords.kyva.exc201711.Rdata")

ran <- get(load("/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/ky_swd/4.Overwintering_Analysis/pca.ran.kyva.exc201711.Rdata"))


##################################################################################
################################################################################## 
### Graphs ###
##################################################################################

###### Global PCA
load("/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/ky_swd/4.Overwintering_Analysis/pca.all.Rdata")
#> all.pca$perc.var
#[1] 6.620502 3.177675 2.309744 2.102024 2.065136 2.005893 2.001910 1.953411
#[9] 1.948237 1.939104 1.932917 1.926837 1.923702 1.912197 1.908626 1.904652
#[17] 1.897711 1.894595 1.891876 1.889637 1.886755 1.879150 1.873646 1.870153
#[25] 1.865852 1.859608 1.857100 1.855714 1.848746 1.846287 1.844100 1.835441
#[33] 1.830291 1.823591 1.806633 1.792963 1.791634 1.789299 1.771471 1.764768
#[41] 1.747103 1.740598 1.711112 1.660656 1.649559 1.630631 1.630106 1.602117
#[49] 1.587137 1.573638 1.567755
load("/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/ky_swd/4.Overwintering_Analysis/pca.Name.Rdata")
#Name.pca$perc.var
#[1] 4.888268 3.770809 3.586419 3.538159 3.516198 3.509006 3.487864 3.467869
#[9] 3.459886 3.450258 3.443352 3.425772 3.423039 3.413686 3.401709 3.399235
#[17] 3.387287 3.381214 3.374449 3.359205 3.350126 3.340085 3.331495 3.329094
#[25] 3.324156 3.293411 3.272662 3.168484 2.906803

separate(pca.all, Range, into = c("range","cont"), sep = "_") %>% 
  ggplot(aes(
    x=V1,
    y=V2,
    shape = range,
    fill= continent
  )) + geom_point(size = 2.9) +
  scale_shape_manual(values = 21:22) +
  theme_bw() 

#ggsave(glob.pca12, file = "glob.pca12.pdf",
#       w=5, h = 4)

ggplot(data = world) +
  geom_sf(fill= "grey70") +
    coord_sf(xlim = c(-160, 160), ylim = c(-60, 60), expand = FALSE) + 
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.2), 
        panel.background = element_rect(fill = "white")) +
  geom_jitter(data = separate(pca.all, Range, into = c("range","cont"), sep = "_"),
              color = "black",
              aes(
                x=long,
                y=lat,
                fill = V1,
                shape = range
              ), alpha = 0.9, size = 3.5) + 
  scale_fill_gradient2(low = "steelblue", ,high = "firebrick",
                       midpoint = 0.1) +
  scale_shape_manual(values = 21:22)



pca.Name %>% 
  ggplot(aes(
    x=V1,
    y=V2,
    fill= province
  )) + geom_point(size = 2.9, shape = 21, alpha = 0.7) +
  scale_fill_manual(values = c(
    "orange4","brown","purple","navy","red","springgreen","gold","cyan"
  )) + 
  theme_bw() 



ggplot(data = world) +
  geom_sf(fill= "grey70") +
  coord_sf(xlim = c(-160, -50), ylim = c(15, 50), expand = FALSE) + 
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.2), 
        panel.background = element_rect(fill = "white")) +
  geom_jitter(data = pca.Name,
              color = "black",
              aes(
                x=long,
                y=lat,
                fill = V1,
              ), alpha = 0.9, size = 3.5, shape = 21) + 
  scale_fill_viridis_c()

pca.Name
lm(V1 ~ long*lat, data = pca.Name) %>% anova

##################################################################################
### Figure a ###
##################################################################################
		   
  fast_matrix.meta.ky  %>%
  mutate(comp = case_when(comp == "Lexington_Berea" ~ "Berea_Lexington",
                          TRUE ~ comp)) ->
  fst_ow_dat

fst_ow_dat %>% 
  filter(comp == "Berea_Lexington") %>%
  summarise(m.FST = mean(fst))

fst_ow_dat %>% 
  mutate(M1 = month(as.Date(date1))) %>%
  mutate(M2 = month(as.Date(date2))) %>%
  filter(comp == "Berea_Lexington") %>%
  group_by(year1==year2) %>%
  summarise(m.FST = mean(fst))

fst_ow_dat %>% 
  mutate(M1 = month(as.Date(date1))) %>%
  mutate(M2 = month(as.Date(date2))) %>%
  filter(comp == "Berea_Lexington") %>%
  group_by(year1==year2, M1==M2) %>%
  summarise(m.FST = mean(fst))

fst_ow_dat %>% 
  mutate(M1 = month(as.Date(date1))) %>%
  mutate(M2 = month(as.Date(date2))) %>%
  filter(comp == "Lexington_Lexington") %>%
  group_by(year1,year2) %>%
  summarise(m.FST = mean(fst, na.rm = T))


fst_slope.ky <- fst_ow_dat %>%
  ggplot(aes(
    x=delta,
    y=(fst),
    color = comp
  )) + 
  geom_point(shape = 21, fill = "grey", size = 3.5) +
  geom_smooth(method = "lm") +
  theme_bw() + facet_wrap(.~comp) +
  scale_color_manual(values = c("steelblue","grey2","navyblue")) +
  theme(legend.position = "none")
fst_slope.ky

lm((fst)/(1-(fst)) ~ delta, data = filter(fast_matrix.meta.ky, comp == "Lexington_Berea")) %>% summary
lm((fst)/(1-(fst)) ~ delta, data = filter(fast_matrix.meta.ky, comp == "Berea_Berea")) %>% summary
lm((fst)/(1-(fst)) ~ delta, data = filter(fast_matrix.meta.ky, comp == "Lexington_Lexington")) %>% summary


ggsave(fst_slope.ky, 
       file = paste(root, "fst_slope_ky.pdf", sep = "/" ), 
       w = 8, h = 2.3)


plot.fst_ow_dat.M <- fst_ow_dat %>% 
  mutate(M1 = month(as.Date(date1))) %>%
  mutate(M2 = month(as.Date(date2))) %>%
  mutate(deltaM = abs(M1-M2)) %>%
  mutate(deltaY = abs(year1-year2)) %>%
  filter(deltaY == 0) %>%
  ggplot(aes(
    x=deltaM,
    y=(fst),
    color = comp
  )) + 
  geom_point(shape = 21, fill = "grey", size = 3.5) +
  geom_smooth(method = "lm") +
  theme_bw() + facet_wrap(.~comp) +
  scale_color_manual(values = c("steelblue","grey2","navyblue")) +
  theme(legend.position = "none")
plot.fst_ow_dat.M


#####
##################################################################################
### Figure b ###
##################################################################################
     
 fast_matrix.meta.ky %>%
  mutate(samey = year2==year1) |> mutate(
  samey = if_else(condition = samey == "TRUE",
                  true = "Within Growing Season", 
                  false = "Overwinter")) -> dat_ow

fst.box.ky <- dat_ow%>%
  filter(comp %in% c("Berea_Berea","Lexington_Lexington")) %>%
    ggplot(aes(
    x=comp ,
    y=fst,
    fill=samey)) + 
  geom_boxplot(alpha = 0.7) + 
  scale_fill_manual(values = c("darkslateblue","firebrick"))+
  theme_bw() +
  labs(
	x = NULL)+
  theme(legend.position = "none")
fst.box.ky

ggsave(fst.box.ky, 
       file = paste(root, "fst.box_ky.pdf", sep = "/" ), 
       w= 3, h = 4)

wilcox.test(fst~samey, data = filter(dat_ow,comp == "Berea_Berea"))
wilcox.test(fst~samey, data = filter(dat_ow,comp == "Lexington_Lexington"))

	#ggsave(fst.box.ky, file = "fst.box_ky.png",
	#	   w= 6, h = 5)		   
	#	   
##################################################################################	   
### Figure c ###
##################################################################################

fst.box.kyva <- fast_matrix.meta.kyva |> filter(
	states != "Virginia") |> mutate(
  states = if_else(condition = states == "Kentucky", 
                   true = "Within Kentucky", false = "Between Kentucky and Virginia")) |>
    ggplot(aes(
    x=states,
    y=fst
  )) + 
  geom_boxplot() +
  labs(x = NULL) + 
	ggtitle(expression("c.) F"[ST]~" of Within State vs Different State Populations")) +
	ylab(bquote(F[ST])) +
  theme_bw()+
  theme(legend.position = "none")
fst.box.kyva

	ggsave(fst.box.kyva, 
	       file = paste(root, "fst.box.kyva.pdf", sep = "/" ), 
	       w= 3, h = 4)
	#	ggsave(fst.box.kyva, file = "fst.box.kyva.png",
#		   w= 6, h = 5)

	fast_matrix.meta.kyva %>%
	  filter(states == "DifferentStates") %>%
	  summarise(m.FST = mean(fst))
	
	fast_matrix.meta.kyva %>%
	  filter(states == "Kentucky") %>%
	  summarise(m.FST = mean(fst, na.rm = T))
	
##################################################################################
### Figure d ###	
##################################################################################
names <- c("Berea (KY)", "Charlottesville (VA)", "Lexington (KY)")
PCA12.kyva <- PCA.coords.kyva %>%
  ggplot(aes(
    x=V1,
    y=V2,
    fill= as.character(city)
  )) + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 5, shape = 21) +
  labs(
	x = "Component 1",
	y = "Component 2",
	title = "d.) PCA - Sampling Sites",
	color = NULL) + 
  scale_color_hue(labels = names) +
  theme_bw() +
  scale_fill_manual(values = c("steelblue","orange","navyblue"))
PCA12.kyva
  
  
ggsave(PCA12.kyva, 
       file = paste(root, "PCA12_kyva.pdf", sep = "/" ), 
       w= 8, h = 4)
#ggsave(PCA12.kyva, file = "PCA12_kyva.png",
#       w= 6, h = 5)


PCA.year.kyva <- PCA.coords.kyva %>%
  filter(province == "Kentucky") %>%
  ggplot(aes(
    x=as.Date(Collection_date),
    y=V7,
  )) + 
  geom_point(size = 2.9, shape = 21, aes(fill= as.character(city))) +
  geom_smooth(method = "lm", se = F, linetype = "dashed", color = "black") +
  theme_bw() +
  scale_fill_manual(values = c("steelblue","navyblue")) 
PCA.year.kyva

###
PCA.coords.kyva %>%
  filter(province == "Kentucky") ->
  ky_data.pca

lm(V7 ~ as.factor(year(ky_data.pca$Collection_date))*city, data = ky_data.pca) %>% anova
#Response: V7
#Df   Sum Sq  Mean Sq F value   Pr(>F)
#as.factor(year(ky_data.pca$Collection_date))       3 0.158400 0.052800 10.0978 0.003071
#city                                               1 0.039793 0.039793  7.6103 0.022160
#as.factor(year(ky_data.pca$Collection_date)):city  3 0.033381 0.011127  2.1280 0.166787
#Residuals                                          9 0.047060 0.005229                 
#
#as.factor(year(ky_data.pca$Collection_date))      **
#  city                                              * 
#  as.factor(year(ky_data.pca$Collection_date)):city   
#Residuals                  
#
