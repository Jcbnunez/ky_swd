#########################
### Loading Libraries ###
#########################

library(tidyverse)
library(magrittr)
library(data.table)
library(reshape2)
library(poolfstat)
library(FactoMineR)
require(gtools)
require(foreach)
library(ggpmisc)
library(gt)

##################################################################################
################################################################################## 
### Graphs ###
##################################################################################
##################################################################################
### Figure a ###
##################################################################################
		   
fst_slope.ky <- fast_matrix.meta.ky  %>%
# mutate(same_city = city1==city2,
#		same_fruit = fruit1==fruit2,
#		fruit_comp = paste(fruit1, fruit2, sep="_")) %>%
  ggplot(aes(
    x=delta,
    y=(fst)
  )) + 
  geom_point() +
  geom_smooth(method = "lm") +
  labs(
	title = "a.) Fst Slope Over Time",
	x = "Difference in Days between Pairwise Samples") + 
  ylab(bquote(F[ST])) +
	stat_poly_eq(formula = my.formula, 
		   aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
		   parse = TRUE)
  
ggsave(fst_slope.ky, file = "fst_slope_ky.pdf", w = 6, h = 5)
ggsave(fst_slope.ky, file = "fst_slope_ky.png", w = 6, h = 5)		   

##################################################################################
### Figure b ###
##################################################################################
     
fst.box.ky <- fast_matrix.meta.ky %>%
  mutate(samey = year2==year1) |> mutate(
  samey = if_else(condition = samey == "TRUE", true = "Within Growing Season", false = "Overwinter"))%>%
    ggplot(aes(
    x=samey,
    y=fst )) + geom_boxplot() + labs(
	x = NULL) + 
	ggtitle(expression("b.) F"[ST]~" of Within Season vs Overwintering Populations")) +
	ylab(bquote(F[ST]))
	
	
	ggsave(fst.box.ky, file = "fst.box_ky.pdf",
		   w= 6, h = 5)
	ggsave(fst.box.ky, file = "fst.box_ky.png",
		   w= 6, h = 5)		   
		   
##################################################################################	   
### Figure c ###
##################################################################################

fst.box.kyva <- fast_matrix.meta.kyva |> filter(
	states != "Virginia") |> mutate(
  states = if_else(condition = states == "Kentucky", true = "Within Kentucky", false = "Between Kentucky and Virginia")) |>
    ggplot(aes(
    x=states,
    y=fst
  )) + 
  geom_boxplot() +
  labs(x = NULL) + 
	ggtitle(expression("c.) F"[ST]~" of Within State vs Different State Populations")) +
	ylab(bquote(F[ST]))
  
	ggsave(fst.box.kyva, file = "fst.box.kyva.pdf",
		   w= 6, h = 5)		
	ggsave(fst.box.kyva, file = "fst.box.kyva.png",
		   w= 6, h = 5)

##################################################################################
### Figure d ###	
##################################################################################
names <- c("Berea (KY)", "Charlottesville (VA)", "Lexington (KY)")
PCA12.kyva <- PCA.coords.kyva %>%
  ggplot(aes(
    x=V1,
    y=V2,
    color= as.character(city)
  )) + 
  geom_point(size = 5) +
  labs(
	x = "Component 1",
	y = "Component 2",
	title = "d.) PCA - Sampling Sites",
	color = NULL) + 
  scale_color_hue(labels = names)
  
  
ggsave(PCA12.kyva, file = "PCA12_kyva.pdf",
       w= 6, h = 5)
ggsave(PCA12.kyva, file = "PCA12_kyva.png",
       w= 6, h = 5)
	   
	   