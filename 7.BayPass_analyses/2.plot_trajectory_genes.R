### Plot the gene trajectories.

library(knitr)
library(zoo)
library(RColorBrewer)
library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
require(readODS)
require(poolfstat)
require(vioplot)
require(data.table)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(foreach)
library(forcats)
library(magrittr)


####
##### Tracking analysis
samps <- fread("/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023/2024_Joaquin_KY/swd.metadata.txt")

dat_f <- "all.pooldata.new.rds"    
all_dat <- readRDS(dat_f)

sample.names <- all_dat@poolnames
exclude <- c("KY20")#,"KY17", "KY11")

ky.tag <- grep("KY", sample.names)
indexer.ky <- data.frame(name = sample.names[grep("KY", sample.names)],
           index = ky.tag)   
ky.tag.flt <- indexer.ky$index[which(!indexer.ky$name %in% exclude)]

### load snp info
load("snp_info.Rdata")
snp_info %>% 
filter(snp_id %in% 
c(
"chr3_63316068",
"chr3_67001738",
"chr3_65836197",
"chrX_11184269",
"chr2R_20147914")
) %>%
mutate(rs.id = rownames(.)) %>%
separate(remove = F, rs.id, into = c("feat", "index"),
sep = "s") ->
snp_info.df

## true seasonal
#2 XM_017070295.3 (kra) (3) 63316068 0.006262193 0.003147274
#3 XM_036816734.2 (SNF4Aγ) (3) 67001738 0.004902951 0.044552508
#4 XM_036819182.2 (Nlg1) (3) 65836197 0.033690047 0.022154055
#5 XM_036820755.2 (mamo) (X) 11184269 0.007601839 0.013188811

## spatially discordant seasonal
# 2R 20147914 -- Dg

Trajectory_genes <-
	pooldata.subset(
	all_dat,
	pool.index = ky.tag.flt,
	snp.index = as.numeric(snp_info.df$index),
	min.cov.per.pool = 10,
	max.cov.per.pool = 150,
	min.maf = 0.05,
	return.snp.idx = TRUE,
	verbose = TRUE
	)

ref_count <- Trajectory_genes@refallele.readcount
coverage <- Trajectory_genes@readcoverage
afs <- ref_count/coverage

afs %>%
as.data.frame %>%
mutate(snp_id = snp_info.df$snp_id ) ->
afs.id

names(afs.id) = c(Trajectory_genes@poolnames, "snp_id")

traits <- fread("Means_CT_min_FINAL.csv")

afs.id %>%
melt(id = "snp_id", variable.name = "sampleId_orig") %>%
left_join(samps) %>%
separate(remove = F, snp_id,
into = c("chr", "pos"),
sep = "_") %>%
left_join(traits) ->
afs.id.annot

afs.id.annot$pos = as.numeric(afs.id.annot$pos)

afs.id.annot %>%
group_by(snp_id) %>%
ggplot(aes(
x=`Juliean Date`,
y=value,
color = Time.point,
shape=as.factor(year),,
group=snp_id
)) + geom_point()+
geom_smooth(method = "lm", se = F, 
color = "black" 
) +
facet_grid(site~factor(snp_id, levels = c(
"chr3_63316068" ,
"chr3_65836197",
"chr3_67001738",
"chrX_11184269",
"chr2R_20147914"  
))) + theme_bw() ->
af_trajectories.yday

ggsave(af_trajectories.yday, file = "af_trajectories.yday.pdf",
w=6, h = 2.5)

