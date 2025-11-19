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
samps <- fread("/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023_BAYPASS_MG/2024_Joaquin_KY/swd.metadata.txt")

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
##Chromosome Position RefAllele AltAllele    snp_id
##rs2        chrX     4133         C         T chrX_4133
##rs3        chrX     4139       ATA        AA chrX_4139
##rs4        chrX     4180         A         T chrX_4180
##rs27       chrX     4577         A         G chrX_4577
##rs29       chrX     4633         G         A chrX_4633
##rs39       chrX     4897         A         G chrX_4897
traits <- fread("Means_CT_min_FINAL.csv")

#### part 1. Top C2 alleles

c2_top_hits_annotated <- fread("Top Genes in C2_SNPs.tsv")
c2_top_hits_annotated %<>% mutate(SNP_id = paste(chr, pos, sep ="_")) #%>%
#  filter(!Gene_Name %in% c("No current annotation"))  %>%
#  filter(!is.na(Gene_Name))

snp_info %>% 
filter(snp_id %in% 
c(
  c2_top_hits_annotated$SNP_id)
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

### PCA bit
pca_result <- randomallele.pca(Trajectory_genes, scale = TRUE)

pca_result$pop.loadings %>% as.data.frame() %>%
  mutate(sampleId_orig=row.names(.)) %>%
  left_join(samps) %>%
  left_join(traits)-> pc_loadings

pc_loadings %>%
  ggplot(aes(
    x=V1,
    y=V2,
    label = sampleId_orig,
    fill=Time.point,
  )) + geom_point(size = 4, shape = 21) +
  geom_text() +
  theme_bw()->
  PCA_plot
ggsave(PCA_plot, file = "PCA12.dim.png", w = 5, h =4)

## Trajectory
ref_count <- Trajectory_genes@refallele.readcount
coverage <- Trajectory_genes@readcoverage
afs <- ref_count/coverage

afs %>%
as.data.frame %>%
mutate(snp_id = snp_info.df$snp_id ) ->
afs.id

names(afs.id) = c(Trajectory_genes@poolnames, "snp_id")


afs.id %>%
reshape2::melt(id = "snp_id", variable.name = "sampleId_orig") %>%
left_join(samps) %>%
separate(remove = F, snp_id,
into = c("chr", "pos"),
sep = "_") %>%
left_join(traits) ->
afs.id.annot

afs.id.annot$pos = as.numeric(afs.id.annot$pos)

afs.id.annot %>%
group_by(snp_id, Time.point, year) %>%
  summarize(AFm = mean(value)) %>%
  filter(AFm > 0 & AFm < 1) %>%
  reshape2::dcast(snp_id~Time.point+year) %>%
  .[complete.cases(.),] %>%
  reshape2::melt(id = c("snp_id"), value.var = "AFm" ) ->
  plot_freq_data

ggplot() +
geom_line( data =plot_freq_data,
           aes(
  x=factor(variable, levels = c("First_2020","Second_2020","Third_2020",
                                "First_2021","Second_2021","Third_2021",
                                "First_2022","Second_2022","Third_2022",
                                "First_2023","Second_2023","Third_2023"
  )),
  y=value,
  group=snp_id
),
alpha = 0.3) +
  geom_line( data =filter(plot_freq_data, snp_id == "chr2R_9856476"),
             aes(
               x=factor(variable, levels = c("First_2020","Second_2020","Third_2020",
                                             "First_2021","Second_2021","Third_2021",
                                             "First_2022","Second_2022","Third_2022",
                                             "First_2023","Second_2023","Third_2023"
               )),
               y=value,
               group=snp_id
             ),
             alpha = 0.9, color = "blue", size = 1.3) +
  geom_line( data =filter(plot_freq_data, snp_id == "chr3_63767211"),
             aes(
               x=factor(variable, levels = c("First_2020","Second_2020","Third_2020",
                                             "First_2021","Second_2021","Third_2021",
                                             "First_2022","Second_2022","Third_2022",
                                             "First_2023","Second_2023","Third_2023"
               )),
               y=value,
               group=snp_id
             ),
             alpha = 0.9, color = "red", size = 1.3) +
theme_bw() ->
af_trajectories.yday

ggsave(af_trajectories.yday, file = "af_trajectories.yday.pdf",
w=6, h = 2.5)

#chr2R_9856476
#chr2R_9856487

###
afs.id.annot %>%
  group_by(snp_id) %>%
  ggplot(aes(
    x=as.Date(Date, format = "%m/%d/%Y"),
    y=value,
    color=site,
    #shape=as.factor(year),,
  )) + geom_line()+
  #geom_smooth(method = "lm", se = F, 
  #            color = "black" ) +
  facet_grid(year~factor(snp_id), scales = "free_x") + theme_bw() ->
  af_trajectories.yday

ggsave(af_trajectories.yday, file = "af_trajectories.yday.pdf",
       w=6, h = 2.5)


#### Part 2. 
#### Part 2. Plot T32 and CTmin alleles
#### Part 2. 
ctmin_t32_BF <- fread("ctmin.t32.annot.snps.txt")
Weather_Data <- fread("/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023_BAYPASS_MG/Ellie_MS/Weather_data/SlicedWeatherData.mgarvin.Feb25.csv")

snp_info %>% 
  filter(snp_id %in% 
           c(
             ctmin_t32_BF$SNP_id)
  ) %>%
  mutate(rs.id = rownames(.)) %>%
  separate(remove = F, rs.id, into = c("feat", "index"),
           sep = "s") ->
  snp_info.df2

Trajectory_genes2 <-
  pooldata.subset(
    all_dat,
    pool.index = ky.tag.flt,
    snp.index = as.numeric(snp_info.df2$index),
    min.cov.per.pool = 10,
    max.cov.per.pool = 150,
    min.maf = 0.05,
    return.snp.idx = TRUE,
    verbose = TRUE
  )

### Trajectories
ref_count2 <- Trajectory_genes2@refallele.readcount
coverage2 <- Trajectory_genes2@readcoverage
afs2 <- ref_count2/coverage2

afs2 %>%
  as.data.frame %>%
  mutate(snp_id = snp_info.df2$snp_id ) ->
  afs.id2

names(afs.id2) = c(Trajectory_genes2@poolnames, "snp_id")

afs.id2 %>%
  reshape2::melt(id = "snp_id", variable.name = "sampleId_orig") %>%
  left_join(samps) %>%
  separate(remove = F, snp_id,
           into = c("chr", "pos"),
           sep = "_") %>%
  left_join(traits) ->
  afs.id.annot2

afs.id.annot2$pos = as.numeric(afs.id.annot2$pos)

Weather_Data %>%
  filter(stat == "prop. max") %>%
  filter(win == "c(0, 15)") %>%
  mutate(sampleId_orig = Sample)->
  T32_data

afs.id.annot2 %>% 
  left_join(T32_data, by = "sampleId_orig") ->
  afs.id.annot2.weath

unique(afs.id.annot2.weath$snp_id) -> snps2

corrs = foreach(i= snps2, .combine = "rbind")%do%{
  
  afs.id.annot2.weath %>%
    filter(snp_id == i) -> tmp
  cor.test(~value.x+ctmin, data = tmp) -> ctmin
  cor.test(~value.x+value.y, data = tmp) -> t32
  data.frame(
    snp_id = i,
    ctminp = ctmin$p.value,
    ctminr = ctmin$estimate,
    t32p = t32$p.value,
    t32r = t32$estimate
    
  )
}


corrs %>%
  ggplot(aes(
    x=ctminr,
    y=t32r,
  )) +
  geom_point()  ->
  ctmin.aft32.corr

ggsave(ctmin.aft32.corr, file = "ctmin.aft32.corr.pdf")

corrs %>%
  filter(ctminr > 0.65 & t32r < 0.65) -> topR1
corrs %>%
  filter(ctminr < 0.65 & t32r > 0.65) -> topR2

afs.id.annot2.weath %>%
  filter(snp_id %in% "chr3_80328254") %>%
  filter(value.x > 0 & value.x < 1) %>%
  ggplot(aes(
    x= value.x,
    y= ctmin,
    fill = value.x,
    #linetype = city
  )) + geom_point(shape = 21, size = 3)+ geom_smooth(method = "lm", se = F) +
  scale_fill_gradient2(midpoint = 0.1, low = "steelblue", high = "firebrick") +
  xlab("AF") + ylab("CTmin") + theme_bw()->
  ctmin.cor.plot

ggsave(ctmin.cor.plot, file = "ctmin.cor.plot.pdf", h = 3, w = 4)
