### Plot the gene trajectories.

library(knitr)
library(zoo)
library(viridis)
library(RColorBrewer)
library(data.table)
library(ggplot2)
library(grid)
library(gmodels)
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
source("/netfiles/nunezlab/Shared_Resources/Software/baypass_public/utils/baypass_utils.R")

####
##### Tracking analysis
samps <- fread("/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023_BAYPASS_MG/2024_Joaquin_KY/swd.metadata.txt")

dat_f <- "/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023_BAYPASS_MG/2024_Joaquin_KY/all.pooldata.new.rds"    
all_dat <- readRDS(dat_f)

sample.names <- all_dat@poolnames
exclude <- c("KY20")#,"KY17", "KY11")

ky.tag <- grep("KY", sample.names)
indexer.ky <- data.frame(name = sample.names[grep("KY", sample.names)],
           index = ky.tag)   
ky.tag.flt <- indexer.ky$index[which(!indexer.ky$name %in% exclude)]

### load snp info
load("/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023_BAYPASS_MG/2024_Joaquin_KY/snp_info.Rdata")
##Chromosome Position RefAllele AltAllele    snp_id
##rs2        chrX     4133         C         T chrX_4133
##rs3        chrX     4139       ATA        AA chrX_4139
##rs4        chrX     4180         A         T chrX_4180
##rs27       chrX     4577         A         G chrX_4577
##rs29       chrX     4633         G         A chrX_4633
##rs39       chrX     4897         A         G chrX_4897
traits <- fread("/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023_BAYPASS_MG/2024_Joaquin_KY/Means_CT_min_FINAL.csv")
Weather_Data <- fread("/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023_BAYPASS_MG/Ellie_MS/Weather_data/SlicedWeatherData.mgarvin.Feb25.csv")
names(Weather_Data)[1] = "sampleId_orig"

Weather_Data %>% 
  filter(win == "c(7, 15)") %>%
  filter(variable == "T2M") %>%
  filter(stat == "prop. max") -> t32dat

#### import data
annots <- "/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023_BAYPASS_MG/Final_Datasets/annots.flt.Rdata"
load(annots)
annots.flt$pos = as.numeric(annots.flt$pos)
#annots.flt %>% filter(chr == "chrX") %>% filter(pos == 16121948)

adapt_dat="/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023_BAYPASS_MG/Final_Datasets/C2_BFct_BFT_df.allAnalyses.Rdata"
adapt=get(load(adapt_dat))

### Left join
adapt %>% left_join(annots.flt) -> adapt.annot


### Lindley Process
### Adaptive Tracking Enrichment
seas.c2.ls.xi2=compute.local.scores(data.frame(adapt.annot[,1:2]),
                                    snp.pi=adapt.annot$freqC2,
                                    snp.pvalue = adapt.annot$P_C2 ,
                                    xi=2)

seas.wins = seas.c2.ls.xi2$significant.windows 

###### find overlaps
setDT(adapt.annot)
setDT(seas.wins)

# SNPs become 1 bp intervals
adapt.annot[, `:=`(beg = pos, end = pos)]
setkey(adapt.annot, chr, beg, end)
setkey(seas.wins, chr, beg, end)

### seasonality
hits.seas <- foverlaps(
  x = adapt.annot,
  y = seas.wins,
  by.x = c("chr", "beg", "end"),
  by.y = c("chr", "beg", "end"),
  type = "within",
  nomatch = NULL
)
seasonal.hits = hits.seas$SNP_id

####
adapt.annot %<>%
  mutate(KYseasonal_outlier = case_when(SNP_id %in% seasonal.hits ~ TRUE,
                                        TRUE ~ FALSE) ) 

top_rna_genes =
c(
  "108021406",
  "108009498",
  "108018373"
)

adapt.annot %>%
  filter(Gene %in% top_rna_genes) %>%
  filter(KYseasonal_outlier == TRUE) %>%
  filter(P_C2 > 3)->
  top_hits_annotated

#### part 1. Top C2 alleles

snp_info %>% 
filter(snp_id %in% 
c(
  top_hits_annotated$SNP_id)
) %>%
mutate(rs.id = rownames(.)) %>%
separate(remove = F, rs.id, into = c("feat", "index"),
sep = "s") ->
snp_info.df

Trajectory_genes <-
	pooldata.subset(
	all_dat,
	pool.index = ky.tag.flt,
	snp.index = as.numeric(snp_info.df$index),
	return.snp.idx = TRUE,
	verbose = TRUE
	)

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
left_join(traits) %>%
left_join(select(t32dat, sampleId_orig, T32= value), by = "sampleId_orig" )  ->
afs.id.annot

afs.id.annot$pos = as.numeric(afs.id.annot$pos)

afs.id.annot %>%
group_by(snp_id, Time.point) %>%
  summarize(AFm = ci(value)[1],
            AFlci=ci(value)[2],
            AFuci=ci(value)[3],
            #mT32=mean(T32), 
            mCTmin = ci(ctmin)[1],
            lciCTmin = ci(ctmin)[2],
            uciCTmin = ci(ctmin)[3],
            ) ->
  plot_freq_data


### cors
cor_results <- plot_freq_data %>%
  group_by(snp_id) %>%
  filter(!any(AFm > 0.95 |  AFm < 0.05)) %>%
  ungroup() %>%
  filter(!is.na(AFm), !is.na(mCTmin)) %>%
  group_by(snp_id) %>%
  summarise(
    n = n(),
    r = cor(AFm, mCTmin, method = "pearson"),
    p_value = cor.test(AFm, mCTmin, method = "pearson")$p.value,
    .groups = "drop"
  )

cor_results %>%
  separate(snp_id, remove = F, into = c("chr","pos"), sep = "_") %>%
  group_by(chr) %>%
  slice_min(p_value, n=1) %>%
  print() %>%
  .$snp_id -> select_snps

adapt.annot %>%
  filter(SNP_id %in% select_snps)

plot_freq_data %>%
  filter(snp_id %in% select_snps) %>%
  mutate(numTim = case_when(
    Time.point == "First" ~ 1,#+(year-2020)* 3,
    Time.point == "Second" ~ 2,#+(year-2020)* 3,
    Time.point == "Third" ~ 3,#+(year-2020)* 3
  )) %>%
  ggplot(aes(
    x=numTim,
    y=AFm,
    ymin=AFlci , ymax=AFuci
  )) + 
  geom_errorbar(width = 0.05) +
  geom_line() +
  geom_point(aes(fill = mCTmin), shape  = 21, size = 3.5) + 
  scale_fill_gradient2(midpoint = 4.4, low = "steelblue",
                       high="firebrick")+
  scale_x_continuous(breaks = 1:3) +
  facet_wrap(~snp_id, nrow = 1) + theme_bw() ->
  ctmin_traj
ggsave(ctmin_traj, file = "ctmin_traj.pdf", w = 8.2, h = 3)


####extra-stuff

afs.id.annot %>%
  group_by(snp_id, Time.point) %>%
  summarize(AFm = ci(value)[1],
            AFlci=ci(value)[2],
            AFuci=ci(value)[3],
            #mT32=mean(T32), 
            mCTmin = ci(ctmin)[1],
            lciCTmin = ci(ctmin)[2],
            uciCTmin = ci(ctmin)[3],
  ) ->
  plot_freq_data


