### spatio temporal FST analysis
library(tidyverse)
library(magrittr)
library(foreach)
library(vroom)
library(forcats)
library(data.table)
library(gmodels)
require(poolfstat)
library(SeqArray)
library(gdsfmt)
library(SNPRelate)
library(fastglm)
library(rnaturalearth)
library(rnaturalearthdata)
library(factoextra)
library(FactoMineR)
library(viridis)
library(geosphere)
library(patchwork)
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
traits <- fread("/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023_BAYPASS_MG/2024_Joaquin_KY/Means_CT_min_FINAL.csv")
Weather_Data <- fread("/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023_BAYPASS_MG/Ellie_MS/Weather_data/SlicedWeatherData.mgarvin.Feb25.csv")
names(Weather_Data)[1] = "sampleId_orig"
load("/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023_BAYPASS_MG/2024_Joaquin_KY/snp_info.Rdata")

Weather_Data %>% 
  filter(win == "c(7, 15)") %>%
  filter(variable == "T2M") %>%
  filter(stat == "prop. max") -> t32dat

### sample some SNPs...
set.seed(123)

snp_info %>%
  filter(Chromosome != "chrX") %>%
  slice_sample(n = 50000) ->
  snp_info.slice

snp_info %>% 
  filter(snp_id %in% 
           c(
             snp_info.slice$snp_id)
  ) %>%
  mutate(rs.id = rownames(.)) %>%
  separate(remove = F, rs.id, into = c("feat", "index"),
           sep = "s") ->
  snp_info.df.sampled ###<--- SNPs be here!

#####
samps %>%
  separate(Collection_date, remove = F,
           into = c("y","m","d"),
           sep = "\\-") %>%
  mutate(m = as.numeric(m)) ->
  samps.mod
samps.mod %>% filter(province == "Kentucky")

samps.mod %>% filter(continent %in%
                       c("North_America")
) %>% 
  .$sampleId_orig -> select_samps

###
FST_loci <-
  pooldata.subset(
    all_dat,
    pool.index = which(samps$sampleId_orig %in% select_samps),
    snp.index = as.numeric(snp_info.df.sampled$index),
    return.snp.idx = TRUE,
    verbose = TRUE
  )

fst.NAme <- compute.pairwiseFST(FST_loci,  method = "Anova")		   

pairwise_df <- fst.NAme@PairwiseFSTmatrix  %>%
  as.data.frame() %>%
  rownames_to_column("pop1") %>%
  pivot_longer(
    cols = -pop1,
    names_to = "pop2",
    values_to = "FST"
  ) %>%
  filter(pop1 < pop2)

####
samps.mod %>%
  select(
    pop1=sampleId_orig,
    lat1=lat,
    long1=long,
    city1=city, 
    date1=Collection_date,
  ) -> samps1
samps.mod %>%
  select(
    pop2=sampleId_orig,
    lat2=lat,
    long2=long,
    city2=city, 
    date2=Collection_date,
  ) -> samps2
####

pairwise_df %>%
  left_join(samps1) %>%
  left_join(samps2) %>%
  mutate(
    distance_km = distHaversine(
      cbind(long1, lat1),
      cbind(long2, lat2)
    ) / 1000,
    days_apart = abs(
      as.numeric(as.Date(date1) - as.Date(date2))
    ),
    pairwise_df %>%
  left_join(samps1) %>%
  left_join(samps2) %>%
  mutate(
    distance_km = distHaversine(
      cbind(long1, lat1),
      cbind(long2, lat2)
    ) / 1000,
    days_apart = abs(
      as.numeric(as.Date(date1) - as.Date(date2))
    ),
    years_apart = abs(
      as.numeric(format(as.Date(date1), "%Y")) -
        as.numeric(format(as.Date(date2), "%Y"))
    )
  )) -> pairwise_df.annot

####
pairwise_df.annot %>% 
  filter(
    city1 %in% c("Lexington", "Berea"),
    city2 %in% c("Lexington", "Berea")
  ) %>%
  filter(city1 == city2) %>%
  filter(years_apart %in% 0:1) %>%
  group_by(city1, years_apart) %>%
  summarise(medFST = median((FST / (1 - FST))))


#### plots -- time
pairwise_df.annot %>% 
  filter(
    city1 %in% c("Lexington", "Berea"),
    city2 %in% c("Lexington", "Berea")
  ) %>%
  filter(city1 == city2) %>%
  filter(years_apart %in% 0:1) %>%
  ggplot(
    aes(
      fill=as.factor(years_apart),
      y=(FST / (1 - FST)),
      x = city1
    )
  ) + geom_boxplot() +
  theme_bw() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.8, 0.8)
  ) + ylim(0,0.22) ->
  time_fst
#ggsave(time_fst, file = "time_fst.pdf")

###
pairwise_df.annot %>% 
  filter(
    (city1 == "Lexington" & city2 == "Charlottesville") |
    (city1 == "Berea" & city2 == "Charlottesville")
  ) %>%
  filter(city1 != city2) %>%
  mutate(comp = "VA vs KY") %>%
  ggplot(
    aes(
      y=(FST / (1 - FST)),
      x = paste(city1,city2)
    )
  ) + geom_boxplot() +
  theme_bw() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.8, 0.8)
  ) + ylim(0,0.22) ->
  spac_fst_short

#### plots -- space
pairwise_df.annot %>% 
  filter(city1 != city2) %>%
  group_by(paste(city1, city2)) %>%
  slice_head() %>%
  ggplot(
    aes(
      x=(distance_km),
      y=(FST / (1 - FST))
    )
  ) + 
  geom_hline(yintercept = 0.0243) +
  geom_hline(yintercept = 0.0424) +
  geom_smooth(method = "lm", se = F) +
  geom_point(size = 3, shape = 21, fill = "grey") +
  theme_bw() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.8, 0.8)
  ) + ylim(0,0.22) ->
  spac_fst_long


ggsave(time_fst + spac_fst_short + spac_fst_long +
       plot_layout(widths = c(0.6, 0.6, 1.6)),
       file = "spac_fst.pdf",
       w = 11.5, h = 3)

