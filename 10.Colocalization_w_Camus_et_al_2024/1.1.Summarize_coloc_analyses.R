### Co-localization of invasiveness and adaptive tracking SWD
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
source("/netfiles/nunezlab/Shared_Resources/Software/baypass_public/utils/baypass_utils.R")

####
load("/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023_BAYPASS_MG/Final_Datasets/C2_BFct_BFT_df.Invsasion_outliers.allAnalyses.Rdata")
#adapt_inva.final

####
adapt_inva.final %<>%
  filter(!is.na(KYseasonal_outlier)) %>%
  filter(!is.na(Europe_invasion_outlier)) %>%
  filter(!is.na(NAme_invasion_outlier)) 

all_N = dim(adapt_inva.final)[1]
### define classes

### First question are adaptive tracking outliers enriched in the top %1 of invasiveness outliers
adapt_inva.final %>%
  filter(KYseasonal_outlier == TRUE) %>%
  filter(NAme_invasion_outlier == TRUE) %>%
  dim(.) %>% .[1] -> CLASS_OUTLIER

adapt_inva.final %>%
  filter(KYseasonal_outlier == FALSE) %>%
  filter(NAme_invasion_outlier == TRUE) %>%
  dim(.) %>% .[1] -> noCLASS_OUTLIER

adapt_inva.final %>%
  filter(KYseasonal_outlier == TRUE) %>%
  filter(NAme_invasion_outlier == FALSE) %>%
  dim(.) %>% .[1] -> CLASS_noOUTLIER

adapt_inva.final %>%
  filter(KYseasonal_outlier == FALSE) %>%
  filter(NAme_invasion_outlier == FALSE) %>%
  dim(.) %>% .[1] -> noCLASS_noOUTLIER

CLASS_OUTLIER+noCLASS_OUTLIER+CLASS_noOUTLIER+noCLASS_noOUTLIER -> sanChe1
sanChe1==all_N

Matrix <-
  matrix(c(CLASS_OUTLIER, noCLASS_OUTLIER, 
           CLASS_noOUTLIER, noCLASS_noOUTLIER),
         nrow = 2)

fisher.test(Matrix) -> fet.AdapTrack.InvasiAME

#### Inasiveness and CTmin?
adapt_inva.final %>%
  filter(BF_ctmin >= 15) %>%
  filter(NAme_invasion_outlier == TRUE) %>%
  dim(.) %>% .[1] -> CLASS_OUTLIER

adapt_inva.final %>%
  filter(BF_ctmin < 15) %>%
  filter(NAme_invasion_outlier == TRUE) %>%
  dim(.) %>% .[1] -> noCLASS_OUTLIER

adapt_inva.final %>%
  filter(BF_ctmin >= 15) %>%
  filter(NAme_invasion_outlier == FALSE) %>%
  dim(.) %>% .[1] -> CLASS_noOUTLIER

adapt_inva.final %>%
  filter(BF_ctmin < 15) %>%
  filter(NAme_invasion_outlier == FALSE) %>%
  dim(.) %>% .[1] -> noCLASS_noOUTLIER

CLASS_OUTLIER+noCLASS_OUTLIER+CLASS_noOUTLIER+noCLASS_noOUTLIER -> sanChe1
sanChe1==all_N

Matrix <-
  matrix(c(CLASS_OUTLIER, noCLASS_OUTLIER, 
           CLASS_noOUTLIER, noCLASS_noOUTLIER),
         nrow = 2)

fisher.test(Matrix) -> fet.CTmin.InvasiAME

### invasiveness and t32
adapt_inva.final %>%
  filter(BF_T32 >= 15) %>%
  filter(NAme_invasion_outlier == TRUE) %>%
  dim(.) %>% .[1] -> CLASS_OUTLIER

adapt_inva.final %>%
  filter(BF_T32 < 15) %>%
  filter(NAme_invasion_outlier == TRUE) %>%
  dim(.) %>% .[1] -> noCLASS_OUTLIER

adapt_inva.final %>%
  filter(BF_T32 >= 15) %>%
  filter(NAme_invasion_outlier == FALSE) %>%
  dim(.) %>% .[1] -> CLASS_noOUTLIER

adapt_inva.final %>%
  filter(BF_T32 < 15) %>%
  filter(NAme_invasion_outlier == FALSE) %>%
  dim(.) %>% .[1] -> noCLASS_noOUTLIER

CLASS_OUTLIER+noCLASS_OUTLIER+CLASS_noOUTLIER+noCLASS_noOUTLIER -> sanChe1
sanChe1==all_N

Matrix <-
  matrix(c(CLASS_OUTLIER, noCLASS_OUTLIER, 
           CLASS_noOUTLIER, noCLASS_noOUTLIER),
         nrow = 2)

fisher.test(Matrix) -> fet.T32.InvasiAME

#### Europe?
#### Europe?
#### Europe?
#### Europe?
#### Europe?
adapt_inva.final %>%
  filter(KYseasonal_outlier == TRUE) %>%
  filter(Europe_invasion_outlier == TRUE) %>%
  dim(.) %>% .[1] -> CLASS_OUTLIER

adapt_inva.final %>%
  filter(KYseasonal_outlier == FALSE) %>%
  filter(Europe_invasion_outlier == TRUE) %>%
  dim(.) %>% .[1] -> noCLASS_OUTLIER

adapt_inva.final %>%
  filter(KYseasonal_outlier == TRUE) %>%
  filter(Europe_invasion_outlier == FALSE) %>%
  dim(.) %>% .[1] -> CLASS_noOUTLIER

adapt_inva.final %>%
  filter(KYseasonal_outlier == FALSE) %>%
  filter(Europe_invasion_outlier == FALSE) %>%
  dim(.) %>% .[1] -> noCLASS_noOUTLIER

CLASS_OUTLIER+noCLASS_OUTLIER+CLASS_noOUTLIER+noCLASS_noOUTLIER -> sanChe1
sanChe1==all_N

Matrix <-
  matrix(c(CLASS_OUTLIER, noCLASS_OUTLIER, 
           CLASS_noOUTLIER, noCLASS_noOUTLIER),
         nrow = 2)

fisher.test(Matrix) -> fet.AdapTrack.InvasiEU

#### Inasiveness and CTmin?
adapt_inva.final %>%
  filter(BF_ctmin >= 15) %>%
  filter(Europe_invasion_outlier == TRUE) %>%
  dim(.) %>% .[1] -> CLASS_OUTLIER

adapt_inva.final %>%
  filter(BF_ctmin < 15) %>%
  filter(Europe_invasion_outlier == TRUE) %>%
  dim(.) %>% .[1] -> noCLASS_OUTLIER

adapt_inva.final %>%
  filter(BF_ctmin >= 15) %>%
  filter(Europe_invasion_outlier == FALSE) %>%
  dim(.) %>% .[1] -> CLASS_noOUTLIER

adapt_inva.final %>%
  filter(BF_ctmin < 15) %>%
  filter(Europe_invasion_outlier == FALSE) %>%
  dim(.) %>% .[1] -> noCLASS_noOUTLIER

CLASS_OUTLIER+noCLASS_OUTLIER+CLASS_noOUTLIER+noCLASS_noOUTLIER -> sanChe1
sanChe1==all_N

Matrix <-
  matrix(c(CLASS_OUTLIER, noCLASS_OUTLIER, 
           CLASS_noOUTLIER, noCLASS_noOUTLIER),
         nrow = 2)

fisher.test(Matrix) -> fet.CTmin.InvasiEU

### invasiveness and t32
adapt_inva.final %>%
  filter(BF_T32 >= 15) %>%
  filter(Europe_invasion_outlier == TRUE) %>%
  dim(.) %>% .[1] -> CLASS_OUTLIER

adapt_inva.final %>%
  filter(BF_T32 < 15) %>%
  filter(Europe_invasion_outlier == TRUE) %>%
  dim(.) %>% .[1] -> noCLASS_OUTLIER

adapt_inva.final %>%
  filter(BF_T32 >= 15) %>%
  filter(Europe_invasion_outlier == FALSE) %>%
  dim(.) %>% .[1] -> CLASS_noOUTLIER

adapt_inva.final %>%
  filter(BF_T32 < 15) %>%
  filter(Europe_invasion_outlier == FALSE) %>%
  dim(.) %>% .[1] -> noCLASS_noOUTLIER

CLASS_OUTLIER+noCLASS_OUTLIER+CLASS_noOUTLIER+noCLASS_noOUTLIER -> sanChe1
sanChe1==all_N

Matrix <-
  matrix(c(CLASS_OUTLIER, noCLASS_OUTLIER, 
           CLASS_noOUTLIER, noCLASS_noOUTLIER),
         nrow = 2)

fisher.test(Matrix) -> fet.T32.InvasiEU


#### MErge!!!
#### MErge!!!
#### MErge!!!
#### MErge!!!

All_FETs =
  rbind(
    data.frame(
      Continent = "N.Ame",
      Test = "Seasonality",
      p=fet.AdapTrack.InvasiAME$p.value,
      OR=fet.AdapTrack.InvasiAME$estimate,
      lci=fet.AdapTrack.InvasiAME$conf.int[1],
      uci=fet.AdapTrack.InvasiAME$conf.int[2]
    ),
    data.frame(
      Continent = "N.Ame",
      Test = "CTmin",
      p=fet.CTmin.InvasiAME$p.value,
      OR=fet.CTmin.InvasiAME$estimate,
      lci=fet.CTmin.InvasiAME$conf.int[1],
      uci=fet.CTmin.InvasiAME$conf.int[2]
    ),
    data.frame(
      Continent = "N.Ame",
      Test = "T32",
      p=fet.T32.InvasiAME$p.value,
      OR=fet.T32.InvasiAME$estimate,
      lci=fet.T32.InvasiAME$conf.int[1],
      uci=fet.T32.InvasiAME$conf.int[2]
    ),
    data.frame(
      Continent = "EU",
      Test = "Seasonality",
      p=fet.AdapTrack.InvasiEU$p.value,
      OR=fet.AdapTrack.InvasiEU$estimate,
      lci=fet.AdapTrack.InvasiEU$conf.int[1],
      uci=fet.AdapTrack.InvasiEU$conf.int[2]
    ),
    data.frame(
      Continent = "EU",
      Test = "CTmin",
      p=fet.CTmin.InvasiEU$p.value,
      OR=fet.CTmin.InvasiEU$estimate,
      lci=fet.CTmin.InvasiEU$conf.int[1],
      uci=fet.CTmin.InvasiEU$conf.int[2]
    ),
    data.frame(
      Continent = "EU",
      Test = "T32",
      p=fet.T32.InvasiEU$p.value,
      OR=fet.T32.InvasiEU$estimate,
      lci=fet.T32.InvasiEU$conf.int[1],
      uci=fet.T32.InvasiEU$conf.int[2]
    )
  )

### object ***. All_FETs

##### Part2
##### FINAL PCA
##### FINAL PCA
##### FINAL PCA
##### FINAL PCA
##### FINAL PCA

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

#(seas_lind == TRUE) 
#(inv_lind == TRUE) 

adapt_inva.final %>%
  filter(KYseasonal_outlier == TRUE & NAme_invasion_outlier == TRUE ) ->
  Seas_Invasion_invas_SNPs

snp_info %>% 
  filter(snp_id %in% 
           c(
             Seas_Invasion_invas_SNPs$SNP_id)
  ) %>%
  mutate(rs.id = rownames(.)) %>%
  separate(remove = F, rs.id, into = c("feat", "index"),
           sep = "s") ->
  snp_info.df.finalPCA

samps %>%
  separate(Collection_date, remove = F,
           into = c("y","m","d"),
           sep = "\\-") %>%
  mutate(m = as.numeric(m)) ->
  samps.mod

samps.mod %>% filter(continent %in%
                       c("Asia","Europe","North_America")
) %>% 
  .$sampleId_orig -> select_samps

# 17514 SNPs ...
Trajectory_genes_FINAL <-
  pooldata.subset(
    all_dat,
    pool.index = which(samps$sampleId_orig %in% select_samps),
    snp.index = as.numeric(snp_info.df.finalPCA$index),
    return.snp.idx = TRUE,
    verbose = TRUE
  )

set.seed(1000)
pca_result.final <- randomallele.pca(Trajectory_genes_FINAL, scale = TRUE)
#pca_result.final$perc.var
#[1] 4.7813663 
#[2] 3.8720960

pca_result.final$pop.loadings %>% as.data.frame() %>%
  mutate(sampleId_orig=row.names(.)) %>%
  left_join(samps.mod) -> pc_loadings

### Part 3 --- candidates
### Part 3 --- candidates
### Part 3 --- candidates
### Part 3 --- candidates
### Part 3 --- candidates
### Part 3 --- candidates
### Part 3 --- candidates
### Part 3 --- candidates

### find a candidate
#(seas_lind == TRUE) 
#(inv_lind == TRUE) 

adapt_inva.final %>%
  filter(C2_AM > 4 & P_C2 > 4) %>%
  filter(KYseasonal_outlier == TRUE) %>%
  filter(NAme_invasion_outlier == TRUE) %>%
  mutate(joint_score = C2_AM+P_C2 ) %>%
  arrange(-joint_score)

adapt_inva.final %>%
  filter(C2_AM > 4 & P_C2 > 4) %>%
  filter(KYseasonal_outlier == TRUE) %>%
  filter(NAme_invasion_outlier == TRUE) %>%
  filter(SNP_id %in% c("chrX_16121948")) -> SNP_outlier_info

snp_info %>% 
  filter(snp_id %in% 
           c(
             "chrX_16121948")
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
    #min.cov.per.pool = 10,
    #max.cov.per.pool = 150,
    #min.maf = 0.05,
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


### worldwide
### worldwide
### worldwide
### worldwide
### worldwide

Trajectory_genes_all <-
  pooldata.subset(
    all_dat,
    #pool.index = ky.tag.flt,
    snp.index = as.numeric(snp_info.df$index),
    #min.cov.per.pool = 10,
    #max.cov.per.pool = 150,
    #min.maf = 0.05,
    return.snp.idx = TRUE,
    verbose = TRUE
  )

ref_count_a <- Trajectory_genes_all@refallele.readcount
coverage_a <- Trajectory_genes_all@readcoverage
afs_a <- ref_count_a/coverage_a

afs_a %>%
  as.data.frame %>%
  mutate(snp_id = snp_info.df$snp_id ) ->
  afs.all.id
names(afs.all.id) = c(Trajectory_genes_all@poolnames, "snp_id")

afs.all.id %>%
  reshape2::melt(id = "snp_id", variable.name = "sampleId_orig") %>%
  left_join(samps) %>%
  separate(remove = F, snp_id,
           into = c("chr", "pos"),
           sep = "_") %>%
  group_by(province, Range) %>%
  summarize(AFm = mean(value, na.rm = T),
            latm = mean(lat, na.rm = T),
            longm =  mean(long, na.rm = T),
  ) ->
  plot_freq_data.ALL

#### save object
save(All_FETs,
     pca_result.final,
     SNP_outlier_info,
     pc_loadings,
     afs.id.annot,
     plot_freq_data.ALL,
     file = "plots_invasion_seasonality_figure.Rdata")





