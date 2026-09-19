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

### object All_FETs

### Part 2 --- candidates

### find a candidate
#(seas_lind == TRUE) 
#(inv_lind == TRUE) 

adapt_inva.final %>%
  filter(seas_lind == TRUE) %>%
  filter(inv_lind == TRUE) %>%
  filter(Gene %in% c(108005504))
#chrX_17192750
#chrX_17192811

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

snp_info %>% 
  filter(snp_id %in% 
           c(
             "chrX_17192811", "chrX_17192750",
             "chrX_17192776","chrX_17192781",
             "chrX_17192803")
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

afs.id.annot %>%
  group_by(snp_id, 
           city, 
           year, 
           Time.point,
           #Collection_date
  ) %>%
  summarize(AFm = ci(value, na.rm = T)[1],
            AFlci= ci(value, na.rm = T)[2],
            AFuci= ci(value, na.rm = T)[3],
            mT32=mean(T32), mCTmin = mean(ctmin) ) ->
  plot_freq_data

plot_freq_data %>%
  mutate(numTim = case_when(
    Time.point == "First" ~ 1+(year-2020)* 3,
    Time.point == "Second" ~ 2+(year-2020)* 3,
    Time.point == "Third" ~ 3+(year-2020)* 3)
  ) %>%
  ggplot(aes(
    x=as.numeric(numTim),
    y=AFm,
    group=city,
    #ymin = AFlci, ymax = AFuci,
    #shape = as.character(year)
  )) + 
  geom_line() +
  #geom_errorbar(width = 0.2) +
  geom_point(aes(color = mT32, 
                 #shape=as.character(year)
  ),  size = 3, alpha = 0.95
  ) + 
  scale_color_gradient2(low = "blue", 
                        high = "red", 
                        midpoint = 0.25) +
  scale_shape_manual(values = 21:24) +
  scale_x_continuous(breaks = 1:12) + 
  ylim(0.00,1.0) + facet_grid(city~snp_id) +
  theme_bw()->
  ctmin_traj
ggsave(ctmin_traj, file = "ctmin_traj.pdf", w = 6.5, h = 6.0)




### with ctmin as y axis
afs.id.annot %>%
  group_by(snp_id, 
           city, 
           #year, 
           Time.point,
           #Collection_date
  ) %>%
  summarize(ctminm = ci(ctmin, na.rm = T)[1],
            ctminlci= ci(ctmin, na.rm = T)[2],
            ctminuci= ci(ctmin, na.rm = T)[3],
            mT32=mean(T32), AFm = mean(value) ) ->
  plot_freq_data_ctmin

afs.id.annot %>%
  mutate(numTim = case_when(
    Time.point == "First" ~ 1,
    Time.point == "Second" ~ 2,
    Time.point == "Third" ~ 3)
  ) %>%
  ggplot(aes(
    x=as.factor(numTim),
    y=ctmin,
    #ymin = ctminlci, ymax = ctminuci,
    #shape = city
    #group = year
  )) + 
  #geom_line() +
  #geom_errorbar(width = 0.2) +
  #geom_point(aes(color = value, 
  #               #shape=as.character(year)
  #),  size = 3
  #) +
  geom_boxplot() +
  scale_color_viridis(option = "A"#, midpoint = 0.7
  ) +
  scale_shape_manual(values = 21:24) +
  #scale_x_continuous(breaks = c(1,2,3)) + 
  #ylim(0.00,1.0) + 
  #facet_grid(city~.) +
  theme_bw()->
  ctmin_traj2
ggsave(ctmin_traj2, file = "ctmin_traj2.pdf", w = 3.5, h = 3.0)


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

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

ggplot(data = world) +
  geom_sf(fill= "white") +
  coord_sf(xlim = c(-160.15, 160.99), ylim = c(-55.00, 69.00), expand = FALSE) + 
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.2), 
        panel.background = element_rect(fill = "aliceblue")) +
  geom_point(data = plot_freq_data.ALL, 
             aes(x= longm, y = latm, fill = AFm, shape = Range), 
             size = 3) +
  scale_shape_manual(values = 21:23) +
  scale_fill_viridis(option = "A"#, midpoint = 0.7
  ) -> SAMP.AFs.map
ggsave(SAMP.AFs.map, file = "SAMP.AFs.map.pdf", h = 4, w = 5)

### extra
### extra
### extra
### extra
### extra
### extra
### extra
### extra
### extra
### extra
### extra
### extra
### extra
### extra
### extra

adapt_inva.final %>%
  filter(seas_lind == TRUE ) ->
  seasonal_SNPs

snp_info %>% 
  filter(snp_id %in% 
           c(
             seasonal_SNPs$SNP_id)
  ) %>%
  mutate(rs.id = rownames(.)) %>%
  separate(remove = F, rs.id, into = c("feat", "index"),
           sep = "s") ->
  KY.PCA

samps.mod %>% filter(province %in%
                       c("Kentucky")
) %>% filter(fly_type=="wild") %>%
  .$sampleId_orig -> select_sampsKY


# 17514 SNPs ...
Trajectory_genes_KY <-
  pooldata.subset(
    all_dat,
    pool.index = which(samps$sampleId_orig %in% select_sampsKY),
    snp.index = as.numeric(KY.PCA$index),
    return.snp.idx = TRUE,
    min.cov.per.pool = 10,
    max.cov.per.pool = 150,
    min.maf = 0.11,
    verbose = TRUE
  )

set.seed(1234)
pca_result.ky <- randomallele.pca(Trajectory_genes_KY, scale = TRUE)

pca_result.ky$pop.loadings %>% as.data.frame() %>%
  mutate(sampleId_orig=row.names(.)) %>%
  left_join(samps.mod) -> ky_pc_loadings
#full_join(traits) 
####
ky_pc_loadings %>%
  filter(!is.na(m)) %>%
  mutate(season = case_when(m <= 8 ~ "spring",
                            m > 8 ~ "fall")) %>%
  ggplot(aes(
    x=V1,
    y=V2,
    fill =m
  )) + geom_point(size = 3, shape = 21) + 
  #scale_shape_manual(values = 21:24) +
  scale_fill_gradient2(low = "springgreen", 
                       high = "firebrick", 
                       midpoint = 8.0) +
  theme_bw() ->
  PCA_plot.KY
ggsave(PCA_plot.KY, file = "PCA.KY.FINAL.pdf", w = 4.9, h =4)

ky_pc_loadings %>%
  filter(!is.na(m)) %>%
  mutate(season = case_when(m <= 8 ~ "spring",
                            m > 8 ~ "fall")) %>%
  ggplot(aes(
    x=fct_reorder(as.character(m),m),
    y=V1,
    fill = m
  )) + geom_boxplot() +
  scale_fill_gradient2(low = "springgreen", 
                       high = "firebrick", 
                       midpoint = 8.0) +
  theme_bw() ->
  V1.ky._plot
ggsave(V1.ky._plot, file = "V1.ky._plot.pdf", w = 4.8, h =4)

#### inversion
#### inversion
#### inversion
#### inversion
#### inversion
#### inversion
#### inversion
#### inversion

adapt_inva.final %>%
  filter(chr == "chr2R")  ->
  adapt_inva_Sanja.2R

adapt_inva_Sanja.2R %<>%
  mutate(inv=case_when(pos > 7412539 & pos < 17594956 ~ "inv",
                       TRUE ~ "notInv"))

adapt_inva_Sanja.2R %>%
  filter(BF_T32 >= 15 & inv_lind == TRUE) %>%
  filter(inv == "inv") %>%
  dim(.) %>% .[1] -> CLASS_OUTLIER

adapt_inva_Sanja.2R %>%
  filter(BF_T32 < 15 | inv_lind != TRUE) %>%
  filter(inv == "inv") %>%
  dim(.) %>% .[1] -> noCLASS_OUTLIER

adapt_inva_Sanja.2R %>%
  filter(BF_T32 >= 15 & inv_lind == TRUE) %>%
  filter(inv != "inv") %>%
  dim(.) %>% .[1] -> CLASS_noOUTLIER

adapt_inva_Sanja.2R %>%
  filter(BF_T32 < 15 | inv_lind != TRUE) %>%
  filter(inv != "inv") %>%
  dim(.) %>% .[1] -> noCLASS_noOUTLIER

#CLASS_OUTLIER+noCLASS_OUTLIER+CLASS_noOUTLIER+noCLASS_noOUTLIER -> sanChe1
#sanChe1==all_N

Matrix.Sanja <-
  matrix(c(CLASS_OUTLIER, noCLASS_OUTLIER, 
           CLASS_noOUTLIER, noCLASS_noOUTLIER),
         nrow = 2)

fisher.test(Matrix.Sanja) -> fet.T32.InvasiAME.Matrix.Sanja
fet.T32.InvasiAME.Matrix.Sanja





