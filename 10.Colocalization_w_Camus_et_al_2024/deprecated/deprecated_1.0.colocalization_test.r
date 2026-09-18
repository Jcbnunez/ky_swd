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

annots <- "/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023_BAYPASS_MG/2024_Joaquin_KY/annots.flt.Rdata"
load(annots)
annots.flt$pos = as.numeric(annots.flt$pos)

###
inva_dat="/netfiles/nunezlab/D_suzukii_resources/Datasets/2024_Camus_et_al_MolEcol_Invasiveness/res.baypass.pi_xtx_c2.rds"
inva=readRDS(inva_dat)
names(inva)[c(1,2)] = c("chr", "pos")

adapt_dat="/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023_BAYPASS_MG/2024_Joaquin_KY/C2_BFct_BFT_df.allAnalyses.Rdata"
adapt=get(load(adapt_dat))

left_join(adapt,
inva[,c("chr", "pos", "C2_AM", "C2_EU", "C2_WW")]
) -> adapt_inva

adapt_inva %>%
  select(
    chr,
    pos,
    SNP_id,
    freq_in_KY=freqC2,
    P_C2,
    BF_T32,
    BF_ctmin,
    C2_AM,
    C2_EU,
    C2_WW,
    Allele,
    Gene,
    Feature,
    Feature_type,
    Consequence,
    cDNA_position,
    CDS_position,
    Protein_position,
    Amino_acids,
    Codons,
    Existing_variation,
    Extra
    ) -> adapt_inva.clean

save(adapt_inva.clean, 
     file = "AdaptiveTrackingSNPs_w_InvasiveSNPs.annot.Rdata")

####
adapt_inva %<>%
  filter(!is.na(P_C2)) %>%
  filter(!is.na(C2_AM)) 

adapt_inva %<>%
  left_join(annots.flt) 

  
all_N = dim(adapt_inva)[1]
### define classes


### First question are adaptive tracking outliers enriched in the top %1 of invasiveness outliers
adapt_inva %>%
  filter(P_C2 >= 2) %>%
  filter(C2_AM >= 2) %>%
  dim(.) %>% .[1] -> CLASS_OUTLIER

adapt_inva %>%
  filter(P_C2 < 2) %>%
  filter(C2_AM >= 2) %>%
  dim(.) %>% .[1] -> noCLASS_OUTLIER

adapt_inva %>%
  filter(P_C2 >= 2) %>%
  filter(C2_AM < 2) %>%
  dim(.) %>% .[1] -> CLASS_noOUTLIER

adapt_inva %>%
  filter(P_C2 < 2) %>%
  filter(C2_AM < 2) %>%
  dim(.) %>% .[1] -> noCLASS_noOUTLIER

CLASS_OUTLIER+noCLASS_OUTLIER+CLASS_noOUTLIER+noCLASS_noOUTLIER -> sanChe1
sanChe1==all_N

Matrix <-
  matrix(c(CLASS_OUTLIER, noCLASS_OUTLIER, 
           CLASS_noOUTLIER, noCLASS_noOUTLIER),
         nrow = 2)

fisher.test(Matrix) -> fet.AdapTrack.InvasiAME

#### Inasiveness and CTmin?
adapt_inva %>%
  filter(BF_ctmin >= 15) %>%
  filter(C2_AM >= 2) %>%
  dim(.) %>% .[1] -> CLASS_OUTLIER

adapt_inva %>%
  filter(BF_ctmin < 15) %>%
  filter(C2_AM >= 2) %>%
  dim(.) %>% .[1] -> noCLASS_OUTLIER

adapt_inva %>%
  filter(BF_ctmin >= 15) %>%
  filter(C2_AM < 2) %>%
  dim(.) %>% .[1] -> CLASS_noOUTLIER

adapt_inva %>%
  filter(BF_ctmin < 15) %>%
  filter(C2_AM < 2) %>%
  dim(.) %>% .[1] -> noCLASS_noOUTLIER

CLASS_OUTLIER+noCLASS_OUTLIER+CLASS_noOUTLIER+noCLASS_noOUTLIER -> sanChe1
sanChe1==all_N

Matrix <-
  matrix(c(CLASS_OUTLIER, noCLASS_OUTLIER, 
           CLASS_noOUTLIER, noCLASS_noOUTLIER),
         nrow = 2)

fisher.test(Matrix) -> fet.CTmin.InvasiAME

### invasiveness and t32
adapt_inva %>%
  filter(BF_T32 >= 15) %>%
  filter(C2_AM >= 2) %>%
  dim(.) %>% .[1] -> CLASS_OUTLIER

adapt_inva %>%
  filter(BF_T32 < 15) %>%
  filter(C2_AM >= 2) %>%
  dim(.) %>% .[1] -> noCLASS_OUTLIER

adapt_inva %>%
  filter(BF_T32 >= 15) %>%
  filter(C2_AM < 2) %>%
  dim(.) %>% .[1] -> CLASS_noOUTLIER

adapt_inva %>%
  filter(BF_T32 < 15) %>%
  filter(C2_AM < 2) %>%
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

adapt_inva %>%
  filter(P_C2 >= 2) %>%
  filter(C2_EU >= 2) %>%
  dim(.) %>% .[1] -> CLASS_OUTLIER

adapt_inva %>%
  filter(P_C2 < 2) %>%
  filter(C2_EU >= 2) %>%
  dim(.) %>% .[1] -> noCLASS_OUTLIER

adapt_inva %>%
  filter(P_C2 >= 2) %>%
  filter(C2_EU < 2) %>%
  dim(.) %>% .[1] -> CLASS_noOUTLIER

adapt_inva %>%
  filter(P_C2 < 2) %>%
  filter(C2_EU < 2) %>%
  dim(.) %>% .[1] -> noCLASS_noOUTLIER

CLASS_OUTLIER+noCLASS_OUTLIER+CLASS_noOUTLIER+noCLASS_noOUTLIER -> sanChe1
sanChe1==all_N

Matrix <-
  matrix(c(CLASS_OUTLIER, noCLASS_OUTLIER, 
           CLASS_noOUTLIER, noCLASS_noOUTLIER),
         nrow = 2)

fisher.test(Matrix) -> fet.AdapTrack.InvasiEU

#### Inasiveness and CTmin?
adapt_inva %>%
  filter(BF_ctmin >= 15) %>%
  filter(C2_EU >= 2) %>%
  dim(.) %>% .[1] -> CLASS_OUTLIER

adapt_inva %>%
  filter(BF_ctmin < 15) %>%
  filter(C2_EU >= 2) %>%
  dim(.) %>% .[1] -> noCLASS_OUTLIER

adapt_inva %>%
  filter(BF_ctmin >= 15) %>%
  filter(C2_EU < 2) %>%
  dim(.) %>% .[1] -> CLASS_noOUTLIER

adapt_inva %>%
  filter(BF_ctmin < 15) %>%
  filter(C2_EU < 2) %>%
  dim(.) %>% .[1] -> noCLASS_noOUTLIER

CLASS_OUTLIER+noCLASS_OUTLIER+CLASS_noOUTLIER+noCLASS_noOUTLIER -> sanChe1
sanChe1==all_N

Matrix <-
  matrix(c(CLASS_OUTLIER, noCLASS_OUTLIER, 
           CLASS_noOUTLIER, noCLASS_noOUTLIER),
         nrow = 2)

fisher.test(Matrix) -> fet.CTmin.InvasiEU

### invasiveness and t32
adapt_inva %>%
  filter(BF_T32 >= 15) %>%
  filter(C2_EU >= 2) %>%
  dim(.) %>% .[1] -> CLASS_OUTLIER

adapt_inva %>%
  filter(BF_T32 < 15) %>%
  filter(C2_EU >= 2) %>%
  dim(.) %>% .[1] -> noCLASS_OUTLIER

adapt_inva %>%
  filter(BF_T32 >= 15) %>%
  filter(C2_EU < 2) %>%
  dim(.) %>% .[1] -> CLASS_noOUTLIER

adapt_inva %>%
  filter(BF_T32 < 15) %>%
  filter(C2_EU < 2) %>%
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
  Continent = "Europe",
  Test = "Seasonality",
  p=fet.AdapTrack.InvasiEU$p.value,
  OR=fet.AdapTrack.InvasiEU$estimate,
  lci=fet.AdapTrack.InvasiEU$conf.int[1],
  uci=fet.AdapTrack.InvasiEU$conf.int[2]
),
data.frame(
  Continent = "Europe",
  Test = "CTmin",
  p=fet.CTmin.InvasiEU$p.value,
  OR=fet.CTmin.InvasiEU$estimate,
  lci=fet.CTmin.InvasiEU$conf.int[1],
  uci=fet.CTmin.InvasiEU$conf.int[2]
),
data.frame(
  Continent = "Europe",
  Test = "T32",
  p=fet.T32.InvasiEU$p.value,
  OR=fet.T32.InvasiEU$estimate,
  lci=fet.T32.InvasiEU$conf.int[1],
  uci=fet.T32.InvasiEU$conf.int[2]
)
)

#> All_FETs
#Continent        Test            p       OR       lci      uci
#odds ratio      N.Ame Seasonality 2.777254e-14 1.083053 1.0610646 1.105438
#odds ratio1     N.Ame       CTmin 6.305897e-01 0.952952 0.7961653 1.135557
#odds ratio2     N.Ame         T32 1.388704e-02 1.156967 1.0289759 1.298645
#odds ratio3    Europe Seasonality 1.655324e-08 1.062023 1.0400914 1.084364
#odds ratio4    Europe       CTmin 3.759236e-01 1.079310 0.9052409 1.281597
#odds ratio5    Europe         T32 8.180982e-08 1.367827 1.2205640 1.530612

All_FETs %>%
  ggplot(
    aes(
      x=Test,
      y=log2(OR),
      ymin=log2(lci),
      ymax=log2(uci),
      shape=Continent
    )
  ) + 
  geom_hline(yintercept = 0, linetype= "dashed") +
  geom_errorbar(width = 0.5, position=position_dodge(width=0.5)) +
  geom_point(size = 4, fill = "grey", position=position_dodge(width=0.5)) +
  theme_classic() + scale_shape_manual(values = 21:22) ->
  ORFET_plot

ggsave(ORFET_plot, file = "ORFET_plot.pdf", w=4.3, h = 4)

### find a candidate
adapt_inva %>%
  filter(BF_T32 >= 15) %>%
  filter(BF_ctmin >= 15) %>%
  filter(P_C2 >= 2) %>%
  filter(C2_AM >= 2) ->
  top_T32_invas

adapt_inva %>%
  filter(BF_T32 >= 15) %>%
  filter(BF_ctmin >= 15) %>%
  filter(P_C2 >= 2) %>%
  filter(C2_AM >= 2) %>% dim
# 1 .. just one

adapt_inva %>%
  filter(BF_T32 >= 15) %>%
  filter(P_C2 >= 2) %>%
  filter(C2_AM >= 2) %>% dim
# 50 .. in total

adapt_inva %>%
  filter(BF_T32 >= 15) %>%
  filter(C2_AM >= 2) %>% dim
# 394 .. in total

adapt_inva %>%
  filter(P_C2 >= 2) %>%
  filter(C2_AM >= 2) %>% dim
# 12358 .. in total

adapt_inva %>%
  filter(BF_ctmin >= 15) %>%
  filter(C2_AM >= 2) %>% dim
# 164 ... in total  
  
### Top hit: chr2R 9965890
### Top hit: chr2R 9965890
### Top hit: chr2R 9965890
### Top hit: chr2R 9965890
### Top hit: chr2R 9965890

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
             "chr2R_9965890")
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
           #city, 
           #year, 
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
         Time.point == "First" ~ 1,
         Time.point == "Second" ~ 2,
         Time.point == "Third" ~ 3)
         ) %>%
ggplot(aes(
  x=numTim,
  y=AFm,
  ymin = AFlci, ymax = AFuci,
  #group = year
)) + 
  geom_line() +
  geom_errorbar(width = 0.2) +
  geom_point(aes(color = mCTmin, 
                    #shape=as.character(year)
                    ),  size = 3, alpha = 0.95
                ) + 
  scale_color_gradient2(low = "blue", 
                        high = "red", 
                        midpoint = 4.4) +
  scale_shape_manual(values = 21:24) +
  scale_x_continuous(breaks = c(1,2,3)) + 
  ylim(0.00,1.0) + 
  #facet_grid(city~.) +
  theme_bw()->
  ctmin_traj
ggsave(ctmin_traj, file = "ctmin_traj.pdf", w = 3.5, h = 3.0)

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

##### FINAL PCA
##### FINAL PCA
##### FINAL PCA
##### FINAL PCA
##### FINAL PCA

adapt_inva %>%
  #filter(BF_T32 >= 15) %>%
  filter(P_C2 >= 2 ) %>%
  filter(C2_EU >= 2 | C2_AM >= 2) ->
  top_T32_invas_any
snp_info %>% 
  filter(snp_id %in% 
           c(
             top_T32_invas_any$SNP_id)
  ) %>%
  mutate(rs.id = rownames(.)) %>%
  separate(remove = F, rs.id, into = c("feat", "index"),
           sep = "s") ->
  snp_info.df.finalPCA

#samps %>% filter(!is.na(Range)) %>% .$sampleId_orig -> select_samps

# 17514 SNPs ...
Trajectory_genes_FINAL <-
  pooldata.subset(
    all_dat,
    #pool.index = which(samps$sampleId_orig %in% select_samps),
    snp.index = as.numeric(snp_info.df.finalPCA$index),
    #min.cov.per.pool = 10,
    #max.cov.per.pool = 150,
    #min.maf = 0.05,
    return.snp.idx = TRUE,
    verbose = TRUE
  )

pca_result.final <- randomallele.pca(Trajectory_genes_FINAL, scale = TRUE)
#pca_result.final$perc.var
#[1] 4.7813663 
#[2] 3.8720960

samps %>%
  separate(Collection_date, remove = F,
                           into = c("y","m","d"),
                           sep = "\\-") %>%
  mutate(m = as.numeric(m)) ->
  samps.mod

pca_result.final$pop.loadings %>% as.data.frame() %>%
  mutate(sampleId_orig=row.names(.)) %>%
  left_join(samps.mod) %>%
  full_join(traits) -> pc_loadings
####
pc_loadings %>%
  filter(!is.na(m)) %>%
  ggplot(aes(
    x=V1,
    y=V2,
    shape=Range, fill =m
  )) + geom_point(size = 3) +
  scale_shape_manual(values = 21:23) +
  scale_fill_gradient2(low = "springgreen", 
                        high = "firebrick", 
                        midpoint = 8.0) +
  theme_bw()->
  PCA_plot
ggsave(PCA_plot, file = "PCA.FINAL.pdf", w = 4.8, h =4)

