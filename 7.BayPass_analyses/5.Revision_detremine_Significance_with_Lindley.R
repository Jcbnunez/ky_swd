#### Assess significance with Lindley Scores
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

## window properties
sum(seas.wins$size)
mean(seas.wins$size)
adapt.annot$KYseasonal_outlier %>% table
####
library(dplyr)
library(purrr)
library(tibble)

consequences = c("3_prime_UTR_variant","5_prime_UTR_variant","downstream_gene_variant",
                 "intergenic_variant","intron_variant","missense_variant","synonymous_variant",
                 "upstream_gene_variant"
                 )

fet_results <- map_dfr(consequences, function(x) {
  
  # Construct 2 x 2 table
  tab <- matrix(c(
    sum(adapt.annot$Consequence == x & adapt.annot$KYseasonal_outlier == TRUE,  na.rm = TRUE),
    sum(adapt.annot$Consequence == x & adapt.annot$KYseasonal_outlier == FALSE, na.rm = TRUE),
    sum(adapt.annot$Consequence != x & adapt.annot$KYseasonal_outlier == TRUE,  na.rm = TRUE),
    sum(adapt.annot$Consequence != x & adapt.annot$KYseasonal_outlier == FALSE, na.rm = TRUE)
  ),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    Consequence = c(x, "Other"),
    KYseasonal_outlier = c("TRUE", "FALSE")
  ))
  
  # Fisher exact test
  fet <- fisher.test(tab)
  
  tibble(
    Consequence = x,
    outlier_x = tab[1, 1],
    nonoutlier_x = tab[1, 2],
    outlier_other = tab[2, 1],
    nonoutlier_other = tab[2, 2],
    odds_ratio = unname(fet$estimate),
    lci=fet$conf.int[1], uci=fet$conf.int[2],
    p_value = fet$p.value
  )
})

fet_results %>% 
  filter(odds_ratio > 1 & p_value < 0.05)

fet_results %>%
  ggplot(aes(
    x=fct_reorder(Consequence, odds_ratio),
    y=log2(odds_ratio),
    ymin=log2(lci),
    ymax=log2(uci)
  )) +
  geom_hline(yintercept = 0) +
  geom_errorbar(width = 0.5, position=position_dodge(width=0.5)) +
  geom_point(size = 4, fill = "grey", position=position_dodge(width=0.5)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) ->
  ORFET_conseq_plot

ggsave(ORFET_conseq_plot, file = "ORFET_conseq_plot.pdf", w=4.3, h = 4)

 
###PCA
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

samps %>%
  separate(Collection_date, remove = F,
           into = c("y","m","d"),
           sep = "\\-") %>%
  mutate(m = as.numeric(m)) ->
  samps.mod

Weather_Data %>% 
  filter(win == "c(7, 15)") %>%
  filter(variable == "T2M") %>%
  filter(stat == "prop. max") -> t32dat

Weather_Data %>% 
  filter(win == "c(0, 7)") %>%
  filter(variable == "T2M") %>%
  filter(stat == "mean") -> tmeandat


adapt.annot %>%
  filter(KYseasonal_outlier == TRUE ) ->
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
#pca_result.ky$perc.var
#[1] 7.390181 6.487919 6.208532 6.000027 5.889121 5.815960 5.785795 5.526194
#[9] 5.491128 5.389777 5.307499 5.256612 5.208455 5.065626 4.999171 4.874284
#[17] 4.754431 4.549288
#
pca_result.ky$pop.loadings %>% as.data.frame() %>%
  mutate(sampleId_orig=row.names(.)) %>%
  left_join(samps.mod) -> ky_pc_loadings
#full_join(traits) 
####
ky_pc_loadings %>%
  left_join(tmeandat) %>% 
  left_join(traits) ->
  PCdata_ecovars

###
mod1 = lm(V1 ~ `Juliean Date` + `Fruit type` + city + y, data=PCdata_ecovars)
anova(mod1)

####
adapt.annot %>% 
  filter(KYseasonal_outlier == TRUE) %>%
  filter(BF_ctmin > 15) 

save(PCdata_ecovars, file = "PCdata_ecovars.Rdata")


