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

annots <- "/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023_BAYPASS_MG/Final_Datasets/annots.flt.Rdata"
load(annots)
annots.flt$pos = as.numeric(annots.flt$pos)

###
inva_dat="/netfiles/nunezlab/D_suzukii_resources/Datasets/2024_Camus_et_al_MolEcol_Invasiveness/res.baypass.pi_xtx_c2.rds"
inva=readRDS(inva_dat)
names(inva)[c(1,2)] = c("chr", "pos")

inva = filter(inva, chr %in% c("chr2L","chr2R","chr3","chrX"))

adapt_dat="/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023_BAYPASS_MG/Final_Datasets/C2_BFct_BFT_df.allAnalyses.Rdata"
adapt=get(load(adapt_dat))

left_join(adapt,
inva[,c("chr", "pos","M_P", "C2_AM", "C2_EU", "C2_WW")]
) -> adapt_inva

left_join(adapt_inva, annots.flt) %>% 
  select(
    chr,
    pos,
    SNP_id,
    freq_in_KY=freqC2,
    freq_in_INV=M_P,
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

load("AdaptiveTrackingSNPs_w_InvasiveSNPs.annot.Rdata")

adapt_inva.clean %<>%
  as.data.frame() %>%
  group_by(chr, pos) %>% arrange(chr, pos)

### Adaptive Tracking Enrichment
png("c2-xi2.png")
seas.c2.ls.xi2=compute.local.scores(data.frame(adapt_inva.clean[,1:2]),
                                  snp.pi=adapt_inva.clean$freq_in_KY,
                                  snp.pvalue = adapt_inva.clean$P_C2 ,
                                  xi=2,manplot = T,main="xi=2; C2")
dev.off()

### Invasive in NAme
png("c2.inva-xi2.png")
inv.c2.ls.xi2=compute.local.scores(data.frame(inva[,1:2]),
                                      snp.pi=inva$M_P,
                                      snp.pvalue = inva$C2_AM ,
                                      xi=2,manplot = T,main="xi=2; AMinva_C2")
dev.off()

### Invasive in EUr
png("c2.inva-xi2.png")
invEU.c2.ls.xi2=compute.local.scores(data.frame(inva[,1:2]),
                                   snp.pi=inva$M_P,
                                   snp.pvalue = inva$C2_EU ,
                                   xi=2,manplot = T,main="xi=2; EUinva_C2")
dev.off()

save(inv.c2.ls.xi2, file = "lindley.InvasionAM.Rdata")
save(invEU.c2.ls.xi2, file = "lindley.InvasionEU.Rdata")
save(seas.c2.ls.xi2, file = "lindley.SeasonalityKY.Rdata")

inv.wins = inv.c2.ls.xi2$significant.windows 
invEU.wins = invEU.c2.ls.xi2$significant.windows 
seas.wins = seas.c2.ls.xi2$significant.windows 

###### find overlaps
setDT(adapt_inva.clean)
setDT(inv.wins)
setDT(invEU.wins)
setDT(seas.wins)

# SNPs become 1 bp intervals
adapt_inva.clean[, `:=`(beg = pos, end = pos)]
setkey(adapt_inva.clean, chr, beg, end)
setkey(inv.wins, chr, beg, end)
setkey(invEU.wins, chr, beg, end)
setkey(seas.wins, chr, beg, end)

### invasiveness
hits.inv <- foverlaps(
  x = adapt_inva.clean,
  y = inv.wins,
  by.x = c("chr", "beg", "end"),
  by.y = c("chr", "beg", "end"),
  type = "within",
  nomatch = NULL
)
inv.hits = hits.inv$SNP_id

hits.invEU <- foverlaps(
  x = adapt_inva.clean,
  y = invEU.wins,
  by.x = c("chr", "beg", "end"),
  by.y = c("chr", "beg", "end"),
  type = "within",
  nomatch = NULL
)
invEU.hits = hits.invEU$SNP_id

### seasonality
hits.seas <- foverlaps(
  x = adapt_inva.clean,
  y = seas.wins,
  by.x = c("chr", "beg", "end"),
  by.y = c("chr", "beg", "end"),
  type = "within",
  nomatch = NULL
)
seasonal.hits = hits.seas$SNP_id

####
adapt_inva.clean %<>%
  mutate(seas_lind = case_when(SNP_id %in% seasonal.hits ~ TRUE,
                               TRUE ~ FALSE) ) 
adapt_inva.clean %<>%
  mutate(inv_lind = case_when(SNP_id %in% inv.hits ~ TRUE,
                               TRUE ~ FALSE) ) 
adapt_inva.clean %<>%
  mutate(invEU_lind = case_when(SNP_id %in% invEU.hits ~ TRUE,
                              TRUE ~ FALSE) ) 

### Save object
adapt_inva.clean %>% 
  select(
    chr,
    pos,
    SNP_id,
    freq_in_KY,
    freq_in_InvasionSet=freq_in_INV,
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
    Extra,
    KYseasonal_outlier=seas_lind,
    NaAme_invasion_outlier=inv_lind,
    Europe_invasion_outlier=invEU_lind
  ) -> adapt_inva.final


save(adapt_inva.final, 
     file = "/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023_BAYPASS_MG/Final_Datasets/C2_BFct_BFT_df.Invsasion_outliers.allAnalyses.Rdata"
)




