### RNa-seq analysis
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

annots <- "/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023_BAYPASS_MG/2024_Joaquin_KY/annots.flt.Rdata"
load(annots)
annots.flt$pos = as.numeric(annots.flt$pos)

adapt_dat="/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023_BAYPASS_MG/Final_Datasets/C2_BFct_BFT_df.allAnalyses.Rdata"
adapt=get(load(adapt_dat))

left_join(adapt, annots.flt) %>% 
  select(
    chr,
    pos,
    SNP_id,
    freq_in_KY=freqC2,
    P_C2,
    BF_T32,
    BF_ctmin,
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
  ) %>%
  as.data.frame() -> adapt.annot


seas.c2.ls.xi2=compute.local.scores(data.frame(adapt.annot[,1:2]),
                                    snp.pi=adapt.annot$freq_in_KY,
                                    snp.pvalue = adapt.annot$P_C2 ,
                                    xi=2,manplot = T,main="xi=2; C2")
seas.wins = seas.c2.ls.xi2$significant.windows 
###### find overlaps
setDT(adapt.annot)
setDT(seas.wins)
adapt.annot[, `:=`(beg = pos, end = pos)]
setkey(adapt.annot, chr, beg, end)
setkey(seas.wins, chr, beg, end)

hits.seas <- foverlaps(
  x = adapt.annot,
  y = seas.wins,
  by.x = c("chr", "beg", "end"),
  by.y = c("chr", "beg", "end"),
  type = "within",
  nomatch = NULL
)
seasonal.hits = hits.seas$SNP_id

adapt.annot %<>%
  mutate(seas_lind = case_when(SNP_id %in% seasonal.hits ~ TRUE,
                               TRUE ~ FALSE) ) 

### Run enrichment

FETs.chrs <- foreach(chr.i=c("chr2L","chr2R", "chr3", "chr4" , "chrX"),
                     .combine="rbind", 
                     .errorhandling="remove")%do%{
                       
                       win.tmp <- adapt.annot %>%
                         filter(chr == chr.i)
                       
                       all_N = dim(win.tmp)[1]
                       
                       win.tmp %>%
                         filter(seas_lind == TRUE) %>%
                         filter(BF_T32 >= 15) %>%
                         dim(.) %>% .[1] -> CLASS_OUTLIER
                       
                       win.tmp %>%
                         filter(seas_lind == FALSE) %>%
                         filter(BF_T32 >= 15) %>%
                         dim(.) %>% .[1] -> noCLASS_OUTLIER
                       
                       win.tmp %>%
                         filter(seas_lind == TRUE) %>%
                         filter(BF_T32 < 15) %>%
                         dim(.) %>% .[1] -> CLASS_noOUTLIER
                       
                       win.tmp %>%
                         filter(seas_lind == FALSE) %>%
                         filter(BF_T32 < 15) %>%
                         dim(.) %>% .[1] -> noCLASS_noOUTLIER
                       
                       CLASS_OUTLIER+noCLASS_OUTLIER+CLASS_noOUTLIER+noCLASS_noOUTLIER -> sanChe1
                       sanChe1==all_N
                       
                       Matrix <-
                         matrix(c(CLASS_OUTLIER, noCLASS_OUTLIER, 
                                  CLASS_noOUTLIER, noCLASS_noOUTLIER),
                                nrow = 2)
                       
                       fisher.test(Matrix) -> fet.T32
                       
                       ####
                       win.tmp %>%
                         filter(seas_lind == TRUE) %>%
                         filter(BF_ctmin >= 15) %>%
                         dim(.) %>% .[1] -> CLASS_OUTLIER
                       
                       win.tmp %>%
                         filter(seas_lind == FALSE) %>%
                         filter(BF_ctmin >= 15) %>%
                         dim(.) %>% .[1] -> noCLASS_OUTLIER
                       
                       win.tmp %>%
                         filter(seas_lind == TRUE) %>%
                         filter(BF_ctmin < 15) %>%
                         dim(.) %>% .[1] -> CLASS_noOUTLIER
                       
                       win.tmp %>%
                         filter(seas_lind == FALSE) %>%
                         filter(BF_ctmin < 15) %>%
                         dim(.) %>% .[1] -> noCLASS_noOUTLIER
                       
                       CLASS_OUTLIER+noCLASS_OUTLIER+CLASS_noOUTLIER+noCLASS_noOUTLIER -> sanChe1
                       sanChe1==all_N
                       
                       Matrix <-
                         matrix(c(CLASS_OUTLIER, noCLASS_OUTLIER, 
                                  CLASS_noOUTLIER, noCLASS_noOUTLIER),
                                nrow = 2)
                       
                       fisher.test(Matrix) -> fet.ct
                       ###
                       rbind(
                       data.frame(
                         chr = chr.i,
                         var = "t32",
                         OR = fet.T32$estimate,
                         LCI = fet.T32$conf.int[1],
                         UCI = fet.T32$conf.int[2],
                         P = fet.T32$p.value
                       ),
                       data.frame(
                         chr = chr.i,
                         var = "ctmin",
                         OR = fet.ct$estimate,
                         LCI = fet.ct$conf.int[1],
                         UCI = fet.ct$conf.int[2],
                         P = fet.ct$p.value
                       ))
                     }

FETs.chrs %>%
  filter(chr != "chr4") %>%
  ggplot(aes(
    x=chr,
    y=log2(OR), 
    ymin=log2(LCI), 
    ymax=log2(UCI),
    fill = var
  )) + 
  geom_errorbar(
    width = 0.5,
    position=position_dodge(width=0.5)) +
  geom_point(size = 2.2, 
             shape=21,
             position=position_dodge(width=0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() + scale_fill_manual(values = c("steelblue","gold"))->
  FET_chrsplot
ggsave(FET_chrsplot, file = "FET_chrsplot.pdf",
       h=2.5, w = 3.5)


