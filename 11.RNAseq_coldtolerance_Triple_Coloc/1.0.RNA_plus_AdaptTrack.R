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
library(ggrepel)
source("/netfiles/nunezlab/Shared_Resources/Software/baypass_public/utils/baypass_utils.R")

annots <- "/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023_BAYPASS_MG/Final_Datasets/annots.flt.Rdata"
load(annots)
annots.flt$pos = as.numeric(annots.flt$pos)

adapt_dat="/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023_BAYPASS_MG/Final_Datasets/C2_BFct_BFT_df.allAnalyses.Rdata"
adapt=get(load(adapt_dat))

#####
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

####
####

LX28_full <- get(load("/netfiles/nunezlab/D_suzukii_resources/Datasets/Cold_Tolerance_NickTeets/Differential_Expression_Tables/LX28_full.Rdata"))
#LX13_full <- get(load("/netfiles/nunezlab/D_suzukii_resources/Datasets/Cold_Tolerance_NickTeets/Differential_Expression_Tables/LX13_full.Rdata"))

LX28_full %>%
  dplyr::select(gene.name,
                FC =log2FoldChange,
                padj = padj) %>%
  mutate(line = "lx28")->
  lx28.diff

#rbind(lx28.diff, lx13.diff) -> joint.diff.exp
gene_outliers <- adapt.annot %>%
  filter(!is.na(Gene)) %>%
  filter(Gene != "-") %>%
  group_by(gene.name=Gene) %>%
  summarise(
    KYseasonal_outlier = any(KYseasonal_outlier == TRUE, na.rm = TRUE),
    Top1percentC = any(P_C2 > 3, na.rm = TRUE),
    .groups = "drop"
  )

left_join(gene_outliers, lx28.diff) %>%
  filter(!is.na(FC)) ->
  lx28.diff.outliers

save(lx28.diff.outliers, file = "lx28.diff.outliers.Rdata")

