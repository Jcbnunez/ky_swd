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
