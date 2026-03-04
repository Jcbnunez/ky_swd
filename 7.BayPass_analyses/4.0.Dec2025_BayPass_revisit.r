## 4. Revision analysis for BayPass
# -- module load Rgeospatial
library(data.table)
library(tidyverse)
library(foreach)
library(forcats)
library(magrittr)

## load data... this has now been centralized to Nunez lab folder.
datBay <- "/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023_BAYPASS_MG/2024_Joaquin_KY/C2_BFct_BFT_df.allAnalyses.Rdata"

load(datBay)
## loads --> C2_BFct_BFT_df

### load annotations
annots <- "/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023_BAYPASS_MG/2024_Joaquin_KY/annots.flt.Rdata"
load(annots)
annots.flt$pos = as.numeric(annots.flt$pos)
C2_BFct_BFT_df %>%
  left_join(annots.flt) -> C2_BFct_BFT_df



## filter to BF > 15
C2_BFct_BFT_df %>% 
  filter(BF_ctmin >= 15) -> CTMinSNPs
quantile(C2_BFct_BFT_df$BF_ctmin, 0.01)

write.table(CTMinSNPs, 
            file = "CTMinSNPs.txt", 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")


## filter to BF > 15 for T32
C2_BFct_BFT_df %>% 
  filter(BF_T32 >= 15) -> T32_SNPs
quantile(C2_BFct_BFT_df$BF_T32, 0.01)

write.table(T32_SNPs, 
            file = "T32_SNPs.txt", 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")


###
C2_BFct_BFT_df %>% 
  filter(BF_T32 >= 15 & BF_ctmin >= 15) ->
  Both_outliers

write.table(Both_outliers, 
            file = "Both_outliers.txt", 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

### 12 SNPs are in both datasets

### plot dataset
C2_BFct_BFT_df %>%
  filter(BF_T32 > 0) %>%
  .[sample(dim(.)[1],10000),] ->
  setT32
C2_BFct_BFT_df %>%
  filter(BF_ctmin > 0) %>%
  .[sample(dim(.)[1],10000),] ->
  setCtmin

C2_BFct_BFT_df %>%
  filter(BF_T32 >= 15 & BF_ctmin >= 15) ->
  trueT32CT

full_join(setT32, trueT32CT) -> random_Set

ggplot() +
  geom_point(
    data=random_Set,
    aes(
      x=BF_T32,
      y=BF_ctmin,
    ),
    alpha=0.8, size=0.5) +
  geom_point(
    data=trueT32CT,
    aes(
      x=BF_T32,
      y=BF_ctmin,
    ),
    alpha=0.8, size=1.2, color= "red") +
  geom_hline(yintercept = 15, linetype = "dashed") +
  geom_vline(xintercept = 15, linetype = "dashed") +
  geom_text(
    data = filter(C2_BFct_BFT_df, BF_T32 >= 15 & BF_ctmin >= 15 ), 
    aes(label = Gene, x=BF_T32, y=BF_ctmin )) +
  theme_bw() ->
  BF_plots
ggsave(BF_plots, file = "BF_plots.pdf", w = 4, h = 4)

#### Enrichment analyses. 
#####
C2_BFct_BFT_df %>%
  filter(BF_T32 >= 15) %>%
  filter(BF_ctmin >= 15) -> tycy
C2_BFct_BFT_df %>%
  filter(BF_T32 >= 15) %>%
  filter(BF_ctmin < 15) -> tycn
C2_BFct_BFT_df %>%
  filter(BF_T32 < 15) %>%
  filter(BF_ctmin >= 15) -> tncy
C2_BFct_BFT_df %>%
  filter(BF_T32 < 15) %>%
  filter(BF_ctmin < 15) -> tncn

###
FETdat.tc <-
  matrix(c(dim(tycy)[1], dim(tycn)[1], 
           dim(tncy)[1], dim(tncn)[1]),
         nrow = 2)

fisher.test(FETdat.tc, alternative = "greater")


#### Look at C2
quantile(C2_BFct_BFT_df$P_C2)
C2_BFct_BFT_df %>% 
  filter(BF_T32 >= 15 & BF_ctmin >= 15) %>%
  filter(P_C2 >= 2.071122)
  

