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

#LX13_full %>%
#  dplyr::select(gene.name,
#                FC =log2FoldChange,
#                padj = padj)%>%
#  mutate(line = "lx13") ->
#  lx13.diff

#rbind(lx28.diff, lx13.diff) %>%
#  filter(gene.name == 108010308)
#rbind(lx28.diff, lx13.diff) %>%
#  filter(gene.name == 108009231)


#rbind(lx28.diff, lx13.diff) -> joint.diff.exp
adapt.annot %>%
  mutate(
    category = case_when(
      P_C2 >= 2 & BF_T32 >= 15 & BF_ctmin >= 15 ~ "all3",
      TRUE ~ NA
    )
  ) %>% 
  group_by(Gene) %>%
  count(category) ->
  genes_all3

genes_all3$category %>% table

adapt.annot %>%
  mutate(
    category = case_when(
      BF_ctmin >= 15 ~ "ctmin",
      TRUE ~ NA
    )
  ) %>% 
  group_by(Gene) %>%
  count(category) ->
  genes_ctmin

genes_ctmin$category %>% table


adapt.annot %>%
  mutate(
    category = case_when(
      P_C2 >= 2 ~ "C2",
      TRUE ~ NA
    )
  ) %>% 
  group_by(Gene, category) %>%
  summarise(N = n()) %>%
  filter(category == "C2")->
  genes_C2

genes_C2$category %>% table

####
adapt.annot %>%
  mutate(
    category = case_when(
      P_C2 >= 2 & BF_T32 >= 15 & BF_ctmin >= 15 ~ "all3",
      P_C2 >= 2 & BF_T32 >= 15 ~ "C2_t32",
      P_C2 >= 2 & BF_ctmin >= 15 ~ "C2_ctmin",
      BF_ctmin >= 15 ~ "ctmin",
      TRUE ~ NA
    )
  ) %>% 
  group_by(Gene) %>%
  count(category) ->
  genes_t32_c2

genes_t32_c2 %>%
  count(category)

####
left_join(lx28.diff, select(genes_C2,
                            gene.name = Gene, category),
           by = "gene.name") %>%
  group_by(category, padj < 0.01) ->
  lx28.diff.annot

save(lx28.diff.annot, file = "lx28.diff.annot.Rdata")

lx28.diff.annot %>% group_by(category, padj < 0.01) %>%
  summarise(N = n())

#Matrix <-
#  matrix(c(CLASS_OUTLIER, noCLASS_OUTLIER, 
#           CLASS_noOUTLIER, noCLASS_noOUTLIER),
#         nrow = 2)
Matrix <-
  matrix(c(572, 303, 
           6205, 4533),
         nrow = 2)

fisher.test(Matrix) -> fet.gene_expression
fet.gene_expression

lx28.diff %>%
  filter(gene.name == 108010308)

lx28.diff %>%
  filter(gene.name %in% filter(genes_t32_c2, category == "all3")$Gene )

left_join(lx28.diff, select(genes_C2,
                            gene.name = Gene, category)) -> 
  lx28.diff.genesc2t32ct

### plot

ggplot() +
  geom_point(dat = filter(lx28.diff.genesc2t32ct, is.na(category)),
             aes(x=FC, y = -log10(padj)),
             color = "grey", alpha = 0.8
  ) +
  geom_point(dat = filter(lx28.diff.genesc2t32ct, 
                          category == "C2"),
             aes(x=FC, y = -log10(padj)),
             color = "red", alpha = 0.8
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") + 
  geom_point(dat = filter(lx28.diff.genesc2t32ct, 
                          gene.name %in% c(
                            "108007264","108009404","108009449",
                            "108009614","108009718","108009919",
                            "108010308","108012141","108014002",
                            "108016014","108016422","108016967"
                          )),
             aes(x=FC, y = -log10(padj)),
             color = "gold", shape = 18, size = 3 ) +
  geom_text_repel(dat = filter(lx28.diff.genesc2t32ct, 
                               gene.name %in% c(
                                 "108007264","108009404","108009449",
                                 "108009614","108009718","108009919",
                                 "108010308","108012141","108014002",
                                 "108016014","108016422","108016967"
                               )),
                  aes(x=FC, y = -log10(padj), label = gene.name),
                  size = 1.5, max.overlaps = 100 ) +
  theme_bw()->
  diffexp.lx.C2_t32
ggsave(diffexp.lx.C2_t32, 
       file = "diffexp.lx.C2_t32.pdf", 
       w = 3, h = 3)

### Do a ven diagram
adapt.annot %>%
  filter(P_C2 >= 2) %>% 
  group_by(Gene) %>%
  summarise(N = n()) %>%
  filter(!is.na(Gene) | Gene != "-") %>%
  .$Gene->
  genes_C2_names

adapt.annot %>%
  filter(BF_ctmin >= 15) %>% 
  group_by(Gene) %>%
  summarise(N = n()) %>%
  filter(!is.na(Gene) | Gene != "-") %>%
  .$Gene->
  genes_CTmin_names

adapt.annot %>%
  filter(BF_T32 >= 15) %>% 
  group_by(Gene) %>%
  summarise(N = n()) %>%
  filter(!is.na(Gene) | Gene != "-") %>%
  .$Gene->
  genes_T32_names

lx28.diff %>%
  filter(padj <= 0.01) %>%
  .$gene.name ->
  genes_RNA_names
  

library(ggVennDiagram)
library(VennDiagram)

x <- list(
  C2 = genes_C2_names,
  CTmin = genes_CTmin_names,
  T32 = genes_T32_names,
  RNA = genes_RNA_names
)


olap <- calculate.overlap(x)
lapply(olap, length)
names(olap)

ggVennDiagram(x, label = "count") +
scale_fill_gradient(low = "white", high = "white") +
theme_void()-> venn.g

ggsave(venn.g, file = "venn.pdf")


