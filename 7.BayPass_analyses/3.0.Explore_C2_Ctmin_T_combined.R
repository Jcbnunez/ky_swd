### Explore the outliers from different tests

library(data.table)
library(tidyverse)
library(foreach)
library(forcats)
library(magrittr)

load("C2_BFct_BFT_df.allAnalyses.Rdata")
load("annots.flt.Rdata")
## "annots.flt"     "C2_BFct_BFT_df"
#####
#### Lets explore the data a bit. these the top hits for C2. <-- jcbn here
c2_top_hits_annotated <- fread("Top Genes in C2_SNPs.tsv")
c2_top_hits_annotated %<>% mutate(SNP_id = paste(chr, pos, sep = "_"))
C2_BFct_BFT_df %>%
  mutate(top_hit = case_when(SNP_id %in% c2_top_hits_annotated$SNP_id ~ "yes",
                             TRUE ~ "no"))-> 
  C2_BFct_BFT_df


#####
C2_BFct_BFT_df %>%
  filter(BF_T32 >=5)
#27432

C2_BFct_BFT_df %>%
  filter(BF_ctmin >=5)
#11522


C2_BFct_BFT_df %>%
  filter(BF_T32 >=5) %>%
  group_by(top_hit) %>%
  summarize(N=n())
#top_hit     N
#<chr>   <int>
# 1 no      27410
# 2 yes        22

C2_BFct_BFT_df %>%
  filter(BF_ctmin >=5) %>%
  group_by(top_hit) %>%
  summarize(N=n())
#top_hit     N
#<chr>   <int>
#  1 no      11522



###
C2_BFct_BFT_df %>%
  filter(P_C2 > 2) %>%
  ggplot(aes(
    x=P_C2,
    y=BF_ctmin,
    shape = top_hit,
    size = top_hit,
    fill = BF_ctmin > 5
  )) + geom_point() +
  geom_hline(yintercept = 5) +
  scale_size_manual(values = c(1,3)) +
  scale_shape_manual(values = 23:24) ->
  c2_bfctmin.plot

ggsave(c2_bfctmin.plot, file = "c2_bfctmin.plot.pdf")

C2_BFct_BFT_df %>%
  filter(P_C2 > 2) %>%
  ggplot(aes(
    x=P_C2,
    y=BF_T32,
    shape = top_hit,
    size = top_hit,
    fill = BF_T32 > 5
  )) + geom_point() +
  geom_hline(yintercept = 5) +
  scale_size_manual(values = c(1,3)) +
  scale_shape_manual(values = 23:24) ->
  c2_bfT32.plot

ggsave(c2_bfT32.plot, file = "c2_bfT32.plot.pdf")
  
  
C2_BFct_BFT_df %>%
  filter(BF_T32 > 5 & BF_ctmin > 5) %>%
  ggplot(aes(
    x=BF_ctmin,
    y=BF_T32,
    shape = top_hit,
    size = top_hit,
  )) + geom_point() +
  geom_hline(yintercept = 5) +
  scale_size_manual(values = c(1,3)) +
  scale_shape_manual(values = 23:24) ->
  bfCTmin_bfT32.plot

ggsave(bfCTmin_bfT32.plot, file = "bfCTmin_bfT32.plot.pdf")

#####
##data:  FETdat with Ctmin
##p-value < 2.2e-16
##alternative hypothesis: true odds ratio is greater than 1
##95 percent confidence interval:
##  1.702682      Inf
##sample estimates:
##  odds ratio 
##1.897684 
##
##data:  FETdatT32
##p-value < 2.2e-16
##alternative hypothesis: true odds ratio is greater than 1
##95 percent confidence interval:
##  8.671068      Inf
##sample estimates:
##  odds ratio 
##8.981149 
##


#data:  FETdatT32 with only top outliers
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is greater than 1
#95 percent confidence interval:
#  10.22897      Inf
#sample estimates:
#  odds ratio 
#15.3168 
#> FETdatT32
#[,1]    [,2]
#[1,]   22   27410
#[2,]  290 5533884
#

### Enrichment Analyses
C2_BFct_BFT_df %>%
  filter(top_hit == "yes") %>%
  filter(BF_ctmin > 5) -> c2yesBfyes
C2_BFct_BFT_df %>%
  filter(top_hit == "yes") %>%
  filter(BF_ctmin < 5) -> c2yesBfno
C2_BFct_BFT_df %>%
  filter(top_hit == "no") %>%
  filter(BF_ctmin > 5) -> c2noBfyes
C2_BFct_BFT_df %>%
  filter(top_hit == "no") %>%
  filter(BF_ctmin < 5) -> c2noBfno

###
dim(c2yesBfyes)[1] +
dim(c2yesBfno)[1] +
dim(c2noBfyes)[1] +
dim(c2noBfno)[1] 

FETdat <-
matrix(c(dim(c2yesBfyes)[1], dim(c2yesBfno)[1], 
			dim(c2noBfyes)[1], dim(c2noBfno)[1]),
       nrow = 2)
       
fisher.test(FETdat, alternative = "greater")

#####
### Enrichment Analyses
C2_BFct_BFT_df %>%
  filter(top_hit == "yes") %>%
  filter(BF_T32 > 5) -> c2yesBfTyes
C2_BFct_BFT_df %>%
  filter(top_hit == "yes") %>%
  filter(BF_T32 < 5) -> c2yesBfTno
C2_BFct_BFT_df %>%
  filter(top_hit == "no") %>%
  filter(BF_T32 > 5) -> c2noBfTyes
C2_BFct_BFT_df %>%
  filter(top_hit == "no") %>%
  filter(BF_T32 < 5) -> c2noBfTno

###
dim(c2yesBfTyes)[1] +
dim(c2yesBfTno)[1] +
dim(c2noBfTyes)[1] +
dim(c2noBfTno)[1] 

FETdatT32 <-
matrix(c(dim(c2yesBfTyes)[1], dim(c2yesBfTno)[1], 
			dim(c2noBfTyes)[1], dim(c2noBfTno)[1]),
       nrow = 2)
       
fisher.test(FETdatT32, alternative = "greater")

####
C2_BFct_BFT_df %>%
  filter(P_C2 > 2) %>%
  filter(BF_T32 > 5) %>%
  filter(BF_ctmin > 5) -> c2yesBfTyes

### GO Enrichment sets
### Lists
C2_BFct_BFT_df %>%
  filter(BF_T32 > 5) %>%
  select(chr, pos) -> T32_outliers
  
C2_BFct_BFT_df %>%
  filter(BF_ctmin > 5) %>%
  select(chr, pos) -> ctmin_outliers

write.table(T32_outliers, 
			file = "T32_outliers.txt", 
			append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")
            
write.table(ctmin_outliers, 
			file = "ctmin_outliers.txt", 
			append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")  
            
