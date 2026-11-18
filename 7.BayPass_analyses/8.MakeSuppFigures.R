### Extra plots for CTmin and T32

library(data.table)
library(tidyverse)
library(foreach)
library(forcats)
library(magrittr)

baypassdf <- get(load("/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023_BAYPASS_MG/Final_Datasets/C2_BFct_BFT_df.Invsasion_outliers.allAnalyses.Rdata"))

baypassdf %>%
  filter(chr %in% c("chr2L", "chr2R", "chr3", "chrX" )) %>%
  filter(BF_T32 > 0 & BF_T32 < 10) %>%
  slice_sample(n = 10000) ->
  bay_sub
baypassdf %>%
  filter(chr %in% c("chr2L", "chr2R", "chr3", "chrX" )) %>%
  filter(BF_T32 > 5) ->
  bay_20T32

baypassdf %>%
  filter(chr %in% c("chr2L", "chr2R", "chr3", "chrX" )) %>%
  filter(BF_T32 > 40) %>%
  group_by(Gene) %>%
  slice_max(BF_T32)->
  genest32

ggplot() +
  geom_point(data = bay_20T32,
             aes(x= pos/1e6, y = BF_T32,
                 color =  BF_T32 >= 20)) +
  scale_color_manual(values = c("black", "red")) +
  geom_point(data = bay_sub,
             aes(x= pos/1e6, y = BF_T32)) +
  geom_text(data = genest32,
             aes(x= pos/1e6, y = BF_T32, label = Gene )) +
  geom_hline(yintercept = 20, linetype = "dashed") +
  theme_classic() + theme(legend.position = "none") +
  facet_grid(~chr, scales = "free_x") +
  xlab("Genomic Position (Mb)")+
  ylab("Bayes Factor")->
  bf32_plot

ggsave(bf32_plot, file= "bf32_plot.pdf",
       w= 8, h = 3)


###


baypassdf %>%
  filter(chr %in% c("chr2L", "chr2R", "chr3", "chrX" )) %>%
  filter(BF_ctmin > 0 & BF_ctmin < 10) %>%
  slice_sample(n = 10000) ->
  bay_subct

baypassdf %>%
  filter(chr %in% c("chr2L", "chr2R", "chr3", "chrX" )) %>%
  filter(BF_ctmin > 5) ->
  bay_20ct

baypassdf %>%
  filter(chr %in% c("chr2L", "chr2R", "chr3", "chrX" )) %>%
  filter(BF_ctmin > 35) %>%
  group_by(Gene) %>%
  slice_max(BF_ctmin)->
  genesct

ggplot() +
  geom_point(data = bay_20ct,
             aes(x= pos/1e6, y = BF_ctmin,
                 color =  BF_ctmin >= 20)) +
  scale_color_manual(values = c("black", "red")) +
  geom_point(data = bay_subct,
             aes(x= pos/1e6, y = BF_ctmin)) +
  geom_text(data = genesct,
            aes(x= pos/1e6, y = BF_ctmin, label = Gene )) +
  geom_hline(yintercept = 20, linetype = "dashed") +
  theme_classic() + theme(legend.position = "none") +
  facet_grid(~chr, scales = "free_x") +
  xlab("Genomic Position (Mb)")+
  ylab("Bayes Factor")->
  bfct_plot

ggsave(bfct_plot, file= "bfct_plot.pdf",
       w= 8, h = 3)




