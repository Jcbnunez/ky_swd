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
library(patchwork)
source("/netfiles/nunezlab/Shared_Resources/Software/baypass_public/utils/baypass_utils.R")

### Load plots

load("PCdata_ecovars.Rdata")


ky_pc_loadings %>%
  mutate(julian_date = yday(Collection_date)) %>%
  filter(!is.na(m)) %>%
  mutate(season = case_when(m <= 8 ~ "spring",
                            m > 8 ~ "fall")) %>%
  ggplot(aes(
    x=V1,
    y=V2,
    fill =julian_date
  )) + geom_point(size = 3, shape = 21) + 
  scale_fill_gradient2(low = "springgreen", 
                       high = "firebrick", 
                       midpoint = 250) +
  theme_bw()  ->
  PCA_plot.KY
#ggsave(PCA_plot.KY, file = "PCA.KY.FINAL.pdf", w = 4.9, h =4)

PCdata_ecovars %>%
  ggplot(aes(
    x=`Juliean Date`,
    y=V1,
    fill =value
  )) + geom_smooth(method = "lm", color = "black") + 
  geom_point(size = 3, shape = 21) + 
  scale_fill_gradient2(low = "steelblue", 
                       high = "firebrick", 
                       midpoint = 16) +
  theme_bw() ->
  V1.ky._plot
#ggsave(V1.ky._plot, file = "V1.ky._plot.pdf", w = 4.8, h =4)


###
### plot
load("lx28.diff.outliers.Rdata")
filter(lx28.diff.outliers, 
       KYseasonal_outlier == TRUE &
         padj < 0.05)

filter(lx28.diff.outliers, 
       Top1percentC == TRUE &
         padj < 0.05)

# Define significant expression
lx28.diff.outliers$DEG <- lx28.diff.outliers$padj < 0.01

# Construct 2 x 2 contingency table
tab <- table(
  Top1percentC = lx28.diff.outliers$Top1percentC,
  DEG = lx28.diff.outliers$DEG
)

tab

# Fisher's exact test
fet <- fisher.test(tab)

fet

ggplot() +
  geom_point(dat = filter(lx28.diff.outliers, 
                          KYseasonal_outlier == FALSE &
                            Top1percentC == FALSE
                            ),
             aes(x=FC, y = -log10(padj)),
             color = "grey", alpha = 0.6
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") + 
  geom_point(dat = filter(lx28.diff.outliers, 
                            Top1percentC == TRUE,
                          padj < 0.05
  ),
  aes(x=FC, y = -log10(padj)),
  color = "#8FBC8F", alpha = 0.7
  ) +
  geom_point(dat = filter(lx28.diff.outliers, 
                          KYseasonal_outlier == TRUE,
                          padj < 0.05),
             aes(x=FC, y = -log10(padj)),
             fill = "gold", color = "black", 
             alpha = 1.0, shape = 23, size = 3
  ) +
  theme_bw()->
  diffexp.lx.C2_t32
ggsave(diffexp.lx.C2_t32, 
       file = "diffexp.lx.C2_t32.pdf", 
       w = 3, h = 3)


plots = (
  PCA_plot.KY + V1.ky._plot + diffexp.lx.C2_t32
)

ggsave(plots, file = "fig3downpanles.pdf", w = 10, h = 3)

