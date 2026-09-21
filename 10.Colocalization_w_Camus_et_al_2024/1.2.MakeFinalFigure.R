### Co-localization of invasiveness and adaptive tracking SWD
### Plot
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
library(patchwork)
source("/netfiles/nunezlab/Shared_Resources/Software/baypass_public/utils/baypass_utils.R")

###All data
load("/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023_BAYPASS_MG/Final_Datasets/C2_BFct_BFT_df.Invsasion_outliers.allAnalyses.Rdata")
#adapt_inva.final
adapt_inva.final %>%
  filter(KYseasonal_outlier == TRUE & NAme_invasion_outlier == TRUE ) %>%
  dim

### load plot data
load("plots_invasion_seasonality_figure.Rdata")


All_FETs %>%
  filter(Test %in% c("Seasonality","CTmin")) %>%
  ggplot(
    aes(
      x=Test,
      y=log2(OR),
      ymin=log2(lci),
      ymax=log2(uci),
      shape=Continent
    )
  ) + 
  geom_hline(yintercept = 0, linetype= "dashed") +
  geom_errorbar(width = 0.5, position=position_dodge(width=0.5)) +
  geom_point(size = 4, fill = "grey", position=position_dodge(width=0.5)) +
  theme_classic() + scale_shape_manual(values = 21:22) + theme(
    legend.position = "inside",
    legend.position.inside = c(0.03, 0.97),
    legend.justification = c(0, 1),
    legend.background = element_rect(
      fill = scales::alpha("white", 0.8),
      color = NA
    )) ->
  ORFET_plot

#ggsave(ORFET_plot, file = "ORFET_plot.pdf", w=4.3, h = 4)

### PCA
pca_result.final$perc.var

####
pc_loadings %>%
  filter(!is.na(m)) %>%
  mutate(season = case_when(m <= 8 ~ "spring",
                            m > 8 ~ "fall")) %>%
  ggplot(aes(
    x=V1,
    y=V2,
    shape=Range, fill =m, #label = sampleId_orig
    #shape=Range, fill =season, label = sampleId_orig
  )) + geom_point(size = 3) + #geom_text(size = 0.5) +
  scale_shape_manual(values = 21:23) +
  scale_fill_gradient2(low = "springgreen", 
                       high = "firebrick", 
                       midpoint = 8.0) +
  theme_bw() + theme(
    legend.position = "inside",
    legend.position.inside = c(0.03, 0.97),
    legend.justification = c(0, 1),
    legend.background = element_rect(
      fill = scales::alpha("white", 0.8),
      color = NA
    )) ->
  PCA_plot
ggsave(PCA_plot, file = "PCA.FINAL.pdf", w = 4.8, h =4)

###
### allele freq plot

afs.id.annot %>%
  group_by(snp_id, 
           #city, 
           #year, 
           Time.point,
           #Collection_date
  ) %>%
  summarize(AFm = ci(value, na.rm = T)[1],
            AFlci= ci(value, na.rm = T)[2],
            AFuci= ci(value, na.rm = T)[3],
            mT32=mean(T32), mCTmin = mean(ctmin) ) ->
  plot_freq_data

plot_freq_data %>%
  mutate(numTim = case_when(
    Time.point == "First" ~ 1,#+(year-2020)* 3,
    Time.point == "Second" ~ 2,#+(year-2020)* 3,
    Time.point == "Third" ~ 3,#+(year-2020)* 3
  )
  ) %>%
  ggplot(aes(
    x=as.numeric(numTim),
    y=AFm,
    #group=city,
    ymin = AFlci, ymax = AFuci,
    #shape = as.character(year)
  )) + 
  geom_line() +
  geom_errorbar(width = 0.1) +
  geom_point(aes(color = mCTmin, 
                 #shape=as.character(year)
  ),  size = 3, alpha = 0.95
  ) + 
  scale_color_gradient2(low = "blue", 
                        high = "red", 
                        midpoint = 4.35) +
  scale_shape_manual(values = 21:24) +
  scale_x_continuous(breaks = 1:3) + 
  ylim(0.00,1.0) + facet_grid(~snp_id) +
  theme_bw() + theme(
    legend.position = "inside",
    legend.position.inside = c(0.03, 0.97),
    legend.justification = c(0, 1),
    legend.background = element_rect(
      fill = scales::alpha("white", 0.8),
      color = NA
    ))->
  ctmin_traj
#ggsave(ctmin_traj, file = "ctmin_traj.pdf", w = 6.5, h = 6.0)


### worldwide plot

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

ggplot(data = world) +
  geom_sf(fill= "white") +
  coord_sf(xlim = c(-160.15, 160.99), ylim = c(-55.00, 69.00), expand = FALSE) + 
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.2), 
        panel.background = element_rect(fill = "aliceblue")) +
  geom_point(data = plot_freq_data.ALL, 
             aes(x= longm, y = latm, fill = AFm, shape = Range), 
             size = 3) +
  scale_shape_manual(values = 21:23) +
  scale_fill_viridis(option = "A"#, midpoint = 0.7
  ) + theme(
    legend.position = "inside",
    legend.position.inside = c(0.03, 0.97),
    legend.justification = c(0, 1),
    legend.background = element_rect(
      fill = scales::alpha("white", 0.8),
      color = NA
    ))-> SAMP.AFs.map
#ggsave(SAMP.AFs.map, file = "SAMP.AFs.map.pdf", h = 4, w = 5)


####
top <- ORFET_plot + PCA_plot +
  plot_layout(ncol = 2)

bottom <- ctmin_traj + SAMP.AFs.map +
  plot_layout(
    ncol = 2,
    widths = c(0.5, 1.5)
  )

ggsave(top / bottom,
       file = "Figure4.allpnales.pdf", h = 5, w = 5 )




