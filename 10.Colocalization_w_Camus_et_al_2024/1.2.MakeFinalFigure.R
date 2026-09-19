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
source("/netfiles/nunezlab/Shared_Resources/Software/baypass_public/utils/baypass_utils.R")





All_FETs %>%
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
  theme_classic() + scale_shape_manual(values = 21:22) ->
  ORFET_plot

ggsave(ORFET_plot, file = "ORFET_plot.pdf", w=4.3, h = 4)
