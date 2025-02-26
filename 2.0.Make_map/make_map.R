### Make Map

library(tidyverse)
library(magrittr)
library(reshape2)
library(data.table)
library(foreach)

library(rnaturalearth)
library(rnaturalearthdata)

###
meta <- "/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023/2024_Joaquin_KY/swd.metadata.txt"

metdat <- fread(meta)

####
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)


ggplot(data = world) +
  geom_sf(fill= "antiquewhite") +
  coord_sf(xlim = c(-135.15, 160.99), ylim = c(-55.00, 69.00), expand = FALSE) + 
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.2), 
        panel.background = element_rect(fill = "aliceblue")) +
  geom_point(data = metdat, 
             aes(x=long, y = lat, fill = Range), size = 2.2, shape = 21) -> SAMP.map
ggsave(SAMP.map, file = "SAMP.map.png", h = 4, w = 6)

###
metdat %>%
  mutate(semidate = as.Date(paste(1,min_month, year, sep = "/" ), 
                            format = "%d/%m/%Y")) %>%
  ggplot(aes(
    x=semidate,
    y=province
  )) + geom_point()  -> semi.dates
ggsave(semi.dates, file = "semi.dates.png", h = 4, w = 6)
