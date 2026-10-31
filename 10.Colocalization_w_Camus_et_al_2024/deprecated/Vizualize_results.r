### Analysis and vizualization for colocalization
### libraries
library(FactoMineR)
library(factoextra)
library(SeqArray)
library(data.table)
library(foreach)
library(tidyverse)
library(magrittr)
library(vroom)
library(poolfstat)
library(rnaturalearth)
library(rnaturalearthdata)

### Metadata
meta <- "/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023_BAYPASS_MG/2024_Joaquin_KY/swd.metadata.txt"
samps <- fread(meta)
setDT(samps)
samps$year = year(as.Date(samps$Collection_date, format = "%Y-%m-%d"))
