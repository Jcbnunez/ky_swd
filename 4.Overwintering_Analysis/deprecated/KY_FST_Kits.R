#### SWD KY PCA
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
library(reshape2)
library(gdsfmt)
library(SNPRelate)
library(gtools)

root="/gpfs2/scratch/kaeller/labtasks/"

### Metadata
meta <- "/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023/2024_Joaquin_KY/swd.metadata.txt"
samps <- fread(meta)
setDT(samps)
samps$year = year(as.Date(samps$Collection_date, format = "%Y-%m-%d"))

### open GDS
geno_auto = readRDS("/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023/2024_Joaquin_KY/auto.pooldata.new.rds")

subset <- pooldata.subset(
  geno_auto,
  min.cov.per.pool = -1,
  max.cov.per.pool = 1e+06,
  #creating a min allele freq of 0.05
  min.maf = 0.05,
  verbose = TRUE
)

ky_subset <- pooldata.subset(
	subset,
	pool.index=16:35
	)
	

##### Begin temporal code

kentucky_samps <- samps %>%
filter(province == "Kentucky")

kentucky_samps %>%
mutate(Collection_date = as.Date(Collection_date, 
format = "%Y-%m-%d")) ->
working.obj_ky

# number rows in working.obj_ky (ie length)
L = dim(working.obj_ky)[1]

## Combining each collection with all other collections
comp_vector = combinations(
  L,
  2, 
  v=1:L,
  set=TRUE, 
  repeats.allowed=FALSE)
  
print("Create combination vector")

comp_vector %<>%
  as.data.frame() %>%
  mutate(day_diff = NA)
  
for(i in 1:dim(comp_vector)[1]) {
  
  #showing difference in days between each collection
  date1=working.obj_ky$Collection_date[comp_vector[i,1]]
  date2=working.obj_ky$Collection_date[comp_vector[i,2]]
  comp_vector$day_diff[i] = abs(as.numeric(date1-date2))
  
  comp_vector$pop1[i] = working.obj_ky$city[comp_vector[i,1]]
  comp_vector$pop2[i] = working.obj_ky$city[comp_vector[i,2]]
  
  comp_vector$samp1[i] = working.obj_ky$sampleId[comp_vector[i,1]]
  comp_vector$samp2[i] = working.obj_ky$sampleId[comp_vector[i,2]]
}


### Calculate FST

#Generate outfile object
outfile = data.frame(
  samp1 = rep(NA, dim(comp_vector)[1]),
  samp2 = rep(NA, dim(comp_vector)[1]),
  FST = rep(NA, dim(comp_vector)[1])
)

for(i in 1:dim(comp_vector)[1]){
  
  print(i/dim(comp_vector)[1] * 100)
  
#  samps_to_compare = c(comp_vector$V1[i], comp_vector$V2[i])
  
   pooli <- pooldata.subset(
	ky_subset,
	pool.index=c(comp_vector$V1[i],comp_vector$V2[i])
	)
	
  fst.out <- computeFST(pooli, method = "Anova")
  
  outfile$samp1[i] = comp_vector$samp1[i]
  outfile$samp2[i] = comp_vector$samp2[i]
  outfile$FST[i] = fst.out$Fst
  
}

fst.out <- computeFST(pooli, method = "Anova")

left_join(comp_vector, outfile) -> Out_comp_vector_samepops

save(Out_comp_vector_samepops,
file = paste(root, "kentucky_fst_samples", ".Rdata", sep = ""))