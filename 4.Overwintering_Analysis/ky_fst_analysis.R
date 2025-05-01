#########################
### Loading Libraries ###
#########################

library(tidyverse)
library(magrittr)
library(data.table)
library(reshape2)
library(poolfstat)
library(FactoMineR)
require(gtools)
require(foreach)
#library(ggpmisc)
#library(gt)

##################################################################################
######################
### Metadata - All ###
######################
##################################################################################

samps <- fread("/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023/2024_Joaquin_KY/swd.metadata.txt")

dat_f <- "/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023/2024_Joaquin_KY/auto.pooldata.new.rds"    
auto_dat <- readRDS(dat_f)


#### Global PCA
sample.names <- auto_dat@poolnames
exclude <- c("KY20")


all.set <-
  pooldata.subset(
    auto_dat,
    pool.index = which(sample.names != exclude),
    min.cov.per.pool = 10,
    max.cov.per.pool = 150,
    min.maf = 0.05,
    return.snp.idx = TRUE,
    verbose = TRUE
  )

all.pca =
  randomallele.pca(
    all.set,
    scale = FALSE
  )

all.pca$pop.loadings %>%
  as.data.frame() %>%
  mutate(sampleId_orig = rownames(.)) %>%
  left_join(samps)-> pca.all

save(pca.all, file = "pca.all.Rdata")

### N_Ame PCA
samps %>% 
  filter(continent == "North_America") %>%
  filter(!sampleId_orig %in% c("KY20","KY17", "KY11")) %>%
  .$sampleId_orig -> keep_NAM

Name.set <-
  pooldata.subset(
    auto_dat,
    pool.index = which(sample.names %in% keep_NAM),
    min.cov.per.pool = 10,
    max.cov.per.pool = 150,
    min.maf = 0.05,
    return.snp.idx = TRUE,
    verbose = TRUE
  )

Name.pca =
  randomallele.pca(
    Name.set,
    scale = FALSE
  )

Name.pca$pop.loadings %>%
  as.data.frame() %>%
  mutate(sampleId_orig = rownames(.)) %>%
  left_join(samps)-> pca.Name

save(pca.Name, file = "pca.Name.Rdata")


#####
sample.names <- auto_dat@poolnames
exclude <- c("KY20","KY17", "KY11")

ky.tag <- grep("KY", sample.names)
va.tag <- grep("VA", sample.names)

indexer.ky <- data.frame(name = sample.names[grep("KY", sample.names)],
           index = ky.tag)
		   
ky.tag.flt <- indexer.ky$index[which(!indexer.ky$name %in% exclude)]

info1 <- samps %>%
  dplyr::select(pop1=sampleId_orig,
                city1=city,
                date1=Collection_date,
				fruit1=fruit_type)
info2 <- samps %>%
  dplyr::select(pop2=sampleId_orig,
                city2=city,
                date2=Collection_date,
				fruit2=fruit_type)
				
my.formula <- y ~ x				
				
##################################################################################
############################
######### KY - All #########
##################################################################################

###############
## KY - Pool ##
###############

ky.set <-
	pooldata.subset(
	auto_dat,
	pool.index = ky.tag.flt,
	min.cov.per.pool = 10,
	max.cov.per.pool = 150,
	min.maf = 0.05,
	return.snp.idx = TRUE,
	verbose = TRUE
	)

ky.set@refallele.readcount %>% colMeans()

snp.dt.ky <- 
	ky.set@snp.info %>%
	as.data.frame() %>%
	mutate(
		snp.id = rownames(.))
  
af.ky <- ky.set@refallele.readcount/ky.set@readcoverage
rownames(af.ky) = snp.dt.ky$snp.id

##################################################################################
#####################
### KY - Analysis ###
#####################

pca.ran.ky =
randomallele.pca(
  ky.set,
  scale = FALSE)

PCA.coords.ky <-
	pca.ran.ky$pop.loadings %>%
	as.data.frame() %>%
	mutate(sampleId_orig = rownames(.) ) %>%
	left_join(samps) %>%
	mutate(Year = year(Collection_date))
	
fst.ky <- compute.pairwiseFST(ky.set,  method = "Anova")
	
##################################################################################
#####################
### KY - Matrixes ###
#####################

fast_matrix.ky <-
	fst.ky@PairwiseFSTmatrix %>%
	as.data.frame() %>%
	mutate(pop1 = colnames(.)) %>%
	melt(id = "pop1", value.name = "fst", 
    variable.name = "pop2")
				
fast_matrix.meta.ky <- fast_matrix.ky %>%
  left_join(info1) %>%
  left_join(info2) %>%
  mutate(delta = abs(as.Date(date1)-as.Date(date2))) %>%
  mutate(year1 = year(as.Date(date1))) %>%
  mutate(year2 = year(as.Date(date2))) %>%
  mutate(comp = paste(city1, city2, sep= "_")) |>
  mutate(dataset = "Kentucky")
  
save(fast_matrix.meta.ky, file = "fst.matrix.meta.ky.Rdata")
  
##################################################################################
###############
## VA - Pool ##
###############

va.set <-
	pooldata.subset(
	auto_dat,
	pool.index = va.tag,
	min.cov.per.pool = 10,
	max.cov.per.pool = 150,
	min.maf = 0.05,
	return.snp.idx = TRUE,
	verbose = TRUE
	)		   
		   
va.set@refallele.readcount %>% colMeans()

snp.dt.va <- 
	va.set@snp.info %>%
	as.data.frame() %>%
	mutate(
		snp.id = rownames(.))
  
af.va <- va.set@refallele.readcount/va.set@readcoverage
rownames(af.va) = snp.dt.va$snp.id

##################################################################################		   
#####################
### VA - Analysis ###
#####################

pca.ran.va =
randomallele.pca(
  va.set,
  scale = FALSE)

PCA.coords.va <-
	pca.ran.va$pop.loadings %>%
	as.data.frame() %>%
	mutate(sampleId_orig = rownames(.) ) %>%
	left_join(samps) %>%
	mutate(Year = year(Collection_date))

fst.va <- compute.pairwiseFST(va.set,  method = "Anova")		   

##################################################################################		   	   
#####################
### VA - Matrixes ###
#####################

fast_matrix.va <-
	fst.va@PairwiseFSTmatrix %>%
	as.data.frame() %>%
	mutate(pop1 = colnames(.)) %>%
	melt(id = "pop1", value.name = "fst", 
    variable.name = "pop2")
				
fast_matrix.meta.va <- fast_matrix.va %>%
  left_join(info1) %>%
  left_join(info2) %>%
  mutate(delta = abs(as.Date(date1)-as.Date(date2))) %>%
  mutate(year1 = year(as.Date(date1))) %>%
  mutate(year2 = year(as.Date(date2))) %>%
  mutate(comp = paste(city1, city2, sep= "_")) |>
  mutate(dataset = "Virginia")
  
  
################
## KY-VA Pool ##
################

kyva.set <- pooldata.subset(
  auto_dat,
  pool.index = c(ky.tag.flt,va.tag),
  min.cov.per.pool = 10,
  max.cov.per.pool = 150,
  min.maf = 0.05,
  return.snp.idx = TRUE,
  verbose = TRUE
)

kyva.set@refallele.readcount %>% colMeans()

snp.dt.kyva <- kyva.set@snp.info %>%
  as.data.frame() %>%
  mutate(snp.id = rownames(.))

##################################################################################
####################
## KY VA Analysis ##
####################

pca.ran.kyva =
  randomallele.pca(
    kyva.set,
    scale = FALSE
  )

PCA.coords.kyva <- pca.ran.kyva$pop.loadings %>%
  as.data.frame() %>%
  mutate(sampleId_orig = rownames(.) ) %>%
  left_join(samps) %>%
  mutate(Year = year(Collection_date))
  
fst.kyva <- compute.pairwiseFST(kyva.set,  method = "Anova")	
	
##################################################################################	  
########################
### KY VA - Matrixes ###
########################

fast_matrix.kyva <-
	fst.kyva@PairwiseFSTmatrix %>%
	as.data.frame() %>%
	mutate(pop1 = colnames(.)) %>%
	melt(id = "pop1", value.name = "fst", 
    variable.name = "pop2")  

fast_matrix.meta.kyva <- fast_matrix.kyva %>%
  left_join(info1) %>%
  left_join(info2) %>%
  mutate(delta = abs(as.Date(date1)-as.Date(date2))) %>%
  mutate(year1 = year(as.Date(date1))) %>%
  mutate(year2 = year(as.Date(date2))) %>%
  mutate(comp = paste(city1, city2, sep= "_")) |>
  mutate(states = if_else(
	condition = comp == "Berea_Berea",
	true = "Kentucky",
	false = "DifferentStates"),
	states = if_else(
	condition = comp == "Lexington_Berea",
	true = "Kentucky",
	false = states),
	states = if_else(
	condition = comp == "Berea_Lexington",
	true = "Kentucky",
	false = states),
	states = if_else(
	condition = comp == "Lexington_Lexington",
	true = "Kentucky",
	false = states),
	states = if_else(
	condition = comp == "Charlottesville_Charlottesville",
	true = "Virginia",
	false = states)
	)
	
fast_matrix.kyva.separated <-
	rbind(fast_matrix.meta.va, fast_matrix.meta.ky)	 
	
save(fast_matrix.meta.kyva, file = "fst.matrix.meta.kyva.Rdata")
	