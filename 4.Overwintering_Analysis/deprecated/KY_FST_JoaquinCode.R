### Explore SWD FST

library(tidyverse)
library(magrittr)
library(data.table)
library(reshape2)
library(poolfstat)
library(FactoMineR)
require(gtools)
require(foreach)

samps <- fread("/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023/2024_Joaquin_KY/swd.metadata.txt")

dat_f <- "/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023/2024_Joaquin_KY/auto.pooldata.new.rds"    
auto_dat <- readRDS(dat_f)

sample.names <- auto_dat@poolnames
exclude <- c("KY20")

ky.tag <- grep("KY", auto_dat@poolnames)

data.frame(name = auto_dat@poolnames[grep("KY", auto_dat@poolnames)],
           index = ky.tag) -> indexer
ky.tag.flt <- indexer$index[which(!indexer$name %in% exclude)]
#ky.tag.flt <- indexer$index[which(indexer$name %in% pops)]

####
ky.set <- pooldata.subset(
  auto_dat,
  pool.index = ky.tag.flt,
  min.cov.per.pool = 10,
  max.cov.per.pool = 150,
  min.maf = 0.05,
  return.snp.idx = TRUE,
  verbose = TRUE
)

ky.set@refallele.readcount %>% colMeans()

ky.set@snp.info %>%
  as.data.frame() %>%
  mutate(snp.id = rownames(.)) ->
  snp.dt

af.ky <- ky.set@refallele.readcount/ky.set@readcoverage
rownames(af.ky) = snp.dt$snp.id

pca.ran =
randomallele.pca(
  ky.set,
  scale = FALSE
)

#save(pca.ran, file = "pca.ran.min40.max120.maf1.Rdata")

pca.ran$pop.loadings %>%
  as.data.frame() %>%
  mutate(sampleId_orig = rownames(.) ) %>%
  left_join(samps) %>%
  mutate(Year = year(Collection_date)) ->
  PCA.coords

anova(lm(V1~Year+city+fruit_type, data = PCA.coords))
anova(lm(V2~Year+city+fruit_type, data = PCA.coords))
anova(lm(V3~Year+city+fruit_type, data = PCA.coords))
anova(lm(V4~Year+city+fruit_type, data = PCA.coords))
anova(lm(V5~Year+city+fruit_type, data = PCA.coords))
anova(lm(V6~Year+city+fruit_type, data = PCA.coords))
anova(lm(V7~Year+city+fruit_type, data = PCA.coords))
anova(lm(V8~Year+city+fruit_type, data = PCA.coords))
anova(lm(V9~Year+city+fruit_type, data = PCA.coords))
anova(lm(V10~Year+city+fruit_type, data = PCA.coords))

PCA.coords %>%
  ggplot(aes(
    x=V1,
    y=V6,
    color= as.character(Year),
    shape=fruit_type
  )) + 
  geom_point(size = 5) +
  facet_grid(~city)->
  PCA12
ggsave(PCA12, file = "PCA12.pdf",
       w= 6, h = 5)

fst <- compute.pairwiseFST(ky.set,  method = "Anova")

fst@PairwiseFSTmatrix %>%
  as.data.frame() %>%
  mutate(pop1 = colnames(.)) %>%
  melt(id = "pop1", value.name = "fst", 
       variable.name = "pop2") ->
  fast_matrix
###
samps %>%
  dplyr::select(pop1=sampleId_orig,
                city1=city,
                date1=Collection_date) -> info1
samps %>%
  dplyr::select(pop2=sampleId_orig,
                city2=city,
                date2=Collection_date) -> info2

fast_matrix %>%
  left_join(info1) %>%
  left_join(info2) %>%
  mutate(delta = abs(as.Date(date1)-as.Date(date2))) %>%
  mutate(year1 = year(as.Date(date1))) %>%
  mutate(year2 = year(as.Date(date2))) %>%
  mutate(comp = paste(city1, city2, sep= "_")) ->
  fast_matrix.meta



###

fast_matrix.meta %>%
  ggplot(aes(
    x=delta,
    y=(fst)
  )) + 
  geom_point() +
  geom_smooth(method = "lm")  ->
  fst
ggsave(fst, file = "fst.pdf",
       w= 6, h = 5)



fast_matrix.meta %>%
  mutate(samey = year2==year1) %>%
    ggplot(aes(
    x=samey,
    y=fst
  )) + 
  geom_boxplot()  ->
  fst.box
ggsave(fst.box, file = "fst.box.pdf",
       w= 6, h = 5)

######


ky.tag <- grep("KY", auto_dat@poolnames)
va.tag <- grep("VA", auto_dat@poolnames)


####
ky.va.set <- pooldata.subset(
  auto_dat,
  pool.index = c(ky.tag.flt,va.tag),
  min.cov.per.pool = 10,
  max.cov.per.pool = 150,
  min.maf = 0.05,
  return.snp.idx = TRUE,
  verbose = TRUE
)

ky.va.set@refallele.readcount %>% colMeans()

ky.va.set@snp.info %>%
  as.data.frame() %>%
  mutate(snp.id = rownames(.)) ->
  snp.dt

pca.ran.va.ky =
  randomallele.pca(
    ky.va.set,
    scale = FALSE
  )

pca.ran.va.ky$pop.loadings %>%
  as.data.frame() %>%
  mutate(sampleId_orig = rownames(.) ) %>%
  left_join(samps) %>%
  mutate(Year = year(Collection_date)) ->
  PCA.vaky.coords


PCA.vaky.coords %>%
  ggplot(aes(
    x=V1,
    y=V2,
    color= as.character(city),
    #shape=city
  )) + 
  geom_point(size = 5) ->
  PCA12.va
ggsave(PCA12.va, file = "PCA12.va.pdf",
       w= 6, h = 5)
       
       
       
save(fast_matrix.meta,
file = paste(root, "fast_matrix_fst", ".Rdata", sep = ""))