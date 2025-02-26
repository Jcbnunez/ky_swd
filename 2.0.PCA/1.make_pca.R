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
  min.maf = 0.05,
  verbose = TRUE
)

pca_result <- randomallele.pca(subset, scale = TRUE)

pca_result$pop.loadings %>% as.data.frame() %>%
  mutate(sampleId_orig=row.names(.)) %>%
  left_join(samps)-> pc_loadings
####
pc_loadings %>%
  ggplot(aes(
    x=V1,
    y=V2,
    fill=continent
  )) + geom_point(size = 4, shape = 21) +
  theme_bw()->
  PCA_plot
ggsave(PCA_plot, file = "PCA12.dim.png", w = 5, h =4)

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

ggplot(data = world) +
  geom_sf(fill= "antiquewhite") +
  coord_sf(xlim = c(-135.15, 160.99), ylim = c(-55.00, 69.00), expand = FALSE) + 
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.2), 
        panel.background = element_rect(fill = "aliceblue")) +
  geom_point(data = pc_loadings, 
             aes(x=long, y = lat, fill = V1), size = 2.2, shape = 21) +
  scale_fill_gradient2()-> SAMP.PC1.map
ggsave(SAMP.PC1.map, file = "SAMP.PC1.map.png", h = 4, w = 6)

###
samps %>%
  filter(continent == "North_America") %>%
  filter(sampleId_orig != "KY20") %>%
  filter(sampleId_orig != "KY11") %>%
  .$sampleId_orig -> NAM_samps

NAME_ids = which(geno_auto@poolnames %in% NAM_samps)

pool.index = NAME_ids

subset_NAme <- pooldata.subset(
  geno_auto,
  pool.index = NAME_ids,
  min.cov.per.pool = -1,
  max.cov.per.pool = 1e+06,
  min.maf = 0.05,
  verbose = TRUE
)
pca_result_NAme <- randomallele.pca(subset_NAme, scale = TRUE)

pca_result_NAme$pop.loadings %>% as.data.frame() %>%
  mutate(sampleId_orig=row.names(.)) %>%
  left_join(samps)-> pc_loadings_NAme

###
pc_loadings_NAme %>%
  ggplot(aes(
    x=V1,
    y=V2,
    fill=province
  )) + geom_point(size = 4, shape = 21) +
  theme_bw()->
  PCA_plot_name
ggsave(PCA_plot_name, file = "PCA_plot_name.png", w = 5, h =4)

samps %>%
  filter(continent == "North_America") %>%
  filter(sampleId_orig != "KY20") %>%
  filter(sampleId_orig != "KY11") %>%
  .$sampleId_orig -> NAM_samps

NAME_ids = which(geno_auto@poolnames %in% NAM_samps)

pool.index = NAME_ids

subset_NAme <- pooldata.subset(
  geno_auto,
  pool.index = NAME_ids,
  min.cov.per.pool = -1,
  max.cov.per.pool = 1e+06,
  min.maf = 0.05,
  verbose = TRUE
)
pca_result_NAme <- randomallele.pca(subset_NAme, scale = TRUE)

pca_result_NAme$pop.loadings %>% as.data.frame() %>%
  mutate(sampleId_orig=row.names(.)) %>%
  left_join(samps)-> pc_loadings_NAme

###
pc_loadings_NAme %>%
  ggplot(aes(
    x=V3,
    y=V4,
    fill=province
  )) + geom_point(size = 4, shape = 21) +
  theme_bw()->
  PCA_plot_name23
ggsave(PCA_plot_name23, file = "PCA23_plot_name.png", w = 5, h =4)

###
samps %>%
  filter(province == "Kentucky") %>%
  filter(sampleId_orig != "KY20") %>%
  filter(sampleId_orig != "KY11") %>% ##samp 11 behaves weirdly
  .$sampleId_orig -> KY_samps

KY_ids = which(geno_auto@poolnames %in% KY_samps)

subset_KY <- pooldata.subset(
  geno_auto,
  pool.index = KY_ids,
  min.cov.per.pool = -1,
  max.cov.per.pool = 1e+06,
  min.maf = 0.05,
  verbose = TRUE
)
pca_result_KY <- randomallele.pca(subset_KY, scale = TRUE)

pca_result_KY$pop.loadings %>% as.data.frame() %>%
  mutate(sampleId_orig=row.names(.)) %>%
  left_join(samps)-> pc_loadings_KY

pc_loadings_KY %>%
  ggplot(aes(
    x=year,
    y=V6,
    color = city
  )) + geom_smooth(method = "lm", se = F) + 
  geom_point(size = 4, shape = 21) + 
  theme_bw()->
  KY_PCA_year
ggsave(KY_PCA_year, file = "KY_PCA_year.png", w = 5, h =4)

pc_loadings_KY %>%
  ggplot(aes(
    x=year,
    y=V7,
    color = city
  )) + geom_smooth(method = "lm", se = F) + 
  geom_point(size = 4, shape = 21) + 
  theme_bw()->
  KY_PCA_year
ggsave(KY_PCA_year, file = "KY_PCA_year.png", w = 5, h =4)


cor.test(~V5+year, data = pc_loadings_KY)
cor.test(~V6+year, data = pc_loadings_KY)
cor.test(~V7+year, data = pc_loadings_KY)




####
genofile <- seqOpen("/netfiles/nunezlab/D_suzukii_resources/vcfs_gds/SWD_KY.PoolSeq.PoolSNP.001.5.ORCC.gds", allow.duplicate=T)

#####
seqResetFilter(genofile)
seqSetFilter(genofile, sample.id=samps$SequencingId)
snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile, verbose = T))

snps.dt <- snps.dt[nAlleles==2][missing == 0][chr %in% c("2L","2R","3")]
snps.dt %<>% mutate(SNP_id = paste(chr, pos, sep = "_")) 
snps.dt %>% dim


#####
seqResetFilter(genofile)
seqSetFilter(genofile, 
variant.id=snps.dt$variant.id)

ad <- seqGetData(genofile, "annotation/format/AD")$data
dp <- seqGetData(genofile, "annotation/format/DP")

dp[,which(colMeans(dp) > 45)] -> dp.f
ad[,which(colMeans(dp) > 45)] -> ad.f

sampleids <- seqGetData(genofile, "sample.id")
### create the data object
dat = ad.f/dp.f
dim(dat)

dat[,which(colMeans(dat) > 0.1)] -> dat.maf
dim(dat.maf)


#colnames(dat) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") , 
#                       sep="_")
rownames(dat.maf) <- seqGetData(genofile, "sample.id")

samps %>%
  filter(fly_type == "wild") %>%
  .$SequencingId -> samps.select

### PCA
dat.maf %>% 
  as.data.frame() %>%
  filter(rownames(.) %in% samps.select) %>% 
  PCA(graph = F, ncp = 20) ->
  pca.object

#### Plot

pca.object$ind$coord %>%
  as.data.frame() %>% 
  mutate(SequencingId = rownames(.)) %>%
  left_join(samps) ->
  pca.meta.dim

cor.test(pca.meta.dim$Dim.1, pca.meta.dim$year)
cor.test(pca.meta.dim$Dim.2, pca.meta.dim$year)
cor.test(pca.meta.dim$Dim.3, pca.meta.dim$year)
cor.test(pca.meta.dim$Dim.4, pca.meta.dim$year)
cor.test(pca.meta.dim$Dim.5, pca.meta.dim$year)
cor.test(pca.meta.dim$Dim.6, pca.meta.dim$year)
cor.test(pca.meta.dim$Dim.7, pca.meta.dim$year)
cor.test(pca.meta.dim$Dim.8, pca.meta.dim$year)
cor.test(pca.meta.dim$Dim.9, pca.meta.dim$year)

anova(lm(Dim.8 ~ fruit_type, data = pca.meta.dim))
anova(lm(Dim.9 ~ fruit_type, data = pca.meta.dim))
anova(lm(Dim.10 ~ fruit_type, data = pca.meta.dim))

pca.meta.dim %>%
  ggplot(aes(
    x=Dim.7,
    y=Dim.6,
    #shape = as.character(year),
    color = as.character(year)
  )) +
  geom_point(size = 3) +
  theme_bw() ->
  PCA12.dim

ggsave(PCA12.dim, file = "PCA12.dim.pdf", w = 5, h =4)


pca.meta.dim %>%
  dplyr::select(Dim.6, Dim.7, year) %>%
  melt(id="year") %>%
  ggplot(aes(
    x=year,
    y=value,
    #shape = as.character(year),
    color = as.character(variable)
  )) +
  geom_point(size = 3) +
  geom_smooth(method = "lm") + 
  theme_bw() + facet_wrap(~variable) ->
  PCA12.dim.linear

ggsave(PCA12.dim.linear, file = "PCA12.dim.linear.pdf", 
       w = 7, h =3.5)


pca.meta.dim %>%
  ggplot(aes(
    x=Dim.12,
    y=Dim.9,
    #shape = as.character(year),
    color = as.character(fruit_type)
  )) +
  geom_point(size = 3) +
  theme_bw() ->
  PCA89.dim

ggsave(PCA89.dim, file = "PCA89.dim.pdf", w = 5, h =4)
