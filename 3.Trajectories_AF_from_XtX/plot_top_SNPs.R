#### Plot trajectories
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

### Metadata
meta_git <- "/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023/KY_2020_2023.SampleId_metadata.txt"
samps <- fread(meta_git)
setDT(samps)

### open GDS
genofile <- seqOpen("/netfiles/nunezlab/D_suzukii_resources/vcfs_gds/SWD_KY.PoolSeq.PoolSNP.001.5.ORCC.gds", allow.duplicate=T)

top_SNPs = c(
  "3_53887454",
  "X_26883655",
  "2R_1034795",
  "X_27374125",
  "3_63473325",
  "3_65201251",
  "3_66750625",
  "3_58775388",
  "3_2383254",
  "X_28948518"  
)

#####
seqResetFilter(genofile)
seqSetFilter(genofile, sample.id=samps$SequencingId)
snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile, verbose = T))

#snps.dt <- snps.dt[nAlleles==2][missing == 0][chr %in% c("2L","2R","3")]
snps.dt %<>% mutate(SNP_id = paste(chr, pos, sep = "_")) 
snps.dt %>% filter(SNP_id %in% top_SNPs) -> snps.dt.flt

snps.dt.flt %>% dim


#####
seqResetFilter(genofile)
seqSetFilter(genofile, 
             variant.id=snps.dt.flt$variant.id)

ad <- seqGetData(genofile, "annotation/format/AD")$data
dp <- seqGetData(genofile, "annotation/format/DP")

sampleids <- seqGetData(genofile, "sample.id")
### create the data object
dat = ad/dp
dim(dat)


#colnames(dat) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") , 
#                       sep="_")
rownames(dat) <- seqGetData(genofile, "sample.id")

samps %>%
  filter(fly_type == "wild") %>%
  .$SequencingId -> samps.select

####
dat %>% 
  as.data.frame() %>%
  filter(rownames(.) %in% samps.select) %>%
  mutate(SequencingId = rownames(.)) %>%
  left_join(samps) ->
  dat.meta

dat.meta %>% 
  ggplot(aes(
    x=exactDate,
    y=V1
  )) + geom_point() +
  geom_line(linetype = "dashed") +
  ggtitle("chr3:65201251", subtitle = "Downstream variant to Fer1 Gene" )+
  facet_grid(~year, scales = "free_x")->
  plotv1

ggsave(plotv1, file = "plotv1.pdf", w = 8, h=4)
  

dat.meta %>% 
  ggplot(aes(
    x=exactDate,
    y=V2
  )) + geom_point() +
  geom_line(linetype = "dashed") +
  facet_grid(~year, scales = "free_x")->
  plotv2

ggsave(plotv2, file = "plotv2.pdf", w = 8, h=4)
