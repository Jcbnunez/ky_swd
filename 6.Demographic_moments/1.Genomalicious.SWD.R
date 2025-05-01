#### Moments admixture with Nucella -- for trios
# module load Rgeospatial
library(genomalicious)
library(SeqArray)
library(foreach)
library(sp)
library(poolfstat)
##pairsF="VA_KY_pairs.txt"
samps <- fread("/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023/2024_Joaquin_KY/swd.metadata.txt")

args = commandArgs(trailingOnly=TRUE)
i=as.numeric(args[1])
pairsF=args[2]

pairs = fread(pairsF)
#1.Create output the directory
system(paste("mkdir", "dadi_objects_duos", sep = " "))
system(paste("mkdir", "L_meta_objects_duos", sep = " "))

#2. Load the pairs data
SAMPLE_INPUT=c(pairs$pop1[i],pairs$pop2[i])
#

#### Read in GDS file
dat_f <- "/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023/2024_Joaquin_KY/auto.pooldata.new.rds"    
poolobj <- readRDS(dat_f)

sample.names <- poolobj@poolnames

duo.set <-
  pooldata.subset(
    poolobj,
    pool.index = which(sample.names %in% SAMPLE_INPUT),
    min.cov.per.pool = 10,
    max.cov.per.pool = 150,
    min.maf = 0.01,
    return.snp.idx = TRUE,
    verbose = TRUE
  )

### get allele frequency data
ad <- duo.set@refallele.readcount
dp <- duo.set@readcoverage  
dat <- ad/dp

### parse
dim(dat)
colnames(dat) <- duo.set@poolnames

snp.dt = duo.set@snp.info
snp.dt %>% mutate(snp_id = rownames(.)) -> snp.dt
####
### convert to format for genomalicious
tf <- TRUE
dat <- as.data.table(dat)
dat[,locus:=snp.dt$snp_id[tf]]
dat[,ref:=snp.dt$RefAllele[tf]]
dat[,alt:=snp.dt$AltAllele[tf]]

datl <- melt(dat, id.vars=c("locus", "ref", "alt"))
setnames(datl, c("variable", "value"), c("POOL", "FREQ"))
######

### Calculate NEff
dp = data.frame(dp)
names(dp) = duo.set@poolnames
colMeans(dp, na.rm = T) %>%
  data.frame(DP=.) %>%
  mutate( sampleId_orig=rownames(.)) ->
  dp.df
####

samps %>% 
  filter(sampleId_orig %in% SAMPLE_INPUT) ->
  meta

####
dp.df %>%
  left_join(meta) %>%
  group_by(sampleId_orig) %>%
  mutate(ne = floor((nFlies*DP)/(nFlies+DP-1)) )->neff

datl <- merge(datl, neff[,c("sampleId_orig", "ne")], by.x="POOL", by.y="sampleId_orig")
datl[,FREQ:=1-FREQ]
datl <- na.omit(datl)

datl %>%
  filter(POOL == pairs$pop1[i]) %>%
  mutate(POOL = "A_pop1") ->
  datl1

datl %>%
  filter(POOL == pairs$pop2[i]) %>%
  mutate(POOL = "B_pop2") ->
  datl2


rbind(datl1, datl2) -> datl.relab

sfs_method="probs"

dadi <- dadi_inputs(
  datl.relab,
  type = "freqs",
  popCol = "POOL",
  locusCol = "locus",
  refCol = "ref",
  altCol = "alt",
  freqCol = "FREQ",
  indsCol = "ne",
  freqMethod = sfs_method
)
dadi <- na.omit(dadi) 

L=dim(dadi)[1]
####
### Write dadi file

fn <- paste("dadi_objects_duos/",
            sfs_method, ".",
            pairs[i,]$pop1, ".",
            pairs[i,]$pop2, ".",
            "delim", sep="")

write.table(dadi,
            file=fn,
            sep="\t", quote=F, row.names=F)

### Metadata
neff %>%
  filter(sampleId_orig == pairs[i,]$pop1) %>% .$ne -> an1_ne
neff %>%
  filter(sampleId_orig == pairs[i,]$pop2) %>% .$ne -> an2_ne

fn2 <- paste("L_meta_objects_duos/",
            sfs_method, ".",
            pairs[i,]$pop1, ".",
            pairs[i,]$pop2, ".",
            "meta", sep="")

write.table(data.frame(an1_ne, an2_ne, L),
            file=fn2,
            sep="\t", quote=F, row.names=F)

