library(knitr)
library(zoo)
library(RColorBrewer)
library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
require(readODS)
require(poolfstat)
require(vioplot)
require(data.table)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(foreach)
library(forcats)
library(magrittr)


autocor = function(x) {
  abs(cor(x[-1], x[-length(x)]))
}
thresUnif = function(L, cor, xi, alpha = 0.05) {
  a.poly.coef = rbind(c(-5.5, 2.47, 2.04, 0.22), c(6.76, 
                                                   -4.16, -5.76, -4.08), c(-5.66, -1.82, 1.04, 1.16), 
                      c(-2.51, -4.58, -6.95, -9.16))
  b.poly.coef = rbind(c(-1.22, 0.37, 2.55, 3.45), c(3.17, 
                                                    2.14, -0.02, -0.98), c(-1.99, -2.35, -2.31, -2.33))
  if (sum(xi != 1:4) == 4) {
    print("xi must be equal to 1, 2, 3 or 4")
    thres = NULL
  }
  else {
    tmp.cor2 = cor^2
    a = log(L) + a.poly.coef[1, xi] * (cor^3) + a.poly.coef[2, 
                                                            xi] * (cor^2) + a.poly.coef[3, xi] * cor + a.poly.coef[4, 
                                                                                                                   xi]
    b = b.poly.coef[1, xi] * (cor^2) + b.poly.coef[2, 
                                                   xi] * cor + b.poly.coef[3, xi]
    thres = (log(-log(1 - alpha)) - a)/b
  }
  return(thres)
}

#### Annontations
#annots <- fread("/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023/jn.bamlist.new.freebayes.noclump.vep",
#header= T)

#annots %>% 
#filter(Allele %in% c("A","C","T","G")) %>%
#separate(Location, into = c("chr","pos"), sep = ":") ->
#annots.flt
#save(annots.flt, file = "annots.flt.Rdata")

# ---> setwd("2024_Joaquin_KY/")
load("annots.flt.Rdata")

### load cache from MG analyses
tmp=list.files("joaquin_ky_new_cache/html/")
tmp=tmp[grep("rdb",tmp)]
tmp=matrix(unlist(strsplit(tmp,split="\\.")),ncol=2,byrow=T)[,1]
for(i in tmp){
   lazyLoad(filebase = paste0("joaquin_ky_new_cache/html/",i))
}

### Load seasonal data
tmp.all.c2=c(ct.C2.X$log10.1.pval.,ct.C2.A$log10.1.pval.)
tmp.all.pi=c(ct.pixtx.X$M_P,ct.pixtx.A$M_P)
tmp.pos=rbind(read.table(gzfile("baypass_files/X.ky.new.ct.snpdet.gz")),read.table(gzfile("baypass_files/A.ky.new.ct.snpdet.gz")))
ct.c2.ls.xi2=compute.local.scores(tmp.pos[,1:2],snp.pi=tmp.all.pi,snp.pvalue = tmp.all.c2 ,xi=2,manplot = T,main="xi=2")
ct.c2.ls.xi2$significant.windows -> windows_of_interest.C2
ct.c2.ls.xi2$res.local.scores -> C2_genome_scan

setDT(C2_genome_scan)

C2_df = data.frame(tmp.pos[,1:2],snp.pi=tmp.all.pi,snp.pvalue = tmp.all.c2)
names(C2_df)[1:4] = c("chr","pos","freqC2", "P_C2")
setDT(C2_df)

### save objects
save(C2_df, file = "C2_df.Rdata") 

### Describe raw C2s
C2_df %>%
filter(P_C2 > 2) %>% dim
C2_df %>%
filter(P_C2 > 2) %>% 
group_by(chr) %>%
summarize(N=n()) -> sig1
C2_df %>%
group_by(chr) %>%
summarize(all=n()) -> all
full_join(sig1, all) %>%
mutate(frac = N/all*100)

### Add annotations for FET
annots.flt$pos = as.numeric(annots.flt$pos)

left_join(C2_df, annots.flt) ->
C2_df.annot

C2_df.annot %>%
  group_by(Consequence) %>%
  summarize(N=n()) %>%
  filter(N > 100000) %>%
  .$Consequence -> common_cats
        
## filter for common categories
C2_df.annot %>%
filter(
Consequence %in% common_cats
) ->
C2_df.annot.common

### loop to calculate global enrichment
global_enrich =
foreach(k = common_cats,
.combine = "rbind")%do%{

##grand total
#Grand_total = dim(C2_df.annot.common)[1]
sig_class = dim(filter(C2_df.annot.common, Consequence == k & P_C2 >= 2 ))[1]
nonsig_class = dim(filter(C2_df.annot.common, Consequence == k & P_C2 < 2 ))[1]
sig_nonclass = dim(filter(C2_df.annot.common, Consequence != k & P_C2 >= 2 ))[1]
nonsig_nonclass = dim(filter(C2_df.annot.common, Consequence != k & P_C2 < 2 ))[1]
#Grand_total-(sig_class+nonsig_class+sig_nonclass+nonsig_nonclass)

dat <- data.frame(
  "sig" = c(sig_class, nonsig_class),
  "nonsig" = c(sig_nonclass, nonsig_nonclass),
  row.names = c("class", "nonclass"),
  stringsAsFactors = FALSE
)
fet <- fisher.test(dat)

data.frame(
class = k,
test = "global",
OR = fet$estimate,
lci = fet$conf.int[1],
uci = fet$conf.int[2],
p = fet$p.value)

} ### here end global

### chr by chr
chr_enrich =
foreach(k = common_cats,
.combine = "rbind")%do%{
foreach(j = c("chr2R","chr2L","chr3","chr4","chrX"),
.combine = "rbind")%do%{

C2_df.annot.common %>%
filter(chr == j) ->
tmp.chr

##grand total
#Grand_total = dim(tmp.chr)[1]
sig_class = dim(filter(tmp.chr, Consequence == k & P_C2 >= 2 ))[1]
nonsig_class = dim(filter(tmp.chr, Consequence == k & P_C2 < 2 ))[1]
sig_nonclass = dim(filter(tmp.chr, Consequence != k & P_C2 >= 2 ))[1]
nonsig_nonclass = dim(filter(tmp.chr, Consequence != k & P_C2 < 2 ))[1]
#Grand_total-(sig_class+nonsig_class+sig_nonclass+nonsig_nonclass)

dat <- data.frame(
  "sig" = c(sig_class, nonsig_class),
  "nonsig" = c(sig_nonclass, nonsig_nonclass),
  row.names = c("class", "nonclass"),
  stringsAsFactors = FALSE
)
fet <- fisher.test(dat)

data.frame(
class = k,
test = j,
OR = fet$estimate,
lci = fet$conf.int[1],
uci = fet$conf.int[2],
p = fet$p.value)

}} ## here end chr enrich!

#### table
rbind(global_enrich, 
		chr_enrich) ->
		total_enrichment

save(total_enrichment, file = "total_enrichment.Rdata")

#### plot
total_enrichment%>% 
		filter(test == "global") %>%
		ggplot(
		
		aes(
		x=fct_reorder(class,OR),
		y=log2(OR),
		ymin = log2(lci),
		ymax = log2(uci))
		
		) +
		geom_point(size = 3,
		position=position_dodge(width=0.5)) +
		geom_errorbar(width = 0.2,
		position=position_dodge(width=0.5)) +
		geom_hline(yintercept = 0) + theme_bw() +
		scale_x_discrete(guide = guide_axis(angle = 90)) ->
		OR_plot_glob
		
		ggsave(OR_plot_glob, 
		file = "OR_plot_glob.pdf",
		w= 3.0, h= 4)

#####
C2_df.annot.common ->
tmp.chr

chrOnly_enrich =
foreach(j = c("chr2R","chr2L","chr3","chr4","chrX"),
.combine = "rbind")%do%{


##grand total
#Grand_total = dim(tmp.chr)[1]
sig_class = dim(filter(tmp.chr, chr == j & P_C2 >= 2 ))[1]
nonsig_class = dim(filter(tmp.chr, chr == j & P_C2 < 2 ))[1]
sig_nonclass = dim(filter(tmp.chr, chr != j & P_C2 >= 2 ))[1]
nonsig_nonclass = dim(filter(tmp.chr, chr != j & P_C2 < 2 ))[1]
#Grand_total-(sig_class+nonsig_class+sig_nonclass+nonsig_nonclass)

dat <- data.frame(
  "sig" = c(sig_class, nonsig_class),
  "nonsig" = c(sig_nonclass, nonsig_nonclass),
  row.names = c("class", "nonclass"),
  stringsAsFactors = FALSE
)
fet <- fisher.test(dat)

data.frame(
test = j,
OR = fet$estimate,
lci = fet$conf.int[1],
uci = fet$conf.int[2],
p = fet$p.value)

} ## here end chr enrich!


### Describe Lindsay Scores
ct.c2.ls.xi2$significant.windows %>%
group_by(chr) %>%
summarize(N=n()) 

### plot lindley scores
chrInf = C2_df[, .(L = .N, cor = autocor(P_C2)), chr]
chrInf[, `:=`(th01, thresUnif(L, cor, xi = 2, alpha = 0.05)),]

left_join(C2_genome_scan, chrInf) -> C2_genome_scan
save(C2_genome_scan, file = "C2_genome_scan.Rdata")


don <- C2_genome_scan %>% 
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(pos)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(C2_genome_scan, ., by=c("chr"="chr")) %>%
  # Add a cumulative position of each SNP
  arrange(chr, pos) %>%
  mutate( BPcum=pos+tot)

axisdf = don %>%
  group_by(chr) %>%
  summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

### Load BF T>32 data + CTmin
tmp.x.cov.idx=ct.beta.X$COVARIABLE==2 ; tmp.a.cov.idx=ct.beta.A$COVARIABLE==2
tmp.all.bf=c(ct.beta.X$BF.dB.[tmp.x.cov.idx],ct.beta.A$BF.dB.[tmp.a.cov.idx])
#hist(tmp.all.bf,breaks=10,ylab="Pvalue",main="Distribution of p-value associated with XtX* (unilateral and all SNPs)",freq=F)
#abline(h=1,col="red")
tmp.all.pi=c(ct.pixtx.X$M_P,ct.pixtx.A$M_P)
#tmp.pos=rbind(fread("baypass_files/X.ky.new.ct.snpdet.gz",data.table = F),fread("baypass_files/A.ky.new.ct.snpdet.gz",data.table = F))
tmp.pos=rbind(read.table(gzfile("baypass_files/X.ky.new.ct.snpdet.gz")),read.table(gzfile("baypass_files/A.ky.new.ct.snpdet.gz")))

T32BF_df = data.frame(tmp.pos[,1:2],snp.pi=tmp.all.pi,snp.bf=tmp.all.bf)
names(T32BF_df)[1:4] = c("chr","pos", "freqBF", "BF_T32")
setDT(T32BF_df)


left_join(C2_df, T32BF_df, by = c("chr","pos")) -> 
C2_BF_df
setDT(C2_BF_df)

C2_BF_df %>%
  filter(BF_T32 > 1) %>% dim

C2_BF_df %<>% mutate(snp_id = paste(chr, pos, sep = "_"))
C2_BF_df %>%
filter(snp_id %in%
c("chr3_63316068",
"chr3_67001738",
"chr3_65836197",
"chrX_11184269",
"chr2R_20147914"))

#chr      pos    freqC2     P_C2    freqBF     BF_T32         snp_id
#chr3 63316068 0.5664361 4.379690 0.5664361  3.3546233  chr3_63316068 (kra)
#chr3 65836197 0.7491729 3.515994 0.7491729  1.0273694  chr3_65836197 (Nlg1)
#chr3 67001738 0.3994053 4.715222 0.3994053  3.0855904  chr3_67001738 (SNF4Aγ)

#chr2:  chr2R 20147914 0.6753870 2.021395 0.6753870 -7.7791230 chr2R_20147914 (Dg)

##### BFs
C2_BF_df %>%
filter(BF_T32 > 1 & P_C2 > 2) %>% dim
C2_BF_df %>%
filter(BF_T32 > 1 & P_C2 > 2) %>% 
group_by(chr) %>%
summarize(N=n()) -> sigbf
C2_BF_df %>%
group_by(chr) %>%
summarize(all=n()) -> allbg
full_join(sigbf, allbg) %>%
mutate(frac = N/all*100)

annots.flt$pos = as.numeric(annots.flt$pos)

C2_BF_df2=C2_BF_df
names(C2_BF_df2)[1:2] = c("chr","pos")
don %>%
left_join(C2_BF_df2) %>%
left_join(annots.flt) ->
don_BF

ggplot() +
    # Show all points
    geom_line(data = don_BF, 
    aes(x=pos/1e6, y=log10(lindley+1), color=as.factor(chr)), 
    alpha=0.8, size=0.5) +
    geom_line(data = don_BF,aes(x=pos/1e6, y = log10(th01+1))) +
    geom_point(data = slice_max(filter(group_by(don_BF, Feature), BF_T32 > 1 & P_C2 > 2 & lindley >  th01), BF_T32), 
    aes(x=pos/1e6, y=log10(lindley+1)), 
    alpha=0.8, size=1.2, color = "red") +
    scale_color_manual(values = rep(c("skyblue", "grey"), 22 )) +
    # custom X axis:
    scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
    # Custom the theme:
    theme_bw() +
    geom_text(data = slice_max(filter(group_by(don_BF, Feature), BF_T32 > 1 & P_C2 > 2 & lindley >  th01), BF_T32), 
    size=2, aes(label = Feature, x=pos/1e6, y=log10(lindley+1) )) + 
    facet_grid(~chr, scales = "free_x", space = "free") +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) -> coloc_c2BF_manhat
    
ggsave(coloc_c2BF_manhat, file = "coloc_c2BF_manhat.pdf", 
w = 6, h = 3)


filter(don_BF, BF_T32 >= 1 & P_C2 >= 2 & lindley >=  th01) %>%
left_join(annots.flt) %>%
mutate(snp_id = paste(chr, pos, sep = "_")) ->
BF_C2_hits

save(BF_C2_hits, file = "BF_C2hits.Rdata")

BF_C2_hits %>%
group_by(chr) %>%
  summarize(N=n())

BF_C2_hits %>%
  group_by(Consequence) %>%
  summarize(N=n())

BF_C2_hits %>%
  group_by(Feature) %>%
  summarize(N=n())

###conseervative outliers
tmp.chr2=don_BF
TOPHITS_chrOnly_enrich =
  foreach(j = c("chr2R","chr2L","chr3","chr4","chrX"),
          .combine = "rbind")%do%{
            
            
            ##grand total
            #Grand_total = dim(tmp.chr)[1]
            sig_class = dim(filter(tmp.chr2, chr == j & BF_T32 >= 1 & P_C2 >= 2 & lindley >=  th01))[1]
            nonsig_class = dim(filter(tmp.chr2, chr == j & BF_T32 <= 1 & P_C2 >= 2 & lindley >=  th01 ))[1]
            sig_nonclass = dim(filter(tmp.chr2, chr != j & BF_T32 >= 1 & P_C2 >= 2 & lindley >=  th01 ))[1]
            nonsig_nonclass = dim(filter(tmp.chr2, chr != j & BF_T32 <= 1 & P_C2 >= 2 & lindley >=  th01 ))[1]
            #Grand_total-(sig_class+nonsig_class+sig_nonclass+nonsig_nonclass)
            
            dat <- data.frame(
              "sig" = c(sig_class, nonsig_class),
              "nonsig" = c(sig_nonclass, nonsig_nonclass),
              row.names = c("class", "nonclass"),
              stringsAsFactors = FALSE
            )
            fet <- fisher.test(dat)
            
            data.frame(
              test = j,
              OR = fet$estimate,
              lci = fet$conf.int[1],
              uci = fet$conf.int[2],
              p = fet$p.value)
            
          } ## here end chr enrich!
####
##### Tracking analysis
samps <- fread("/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023/2024_Joaquin_KY/swd.metadata.txt")

dat_f <- "all.pooldata.new.rds"    
all_dat <- readRDS(dat_f)

sample.names <- all_dat@poolnames
exclude <- c("KY20")#,"KY17", "KY11")

ky.tag <- grep("KY", sample.names)

indexer.ky <- data.frame(name = sample.names[grep("KY", sample.names)],
           index = ky.tag)   
ky.tag.flt <- indexer.ky$index[which(!indexer.ky$name %in% exclude)]

ky.set <-
	pooldata.subset(
	all_dat,
	pool.index = ky.tag.flt,
	min.cov.per.pool = 10,
	max.cov.per.pool = 150,
	min.maf = 0.05,
	return.snp.idx = TRUE,
	verbose = TRUE
	)

ky.set@snp.info %>%
as.data.frame() %>%
mutate(rs.id = rownames(.)) ->
snp.info.annot
names(snp.info.annot)[1:2] = c("chr","pos")
snp.info.annot %<>%
mutate(snp_id = paste(chr, pos, sep = "_"))

snps_in_C2_BF =
snp.info.annot %>%
filter(snp_id %in% BF_C2_hits$snp_id)

####
snps_in_C2_BF %>%
separate(remove = F, rs.id, into = c("feat", "index"),
sep = "s") ->
snps_in_C2_BF.df


### using now the info to extract SNPs
ky.set.topSNPs.C2BF <-
	pooldata.subset(
	all_dat,
	pool.index = ky.tag.flt,
	snp.index = as.numeric(snps_in_C2_BF.df$index),
	min.cov.per.pool = 10,
	max.cov.per.pool = 150,
	min.maf = 0.05,
	return.snp.idx = TRUE,
	verbose = TRUE
	)

ky.set.topSNPs.C2BF@snp.info %>%
as.data.frame() %>%
mutate(snp_id = paste(Chromosome, Position,sep = "_")) ->
snp_info.C2BF

ref_count.C2BF <- ky.set.topSNPs.C2BF@refallele.readcount
coverage.C2BF <- ky.set.topSNPs.C2BF@readcoverage
afs.C2BF <- ref_count.C2BF/coverage.C2BF

afs.C2BF %>%
as.data.frame %>%
mutate(snp_id = snp_info.C2BF$snp_id ) ->
afs.C2BF.id

names(afs.C2BF.id) = c(ky.set.topSNPs.C2BF@poolnames, "snp_id")

afs.C2BF.id %<>%
melt(id = "snp_id", variable.name = "sampleId_orig") %>%
left_join(samps) %>%
separate(remove = F, snp_id,
into = c("chr", "pos"),
sep = "_") 

afs.C2BF.id$pos = as.numeric(afs.C2BF.id$pos)

traits <- fread("Means_CT_min_FINAL.csv")


afs.C2BF.id %<>%
left_join(BF_CT_hits) %>%
left_join(traits) %>%
mutate(
yday = yday(as.Date(Collection_date, format = "%Y-%m-%d")),
year = as.numeric(year(as.Date(Collection_date, format = "%Y-%m-%d")))
) 

### run correlation with julian day
gens = unique(afs.C2BF.id$Feature)
sites = unique(afs.C2BF.id$site)

yday_cors.C2BF =
foreach(i = gens, .combine = "rbind",
.errorhandling = "remove")%do%{
foreach(j = sites, .combine = "rbind",
.errorhandling = "remove")%do%{

#z=1 ;i = gens[z]; j = sites[z]; k = posi[z]

afs.C2BF.id %>%
filter(Feature == i) %>%
filter(site == j) -> tmp1

posi = unique(tmp1$pos)
foreach(k = posi, .combine = "rbind",
.errorhandling = "remove")%do%{
message(i,j,k)

tmp1 %>%
filter(pos == k) %>%
cor.test(~value+yday, data = .) -> o2

message(o2$estimate)

data.frame(
Feature = i,
site = j,
pos = k,
cor.yday = o2$estimate,
p.yday = o2$p.value)

}}}

yday_cors.C2BF %>%
dcast(Feature+pos~site, value.var = "p.yday") ->
yday_cors.C2BF.dcast
names(yday_cors.C2BF.dcast)[4] = "South"
yday_cors.C2BF.dcast %>%
filter(South < 0.05 & Berea < 0.05)


### True seasonal targets
         #Feature                pos       Berea       South
#2 XM_017070295.3 (kra) (3R) 63316068 0.006262193 0.003147274
#3 XM_036816734.2 (SNF4Aγ) (3R) 67001738 0.004902951 0.044552508
#4 XM_036819182.2 (Nlg1) (3R) 65836197 0.033690047 0.022154055
#5 XM_036820755.2 (mamo) (X) 11184269 0.007601839 0.013188811

yday_cors.C2BF %>%
  dcast(Feature+pos~site, value.var = "cor.ctmin") ->
  yday_cors.C2BF.dcast
names(yday_cors.C2BF.dcast)[4] = "South"
yday_cors.C2BF.dcast %>%
  filter(South > abs(0.5) & Berea > abs(0.5))

snp.info.annot %>%
  filter(snp_id %in% c("chr3_67001738",
                       "chr3_65836197",
                       "chr3_63316068"
                       ))


#### Spatial and Temporal C2 + XtX
#### Spatial and Temporal C2 + XtX#### Spatial and Temporal C2 + XtX#### Spatial and Temporal C2 + XtX
#### Spatial and Temporal C2 + XtX
#### Spatial and Temporal C2 + XtX
#### Spatial and Temporal C2 + XtX
#### Spatial and Temporal C2 + XtX#### Spatial and Temporal
### Load XtX data

tmp.all.xtx=c(ct.pixtx.X$log10.1.pval.,ct.pixtx.A$log10.1.pval.)
tmp.all.xtx[is.infinite(tmp.all.xtx)]=1.1*max(tmp.all.xtx[!is.infinite(tmp.all.xtx)])
tmp.all.pi=c(ct.pixtx.X$M_P,ct.pixtx.A$M_P)
tmp.pos=rbind(read.table(gzfile("baypass_files/X.ky.new.ct.snpdet.gz")),read.table(gzfile("baypass_files/A.ky.new.ct.snpdet.gz")))
ct.xtx.ls.xi2=compute.local.scores(tmp.pos[,1:2],snp.pi=tmp.all.pi,snp.pvalue = tmp.all.xtx,xi=2,manplot = T,main="xi=2")
ct.xtx.ls.xi2$significant.windows -> windows_of_interest.XTX
ct.xtx.ls.xi2$res.local.scores -> XTX_genome_scan


XtX_df = data.frame(tmp.pos[,1:2],snp.pi=tmp.all.pi,snp.pvalue = tmp.all.xtx)
names(XtX_df)[1:4] = c("chr","pos","freqXtX","P_XtX")
setDT(XtX_df)

save(XtX_df, file = "XtX_df.Rdata")

### Describe raw XtXs
XtX_df %>%
filter(snp.pvalue > 2) %>% dim
XtX_df %>%
filter(snp.pvalue > 2) %>% 
group_by(chr) %>%
summarize(N=n()) -> sig2
XtX_df %>%
group_by(chr) %>%
summarize(all=n()) -> all2
full_join(sig2, all2) %>%
mutate(frac = N/all*100)

### Describe Lindsay Scores
ct.xtx.ls.xi2$significant.windows %>%
group_by(chr) %>%
summarize(N=n()) 


### plot
chrInf2 = XtX_df[, .(L = .N, cor = autocor(snp.pvalue)), chr]
chrInf2[, `:=`(th01, thresUnif(L, cor, xi = 2, alpha = 0.05)),]

left_join(XTX_genome_scan, chrInf2) -> XTX_genome_scan

don2 <- XTX_genome_scan %>% 
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(pos)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(XTX_genome_scan, ., by=c("chr"="chr")) %>%
  # Add a cumulative position of each SNP
  arrange(chr, pos) %>%
  mutate( BPcum=pos+tot)

axisdf2 = don %>%
  group_by(chr) %>%
  summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

#now plot
### Co-Loc -- fooverlap
setDT(windows_of_interest.C2)
setDT(windows_of_interest.XTX)

setkey(windows_of_interest.C2, chr, beg, end)
setkey(windows_of_interest.XTX, chr, beg, end)

foverlaps(windows_of_interest.C2[,c("chr","beg","end")],
		  windows_of_interest.XTX[,c("chr","beg","end")]) %>%
		  .[complete.cases(.),] %>%
		  mutate(ymin = -2.8, ymax = 2.2) ->
		  win.overlaps

save(win.overlaps, file = "C2.XTX.win.overlaps.Rdata")

### Co-Loc!
### Co-Loc!
### Co-Loc!
### Co-Loc!
ggplot() +
    # Show all points
    geom_segment(data = win.overlaps,
    aes(x=beg/1e6, xend=end/1e6, 
    y = ymin, yend = ymax), color = "gold") +
    geom_line(data = don, aes(x=pos/1e6, y=log10(lindley+1), color=as.factor(chr)), 
    alpha=0.8, size=0.5) +
    geom_line(data = don,aes(x=pos/1e6, y = log10(th01+1))) +
    geom_line(data = don2, aes(x=pos/1e6, y=-log10(lindley+1), color=as.factor(chr)), 
    alpha=0.8, size=0.5) +
    geom_line(data = don2,aes(x=pos/1e6, y = -log10(th01+1))) +
    scale_color_manual(values = rep(c("skyblue", "grey"), 22 )) +
    # custom X axis:
    scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
    # Custom the theme:
    theme_bw() +
    ylim(-2.8, 2.2) +
    facet_grid(~chr, scales = "free_x", space = "free") +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) -> coloc_manhat
    
ggsave(coloc_manhat, file = "coloc_manhat.pdf", w = 5, h = 3)


#ct.bf.t32.ls.xi1=compute.local.scores(tmp.pos[,1:2],snp.pi=tmp.all.pi,snp.bf=tmp.all.bf,xi=1,manplot = T,main="xi=1")
left_join(C2_df, XtX_df, by = c("chr","pos")) %>%
mutate(beg=pos, end=pos)-> C2_xtx_df
setDT(C2_xtx_df)

setkey(win.overlaps, chr,beg,end)
foverlaps(C2_xtx_df,win.overlaps[,c("chr","beg","end")]) %>%
mutate(win_sig_id = paste(chr,beg,end, sep = "_")) ->
C2_xtx_df_peaksannot
C2_xtx_df_peaksannot$pos = as.character(C2_xtx_df_peaksannot$pos)

###
win.overlaps %>%
mutate(win_sig_id = paste(chr,beg,end, sep = "_"))->
loop_df

C2_xtx_df_peaksannot$pos = as.numeric(C2_xtx_df_peaksannot$pos)

map_genes =
foreach(i = loop_df$win_sig_id,
.combine = "rbind")%do%{

C2_xtx_df_peaksannot %>%
filter(win_sig_id == i) %>%
left_join(annots.flt) -> explore
explore$Feature %>% table -> o
data.frame(o) %>%
mutate(win = i)

}

save(map_genes, file = "map_genes.swd.Rdata")

map_genes$`.` %>% table

names(map_genes)[1] = "Feature"
map_genes %>%
  filter(Feature != "-") ->
  map_genes.xtx.x2.flt

#####
##### Tracking of alleles
ky.set@snp.info %>%
  as.data.frame() %>%
  mutate(rs.id = rownames(.)) ->
  snp.info.annot

names(snp.info.annot)[1:2] = c("chr","pos")
snp.info.annot %<>%
  mutate(snp_id = paste(chr, pos, sep = "_"))

C2_xtx_df %>%
  filter(P_C2 >= 2 & P_XtX >= 2) %>%
  left_join(annots.flt) %>%
  mutate(snp_id = paste(chr, pos, sep = "_") )->
  XtX_C2_hits
  
snps_in_C2_XtX =
  snp.info.annot %>%
  filter(snp_id %in% XtX_C2_hits$snp_id)

snps_in_C2_XtX %>%
  separate(remove = F, rs.id, into = c("feat", "index"),
           sep = "s") ->
  snps_in_C2_XtX.df

### using now the info to extract SNPs
ky.set.topSNPs.C2XtX <-
  pooldata.subset(
    all_dat,
    pool.index = ky.tag.flt,
    snp.index = as.numeric(snps_in_C2_XtX.df$index),
    min.cov.per.pool = 10,
    max.cov.per.pool = 150,
    min.maf = 0.05,
    return.snp.idx = TRUE,
    verbose = TRUE
  )

ky.set.topSNPs.C2XtX@snp.info %>%
  as.data.frame() %>%
  mutate(snp_id = paste(Chromosome, Position,sep = "_")) ->
  snp_info.C2XtX

ref_count.C2XtX <- ky.set.topSNPs.C2XtX@refallele.readcount
coverage.C2XtX <- ky.set.topSNPs.C2XtX@readcoverage
afs.C2XtX <- ref_count.C2XtX/coverage.C2XtX

afs.C2XtX %>%
  as.data.frame %>%
  mutate(snp_id = snp_info.C2XtX$snp_id ) ->
  afs.C2XtX.id

names(afs.C2XtX.id) = c(ky.set.topSNPs.C2XtX@poolnames, "snp_id")

afs.C2XtX.id %<>%
  melt(id = "snp_id", variable.name = "sampleId_orig") %>%
  left_join(samps) %>%
  separate(remove = F, snp_id,
           into = c("chr", "pos"),
           sep = "_") 

afs.C2XtX.id$pos = as.numeric(afs.C2XtX.id$pos)

traits <- fread("Means_CT_min_FINAL.csv")

afs.C2XtX.id %<>%
  left_join(XtX_C2_hits) %>%
  left_join(traits) %>%
  mutate(
    yday = yday(as.Date(Collection_date, format = "%Y-%m-%d")),
    year = as.numeric(year(as.Date(Collection_date, format = "%Y-%m-%d")))
  ) 

### run correlation with julian day for XtX
gens = unique(afs.C2XtX.id$Feature)
sites = unique(afs.C2XtX.id$site)

yday_cors.C2XtX =
  foreach(i = gens, .combine = "rbind",
          .errorhandling = "remove")%do%{
            foreach(j = sites, .combine = "rbind",
                    .errorhandling = "remove")%do%{
                      
                      #z=1 ;i = gens[z]; j = sites[z]; k = posi[z]
                      
                      afs.C2XtX.id %>%
                        filter(Feature == i) %>%
                        filter(site == j) -> tmp1
                      
                      posi = unique(tmp1$pos)
                      foreach(k = posi, .combine = "rbind",
                              .errorhandling = "remove")%do%{
                                message(i,j,k)
                                
                                tmp1 %>%
                                  filter(pos == k) %>%
                                  cor.test(~value+yday, data = .) -> o2
                                
                                message(o2$estimate)
                                
                                data.frame(
                                  chr= unique(tmp1$chr),
                                  pos = k,
                                  Feature = i,
                                  site = j,
                                  cor.yday = o2$estimate,
                                  p.yday = o2$p.value)
                                
                              }}}

yday_cors.C2XtX %>%
  left_join(map_genes.xtx.x2.flt) ->
  yday_cors.C2XtX.win

save(yday_cors.C2XtX.win, file = "XtX_C2_yday_correlations.Rdata")

yday_cors.C2XtX.win %>% filter(Feature == "XM_017084109.3")

yday_cors.C2XtX.win %>%
  dcast(Feature+chr+pos+win~site, value.var = "p.yday") ->
  yday_cors.C2XtX.dcast
names(yday_cors.C2XtX.dcast)[6] = "South"

yday_cors.C2XtX.dcast %>%
  filter(!is.na(win)) %>%
  filter(South < 0.05 & Berea < 0.05) 

yday_cors.C2XtX.dcast %>% filter(Feature == "XM_017073893.3")

yday_cors.C2XtX.win %>%
  dcast(Feature+chr+pos+win~site, value.var = "cor.yday") ->
  yday_cors.C2XtX.CORS.dcast
names(yday_cors.C2XtX.CORS.dcast)[6] = "South"

yday_cors.C2XtX.CORS.dcast %>% 
  filter(chr == "chr2R" & pos == 20147914)
#chr2R 20147914
###
yday_cors.C2XtX.CORS.dcast%>%
  group_by(chr,pos) %>%
  slice_head(n =1) %>%
  left_join(XtX_C2_hits) %>%
  left_join(T32BF_df) %>%
  filter(P_C2 >= 2 & P_XtX >= 2) ->
  yday_cors.C2XtX.CORS.dcast

tot = dim(yday_cors.C2XtX.CORS.dcast)[1]

yday_cors.C2XtX.CORS.dcast %>%
  mutate(sign =Berea*South) %>%
  mutate(sign_test = case_when(sign >= 0 ~ "concordant",
                               sign < 0 ~ "discordant")) ->
  yday_cors.C2XtX.CORS.dcast

yday_cors.C2XtX.CORS.dcast %>%
  group_by(sign_test) %>%
  summarize(N = n())

yday_cors.C2XtX.CORS.dcast %>%
  group_by(sign_test, !is.na(win) ) %>%
  summarize(N = n())

###
yday_cors.C2XtX.CORS.dcast %>%
  ggplot(aes(
    x=Berea, y =South,
    color = !is.na(win),
    size =!is.na(win),
  )) + geom_point() +
  scale_size_manual(values = c(0.5, 3.0)) +
  xlim(-1,1) + ylim(-1,1) + theme_bw() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) ->
  cors.XtXC2

ggsave(cors.XtXC2, file = "cors.XtXC2.pdf")

### with BF T32
yday_cors.C2XtX.CORS.dcast %>%
  ggplot(aes(
    x=Berea, y =South,
    color = BF_T32 > 1,
    size = BF_T32 > 1,
  )) + geom_point() +
  scale_size_manual(values = c(0.5, 3.0)) +
  xlim(-1,1) + ylim(-1,1) + theme_bw() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) ->
  cors.XtXC2.BF_T32

ggsave(cors.XtXC2.BF_T32, file = "cors.XtXC2.BF_T32.pdf")



