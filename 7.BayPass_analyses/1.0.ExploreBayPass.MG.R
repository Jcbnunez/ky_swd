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

### load cache from MG analyses
setwd("2024_Joaquin_KY/")
require(knitr)
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
ct.c2.ls.xi2=compute.local.scores(tmp.pos[,1:2],snp.pi=tmp.all.pi,snp.pvalue = tmp.all.c2,xi=2,manplot = T,main="xi=2")

### load GEA T>32 data
tmp.x.cov.idx=ct.beta.X$COVARIABLE==2 ; tmp.a.cov.idx=ct.beta.A$COVARIABLE==2
tmp.all.bf=c(ct.beta.X$BF.dB.[tmp.x.cov.idx],ct.beta.A$BF.dB.[tmp.a.cov.idx])

tmp.all.pi=c(ct.pixtx.X$M_P,ct.pixtx.A$M_P)
tmp.pos=rbind(read.table(gzfile("baypass_files/X.ky.new.ct.snpdet.gz")),read.table(gzfile("baypass_files/A.ky.new.ct.snpdet.gz")))
ct.bf.t32.ls.xi1=compute.local.scores(tmp.pos[,1:2],snp.pi=tmp.all.pi,snp.bf=tmp.all.bf,xi=1,manplot = T,main="xi=1")
