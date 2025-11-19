
setwd("C:/Users/andre/Desktop/GitHub/suzukii_FST")

library(foreach)
library(data.table)
library(tidyverse)
library(psych)
library(poolfstat)
library(magrittr)

setwd("C:/Users/andre/Desktop/GitHub/suzukii_FST")

files <- list.files(pattern = "^AF\\.bwt.*\\.txt$", full.names = TRUE)
files <- files[1:8]

alldf <- foreach(file = files, .combine = bind_rows) %do% {
  
  tmpdf <- fread(file)
  
  colnames(tmpdf) <- c("rep","cycle","pos","id","freq","wt","season","y","tS","tW","N")
  
  tmpdf$source <- basename(file)
  
  return(tmpdf)  
}




alldf %>%
  filter(cycle != 1) %>%
  group_by(cycle, wt) %>%
  summarize(N = mean(N)) %>%
  ggplot(
    aes(
      x=cycle,
      y=N,
      color=as.character(wt)
    )) + geom_line() ->
  demogrpahic_N

ggsave(demogrpahic_N, file = "coarse_demogrpahic_N.pdf",
       w= 9, h = 4)

alldf %>%
  filter(cycle != 1) %>%
  group_by(wt, season, rep) %>%
  summarize(Ne = harmonic.mean(N)) %>%
  ggplot(
    aes(
      x=as.character(wt),
      y=Ne,
      color=season
    )) + geom_violin() ->
  seasonal_drops

ggsave(seasonal_drops, file = "coarse_seasonal_drops.pdf",
       w= 9, h = 4)

###
alldf %>%
  filter(cycle != 1) %>%
  group_by(wt, rep) %>%
  summarize(NeHarm = harmonic.mean(N)) %>%
  ggplot(
    aes(
      x=as.character(wt),
      y=NeHarm,
    )) + geom_violin() ->
  seasonal_drops_delta

ggsave(seasonal_drops_delta, file = "coarse_seasonal_drops_delta.pdf",
       w= 9, h = 4)

#####
#####
##### --- Create wide-form
#####
#####
alldf[, cycle_col := cycle]
mafdf <- alldf[freq >= 0.05, ]
mafdf <- mafdf[freq <= 0.95, ]

wide_df <- dcast(
  filter(mafdf, cycle_col != 1),
  rep + pos + id + wt ~ season+y+cycle_col,
  value.var = "freq",
  fill = 0
)

####
#r, w
alldf$wt %>% unique -> ws
alldf$rep %>% unique -> rs

outres = 
  foreach(r=rs, .combine="rbind")%do%{
    foreach(w=ws, .combine="rbind")%do%{
      
      wide_df %>%
        filter(wt == w) %>%
        filter(rep == r) %>%
        select(!c(rep, pos, id, wt)) ->
        tmp_slice
      
      COV = 1000
      num_rows=dim(tmp_slice)[1]
      num_cols=dim(tmp_slice)[2]
      pool_sizes =
        
        #### Create pool object
        cov_artificial = matrix(rep(COV, num_rows * num_cols),
                                nrow = num_rows,
                                ncol = num_cols)
      
      ### extract number of invidivials
      alldf %>%
        filter(wt == w) %>%
        filter(rep == r) %>%
        filter(cycle_col != 1) %>%
        group_by(season, cycle_col, y) %>%
        summarize(N = mean(N)) %>%
        .$N -> pool_sizes
      
      pool <- new("pooldata",
                  npools=dim(tmp_slice)[2], #### Rows = Number of pools
                  nsnp=dim(tmp_slice)[1], ### Columns = Number of SNPs
                  refallele.readcount=as.matrix(ceiling(tmp_slice*COV)),
                  readcoverage=cov_artificial,
                  poolsizes=pool_sizes * 2,
                  poolnames = colnames(tmp_slice))
      
      
      fst.out <- compute.pairwiseFST(pool, method = "Anova")
      
      fst.out@values %>%
        as.data.frame() %>%
        mutate(comp = rownames(.)) %>%
        separate(comp, into = c("samp1","samp2"), sep = "\\;") %>%
        separate(samp1, into = c("season1","year1","cycle1"), sep = "_") %>%
        separate(samp2, into = c("season2","year2","cycle2"), sep = "_") %>% 
        mutate(deltay=abs(as.numeric(year1)-as.numeric(year2))) %>%
        mutate(deltaC=abs(as.numeric(cycle1)-as.numeric(cycle2))) %>%
        filter(deltay==0) %>% 
        group_by(paste(season1,season2)) %>%
        summarize(medFST = median(`Fst Estimate`),
                  meanFST = mean(`Fst Estimate`),
                  varFST = var(`Fst Estimate`),
        ) -> seas_fst
      
      fst.out@values %>%
        as.data.frame() %>%
        mutate(comp = rownames(.)) %>%
        separate(comp, into = c("samp1","samp2"), sep = "\\;") %>%
        separate(samp1, into = c("season1","year1","cycle1"), sep = "_") %>%
        separate(samp2, into = c("season2","year2","cycle2"), sep = "_") %>% 
        mutate(deltay=abs(as.numeric(year1)-as.numeric(year2))) %>%
        mutate(deltaC=abs(as.numeric(cycle1)-as.numeric(cycle2))) %>%
        group_by(deltay) %>%
        summarize(medFST = median(`Fst Estimate`),
                  meanFST = mean(`Fst Estimate`),
                  varFST = var(`Fst Estimate`),
        ) -> deltay_fst
      
      reshape2::melt(seas_fst) -> mel1
      names(mel1)[1] = "comparison"
      mel1 %<>% mutate(type = "seasonal", r = r, w =w)
      reshape2::melt(deltay_fst, id = "deltay") -> mel2
      names(mel2)[1] = "comparison"
      mel2 %<>% mutate(type = "yearly", r = r, w =w)
      
      rbind(mel1, mel2) -> o
      
      return(o)
    }}


outres %>%
  filter(variable == "medFST") %>%
  ggplot(
    aes(
      x= comparison,
      y= value,
      color = as.character(w)
    )) + geom_boxplot() +
  facet_wrap(~type, scales = "free") ->
  boxplots_all

ggsave(boxplots_all, file = "coarse_boxplots_all.pdf",
       w= 9, h = 4)
