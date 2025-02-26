### SNP calling job making array -- SWD

library(tidyverse)
library(magrittr)
library(data.table)
library(foreach)

### SWD wins:
gen_Ls =
data.frame(
  chr = c("chrX", "chr2L", "chr2R", "chr3" , "chr4", "chrY"),
  L = c(
    32723399,
    26650000,
    24272511,
    93255888,
    2560000,
    13207119)
)


###
win.bp <- 1.2e5
step.bp <-win.bp+1

## prepare windows
wins <- foreach(chr.i=gen_Ls$chr,
                .combine="rbind", 
                .errorhandling="remove")%do%{
                  
                  tmp <- gen_Ls %>%
                    filter(chr == chr.i)
                  
                  data.table(chr=chr.i,
                             start=seq(from=1, to=tmp$L-win.bp, by=step.bp),
                             end=seq(from=1, to=tmp$L-win.bp, by=step.bp) + win.bp)
                }

wins[,i:=1:dim(wins)[1]]

dim(wins)

write.table(wins, file = "/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/ky_swd/1.2.SNP_calling/SWD.jobs.txt", 
            append = FALSE, quote = FALSE, sep = ",",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")
