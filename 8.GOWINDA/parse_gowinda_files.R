## make GO func file

library(tidyverse)
library(data.table)
library(magrittr)
library(reshape2)
library(foreach)

##this downloaded from NCBI + a small bit of manual curation
orig_go <- fread("/gpfs2/scratch/jcnunez/SWD_GO/GCF_043229965.1-RS_2025_01_gene_ontology.txt")


orig_go$GO_ID %>% unique() %>% sort() -> uniq_GOTERMS

parsed_go = 
foreach(i=uniq_GOTERMS,
        .combine = "rbind", .errorhandling = "remove")%do%{
          
          message(i)
          orig_go %>% filter(GO_ID == i)  -> tmp
          
          GId = i
          GOes = paste(tmp$GeneID, sep="", collapse=" ")
          
          data.frame(GO_ID = GId,
                     GENE_IDS = GOes) -> o
             
          return(o)
        }

GO_identifiers <- fread("basic_GO_terms.txt", 
                             sep = "\t", header = F)

GO_identifiers[,c(1,2)] -> go_terms
names(go_terms) = c("GO_ID", "GO_term")

left_join(parsed_go, go_terms) -> parsed_go_named
parsed_go_named[c("GO_ID","GO_term","GENE_IDS")] -> parsed_go_named

write.table(parsed_go_named, 
            file = "parsed_go_named.SWD.funcassione.txt", 
            append = FALSE, quote = F, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = F, qmethod = c("escape", "double"),
            fileEncoding = "")



