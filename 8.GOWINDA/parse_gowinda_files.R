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


#### step 2 ## Parse GTF file
raw_gtf <- fread("/netfiles/nunezlab/Shared_Resources/in_transit/for_kit/GCF_043229965.1_CBGP_Dsuzu_IsoJpt1.0_genomic.gtf")

raw_gtf %>%
  separate(V9, into = c("gene_id","transcript_id", "db_xref"), sep = ";") %>%
  separate(db_xref, into = c("db","real_gene_id"), sep = "\\:") ->
  raw_gtf.s

raw_gtf.s$real_gene_id = gsub('"', '', raw_gtf.s$real_gene_id )

#2L	CM091707.1	NC_092080.1	26,650,000	40	0	
#2R	CM091708.1	NC_092081.1	24,272,511	41.5	0	
#3	CM091709.1	NC_092082.1	93,255,888	40.5	0	
#4	CM091710.1	NC_092083.1	2,560,000	37	0	
#X	CM091711.1	NC_092084.1	32,723,399	40	1	
#Y	CM091712.1	NC_092085.1	13,207,119	42.5	1	
names(raw_gtf.s)[1] = "V1"
raw_gtf.s %>%
  mutate(chr_real = case_when(
    V1 == "NC_092080.1" ~ "chr2L",
    V1 == "NC_092081.1" ~ "chr2R",
    V1 == "NC_092082.1" ~ "chr3",
    V1 == "NC_092083.1" ~ "chr4",
    V1 == "NC_092084.1" ~ "chrX",
    V1 == "NC_092085.1" ~ "chrY"
  )) %>%
  mutate(id1 = 
           case_when(
             V3 == "gene" ~
           paste("gene_id ", '"' ,real_gene_id,'"',
                     "; gene_symbol ", '"LOC' ,real_gene_id,'"',";",
                     sep = ""),
           V3 != "gene" ~
             paste("gene_id ", '"' ,real_gene_id,'"',
                   "; gene_symbol ", '"LOC' ,real_gene_id,'"',
                   "; transcript_id ", '"LOC' ,real_gene_id,'"',
                   "; transcript_symbol ", '"LOC' ,real_gene_id,'"',";",
                   sep = ""),
           ))-> raw_gtf.s.chr

raw_gtf.s.chr %>%
  select(chr_real, 
         V2,
         V3,
         V4, V5, V6, V7, V8,
         id1
         ) -> gtf_gowind


write.table(gtf_gowind, 
            file = "gtf_gowind.swd.gtf", 
            append = FALSE, quote = F, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = F, qmethod = c("escape", "double"),
            fileEncoding = "")


#### Create sets
load("/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023_BAYPASS_MG/2024_Joaquin_KY/C2_df.Rdata")
write.table(C2_df[,c("chr","pos")], 
            file = "c2_universe_snps.txt", 
            append = FALSE, quote = F, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = F, qmethod = c("escape", "double"),
            fileEncoding = "")


C2_df %>%
  filter(P_C2 > 2) ->
  C2_outs
write.table(C2_outs, 
            file = "c2_top1perc_snps.txt", 
            append = FALSE, quote = F, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = F, qmethod = c("escape", "double"),
            fileEncoding = "")



