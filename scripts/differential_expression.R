setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/"))

library(tidyverse)
library(tximport)
library(DESeq2)



# live protein_coding genes
c_elegans.PRJNA13758.WS276.geneIDs  <- read.csv("../raw_data/c_elegans.PRJNA13758.WS276.geneIDs.txt.gz", header=FALSE, stringsAsFactors=FALSE)

pcoding <- c_elegans.PRJNA13758.WS276.geneIDs  %>% 
  dplyr::filter(V5=="Live" & V6=="protein_coding_gene") %>% 
  dplyr::mutate(ext_gene=ifelse(V3=="",V4,V3)) %>% 
  dplyr::select(ens_gene=V2,ext_gene)



# gene and transcript IDs from WormBase WS276
WS276_t2g_all <- read.delim("../raw_data/WS276_t2g_all.tsv", stringsAsFactors=FALSE)

tx2gene <- WS276_t2g_all %>% 
  dplyr::filter(ens_gene %in% pcoding$ens_gene) %>% 
  dplyr::filter(biotype=="mRNA") %>% 
  dplyr::select(transcript,ens_gene) %>% 
  distinct()

