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


ext_ens <- WS276_t2g_all %>% 
  dplyr::filter(ens_gene %in% pcoding$ens_gene) %>% 
  dplyr::filter(biotype=="mRNA") %>% 
  dplyr::select(ens_gene,ext_gene) %>% 
  distinct()

# sample infor

dir <- ("../raw_data/kallisto")

sample_id <- dir(file.path(dir))


sample_condition <- sub("(kallisto_)(.*)(_NIDDK_)(.*)(_)(.*)(_)(.*)","\\2",sample_id)

sample_sync <- sub("(kallisto_)(.*)(_NIDDK_)(.*)(_)(.*)(_)(.*)","\\8",sample_id)

sample_name <- sub("(kallisto_)(.*)(_)(NIDDK_)(.*)(_)(.*)(_)(.*)","\\2\\3\\5\\7",sample_id)


samples <- data.frame(sample = sample_id,
                      condition = sample_condition,
                      sync_meth = sample_sync,
                      samplename = sample_name,
                      filepath = dir) 





##### Pairwise Wald test  ####



pw_wald_list = list()


comp <- c("BRC20067_vs_N2","CB4856_vs_N2","DL238_vs_N2","BRC20067_vs_CB4856","BRC20067_vs_DL238","CB4856_vs_DL238")

for (cp_name in comp){
  
  #cp_name<- "DL238_vs_N2"
  
  cp <- c(sub("(.*)(_vs_)(.*)","\\1",cp_name),sub("(.*)(_vs_)(.*)","\\3",cp_name))
  
  
  samples_sub <- samples %>%
    dplyr::filter(condition %in% cp)
  
  files_sub <- file.path(samples_sub$filepath, samples_sub$sample,"abundance.h5")
  names(files_sub) <- samples_sub$samplename
  
  txi.kallisto_sub <- tximport(files_sub, type = "kallisto", tx2gene = tx2gene)
  
  
  ################### PAIRWISE wald 
  
  
  
  dds <- DESeqDataSetFromTximport(txi.kallisto_sub,
                                  colData = samples_sub,
                                  design = ~ condition)
  
  
  
  #a minimal pre-filtering to keep only rows that have at least 120 reads total.
  keep <- rowSums(counts(dds)) >= 120
  
  dds <- dds[keep,]
  
  # set reference
  cref=sub("(.*)(_vs_)(.*)","\\3",cp_name)
  
  dds$condition <- relevel(dds$condition, ref = cref)
  
  dds <- DESeq(dds)
  
  res  <- results(dds) 
  
  res_df <- res %>%
    data.frame() %>%
    rownames_to_column(var="ens_gene") %>%
    as_tibble() %>%
    dplyr::left_join(ext_ens) %>%
    dplyr::mutate(comparison=cp_name) 
  
  
  
  
  pw_wald_list[[cp_name]] <- res_df
  
  
}


pw_wald <- bind_rows(pw_wald_list)


### Set thresholds for wald DEGs
padj.cutoff <- 0.05
lfc.cutoff <- log2(1.5)

wald_deg <- pw_wald %>%
  dplyr::mutate(differential_express=ifelse(padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff,"Yes","No")) %>%
  dplyr::select(comparison,everything())

#unique(wald_deg$comparison)

 write.csv(wald_deg,paste("../DEresults/file_2_wald_6pairwise_comparison_DEGs.csv",sep = ""), row.names = F, quote = FALSE)

