
### load R packages####

library(tidyverse)
library(ggtree)
library(ape)
library(ggforce)
library(ggpubr)
library(gridExtra)
library(limma)
library(cowplot)



#####

##########################################
#           Figure 1                     #
#             tree                       #
##########################################

## Figure 1B #####

load("processed_data/S1B_File.RData")

strains_to_plot <- c("N2","CB4856","DL238","BRC20067")

tree_raw <- ggtree(as.phylo(as.hclust(treeExon2010_03)), layout="circular",branch.length="rate",size=0.2) + 
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), 
        legend.position="bottom", 
        legend.text = element_text(size=12),
        legend.title =  element_text(size=12),
        text=element_text(family="Arial"), 
        legend.spacing.x = unit(0.1, 'mm')) + 
  geom_tippoint(aes(subset=(label %in% strains_to_plot),  fill=label), shape=21, size=0.8) + 
  scale_fill_manual(values=c("BRC20067"="hotpink3", "CB4856"="blue","DL238"= "cadetblue3", "N2"="orange")) +
  labs(fill='Strains') 


fig_1B <- tree_raw %>% 
  flip(759, 760) %>% 
  flip(742, 743) %>% 
  flip(716, 134) %>% 
  flip(128, 703) %>% 
  flip(685, 5) %>% 
  flip(175, 659) %>% 
  flip(652, 653) %>% 
  flip(172, 639) %>% 
  flip(595, 596) %>% 
  flip(590, 591) %>% 
  flip(578, 579) %>% 
  flip(560, 561) %>% 
  flip(558, 559) %>% 
  flip(555, 556)


ggsave(fig_1B , filename = paste("figures/Fig_1B.svg",sep = ""),  units = "mm",height = 114, width = 114)



##########################################
##########################################
#           Figure 2                     #
#       Genetic variants                 #
##########################################

## Figure 2A #####

load("processed_data/S2A_File.RData")

## bcd ##
plt_pies_list = list()

for (st in c("BRC20067","CB4856","DL238"))
{
  df <- strain_variants_bcd  %>% 
    dplyr::filter(strain==st) %>% 
    as.data.frame() %>%
    dplyr::mutate(end = 2 * pi * cumsum(n)/sum(n),
                  start = lag(end, default = 0),
                  middle = 0.5 * (start + end),
                  hjust = ifelse(middle > pi, 1, 0),
                  vjust = ifelse(middle < pi/2 | middle > 3 * pi/2, 0, 1)) %>% 
    ungroup() %>% 
    dplyr::mutate(nn=round((n/1000),digits = 0)) %>% 
    dplyr::select(-n) %>% 
    dplyr::mutate(n=paste(nn,"K",sep = ""))
  
  df$Region = factor(df$region, levels = c("exon","intron","intergenic"))
  
  plt_pie <- ggplot(df) + 
    geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = 1,
                     start = start, end = end, fill = Region)) +
    geom_text(aes(x = 0.1 * sin(middle), y = 0.4 * cos(middle), label = n,
                  hjust = hjust, vjust = vjust), size=(5/14) *12,family="Arial") +
    coord_fixed() +
    scale_x_continuous(limits = c(-1, 1),  # Adjust so labels are not cut off
                       name = "", breaks = NULL, labels = NULL) +
    scale_y_continuous(limits = c(-1, 1 ),      # Adjust so labels are not cut off
                       name = "", breaks = NULL, labels = NULL) + 
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size=12),
          legend.text = element_text(size=12),
          legend.title =  element_text(size=12),plot.margin = unit(c(0, 0, 0, 0), "mm"),
          text=element_text(family="Arial")) +
    scale_fill_manual(values = c("gold2", "darkorange1", "lightskyblue2"), labels = c("Protein-coding exon","Protein-coding intron", "Intergenic")) + 
    ggtitle(st) + 
    labs(fill='Genetic regions')
  
  plt_pies_list[[st]] <- plt_pie
  
}

## n2 ##

df_N2 <- region_length  %>% 
  as.data.frame() %>%
  dplyr::rename(n=total_length,on=region) %>%
  arrange(on) %>% 
  dplyr::mutate(end = 2 * pi * cumsum(n)/sum(n),
                start = lag(end, default = 0),
                middle = 0.5 * (start + end),
                hjust = ifelse(middle > pi, 1, 0),
                vjust = ifelse(middle < pi/2 | middle > 3 * pi/2, 0, 1)) %>% 
  ungroup() %>% 
  dplyr::mutate(nn=round((n/1000),digits = 0)) %>% 
  dplyr::select(-n) %>% 
  dplyr::mutate(n=paste(nn,"K",sep = ""))

df_N2$region = factor(df_N2$on, levels = c("exon","intron","intergenic"))

plt_pie_N2 <- ggplot(df_N2) + 
  geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = 1,
                   start = start, end = end, fill = region)) +
  geom_text(aes(x = 0.1 * sin(middle), y = 0.4 * cos(middle), label = n,
                hjust = hjust, vjust = vjust), size=(5/14) *9,family="Arial") +
  coord_fixed() +
  scale_x_continuous(limits = c(-1, 1),  # Adjust so labels are not cut off
                     name = "", breaks = NULL, labels = NULL) +
  scale_y_continuous(limits = c(-1, 1 ),      # Adjust so labels are not cut off
                     name = "", breaks = NULL, labels = NULL) + 
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size=12),
        legend.text = element_text(size=12),plot.margin = unit(c(0,0, 0, 0), "mm"),
        legend.title =  element_blank(),
        text=element_text(family="Arial")) +
  scale_fill_manual(values = c("gold2", "darkorange1", "lightskyblue2"), labels = c("Exon","Intron", "Intergenic")) + 
  ggtitle("N2") + 
  labs(fill='Genetic regions') + 
  guides(fill = guide_legend(ncol = 3))



### cowplot Fig2A

legend <- get_legend(plt_pie_N2)

pies3 <- cowplot::plot_grid(plt_pie_N2+ theme(legend.position="none"),
                            plt_pies_list$BRC20067 + theme(legend.position="none"),
                            plt_pies_list$CB4856 + theme(legend.position="none"),
                            plt_pies_list$DL238 + theme(legend.position="none"), 
                            nrow = 2,ncol = 2)

Fig_2A <- cowplot::plot_grid( pies3,legend,
                                nrow = 2,ncol = 1,
                                rel_heights = c(4,1),
                                axis = "lr")

## Figure 2B #####

data_S2B <- read.csv("processed_data/S2B_File.csv", stringsAsFactors=FALSE)

data_S2B$cate = factor(data_S2B$category, levels = c("Protein-coding","iCEL","Transcription factors","Homologies"))


Fig_2B <- ggplot(data_S2B,aes(x=cate, y=coding_variants_perc*100)) + 
  geom_bar(position="dodge", stat="identity",aes(fill=factor(strain)),color="black") + 
  scale_fill_manual(values=c("hotpink3", "blue", "cadetblue3", "orange")) +
  theme_classic() + 
  theme(axis.text.x =  element_text(size=12, angle = 90, hjust = 1,color = "black"),
        axis.text.y =  element_text(size=12,color = "black"),
        legend.text = element_text(size=12),
        legend.title = element_blank(),
        axis.title =  element_text(size=12,color = "black"),
        axis.title.x = element_blank(),
        strip.text.x = element_text(size=12, vjust = 1),
        strip.background = element_blank(), 
        legend.position="top", 
        legend.spacing.x = unit(1, 'mm'),
        text=element_text(family="Arial")) +  
  ylab("Proteins with \n coding variants (%)")  + 
  xlab("Gene category") + 
  labs(fill='Strains') + 
  guides(fill=guide_legend(ncol = 3,byrow=TRUE)) + 
  scale_x_discrete(labels=c("Protein-coding" = "All proteins", "iCEL" = "iCEL1314",  "Homologies" = "Other metabolic\nenzymes","Transcription factors"="Transcription\nfactors"))


## Figure 2C #####

data_S2C <- read.csv("processed_data/S2C_File.csv", stringsAsFactors=FALSE)

Fig_2C <- ggplot(subset(data_S2C,level=="LEVEL.4" & impact == "High"), aes(y=fct_reorder(pathway_count,total_gene, .desc = TRUE), x=gene_w_variants)) + 
  geom_point(aes(col=strain), size=1, alpha=0.8) +  
  scale_color_manual(values=c("hotpink3",  "blue","cadetblue3","orange"))  + 
  theme_bw() +
  theme( axis.text.x =  element_text(size=12, color = "black"),
    axis.text.y =  element_text(size=8,color = "black"),
    strip.text = element_text(size=12),
    strip.background = element_blank(),
    axis.title=element_text(size=12), 
    axis.title.y = element_blank(),
    strip.text.y = element_text(size = 12),
    legend.position="none" ,
    legend.text = element_text(size=12),
    plot.margin = unit(c(5, 5, 5, 15), "mm"),
    text=element_text(family="Arial")) +  
  facet_grid(.~impact, scales = "free")  +
  labs(color='Strains') +
  xlab("Number of genes with\nhigh impact genetic variants") + 
  scale_y_discrete(labels=c("AMINO SUGAR AND NUCLEOTIDE SUGAR METABOLISM (64)" = "AMINO SUGAR AND NUCLEOTIDE\nSUGAR METABOLISM (64)",
                            "PENTOSE AND GLUCURONATE INTERCONVERSIONS (77)"="PENTOSE AND GLUCURONATE\nINTERCONVERSIONS (77)",
                            "NICOTINATE AND NICOTINAMIDE METABOLISM (25)"="NICOTINATE AND NICOTINAMIDE\nMETABOLISM (25)"))

## Figure 2D #####

data_S2D <- read.csv("processed_data/S2D_File.csv", stringsAsFactors=FALSE)


Fig_2D <- ggplot(subset(data_S2D,level=="LEVEL.4"), aes(y=fct_reorder(pathway_count,total_gene, .desc = TRUE), x=variants_perc*100)) + 
  geom_point(aes(col=strain,shape=enrich),size=2) +  
  scale_color_manual(values=c("hotpink3",  "blue","cadetblue3","orange"))  + 
  theme_bw() + 
  facet_grid(.~impact, scales = "free",space = "free")  +
  theme(axis.text.x =  element_text(size=12,color = "black"),
        axis.text.y =  element_text(size=8,color = "black"),
        strip.text = element_text(size=12),
        strip.background = element_blank(),
        axis.title=element_text(size=12), 
        axis.title.y=element_blank(),
        strip.text.y = element_text(size = 12),
        legend.position=c(0.75,0.3),
        legend.text = element_text(size=12),
        text=element_text(family="Arial"),
        plot.margin = unit(c(0, 2, 2, 10), "mm"))+  
  labs(color='Strains',shape='Enrichment') +
  xlab("Percentage of genes with enriched\nmoderate impact genetic variants (%)") + 
  xlim(20,100)+
  guides(color=FALSE) + 
  scale_y_discrete(labels=c("AMINO SUGAR AND NUCLEOTIDE SUGAR METABOLISM (64)" = "AMINO SUGAR AND NUCLEOTIDE\nSUGAR METABOLISM (64)",
                            "GLYCEROPHOSPHOLIPID METABOLISM (93)"="GLYCEROPHOSPHOLIPID\nMETABOLISM (93)"))


###### cowplot Fig 2 ######

fig_2AB <- cowplot::plot_grid(Fig_2A,Fig_2B,
                              labels = c("A","B"), 
                              nrow = 2,
                              ncol = 1,
                              label_size = 12, 
                              label_fontfamily="Arial")


fig_2CD <- cowplot::plot_grid(Fig_2C,Fig_2D,
                              labels = c("C","D"), 
                              nrow = 2,
                              ncol = 1,
                              label_size = 12, 
                              label_fontfamily="Arial",
                              rel_widths = c(1, 0.5),
                              rel_heights = c(2, 0.8),
                              axis = "lr",
                              align = "lr")





Fig_2 <- cowplot::plot_grid(fig_2AB,fig_2CD,
                            nrow = 1,
                            ncol = 2,
                            label_size = 12, 
                            label_fontfamily="Arial",
                            rel_widths = c(1.2, 2),
                            rel_heights = c(1, 1))

ggsave(Fig_2 , filename = paste("figures/Fig_2.svg",sep = ""),  units = "mm",height = 225, width = 175)

##########################################
##########################################
#           Figure 3                     #
#       Transcriptomics                  #
##########################################

## Figure 3A #####

data_S3A <- read.csv("processed_data/S3A_File.csv", stringsAsFactors=FALSE)

data_S3A$Strains<- factor(data_S3A$strain,levels = c("N2","BRC20067","CB4856","DL238"))

Fig_3A <- ggplot(data_S3A, aes_string(x = "PC1", y = "PC2", color = "Strains")) + 
  geom_point(size = 2,aes(shape=sync_meth)) + 
  xlab(paste0("PC1: 43% variance")) + 
  ylab(paste0("PC2: 18% variance")) + 
  coord_fixed() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size=12),
        axis.title=element_text(size=12),
        text=element_text(family="Arial"),
        legend.text = element_text(size=9),
        legend.title =  element_text(size=9)) + 
  scale_color_manual(values=c("orange","hotpink3","blue", "cadetblue3")) + 
  scale_shape_manual(values=c(16,1),labels=c("Bleach", "Settle")) +
  labs(color = "Strain", shape="Synchronization")+ 
  scale_x_continuous(breaks=seq(-30,30,15))

## Figure 3B #####

data_S3B <- read.csv("processed_data/S3B_File.csv", stringsAsFactors=FALSE)



df_limma_all <- data_S3B %>%
  ungroup() %>%
  dplyr::select(c("BRC20067","DL238","CB4856")) 



## venn data 
vdc_all <- vennCounts(df_limma_all)

class(vdc_all) <- 'matrix'


df.vdc_all <- as.data.frame(vdc_all)[-1,] %>%
  mutate(x = c(0, 1.2, 0.8, -1.2, -0.8, 0, 0),
         y = c(1.2, -0.6, 0.5, -0.6, 0.5, -1, 0))


## circle coordinates

df.venn <- data.frame(x = c(0, 0.866, -0.866),
                      y = c(1, -0.5, -0.5),
                      r0 = c(1.5, 1.5, 1.5),
                      labels = c('CB4856', 'DL238','BRC20067'))

Fig_3B <- ggplot() +
  geom_circle(data=df.venn,aes(x0 = x, y0 = y, r = r0, fill = labels, colour= labels), alpha = 0, size = 1) +
  coord_fixed() +
  annotate("text", x = df.vdc_all$x, y = df.vdc_all$y, label = df.vdc_all$Counts, size = 4) +
  scale_fill_manual(values = c( 'hotpink3','blue', 'cadetblue3', "orange"), guide = FALSE) +
  scale_colour_manual(values = c('hotpink3','blue', 'cadetblue3', "orange")) +
  labs(fill = NULL) +
  theme_void()  +
  theme(legend.position = 'none', 
        legend.text = element_text(size=12),
        legend.title =  element_text(size=12),
        plot.title = element_text(size=12, face = "bold",hjust=0.1),
        text=element_text(family="Arial"),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  labs(color = "Strain")

## Figure 3C #####

data_S3C <- read.csv("processed_data/S3C_File.csv", stringsAsFactors=FALSE)


df_limma_class <- data_S3C %>% 
  dplyr::select(c("BRC20067","DL238","CB4856")) 



## venn data 
vdc_class <- vennCounts(df_limma_class)

class(vdc_class) <- 'matrix'


df.vdc_class <- as.data.frame(vdc_class)[-1,] %>%
  mutate(x = c(0, 1.2, 0.8, -1.2, -0.8, 0, 0),
         y = c(1.2, -0.6, 0.5, -0.6, 0.5, -1, 0))


## circle coordinates


Fig_3C <- ggplot() +
  geom_circle(data=df.venn,aes(x0 = x, y0 = y, r = r0, fill = labels, colour= labels), alpha = 0, size = 1,linetype=1) +
  coord_fixed() +
  annotate("text", x = df.vdc_class$x, y = df.vdc_class$y, label = df.vdc_class$Counts, size = 4) +
  scale_fill_manual(values = c('hotpink3','blue', 'cadetblue3', "orange"), guide = FALSE) +
  scale_colour_manual(values = c('hotpink3','blue', 'cadetblue3', "orange")) +
  labs(fill = NULL) +
  theme_void()  +
  theme(legend.position = 'none', 
        legend.text = element_text(size=12),
        legend.title =  element_text(size=12),
        text=element_text(family="Arial"),
        plot.title = element_text(size=12, face = "bold",hjust=0.1),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  labs(color = "Strain")

## Figure 3DEF #####

data_S3DEF <- read.csv("processed_data/S3DEF_File.csv", stringsAsFactors=FALSE)



Fig_3DEF_list = list()

reference <- c("BRC20067_vs_N2","CB4856_vs_N2","DL238_vs_N2")



for (i in reference){
  
  data_all <- data_S3DEF %>% 
    dplyr::filter(comparison==i)
  
  data_class <-data_all %>% dplyr::filter(iCEL=="Yes")
  
  
  st_color<-ifelse(i=="BRC20067_vs_N2","hotpink3",
                   ifelse(i=="CB4856_vs_N2","blue",
                          ifelse(i=="DL238_vs_N2","cadetblue3","orange")))
  
  
  fc = 1.5
  
  
  n_down<-paste("Down: ",nrow(dplyr::filter(data_all,fold=="down"))," (", nrow(dplyr::filter(data_class,fold=="down")),")",sep = "") 
  
  n_up <-paste("Up: ",nrow(dplyr::filter(data_all,fold=="up"))," (", nrow(dplyr::filter(data_class,fold=="up")),")",sep = "") 
  
  
  #
  
  plt_ma_class <- ggplot()  +
    geom_point(data=data_all, aes(x = log2(baseMean+1), y = log2FoldChange),color = st_color,  size = 0.4) +
    geom_point(data=data_class,  aes(x = log2(baseMean+1), y = log2FoldChange), shape=17,color = "springgreen4",  size = 0.8) + 
    geom_hline(yintercept = c(0, -log2(fc), log2(fc)), linetype = c(1, 2, 2),
               color = c("black", "black", "black")) +
    theme_minimal() + 
    ylim(-15,15) + xlim(0,17)+ 
    theme(axis.text =  element_text(size=12,color = "black"),
          plot.title = element_text(size=12,hjust = 0.5),
          axis.title =  element_text(size=12),          
          legend.text = element_text(size=12),
          text=element_text(family="Arial"),
          legend.spacing.x = unit(0, 'cm'),
          legend.title = element_blank(),
          legend.position = "none")+  
    #guides(shape = guide_legend(reverse=TRUE,nrow = 1)) +  
    xlab("Log2 mean expression")  + 
    ylab("Log2 fold change") + 
    annotate(geom="text", x=4.5, y=15, label=n_up,size=12/14*5)+ 
    annotate(geom="text", x=4.5, y=-15, label=n_down,size=12/14*5)

  
  Fig_3DEF_list[[i]] <- plt_ma_class
  
}


## Figure 3G #####

data_S3G <- read.csv("processed_data/S3G_File.csv", stringsAsFactors=FALSE)

data_S3G$regulation = factor(data_S3G$Regulation, levels = c("Up","Down"))

Fig_3G <- ggplot(subset(data_S3G,level=="LEVEL.4"),
                 aes(y=DEG_perc*100,
                     x=fct_reorder(pathway_count,total_gene, .desc = TRUE))) + 
  geom_bar( position="dodge",
            stat="identity",
            aes(fill=factor(strain),
                color=enrich),
            size=1) + 
  coord_flip() +
  scale_fill_manual(values=c("hotpink3",  "blue","cadetblue3","orange"),guide = 'none') +
  scale_color_manual(values=c("black",  "gold2"))+
  labs(fill='Strains',color='Enrichment') + 
  theme_bw() +
  theme(axis.text.x =  element_text(size=12,color = "black"),
        axis.text.y =  element_text(size=9,color = "black"),
        legend.text = element_text(size=12),
        legend.title =  element_text(size=12),
        legend.position = c(0.5,0.8),
        legend.background = element_rect(fill=alpha(0.4)),
        axis.title.x =  element_text(size=12),
        axis.title.y =   element_blank(),
        strip.text.y = element_text(size=12),
        strip.text.x = element_text(size=12),
        strip.background = element_blank(),
        text=element_text(family="Arial"),
        plot.margin = unit(c(0, 0, 5, 10), "mm"),
        legend.key = element_rect(fill = "transparent", color = "transparent")) +  
  ylab("Percentage of differential expressed genes (%)")  + 
  xlab("Pathways") + 
  facet_grid(regulation~strain,scales = "free_y", space = "free")  + 
  scale_y_continuous(breaks=c(0,50,100)) + 
  scale_x_discrete(labels=c("PENTOSE AND GLUCURONATE INTERCONVERSIONS (77)" = "PENTOSE AND GLUCURONATE\nINTERCONVERSIONS (77)"))




###### cowplot Fig 3 ####

Fig_3ABC <- cowplot::plot_grid(Fig_3A,Fig_3B,Fig_3C,
                               labels = c("A","B","C"), 
                               nrow = 1,
                               ncol = 3,
                               label_size = 12, 
                               label_fontfamily="Arial",
                               rel_widths = c(2, 1,1))


Fig_3DEF <- cowplot::plot_grid(Fig_3DEF_list$BRC20067_vs_N2,Fig_3DEF_list$CB4856_vs_N2,Fig_3DEF_list$DL238_vs_N2,
                               labels = c("D","E","F"), 
                               nrow = 1,
                               ncol = 3,
                               label_size = 12, 
                               label_fontfamily="Arial",
                               rel_widths = c(1, 1,1))



Fig_3 <- cowplot::plot_grid(Fig_3ABC,Fig_3DEF,Fig_3G,
                            labels = c("","","G"), 
                            nrow = 3,
                            ncol = 1,
                            label_size = 12, 
                            label_fontfamily="Arial",
                            rel_heights =  c(0.9,0.9,1.3))

ggsave(Fig_3, filename = paste("figures/Fig_3.svg",sep = ""), units = "mm",height = 225, width = 174)



##########################################
##########################################
#           Supplemental                 #
#              figures                   #
##########################################

## Fig 2C supp1 ####
# High impact, pathway Level 1,2,3
Fig2C_supp1 <- ggplot(subset(data_S2C,(!level=="LEVEL.4") & impact == "High"), aes(y=fct_reorder(pathway_count,total_gene, .desc = TRUE), x=gene_w_variants)) + 
  geom_point(aes(col=strain), size=1,alpha=0.8) +   
  scale_color_manual(values=c("hotpink3",  "blue","cadetblue3","orange"))  + 
  theme_bw() +
  theme(axis.text.x =  element_text(size=8, color = "black"),
    axis.text.y =  element_text(size=6,color = "black"),
    strip.background = element_blank(),
    axis.title=element_text(size=12), 
    axis.title.y = element_text(hjust = 0.8),
    strip.text.y = element_text(size = 8),
    legend.text = element_text(size=12),
    text=element_text(family="Arial"))+  
  facet_grid(level~impact, scales = "free",space = "free")  +
  labs(color='Strains') +
  xlab("Number of genes") + 
  ylab("Pathways")


ggsave(Fig2C_supp1, filename = paste("figures/Fig_2C_supp1.svg",sep = ""), units = "mm",height = 225, width = 175)


## Fig 2C supp2 ####
# Moderate impact, pathway Level 4

Fig2C_supp2 <- ggplot(subset(data_S2C,(level=="LEVEL.4") & impact == "Moderate"), aes(y=fct_reorder(pathway_count,total_gene, .desc = TRUE), x=gene_w_variants)) + 
  geom_point(aes(col=strain), size=1,alpha=0.8) +   # Draw points
  scale_color_manual(values=c("hotpink3",  "blue","cadetblue3","orange"))  + 
  theme_bw() +
  theme(axis.text.x =  element_text(size=12,color = "black"),
        axis.text.y =  element_text(size=7,color = "black"),
        strip.text = element_text(size=12),
        strip.background = element_blank(),
        axis.title=element_text(size=12), 
        axis.title.y = element_blank(),
        strip.text.y = element_text(size = 12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        text=element_text(family="Arial")) +  
  facet_grid(level~impact, scales = "free")  +
  labs(color='Strains') +
  xlab("Number of genes with\ngenetic variants")


ggsave(Fig2C_supp2, filename = paste("figures/Fig_2C_supp2.svg",sep = ""), units = "mm",height = 225, width = 175)

## Fig 3G supp1 ####
# pathway Level 1,2,3


Fig3G_supp1 <- ggplot(subset(data_S3G,(!level=="LEVEL.4")),
                 aes(y=DEG_perc*100,
                     x=fct_reorder(pathway_count,total_gene, .desc = TRUE))) + 
  geom_point(aes(col=strain,shape=enrich),size=2) +   
  coord_flip() +
  scale_color_manual(values=c("hotpink3",  "blue","cadetblue3","orange")) +
  theme_bw() +
  theme(axis.text.x =  element_text(size=12,color = "black"),
        axis.text.y =  element_text(size=8,color = "black"),
        legend.text = element_text(size=12),
        legend.title =  element_text(size=12),
        axis.title.y = element_blank(),
        legend.position = "bottom",
        axis.title =  element_text(size=12),
        strip.text.y = element_text(size=12,angle = 360),
        strip.text.x = element_text(size=12),
        strip.background = element_blank(),legend.direction = "vertical",
        text=element_text(family="Arial")) +
  ylab("Percentage of differential expressed genes (%)")  + 
  labs(shape="Enrichment",color="Strains")+
  facet_grid(level~regulation,scales = "free_y", space = "free")

ggsave(Fig3G_supp1, filename = paste("figures/Fig_3G_supp1.svg",sep = ""), units = "mm",height = 225, width = 175)


## Supplemental Fig 1 ####


data_SS1 <- read.csv("processed_data/SS1_File.csv", stringsAsFactors=FALSE)

data_SS1$Strains<- factor(data_SS1$strain,levels = c("N2","BRC20067","CB4856","DL238"))


S1_Fig <- ggplot(data_SS1,aes(x=log_bleach,y=log_settle,color=Strains))+
  geom_point(size=0.1)  +
  theme_bw() +
  theme(axis.text.x =  element_text(size=12,color = "black"),
        axis.text.y =  element_text(size=12,color = "black"),
        plot.title = element_text(hjust = 0.1, size=12),
        axis.title =  element_text(size=12,color = "black"),
        strip.text = element_text(size=12, vjust = 1,color = "black"),
        text=element_text(family="Arial"),
        legend.position = "none",
        strip.background = element_blank()) + 
  xlab("Bleach: log2(normalized counts)") +
  ylab("Settle: log2(normalized counts)") +
  facet_wrap(.~Strains) + 
  scale_color_manual(values=c("orange","hotpink3","blue", "cadetblue3")) 


ggsave(S1_Fig, filename = paste("figures/Supp_Fig_1.svg",sep = ""),  units = "mm",height = 174, width = 174)






















