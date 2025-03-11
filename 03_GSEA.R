library(ggplot2)
library(ggrepel)
library(dplyr)
library(Seurat)
library(patchwork)
library(viridis)
library(magrittr)
library(pheatmap)
library(clusterProfiler)
library(cowplot) #combine plots

# GSEA 
load("D:/Documents/wqt/ZXJ_data/DEGS_all.RData")
genelist<- degs_E3$avg_log2FC
#genelist<- degs_E3_R$avg_log2FC

names(genelist) <- rownames(degs_E3_R)
genelist <- sort(genelist,decreasing = T)

T_CELL_RECEPTOR_geneset <- read.gmt("D:/Documents/wqt/ZXJ_data/GSEA/KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY.v2024.1.Hs.gmt") 
length(unique(KEGG_geneset$term))
egmt <- GSEA(genelist, TERM2GENE=KEGG_geneset, 
             minGSSize = 1,
             pvalueCutoff = 0.99,
             verbose=FALSE)
gseaplot2(egmt,geneSetID = 'KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY',pvalue_table=T)

wp_tcr_co_geneset <- read.gmt("D:/Documents/wqt/ZXJ_data/GSEA/WP_T_CELL_RECEPTOR_AND_COSTIMULATORY_SIGNALING.v2024.1.Hs.gmt") 
egmt <- GSEA(genelist, TERM2GENE=wp_tcr_co_geneset, 
             minGSSize = 1,
             pvalueCutoff = 0.99,
             verbose=FALSE)
gseaplot2(egmt,geneSetID = 'WP_T_CELL_RECEPTOR_AND_COSTIMULATORY_SIGNALING',pvalue_table=T)

