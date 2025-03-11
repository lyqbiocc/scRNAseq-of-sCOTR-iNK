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

#Figure1 violinplot
iNK <- readRDS("./iNK_all_ref.rds")
table(iNK$orig.ident)
iNK<-subset(iNK,orig.ident %in% c("activated.UCB.NK","ES.iNK"))
gc()
VlnPlot(iNK,c("CD3D","CD3E","CD3G","CD247","CD8A","CD8B","CD28"),stack = T,fill.by = "ident")+NoLegend()

