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

##GO function
GO_inbatch <- function(deg_res,clustern,db){
  # get gene symbol of up-regulation in one cluster
  subset(deg_res, cluster == clustern) %>% # get DEGs information of one cluster
    rownames %>% # get gene symbol
    bitr(.,fromType="SYMBOL",toType="ENTREZID",OrgDb = db) %$%ENTREZID %>%
    enrichGO(gene=.,
             Org = db,
             ont="BP",
             pvalueCutoff=0.05,
             pAdjustMethod="BH",
             readable="TRUE") # GO analysis
}

#data input of E3_iNK
data.seu <- Read10X(data.dir = "./E3_iNK/count/filtered_feature_bc_matrix/")
E3_iNK <- CreateSeuratObject(counts = data.seu, project = "E3_iNK", min.cells = 5, min.features = 100,meta.data = meta);rm(data.seu);gc()
E3_iNK[["percent.mt"]] <- PercentageFeatureSet(E3_iNK, pattern = "^MT-")
P=VlnPlot(E3_iNK,features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 4,pt.size = 0,raster = F)&labs(x="",y="")
P
E3_iNK<- subset(E3_iNK, subset=nFeature_RNA<7500&nCount_RNA<40000&percent.mt<10)

D27=readRDS(file = "./ESC_OL.rds")
D27$orig.ident="ESC_iNK"

TCRiNK_integ_ctrl=merge(x=NormalizeData(E1_iNK,normalization.method = "LogNormalize", scale.factor =median(E1_iNK$nCount_RNA)),
                        y=NormalizeData(D27,normalization.method = "LogNormalize", scale.factor =median(D27$nCount_RNA))
)

#differential expressed genes between E3_iNK and iNK
degs_E3=FindMarkers(TCRiNK_integ_ctrl,ident.1 ="E3_iNK",ident.2 = "ESC_iNK" ,min.pct = 0.25,logfc.threshold = 0.25,group.by = "orig.ident",test.use = "negbinom",slot="counts")
dim(degs_E3)
# remain the significant DEGs
degs_E3=degs_E3[degs_E3$p_val_adj<0.05,]
degs_E3$gene=rownames(degs_E3)
degs_E3$cluster[degs_E3$avg_log2FC>0]="E3_iNK"
degs_E3$cluster[degs_E3$avg_log2FC<0]="ESC_iNK"
table(degs_E3$cluster)
degs_E3_GO <- sapply(unique(degs_E3$cluster), function(clustern){GO_inbatch(degs_E3, clustern,"org.Hs.eg.db")}, USE.NAMES = TRUE)

GO_plots=lapply(names(degs_E3_GO),function(x){
  if(any(degs_E3_GO[[x]]@result$p.adjust<0.05)){
    dotplot(degs_E3_GO[[x]],title = paste(x,"_GOterms",collapse = ""),showCategory=15,font.size = 15)+
      theme(
        axis.text.x=element_text(colour="black",size=8,face="bold"),
        axis.title.x =element_text(colour="black",size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5,size=20),legend.position = "bottom")+guides(size = guide_legend(nrow = 3, byrow = TRUE,order = 1,label.theme = element_text(size=8)),color=guide_colorbar(direction="vertical",barwidth = 0.5,barheight = 3,label.theme = element_text(size=8,angle = 0),order = 2))
    
  }
})
p=CombinePlots(plots=GO_plots,ncol=2)
p
ggsave(p,filename ="GO_terms_E3_iNK_vs_ESC.pdf",width=20,heigh=10)
save(degs_E3,degs_E3_GO,file = "./DEGS_all.RData")