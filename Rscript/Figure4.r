
##################数据处理
GSE161529_rpca.anno<-readRDS("/mnt/nas_002/home/limma/02.project/04.BRCA/00.data/GSE161529_SingleCell/GSE161529.Seurat3.rpca/anno/GSE161529_rpca.Basic.merge_umap.maker.rds")
GSE161529_rpca.anno$Group="Normal"
GSE161529_rpca.anno@meta.data[which(GSE161529_rpca.anno$orig.ident %in% c("GSM4909277_B1-KCF0894","GSM4909278_B1-MH0033")),"Group"]="preneoplastic BRCA1+"

GSE161529_rpca.anno@meta.data[which(GSE161529_rpca.anno$orig.ident %in% c("GSM4909282_TN-MH0135")),"Group"]="TNBC BRCA1-"
GSE161529_rpca.anno@meta.data[which(GSE161529_rpca.anno$orig.ident %in% c("GSM4909287_TN-B1-Tum0554","GSM4909288_TN-B1-MH0177")),"Group"]="TNBC BRCA1+"

GSE161529_rpca.anno@meta.data[which(GSE161529_rpca.anno$orig.ident %in% c("GSM4909290_HER2-PM0337","GSM4909294_HER2-MH0176")),"Group"]="HER2"
GSE161529_rpca.anno@meta.data[which(GSE161529_rpca.anno$orig.ident %in% c("GSM4909302_ER-MH0025","GSM4909315_ER-MH0167-T")),"Group"]="ER"

######################maker基因的表达分析
know_markers=read.table("/mnt/nas_002/home/limma/02.project/02.COAD/01.analysis/02.TCGA/02.Singlecell/04.Figure/maker.list",sep="\t",header=F,check.names=F) 
selected_genes = split(know_markers$V1,know_markers$V2)  
Idents(GSE161529_rpca.anno)=GSE161529_rpca.anno@meta.data$seurat_clusters
selected_genes = split(know_markers$V1,know_markers$V2)        ###V1:gene     V2:celltype
GSE161529_rpca.anno.p2=DotPlot(GSE161529_rpca.anno, features = selected_genes , cols = c("lightgrey", "red"),assay = "RNA",cluster.idents =T) +theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=0.5))+ theme (legend.position = "bottom")
GSE161529_rpca.anno.p2[["theme"]][["strip.text"]]$angle=90
GSE161529_rpca.anno.p2[["theme"]][["strip.text"]]$size=9
GSE161529_rpca.anno.p2
#####################umap图可视化
GSE161529.Dim=DimPlot(GSE161529_rpca.anno,group.by = 'celltype',label = T)+ theme (legend.position = "bottom")+labs(size=10)
GSE161529.Dim
##################不同细胞中NSUNs的表达分析
library(reshape2)
gene<-c("NOP2","NSUN2","NSUN3","NSUN4","NSUN5","NSUN6","NSUN7")
vln.df=as.matrix(GSE161529_rpca.anno[["RNA"]]@data[gene,])

vln.df<-as.data.frame(vln.df)
vln.df$gene=rownames(vln.df)
vln.df=melt(vln.df,id="gene")
colnames(vln.df)[c(2,3)]=c("CB","exp")
#head(vln.df)

library("tidyverse")
anno=GSE161529_rpca.anno[[c("celltype","Group")]]
anno$CB=rownames(anno)
vln.df=inner_join(vln.df,anno,by="CB")
vln.df$gene=factor(vln.df$gene,levels = gene) #为了控制画图的基因顺序
vln.df<-dplyr::filter(vln.df,exp>0)

expre.draw<-ggplot(vln.df,aes(celltype,exp))+geom_violin(aes(fill=celltype),scale = "count")+
  facet_grid(vln.df$gene~.,scales = "free_y")+
  scale_fill_brewer(palette = "Set3",direction = 1)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(text=element_text(size=14,face="plain",color="black"),
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "none"
  )
options(repr.plot.height=8, repr.plot.width=6)
expre.draw
#####################不同分型患者的上皮细胞中NSUNs的 表达分析
Epithelial.expres<-dplyr::filter(vln.df,celltype=="Epithelial")
library(Hmisc)
Epithelial.draw<-ggplot(Epithelial.expres,aes(Group,exp))+ 
#geom_violin(fill="#56B4E9",scale="count") +

geom_violin(aes(fill=Group),scale = "count")+
stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),
               geom="pointrange",color = "red")+  
facet_grid(Epithelial.expres$gene~.,scales = "free_y")+
  scale_fill_brewer(palette = "Set3",direction = 1)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(text=element_text(size=14,face="plain",color="black"),
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "none"
  )
options(repr.plot.height=8, repr.plot.width=6)
Epithelial.draw

##########################KEGG的富集分析
GSE161529_rpca.anno$Group="Normal"
GSE161529_rpca.anno@meta.data[which(GSE161529_rpca.anno$orig.ident %in% c("GSM4909277_B1-KCF0894","GSM4909278_B1-MH0033")),"Group"]="preneoplastic BRCA1+"

GSE161529_rpca.anno@meta.data[which(GSE161529_rpca.anno$orig.ident %in% c("GSM4909282_TN-MH0135")),"Group"]="TNBC BRCA1-"
GSE161529_rpca.anno@meta.data[which(GSE161529_rpca.anno$orig.ident %in% c("GSM4909287_TN-B1-Tum0554","GSM4909288_TN-B1-MH0177")),"Group"]="TNBC BRCA1+"

GSE161529_rpca.anno@meta.data[which(GSE161529_rpca.anno$orig.ident %in% c("GSM4909290_HER2-PM0337","GSM4909294_HER2-MH0176")),"Group"]="HER2"
GSE161529_rpca.anno@meta.data[which(GSE161529_rpca.anno$orig.ident %in% c("GSM4909302_ER-MH0025","GSM4909315_ER-MH0167-T")),"Group"]="ER"
#head(GSE161529_rpca.anno@meta.data)
GSE161529_rpca.anno.Epithelial<-subset(GSE161529_rpca.anno,celltype=="Epithelial")
Idents(GSE161529_rpca.anno.Epithelial)<-GSE161529_rpca.anno.Epithelial$Group
GSE161529_rpca.anno.Epithelial.mymarker<- FindAllMarkers(GSE161529_rpca.anno.Epithelial, logfc.threshold=0.7,only.pos = TRUE)
GSE161529_rpca.anno.Epithelial.mymarker<-GSE161529_rpca.anno.Epithelial.mymarker %>% 
         dplyr::filter(cluster=="TNBC BRCA1-") %>%dplyr::filter(avg_log2FC>1)
library("clusterProfiler")
TNBC.maker.gene<- bitr(GSE161529_rpca.anno.Epithelial.mymarker$gene, fromType="SYMBOL",
           toType="ENTREZID",
           OrgDb="org.Hs.eg.db")
TNBC.maker.kegg <- enrichKEGG(gene  = TNBC.maker.gene$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)

#############
TNBC.maker.kegg1<-read.table("/mnt/nas_002/home/limma/02.project/04.BRCA/00.data/GSE161529_SingleCell/GSE161529.Seurat3.rpca/anno//TNBC.maker.keggenrich1.xls",header=TRUE,sep="\t")
TNBC.maker.kegg.draw<-dplyr::arrange(TNBC.maker.kegg1,p.adjust)%>%head(10)%>%dplyr::arrange(p.adjust)
TNBC.maker.kegg.draw$ONTOLOGY<-"KEGG"
TNBC.maker.draw<-TNBC.maker.kegg.draw
###
TNBC.maker.draw$p.adjust<-as.numeric(TNBC.maker.draw$p.adjust)
options(repr.plot.height=10, repr.plot.width=10)
TNBC.go_draw<-ggplot(TNBC.maker.draw,aes(Count,Description))+geom_point(aes(size=Count,color=-1*log10(p.adjust) ))+ labs(x="GeneRatio",y="TOP5 GO Term")+
  theme_bw()+theme(axis.text = element_text(face="plain",size=12),axis.title = element_text(face="plain",size = 14)) + facet_wrap(~ONTOLOGY,scales='free_y',ncol=1)+
  labs(color=expression(-log[10](p.adjust)),size="Count",x="Count",y="",title="")+scale_color_gradient(low="yellow",high = "red")
TNBC.go_draw
ggsave("/mnt/nas_002/home/limma/02.project/04.BRCA/00.data/GSE161529_SingleCell/GSE161529.Seurat3.rpca/anno//Enrich.KEGG.pdf",TNBC.go_draw,height=8,width=8)
