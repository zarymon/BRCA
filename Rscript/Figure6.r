#########################TCGA的差异表达数据分析
#####数据处理，NSUN.FPKM来自Figure5.r
NSUNS.fpkm.clinc1<-as.data.frame(t(NSUN.FPKM))
rownames(NSUNS.fpkm.clinc1)<-colnames(NSUN.FPKM)
NSUNS.fpkm.clinc1$treat<-substr(rownames(NSUNS.fpkm.clinc1),14,15)
NSUNS.fpkm.clinc1<-dplyr::filter(NSUNS.fpkm.clinc1,treat %in% c("01"))
NSUNS.fpkm.clinc1$ID=substr(rownames(NSUNS.fpkm.clinc1),1,12)
NSUNS.fpkm.clinc1$SampleID=rownames(NSUNS.fpkm.clinc1)
NSUNS.fpkm.clinc1<-merge(patient,NSUNS.fpkm.clinc1,by="ID")
NSUN2.level1<-dplyr::select(NSUNS.fpkm.clinc1,c("OS.Status","OS.time","NSUN2","SampleID","treat"))%>%dplyr::filter(treat %in% c("01"))%>%dplyr::select(-c(treat)) %>%dplyr::mutate(NSUN2.level=ifelse(NSUN2>median(NSUN2),"High NSUN2","Low NSUN2"))
rownames(NSUN2.level1)<-NSUN2.level1$SampleID
NSUN2.design=model.matrix(~0+NSUN2.level1$NSUN2.level)
rownames(NSUN2.design) = rownames(NSUN2.level1)
colnames(NSUN2.design) <- c("high","low")

######TCGA数据的差异表达分析
library(limma) ####如果是RNA-Seq数据使用的是count
library(edgeR)
######
all<-TCGA.BRCA$tpm$mRNA
MT.gene<-c(grep("^MT-",rownames(all),value=T),grep("^RP[SL]",rownames(all),value=T))
NMT.gene<-rownames(all)[!(rownames(all) %in% MT.gene)]
all<-all[NMT.gene,]
OS.exprs<-all[,rownames(NSUN2.design)]
OS.dge <- DGEList(counts=OS.exprs)
OS.keep=filterByExpr(OS.dge,NSUN2.design)
OS.dge<-OS.dge[OS.keep,,keep.lib.sizes=FALSE] ###移除0和低计数值
OS.dge <- calcNormFactors(OS.dge)   ##calcNormFactors并不会标准化数据，只能标准化因子
OS.num<-dim(OS.dge)[1]
OS.v <- voom(OS.dge,NSUN2.design, normalize="quantile")   ##权重计算
OS.fit=lmFit(OS.v,NSUN2.design) ####RNA-seq 数据是把count转化成log2-counts per million (logCPM)；芯片数据直接取log2
OS.cont.matrix=makeContrasts('high-low',levels = NSUN2.design)  ###如果是多组（contrasts=c("control-patient“)
OS.fit2=contrasts.fit(OS.fit,OS.cont.matrix)
OS.fit=eBayes(OS.fit2)
#FC_result<-topTable(fit,adjust='BH',coef=1,number=dim(expre)[1],lfc=1,p.value=0.05, sort.by="B", resort.by="M" )
OS.DEG_result.NSUN2<-topTable(OS.fit,adjust='BH',coef=1,number=OS.num,sort.by="logFC" )
#####
###定义排序
OS.alldiff.NSUN2 <- OS.DEG_result.NSUN2[order(OS.DEG_result.NSUN2$logFC,decreasing = T),]
OS.id.NSUN2 <- OS.alldiff.NSUN2$logFC
names(OS.id.NSUN2) <- rownames(OS.alldiff.NSUN2)
OS.DEG.gene.NSUN2<- dplyr::filter(OS.alldiff.NSUN2,adj.P.Val<0.05,abs(logFC)>=1)
dim(OS.DEG.gene.NSUN2)


#################################Metabric的差异表达分析
#####数据的处理
Metabric.BRCA.fpkm<-read.table("/mnt/nas_002/home/limma/02.project/04.BRCA/00.data/00.TCGA/Metabric/brca_metabric/Metabric.data_mrna_illumina_microarray.txt",header=TRUE,sep="\t",check.names=FALSE)
Metabric.patient<-read.table("/mnt/nas_002/home/limma/02.project/04.BRCA/00.data/00.TCGA/Metabric/brca_metabric/Metabric.data_clinical_patient.txt",header=TRUE,sep="\t")
Metabric.BRCA.fpkm<-dplyr::select(Metabric.BRCA.fpkm,-c("Entrez_Gene_Id"))
index=order(rowMeans(Metabric.BRCA.fpkm[,-1]),decreasing = T)
#调整表达谱的基因顺序
expr_ordered=Metabric.BRCA.fpkm[index,]
#对于有重复的基因，保留第一次出现的那个，即行平均值大的那个
keep=!duplicated(expr_ordered[,1])
#得到最后处理之后的表达谱矩阵
expr_max.mRNA=expr_ordered[keep,]
Metabric.BRCA.fpkm<-expr_max.mRNA
rownames(Metabric.BRCA.fpkm)=Metabric.BRCA.fpkm$Hugo_Symbol
Metabric.BRCA.fpkm<-dplyr::select(Metabric.BRCA.fpkm,-c(Hugo_Symbol))
########
Metabric.clinical<-read.table("/mnt/nas_002/home/limma/02.project/04.BRCA/00.data/00.TCGA/Metabric/brca_metabric/Metabric.data_clinical_sample.txt",header=TRUE,sep="\t")
Metabric.clinical<-merge(Metabric.patient,Metabric.clinical,by="PATIENT_ID") %>% dplyr::select(c(PATIENT_ID,AGE_AT_DIAGNOSIS,OS_MONTHS,OS_STATUS,CLAUDIN_SUBTYPE,ER_STATUS,HER2_STATUS,PR_STATUS,GRADE,TUMOR_SIZE,TUMOR_STAGE))
Metabric.NSUN2.level<-as.data.frame(t(as.matrix(Metabric.BRCA.fpkm["NSUN2",])))
colnames(Metabric.clinical)=c("Patient_ID","age","OS.time","OS.Status","Subtype","ER.status","HER2.status","PR.status","Grade","TUMOE_SIZE","TUMOR_Stage")
Metabric.clinical$OS.time=Metabric.clinical$OS.time/12
Metabric.clinical$OS.Status<-gsub(":LIVING","",Metabric.clinical$OS.Status)
Metabric.clinical$OS.Status<-gsub(":DECEASED","",Metabric.clinical$OS.Status)
Metabric.clinical <- dplyr::mutate(Metabric.clinical,Age = ifelse(age>=60,">=60","<60"))                                                                                                # ER_STATUS,HER2_STATUS,PR_STATUS,GRADE,TUMOR_SIZE,TUMOR_STAGE))
Metabric.NSUN2.level<-as.data.frame(t(as.matrix(Metabric.BRCA.fpkm["NSUN2",])))
Metabric.NSUN2.level<- Metabric.NSUN2.level %>% dplyr::mutate(NSUN2.level=ifelse(NSUN2>median(NSUN2),"High NSUN2","Low NSUN2"))
Metabric.NSUN2.level$Patient_ID<-rownames(Metabric.NSUN2.level)

#####差异表达分析
Metabric.NSUN2.design=model.matrix(~0+Metabric.NSUN2.level$NSUN2.level)
rownames(Metabric.NSUN2.design) = rownames(Metabric.NSUN2.level)
colnames(Metabric.NSUN2.design) <- c("high","low")

Metabric.OS.exprs<-Metabric.BRCA.fpkm[,rownames(Metabric.NSUN2.design)]
Metabric.OS.exprs<-na.omit(Metabric.OS.exprs)
Metabric.OS.dge <- DGEList(counts=Metabric.OS.exprs)
dim(Metabric.OS.dge)
Metabric.OS.keep=filterByExpr(Metabric.OS.dge,Metabric.NSUN2.design)
Metabric.OS.dge<-Metabric.OS.dge[Metabric.OS.keep,,keep.lib.sizes=FALSE] ###移除0和低计数值
dim(Metabric.OS.dge)
Metabric.OS.dge <- calcNormFactors(Metabric.OS.dge)   ##calcNormFactors并不会标准化数据，只能标准化因子
dim(Metabric.OS.dge)
Metabric.OS.num<-dim(Metabric.OS.dge)[1]
Metabric.OS.v <- voom(Metabric.OS.dge,Metabric.NSUN2.design, normalize="quantile")   ##权重计算
Metabric.OS.fit=lmFit(Metabric.OS.v,Metabric.NSUN2.design) ####RNA-seq 数据是把count转化成log2-counts per million (logCPM)；芯片数据直接取log2
Metabric.OS.cont.matrix=makeContrasts('high-low',levels = Metabric.NSUN2.design)  ###如果是多组（contrasts=c("control-patient“)
Metabric.OS.fit2=contrasts.fit(Metabric.OS.fit,Metabric.OS.cont.matrix)
Metabric.OS.fit=eBayes(Metabric.OS.fit2)
#FC_result<-topTable(fit,adjust='BH',coef=1,number=dim(expre)[1],lfc=1,p.value=0.05, sort.by="B", resort.by="M" )
Metabric.OS.DEG_result.NSUN2<-topTable(Metabric.OS.fit,adjust='BH',coef=1,number=OS.num,sort.by="logFC" )
###
###定义排序
Metabric.OS.alldiff.NSUN2 <- Metabric.OS.DEG_result.NSUN2[order(Metabric.OS.DEG_result.NSUN2$logFC,decreasing = T),]
Metabric.OS.id.NSUN2 <- Metabric.OS.alldiff.NSUN2$logFC
names(Metabric.OS.id.NSUN2) <- rownames(Metabric.OS.alldiff.NSUN2)
Metabric.OS.DEG.gene.NSUN2<- dplyr::filter(Metabric.OS.alldiff.NSUN2,adj.P.Val<0.05,abs(logFC)>=1)
dim(Metabric.OS.DEG.gene.NSUN2)


###################################两者共有差异表达基因的富集分析
OS.DEG.gene.NSUN2<- dplyr::filter(OS.alldiff.NSUN2,adj.P.Val<0.05)
Metabric.OS.DEG.gene.NSUN2<- dplyr::filter(Metabric.OS.alldiff.NSUN2,adj.P.Val<0.05)
DEG.set<-intersect(rownames(OS.DEG.gene.NSUN2),rownames(Metabric.OS.DEG.gene.NSUN2))
length(DEG.set)
head(DEG.set)
library("clusterProfiler")
NSUN2.DEG.gene<- bitr(DEG.set, fromType="SYMBOL",
           toType="ENTREZID",
           OrgDb="org.Hs.eg.db")
NSUN2.DEG.go <- enrichGO(DEG.set, 
               OrgDb = "org.Hs.eg.db", 
               ont='ALL',
               keyType = "SYMBOL",
               pAdjustMethod = 'BH',
               pvalueCutoff =0.05,
               readable=TRUE
               #keyType = 'ENTREZID'
)
dim(NSUN2.DEG.gene)

NSUN2.DEG.kegg <- enrichKEGG(gene  = NSUN2.DEG.gene$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
                 
######差异表达基因富集分析结果的可视化
NSUN2.DEG.go1<-dplyr::filter(NSUN2.DEG.go,p.adjust<0.01)
NSUN2.DEG.kegg1<-dplyr::filter(NSUN2.DEG.kegg,p.adjust<0.01)
NSUN2.DEG.go.BP.draw<- dplyr::filter(NSUN2.DEG.go1,ONTOLOGY=="BP")%>%dplyr::arrange(p.adjust)%>%head(5)
#NSUN2.DEG.go.BP.draw$Genecount=gsub("/470","",NSUN2.DEG.go.BP.draw$GeneRatio)
NSUN2.DEG.go.CC.draw<- dplyr::filter(NSUN2.DEG.go1,ONTOLOGY=="CC")%>%dplyr::arrange(p.adjust)%>%head(5)
#NSUN2.DEG.go.CC.draw$Genecount=gsub("/474","",NSUN2.DEG.go.CC.draw$GeneRatio)
NSUN2.DEG.go.MF.draw<- dplyr::filter(NSUN2.DEG.go1,ONTOLOGY=="MF")%>%dplyr::arrange(p.adjust)%>%head(5)
NSUN2.DEG.go.MF.draw$ONTOLOGY="CMF"
#NSUN2.DEG.go.MF.draw$Genecount=gsub("/472","",NSUN2.DEG.go.MF.draw$GeneRatio)
NSUN2.DEG.kegg.draw<-dplyr::arrange(NSUN2.DEG.kegg1,p.adjust)%>%head(5)
#NSUN2.DEG.kegg.draw$Genecount=gsub("/299","",NSUN2.DEG.kegg.draw$GeneRatio)
#NSUN2.DEG.go.draw<-rbind(NSUN2.DEG.go.BP.draw,NSUN2.DEG.go.CC.draw,NSUN2.DEG.go.MF.draw)
NSUN2.DEG.kegg.draw$ONTOLOGY<-"KEGG"
NSUN2.DEG.draw<-rbind(NSUN2.DEG.go.BP.draw,NSUN2.DEG.go.CC.draw,NSUN2.DEG.go.MF.draw,NSUN2.DEG.kegg.draw)

NSUN2.DEG.draw$p.adjust<-as.numeric(NSUN2.DEG.draw$p.adjust)
options(repr.plot.height=10, repr.plot.width=10)

NSUN2.go_draw<-ggplot(NSUN2.DEG.draw,aes(Count,Description))+geom_point(aes(size=Count,color=-1*log10(p.adjust) ))+ labs(x="GeneRatio",y="TOP5 GO Term")+
  theme_bw()+theme(axis.text = element_text(face="plain",size=12),axis.title = element_text(face="plain",size = 14)) + facet_wrap(~ONTOLOGY,scales='free_y',ncol=1)+
  labs(color=expression(-log[10](p.adjust)),size="Count",x="Count",y="",title="")+scale_color_gradient(low="yellow",high = "red")
  
NSUN2.go_draw

########################GSEA富集分析
library(enrichplot)
library('GSEABase')
library(fgsea)
library(clusterProfiler)
library(ggsci)
#######MYC
#######加载预定义基因集
MYC.gmtfile <- "/mnt/nas_002/home/limma/02.project/02.COAD/01.analysis/MSigDB/MYC.HALLMARK.symbols.gmt.txt"
MYC.hallmark <- read.gmt(MYC.gmtfile)

##############使用不同的方法得到结果
TCGA.MYC.gsea.re1<-clusterProfiler::GSEA(OS.id.NSUN2 ,TERM2GENE=MYC.hallmark,verbose = T)
TCGA.MYC.g1<-as.data.frame(TCGA.MYC.gsea.re1)
TCGA.MYC.g1<-subset(TCGA.MYC.g1,p.adjust<0.05)
TCGA.MYC.g1<-TCGA.MYC.g1[order(TCGA.MYC.g1$NES,decreasing = T),]

################画图
col_gsea1<-pal_simpsons()(10)
num1=1
TCGA.MYC.clustergsea<-gseaplot2(TCGA.MYC.gsea.re1,geneSetID = rownames(TCGA.MYC.g1)[1:num1],
              title = "",#标题
              #color = col_gsea1[1:num1],#颜色
                            color = "red",#颜色
              base_size = 20,#基础大小
              rel_heights = c(1, 0.2, 0.4),#小图相对高度
              subplots = 1:2,#展示小图
              pvalue_table = TRUE,#p值表格
              ES_geom = "line"#line or dot
          )
TCGA.MYC.clustergsea
##############使用不同的方法得到结果
Metabric.MYC.gsea.re1<-clusterProfiler::GSEA(Metabric.OS.id.NSUN2.test ,TERM2GENE=MYC.hallmark,verbose = T)
Metabric.MYC.g1<-as.data.frame(Metabric.MYC.gsea.re1)
Metabric.MYC.g1<-subset(Metabric.MYC.g1,p.adjust<0.05)
Metabric.MYC.g1<-Metabric.MYC.g1[order(Metabric.MYC.g1$NES,decreasing = T),]

################画图
col_gsea1<-pal_simpsons()(10)
num1=1
Metabric.MYC.clustergsea<-gseaplot2(Metabric.MYC.gsea.re1,geneSetID = rownames(Metabric.MYC.g1)[1:num1],
              title = "",#标题
              #color = col_gsea1[1:num1],#颜色
                            color = "red",#颜色
              base_size = 20,#基础大小
              rel_heights = c(1, 0.2, 0.4),#小图相对高度
              subplots = 1:2,#展示小图
              pvalue_table = TRUE,#p值表格
              ES_geom = "line"#line or dot
          )
Metabric.MYC.clustergsea
########
######################MTORC1
#######加载预定义基因集
MTORC1.gmtfile <- "/mnt/nas_002/home/limma/02.project/02.COAD/01.analysis/MSigDB/MTORC1.HALLMARK.symbols.gmt.txt"
MTORC1.hallmark <- read.gmt(MTORC1.gmtfile)

##############使用不同的方法得到结果
TCGA.MTORC1.gsea.re1<-clusterProfiler::GSEA(OS.id.NSUN2 ,TERM2GENE=MTORC1.hallmark,verbose = T)
TCGA.MTORC1.g1<-as.data.frame(TCGA.MTORC1.gsea.re1)
TCGA.MTORC1.g1<-subset(TCGA.MTORC1.g1,p.adjust<0.05)
TCGA.MTORC1.g1<-TCGA.MTORC1.g1[order(TCGA.MTORC1.g1$NES,decreasing = T),]

################画图
col_gsea1<-pal_simpsons()(10)
num1=1
TCGA.MTORC1.clustergsea<-gseaplot2(TCGA.MTORC1.gsea.re1,geneSetID = rownames(TCGA.MTORC1.g1)[1:num1],
              title = "",#标题
              #color = col_gsea1[1:num1],#颜色
                            color = "red",#颜色
              base_size = 20,#基础大小
              rel_heights = c(1, 0.2, 0.4),#小图相对高度
              subplots = 1:2,#展示小图
              pvalue_table = TRUE,#p值表格
              ES_geom = "line"#line or dot
          )
TCGA.MTORC1.clustergsea
##############使用不同的方法得到结果
Metabric.MTORC1.gsea.re1<-clusterProfiler::GSEA(Metabric.OS.id.NSUN2.test ,TERM2GENE=MTORC1.hallmark,verbose = T)
Metabric.MTORC1.g1<-as.data.frame(Metabric.MTORC1.gsea.re1)
Metabric.MTORC1.g1<-subset(Metabric.MTORC1.g1,p.adjust<0.05)
Metabric.MTORC1.g1<-Metabric.MTORC1.g1[order(Metabric.MTORC1.g1$NES,decreasing = T),]

################画图
col_gsea1<-pal_simpsons()(10)
num1=1
Metabric.MTORC1.clustergsea<-gseaplot2(Metabric.MTORC1.gsea.re1,geneSetID = rownames(Metabric.MTORC1.g1)[1:num1],
              title = "",#标题
              #color = col_gsea1[1:num1],#颜色
                            color = "red",#颜色
              base_size = 20,#基础大小
              rel_heights = c(1, 0.2, 0.4),#小图相对高度
              subplots = 1:2,#展示小图
              pvalue_table = TRUE,#p值表格
              ES_geom = "line"#line or dot
          )
Metabric.MTORC1.clustergsea

######################E2F
#######加载预定义基因集
E2F.gmtfile <- "/mnt/nas_002/home/limma/02.project/02.COAD/01.analysis/MSigDB/E2F.HALLMARK.symbols.gmt.txt"
E2F.hallmark <- read.gmt(E2F.gmtfile)

##############使用不同的方法得到结果
TCGA.E2F.gsea.re1<-clusterProfiler::GSEA(OS.id.NSUN2 ,TERM2GENE=E2F.hallmark,verbose = T)
TCGA.E2F.g1<-as.data.frame(TCGA.E2F.gsea.re1)
TCGA.E2F.g1<-subset(TCGA.E2F.g1,p.adjust<0.05)
TCGA.E2F.g1<-TCGA.E2F.g1[order(TCGA.E2F.g1$NES,decreasing = T),]

################画图
col_gsea1<-pal_simpsons()(10)
num1=1
TCGA.E2F.clustergsea<-gseaplot2(TCGA.E2F.gsea.re1,geneSetID = rownames(TCGA.E2F.g1)[1:num1],
              title = "",#标题
              #color = col_gsea1[1:num1],#颜色
                            color = "red",#颜色
              base_size = 20,#基础大小
              rel_heights = c(1, 0.2, 0.4),#小图相对高度
              subplots = 1:2,#展示小图
              pvalue_table = TRUE,#p值表格
              ES_geom = "line"#line or dot
          )
TCGA.E2F.clustergsea
##############使用不同的方法得到结果
Metabric.E2F.gsea.re1<-clusterProfiler::GSEA(Metabric.OS.id.NSUN2.test ,TERM2GENE=E2F.hallmark,verbose = T)
Metabric.E2F.g1<-as.data.frame(Metabric.E2F.gsea.re1)
Metabric.E2F.g1<-subset(Metabric.E2F.g1,p.adjust<0.05)
Metabric.E2F.g1<-Metabric.E2F.g1[order(Metabric.E2F.g1$NES,decreasing = T),]

################画图
col_gsea1<-pal_simpsons()(10)
num1=1
Metabric.E2F.clustergsea<-gseaplot2(Metabric.E2F.gsea.re1,geneSetID = rownames(Metabric.E2F.g1)[1:num1],
              title = "",#标题
              #color = col_gsea1[1:num1],#颜色
                            color = "red",#颜色
              base_size = 20,#基础大小
              rel_heights = c(1, 0.2, 0.4),#小图相对高度
              subplots = 1:2,#展示小图
              pvalue_table = TRUE,#p值表格
              ES_geom = "line"#line or dot
          )
Metabric.E2F.clustergsea
