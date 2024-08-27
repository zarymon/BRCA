
##################1.建模基因的选取
####所用的数据都来自上面的分析
OS.DEG.gene.NSUN2<- dplyr::filter(OS.alldiff.NSUN2,adj.P.Val<0.05)
Metabric.OS.DEG.gene.NSUN2<- dplyr::filter(Metabric.OS.alldiff.NSUN2,adj.P.Val<0.05)
DEG.set<-intersect(rownames(OS.DEG.gene.NSUN2),rownames(Metabric.OS.DEG.gene.NSUN2))
GSE161529_rpca.anno.Epithelial<-subset(GSE161529_rpca.anno,celltype=="Epithelial")
Idents(GSE161529_rpca.anno.Epithelial)<-GSE161529_rpca.anno.Epithelial$Group
GSE161529_rpca.anno.Epithelial.mymarker<- FindAllMarkers(GSE161529_rpca.anno.Epithelial, logfc.threshold=0.7,only.pos = TRUE)
GSE161529_rpca.anno.Epithelial.mymarker<- GSE161529_rpca.anno.Epithelial.mymarker%>% 
         dplyr::filter(cluster=="TNBC BRCA1-") %>%dplyr::filter(avg_log2FC>1)
library(ggvenn)
library(cowplot)
Venn.data <- list(`DEGS` =DEG.set,
             `TNBC BRCA1-` =GSE161529_rpca.anno.Epithelial.mymarker$gene
             
)
Venn.draw<-ggvenn(Venn.data, c("DEGS", "TNBC BRCA1-"),stroke_color = "white", #c("#E41A1C","orange","red","blue")
             fill_color = c("#E41A1C","orange"),set_name_size=4,text_size=3,show_percentage = FALSE)
options(repr.plot.height=5, repr.plot.width=5)
Venn.draw

###################2.选择建模基因进行Consensus分析
test<-intersect(DEG.set,GSE161529_rpca.anno.Epithelial.mymarker$gene)
length(test)
gene<-as.data.frame(TCGA.BRCA.fpkm)%>% dplyr::select(matches("01A"))

gene = log(gene+1) ###表达数据需要log,使得表达数据符合正态分布
MT.gene<-grep("^MT-",rownames(gene),value=T)
NMT.gene<-rownames(gene)[!(rownames(gene) %in% MT.gene)]
gene<-gene[NMT.gene,]
Con.gene<-gene[test,]

options(repr.plot.height=5, repr.plot.width=5)

##########数据处理，去掉正常样本
Con.gene<-as.matrix(Con.gene)
#head(gene[1:3,1:3])
setwd("/mnt/nas_002/home/limma/02.project/04.BRCA/03.CONs")
workDir="/mnt/nas_002/home/limma/02.project/04.BRCA/03.CONs"
dir.create(workDir, recursive=TRUE)
##########
#对上面这个芯片表达数据我们一般会简单的进行normalization (本次采用中位数中心化),
#然后取在各个样品差异很大的那些gene或者探针的数据来进行聚类分析 mads=apply(d,1,mad)
# mad(x) 绝对中位数差 按行(1)取d数据的
########参数的设置
########参考http://www.360doc.com/content/20/0415/23/46752147_906313676.shtml
#######http://www.bioconductor.org/packages/release/bioc/vignettes/ConsensusClusterPlus/inst/doc/ConsensusClusterPlus.pdf

library(ConsensusClusterPlus)
library(caret)
maxK=6
cluster="km"
distance="euclidean"
#distance="pearson"
########
norma = Con.gene
results=ConsensusClusterPlus(d=norma,maxK=maxK,reps=10,pItem=0.8,pFeature=1,title=workDir,
                             clusterAlg=cluster,distance=distance,#计算距离的方法有spearman' , 'euclidean', 'binary', 'maximum', 'canberra' 
                                                                # 'minkowski" or custom distance function.
                             seed=123456,
                             plot="pdf" 
                            )
#results=ConsensusClusterPlus(d=norma,maxK=maxK,reps=1000,pItem=0.8,pFeature=1,title=workDir,
 #                            clusterAlg=cluster,distance=distance,#计算距离的方法有spearman' , 'euclidean', 'binary', 'maximum', 'canberra' 
  #                                                              # 'minkowski" or custom distance function.
   #                          seed=123456,
    #                         plot="pdf" )

saveRDS(results,"ConsensusClusterPlus.rds")
#########################3.PCA主成分分析
group=as.data.frame(results[[2]][["consensusClass"]])
colnames(group)=c("group")
for (i in 1:length(group$group)){group$group[i] <- paste("cluster",group$group[i],sep="")}
data = t(norma)

###########
library("FactoMineR")
library("factoextra")
pca<- PCA(data, graph = FALSE)
cluster2<-fviz_pca_ind(pca,
             geom.ind = "point", # show points only (nbut not "text")
             pointsize =3,
             pointshape = 21,
             fill.ind = group$group, # color by groups
             palette = c("#FF3030","#0000FF"),#c("#FF3030","#0000FF"),#"lacent",# c("#00AFBB", "#E7B800", "#FC4E07")
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",
             title="")+theme_bw() +theme(text=element_text(size=14,face="plain",color="black"),
        axis.title=element_text(size=16,face="plain",color="black"),
        axis.text = element_text(size=14,face="plain",color="black"),
        legend.title = element_text(size=16,face="plain",color="black"),
        legend.text = element_text(size=14,face="plain",color="black"),
        legend.background = element_blank(),
        legend.position=c(0.9,0.1))+xlim(-15,15)+ylim(-15,15)
cluster2
ggsave("Cluster2.pca.pdf",cluster2,width=6,height=6)

#########################4.生存分析
survival_TCGA.BRCA<-TCGA.BRCA$TCGA.clinc%>%dplyr::select(ID,OS.Status,OS.time)

survival_TCGA.BRCA$neID<-survival_TCGA.BRCA$ID
survival_data<-survival_TCGA.BRCA

 
keep=!duplicated(survival_data$ID)
survival_data=survival_data[keep,]
survival_data=survival_data[keep,] %>%dplyr::filter(OS.time!=0) ###过滤生存时间为0的

#dim(survival_data)
group=as.data.frame(results[[2]][["consensusClass"]])
colnames(group)=c("group")
dim(group)
for (i in 1:length(group$group)){group$group[i] <- paste("cluster",group$group[i],sep="")}
group$ID=substr(rownames(group),1,12)
group=group[!duplicated(group$ID),]

km_data<- merge(survival_data,group,by="ID")
###################KM生存分析
library(survminer) 
library(survival) 
#2. 拟合生存曲线
fit <- survfit(Surv(OS.time,OS.Status) ~ group,  data = km_data) #2. 拟合生存曲线

#summary(fit) #查看生存分析结果
#3. 绘制基础曲线
sur3<-ggsurvplot(fit, # 创建的拟合对象
           data = km_data,  # 指定变量数据来源
           conf.int = TRUE, # 显示置信区间
           pval = TRUE, # 添加P值
           risk.table = TRUE, # 绘制累计风险曲线
           surv.median.line = "hv", # 添加中位生存时间线
           #add.all = TRUE, # 添加总患者生存曲线
            palette=c("#FF3030","#0000FF")) + xlab("Time (year)")
sur3
#ggsave("Cluster.KM.pdf",sur,width=6,height=6)
pdf(file = "/mnt/nas_002/home/limma/02.project/04.BRCA/03.CONs/Fig1.Cluster2.KM.pdf",width =6,height = 6)
print(sur3)
dev.off()

##################################5.分组与临床的相关性分析结果
TCGA.clinc<-TCGA.BRCA$TCGA.clinc
group.new<-group
group.new$ID<-substr(rownames(group.new),1,12)
data<-merge(TCGA.clinc,group.new,by="ID")%>% dplyr::distinct(ID,.keep_all=TRUE)
data$Stage<-as.character(data$Stage)
data[which(data$TNM_M %in% c("MX")),]$TNM_M<-NA
data[which(data$TNM_T %in% c("TX","Ti")),]$TNM_T<-NA
data[which(data$TNM_N %in% c("NX")),]$TNM_N<-NA
#data[which(data$TNM_M %in% c("MX")),]$TNM_M<-NA
data[which(data$Stage %in% c("","Stage X")),]$Stage<-NA
data$Stage<-as.character(data$Stage)
#data$Stage=droplevels(data$Stage)
data[which(data$Stage %in% c("Stage IA","Stage IB")),]$Stage="Stage I"
data[which(data$Stage %in% c("Stage IIA","Stage IIB","Stage IIC")),]$Stage="Stage II"
data[which(data$Stage %in% c("Stage IIIA","Stage IIIB","Stage IIIC")),]$Stage="Stage III"
#data[which(data$Stage %in% c("Stage IVA","Stage IVB")),]$Stage="Stage IV"

TNM_T<-table(data$group,data$TNM_T)
p.value<-chisq.test(TNM_T,simulate.p.value=TRUE)$p.value 
p.value
TNM_T<-as.data.frame(TNM_T)
colnames(TNM_T)=c("Group","clinc","Freq")
TNM_T$sub=paste("TNM_T",p.value)

#############
TNM_M<-table(data$group,data$TNM_M)
p.value<-chisq.test(TNM_M,simulate.p.value=TRUE)$p.value 
p.value
TNM_M<-as.data.frame(TNM_M)
colnames(TNM_M)=c("Group","clinc","Freq")
TNM_M$sub=paste("TNM_M",p.value)
###########

TNM_N<-table(data$group,data$TNM_N)
p.value<-chisq.test(TNM_N,simulate.p.value=TRUE)$p.value 
p.value
TNM_N<-as.data.frame(TNM_N)
colnames(TNM_N)=c("Group","clinc","Freq")
TNM_N$sub=paste("TNM_N",p.value)
##############
Age<-table(data$group,data$Age)
p.value<-chisq.test(Age,simulate.p.value=TRUE)$p.value 
p.value
Age<-as.data.frame(Age)
colnames(Age)=c("Group","clinc","Freq")
Age$sub=paste("Age",p.value)
#################
Stage<-table(data$group,data$Stage)
p.value<-chisq.test(Stage,simulate.p.value=TRUE)$p.value 
p.value
Stage<-as.data.frame(Stage)
colnames(Stage)=c("Group","clinc","Freq")
Stage$sub=paste("Stage",p.value)

#########################
clinc.bar<-rbind(TNM_T,TNM_M,TNM_N,Age,Stage)
p.TNM_T<-ggplot(data=clinc.bar,aes(x=Group,y=Freq,fill=clinc))+
  #geom_bar(stat="identity")+
  geom_col(stat="identity",position="fill")+coord_polar(theta = "x")
p.TNM_T<-p.TNM_T+theme_bw() +theme(text=element_text(size=14,face="plain",color="black"),
        axis.title=element_text(size=16,face="plain",color="black"),
        axis.text = element_text(size=14,face="plain",color="black"),
        legend.title = element_text(size=14,face="plain",color="black"),
        legend.text = element_text(size=10,face="plain",color="black"),#"#9ec417"
        legend.background = element_blank(),
        legend.position="bottom")+ scale_fill_manual(values=c("#cc340c","#e8490f","#f18800","#FF0033","#3f60aa","skyblue","#cc340c","#e8490f","#f18800","#FF0033","#3f60aa","skyblue","#cc340c","#e8490f","#f18800","#FF0033"
))+
#+ scale_fill_manual(values=c("#cc340c","#e8490f","#f18800","#FF00CC","#FF0033","#3f60aa","skyblue","#cc340c","#e8490f","#f18800","#FF00CC","#3f60aa","skyblue","#3f60aa","skyblue","#cc340c","#e8490f","#f18800","#FF00CC","#FF0033"
#))+
  ylab("Percent")+ 
  xlab("Group")+facet_wrap(.~sub,ncol=5)
p.TNM_T


#########################################6. 使用ESTIMATE算法计算不同亚型的Stromal评分、immune评分和ESTIMATE评分，并比较亚型之间的评分的差异
#########################
my_comparisons <- list(c("cluster1","cluster2"))
estimate<-read.table("/mnt/nas_002/home/limma/02.project/04.BRCA/00.data/00.TCGA/breast_cancer_RNAseqV2.txt",header = TRUE, sep="\t")
estimate$ID=substr(estimate$ID,1,15)
group.new<-group
group.new$ID=substr(rownames(group.new),1,15)

estimate<-merge(estimate,group.new,by="ID")
estimate1<-estimate%>% pivot_longer(cols=-c("ID","group"), names_to = "type",values_to = "Score")
esti_p<- ggplot(data=estimate1,aes(x=type,y=Score,fill=group))+ geom_boxplot()+
     theme_bw()+
# stat_summary(mapping=aes(x=type,y=Score,color=group),fun="mean",geom="point",shape=21,size=3, fill="grey",position=position_dodge(width=0.8))+
    stat_compare_means(aes(group=group),label="p.signif", method = "wilcox.test"
                       #label="p.signif",label.y = c(0,0.5)
    )+
    theme(text=element_text(size=14,face="plain",color="black"),
          axis.title=element_text(size=16,face="plain",color="black"),
          axis.text = element_text(size=12,face="plain",color="black",angle=45,hjust=0.9),
          legend.title = element_text(size=16,face="plain",color="black"),
          legend.text = element_text(size=12,face="plain",color="black"),
          legend.background = element_blank(),
          legend.position="bottom")+
    #ylim(-3, 6) +
 ylab("Score")+
    xlab("")+
    scale_fill_manual(values = c("#FF3030","#0000FF"))
  #+labs(title="wilcox.test")
options(repr.plot.height=6, repr.plot.width=6)
esti_p

###################################7.免疫细胞的浸润分析
####################免疫评分
infil<-read.table("/mnt/nas_002/home/limma/02.project/02.COAD/infiltration_estimation_for_tcga.csv",header = TRUE, sep=",",row.names = 1)
infil<-dplyr::select(infil,ends_with("CIBERSORT"))
infil$ID<-row.names(infil)
infil.group<-group
infil.group$ID<-substr(rownames(infil.group),1,15)

Tide_infil<-merge(infil,infil.group,by="ID")%>% 
  pivot_longer(cols=matches("CIBERSORT"), names_to = "Celltype",values_to = "Fraction")
Tide_infil$Celltype<-gsub("_CIBERSORT", "",Tide_infil$Celltype)
tme_p<- ggplot(data=Tide_infil,aes(x=Celltype,y=Fraction,fill=group))+ geom_boxplot()+
     theme_bw()+
 #stat_summary(mapping=aes(x=Celltype,y=Fraction,color=group),fun="mean",geom="point",shape=21,size=3, fill="grey",position=position_dodge(width=0.8))+
    stat_compare_means(aes(group=group),label="p.signif", method = "wilcox.test"
                       #label="p.signif",label.y = c(0,0.5)
    )+
    theme(text=element_text(size=14,face="plain",color="black"),
          axis.title=element_text(size=16,face="plain",color="black"),
          axis.text = element_text(size=12,face="plain",color="black",angle=45,hjust=0.9),
          legend.title = element_text(size=16,face="plain",color="black"),
          legend.text = element_text(size=12,face="plain",color="black"),
          legend.background = element_blank(),
          legend.position="bottom")+
    #ylim(-3, 6) +
 ylab("Score")+
    xlab("")+
    scale_fill_manual(values = c("#FF3030","#0000FF"))
options(repr.plot.height=8, repr.plot.width=12)
tme_p

################################8.免疫检查点基因的表达分析
check.gene<-gene[c("CTLA4","HAVCR2","IL2RA","LAG3","PDCD1","TIGIT"),]
check.gene<-as.data.frame(t(check.gene))
check.gene$ID<-substr(rownames(check.gene),1,15)
check.gene<-merge(check.gene,group.new,by="ID")
check.gene<-check.gene %>% pivot_longer(cols=-c("ID","group"), names_to = "Gene",values_to = "Fpkm")



p.check.gene<-ggplot(data=check.gene,aes(x=Gene,y=log2(Fpkm+0.1),fill=group))+ geom_boxplot()+theme_bw()+stat_compare_means(aes(group=group),label="p.signif", method = "wilcox.test" )

p.check.gene<-p.check.gene+ theme_bw()+
  theme(text=element_text(size=14,face="plain",color="black"),
        axis.title=element_text(size=16,face="plain",color="black"),
        axis.text = element_text(size=14,face="plain",color="black",hjust=0.9),
        axis.text.x = element_text(angle=30),
        legend.title = element_text(size=10,face="plain",color="black"),
        legend.text = element_text(size=10,face="plain",color="black"),
        legend.background = element_blank(),
        legend.position="bottom")+
 # facet_wrap(~Gene,nrow=1)+
  #ylim(0, 0.5) +
  ylab("FPKM")+ 
  xlab("")+
  scale_fill_manual(values = c("#FF3030","#0000FF"))
#+labs(title="wilcox.test")

options(repr.plot.height=5, repr.plot.width=8)
p.check.gene

##############################9.TIDE分析
####TIDE箱线图
group.test<-group
group.test$Patient<-rownames(group)
tide<-read.table("/mnt/nas_002/home/limma/02.project/04.BRCA/00.data/00.TCGA/TIDE.csv",header=TRUE,sep=",")
tide$ID=substr(tide$Patient,14,15)
tide<-dplyr::filter(tide,ID=="01") %>% dplyr::select(Patient,TIDE,Responder)
#tide$ID=substr(tide$Patient,1,12)
tide<-inner_join(tide,group.test,by="Patient")
tide[which(tide$Responder=="False"),]$Responder="No-Response"
tide[which(tide$Responder=="True"),]$Responder="Response"
tide_p<- ggplot(data=tide,aes(x=group,y=TIDE,fill=group))+ geom_boxplot()+
     theme_bw()+
 #stat_summary(mapping=aes(x=group,y=TIDE,color=group),fun="mean",geom="point",shape=21,size=3, fill="grey",position=position_dodge(width=0.8))+
    stat_compare_means(aes(group=group),label="p.signif", method = "wilcox.test"
                       #label="p.signif",label.y = c(0,0.5)
    )+
    theme(text=element_text(size=14,face="plain",color="black"),
          axis.title=element_text(size=16,face="plain",color="black"),
          axis.text = element_text(size=12,face="plain",color="black",angle=45,hjust=0.9),
          legend.title = element_text(size=16,face="plain",color="black"),
          legend.text = element_text(size=12,face="plain",color="black"),
          legend.background = element_blank(),
          legend.position="bottom")+
    #ylim(-3, 6) +
 ylab("Score")+
    xlab("")+
    scale_fill_manual(values = c("#FF3030","#0000FF"))
options(repr.plot.height=5, repr.plot.width=5)
tide_p
ggsave("/mnt/nas_002/home/limma/02.project/04.BRCA/03.CONs/Figure1.tide.pdf",tide_p,height=5,width=5)

#####TIDE柱形图
tide_bar<-table(tide$group,tide$Responder)
p.value<-chisq.test(tide_bar,simulate.p.value=TRUE)$p.value 
paste("P.value=",p.value,sep="")
tide_bar<-as.data.frame(tide_bar)
colnames(tide_bar)=c("Group","Response","Freq")
p.tide<-ggplot(data=tide_bar,aes(x=Group,y=Freq,fill=Response))+
  #geom_bar(stat="identity")+
  geom_bar(stat="identity",position="fill")
p.tide<-p.tide+theme_bw() +theme(text=element_text(size=14,face="plain",color="black"),
        axis.title=element_text(size=16,face="plain",color="black"),
        axis.text = element_text(size=14,face="plain",color="black"),
        legend.title = element_text(size=14,face="plain",color="black"),
        legend.text = element_text(size=10,face="plain",color="black"),
        legend.background = element_blank(),
        legend.position="bottom")+
  ylab("Percent")+ 
  xlab("tide")
p.tide<-p.tide+labs(title=paste("Chisq test P.value=",round(p.value, 4),sep = ""))+
  theme(plot.title = element_text(hjust = 0.5))+scale_fill_manual(values = c("orange","purple"))
options(repr.plot.height=5, repr.plot.width=5)
p.tide



##########################################10.使用GSVA方法分析不同亚型之间的通路富集差异
###################
#########1)数据准备工作
library(GSVA)
library(GSEABase)
c7<- getGmt("/mnt/nas_002/home/limma/02.project/02.COAD/01.analysis/MSigDB/HALLMARK.symbols.gmt.txt")
######这里 kcdf 参数设为"Poisson"是因为使用read count数据，
###如果是使用 log 后的CPM, RPKM, TPM等数据就用默认值"Gaussian"。
###参数 parrallel.sz 设置并行线程数，因为每个基因集的计算是独立的。
TCGA.BRCA<-readRDS("/mnt/nas_002/home/limma/02.project/04.BRCA/00.data/00.TCGA//TCGA.BRCA.rds")
all<-TCGA.BRCA$tpm$mRNA
all<-all[,unique(colnames(all))]

MT.gene<-c(grep("^MT-",rownames(all),value=T),grep("^RP[SL]",rownames(all),value=T))
NMT.gene<-rownames(all)[!(rownames(all) %in% MT.gene)]
all<-all[NMT.gene,]
MT.gene<-c(grep("^MT-",rownames(all),value=T),grep("^RP[SL]",rownames(all),value=T))
NMT.gene<-rownames(all)[!(rownames(all) %in% MT.gene)]
all<-all[NMT.gene,]


gsva_data<- gsva(expr=log(as.matrix(all)+1), gset.idx.list=c7, kcdf="Gaussian", parallel.sz=6, method="ssgsea")
gsva_data<-as.data.frame(gsva_data)

###########
group<-dplyr::arrange(group,group)
gsva_data<-gsva_data[,rownames(group)]
rownames(gsva_data)<-gsub("HALLMARK_","",rownames(gsva_data))
bk <- c(seq(-0.2,0.2,by=0.01))
anno_col<- dplyr::select(group,group)
#c("skyblue","tomato","green3")
#"#FF3030","#0000FF"
ann_colors <- list("group" =c("cluster1"="#FF3030","cluster2"="#0000FF"))
p<-pheatmap::pheatmap(gsva_data, scale="row",
                breaks=bk,fontsize_row = 8, annotation_col = anno_col,annotation_colors=ann_colors,show_colnames = F, show_rownames = T,cluster_rows = TRUE,cluster_cols = F,
                color = colorRampPalette(c("blue","white","red"))(50),treeheight_row=0)
options(repr.plot.height=10, repr.plot.width=15)
p

##########################################11."NFKBIA","PTEN","TPI1","GAPDH"等基因在不同亚型中的差异表达分析
group<-dplyr::arrange(group,group)
gsva_data<-gsva_data[,rownames(group)]
rownames(gsva_data)<-gsub("HALLMARK_","",rownames(gsva_data))
bk <- c(seq(-0.2,0.2,by=0.01))
anno_col<- dplyr::select(group,group)
#c("skyblue","tomato","green3")
#"#FF3030","#0000FF"
ann_colors <- list("group" =c("cluster1"="#FF3030","cluster2"="#0000FF"))
p<-pheatmap::pheatmap(gsva_data, scale="row",
                breaks=bk,fontsize_row = 8, annotation_col = anno_col,annotation_colors=ann_colors,show_colnames = F, show_rownames = T,cluster_rows = TRUE,cluster_cols = F,
                color = colorRampPalette(c("blue","white","red"))(50),treeheight_row=0)
options(repr.plot.height=10, repr.plot.width=15)
p
########################################12."NFKBIA","PTEN","TPI1","GAPDH"与NSUN2基因表达的相关性分析
gene.TCGA<-gene[c("NFKBIA","PTEN","TPI1","GAPDH","NSUN2"),]
gene.TCGA<-as.data.frame(t(gene.TCGA))
#gene.TCGA$ID<-substr(rownames(gene.TCGA),1,15)
cor<-cor.test(as.numeric(gene.TCGA$NSUN2),as.numeric(gene.TCGA$NFKBIA),method = "pearson")
cor
NFKBIA.data<-gene.TCGA[,c("NFKBIA","NSUN2")]
NFKBIA.data$cor=paste0("NFKBIA",signif(cor$p.value,3),sep="")
colnames(NFKBIA.data)<-c("gene","NSUN2","cor")

cor<-cor.test(as.numeric(gene.TCGA$NSUN2),as.numeric(gene.TCGA$PTEN),method = "pearson")
cor
PTEN.data<-gene.TCGA[,c("PTEN","NSUN2")]
PTEN.data$cor=paste0("PTEN",signif(cor$p.value,3),sep="")
colnames(PTEN.data)<-c("gene","NSUN2","cor")
cor<-cor.test(as.numeric(gene.TCGA$NSUN2),as.numeric(gene.TCGA$TPI1),method = "pearson")
cor
TPI1.data<-gene.TCGA[,c("TPI1","NSUN2")]
TPI1.data$cor=paste0("TPI1",signif(cor$p.value,3),sep="")
colnames(TPI1.data)<-c("gene","NSUN2","cor")
cor<-cor.test(as.numeric(gene.TCGA$NSUN2),as.numeric(gene.TCGA$GAPDH),method = "pearson")
cor
GAPDH.data<-gene.TCGA[,c("GAPDH","NSUN2")]
GAPDH.data$cor=paste0("GAPDH",signif(cor$p.value,3),sep="")
colnames(GAPDH.data)<-c("gene","NSUN2","cor")
cor.data<-rbind(NFKBIA.data,PTEN.data,TPI1.data,GAPDH.data)
TCGA.cor<-ggplot(cor.data,aes(cor.data$NSUN2,cor.data$gene))+geom_point(size=3,colour = 'tomato3')+geom_smooth(method="lm")+#geom_vline(xintercept = 40)+
 # geom_text(aes(x=3,y=5),label=paste0("TCGA\n",'Cor:',signif(cor$estimate[[1]],3),"\n",'p.value:',
  #                                   signif(cor$p.value,3) ),size=5,,hjust=0.3,color='red' )+
    theme_bw()+theme(
      axis.text.x=element_text(colour="black",size=15,hjust = 0.5, vjust = 0.5),
      panel.background = element_blank(),axis.ticks.x = element_blank(),
     axis.text.y=element_text(colour="black",size=14),
      #x,y轴标题字体
      axis.title.x=element_text(size=14,colour="black"),
      axis.title.y=element_text(size=14,colour="black"),
      legend.title = element_text(colour="black", size=14),
      legend.text = element_text(colour="black", size = 14),
      legend.background = element_rect(fill="white", size=2, linetype="dotted")
    )+labs(x ="NSUN2 expression" ,y="expression")+facet_wrap(.~cor,ncol=4,scales = "free_y")
options(repr.plot.height=4, repr.plot.width=12)
TCGA.cor

