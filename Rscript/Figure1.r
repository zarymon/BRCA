library(Seurat)
library(ggplot2)
library(cowplot)
library(tidyverse)
library("clusterProfiler")
library(data.table)
library(openxlsx)
library("org.Hs.eg.db")
source("/mnt/nas_002/home/limma/00.Test.Script/Soure/SC.r")
library(ggrepel)
library(patchwork)
library(ggrepel)
library(ggpubr)
library(ggsignif)

#####1.数据处理
TCGA.BRCA<-readRDS("/mnt/nas_002/home/limma/02.project/04.BRCA/00.data/00.TCGA//TCGA.BRCA.rds")
TCGA.BRCA.fpkm<-TCGA.BRCA$fpkm$mRNA  ########表达数据
TCGA.BRCA.fpkm<-TCGA.BRCA.fpkm[,unique(colnames(TCGA.BRCA.fpkm))]
NSUN.FPKM<-TCGA.BRCA.fpkm[c("NOP2","NSUN2","NSUN3","NSUN4","NSUN5","NSUN6","NSUN7"),]
rownames(NSUN.FPKM)
rownames(NSUN.FPKM)<-c("NSUN1","NSUN2","NSUN3","NSUN4","NSUN5","NSUN6","NSUN7")

patient<-TCGA.BRCA$TCGA.clinc    ########临床信息数据
NSUNS.fpkm.clinc<-as.data.frame(t(NSUN.FPKM))
rownames(NSUNS.fpkm.clinc)<-colnames(NSUN.FPKM)
NSUNS.fpkm.clinc$treat<-substr(rownames(NSUNS.fpkm.clinc),14,15)
NSUNS.fpkm.clinc<-dplyr::filter(NSUNS.fpkm.clinc,treat %in% c("01"))
NSUNS.fpkm.clinc$ID=substr(rownames(NSUNS.fpkm.clinc),1,12)
NSUNS.fpkm.clinc$SampleID=substr(rownames(NSUNS.fpkm.clinc),1,15)
NSUNS.fpkm.clinc<-merge(patient,NSUNS.fpkm.clinc,by="ID")
##########TNM_T的临床信息画图
NSUNS.fpkm.clinc.TNM_T<-dplyr::select(NSUNS.fpkm.clinc,c("TNM_T","NSUN1","NSUN2","NSUN3","NSUN4", "NSUN5","NSUN6","NSUN7")) %>%dplyr::filter(TNM_T %in% c("T1","T2","T3","T4"))%>% dplyr::mutate(TNM_T1= ifelse(TNM_T %in% c("T1"),"T1","T2-T4"))%>% pivot_longer(cols=-c("TNM_T","TNM_T1"), names_to = "Gene",values_to = "FPKM")
NSUNS.fpkm.clinc.TNM_T$TNM_T<-NSUNS.fpkm.clinc.TNM_T$TNM_T1
table(NSUNS.fpkm.clinc.TNM_T$TNM_T)
NSUNS.TNM_T.draw<- ggplot(data=NSUNS.fpkm.clinc.TNM_T,aes(x=Gene,y=log2(FPKM+0.01),fill=TNM_T))+ geom_boxplot()+
     theme_bw()+
 stat_summary(mapping=aes(x=Gene,y=log2(FPKM+0.01),color=TNM_T),fun="mean",geom="point",shape=21,size=3, fill="red",position=position_dodge(width=0.8))+
    stat_compare_means(aes(group=TNM_T),label="p.signif", method = "wilcox.test"
                       #label="p.signif",label.y = c(0,0.5)
    )+theme(text=element_text(size=14,face="plain",color="black"),
          axis.title=element_text(size=16,face="plain",color="black"),
          axis.text = element_text(size=12,face="plain",color="black",angle=45,hjust=0.9),
          legend.title = element_text(size=16,face="plain",color="black"),
          legend.text = element_text(size=12,face="plain",color="black"),
          legend.background = element_blank(),
          legend.position="bottom")+
    
    #ylim(-3, 6) +
 ylab("Log2(FPKM+0.01)")+
    xlab("")+
    scale_fill_manual(values = c("skyblue","tomato","green3","purple"))
  #+labs(title="wilcox.test")
options(repr.plot.height=5, repr.plot.width=5)
NSUNS.TNM_T.draw

#############TNM_M的临床信息
NSUNS.fpkm.clinc.TNM_M<-dplyr::select(NSUNS.fpkm.clinc,c("TNM_M","NSUN1","NSUN2","NSUN3","NSUN4", "NSUN5","NSUN6","NSUN7"))%>%dplyr::filter(TNM_M %in% c("M0","M1"))%>% pivot_longer(cols=-c("TNM_M"), names_to = "Gene",values_to = "FPKM")
table(NSUNS.fpkm.clinc.TNM_M$TNM_M)
NSUNS.TNM_M.draw<- ggplot(data=NSUNS.fpkm.clinc.TNM_M,aes(x=Gene,y=log2(FPKM+0.01),fill=TNM_M))+ geom_boxplot()+
     theme_bw()+
 stat_summary(mapping=aes(x=Gene,y=log2(FPKM+0.01),color=TNM_M),fun="mean",geom="point",shape=21,size=3, fill="red",position=position_dodge(width=0.8))+
    stat_compare_means(aes(group=TNM_M),label="p.signif", method = "wilcox.test"
                       #label="p.signif",label.y = c(0,0.5)
    )+theme(text=element_text(size=14,face="plain",color="black"),
          axis.title=element_text(size=16,face="plain",color="black"),
          axis.text = element_text(size=12,face="plain",color="black",angle=45,hjust=0.9),
          legend.title = element_text(size=16,face="plain",color="black"),
          legend.text = element_text(size=12,face="plain",color="black"),
          legend.background = element_blank(),
          legend.position="bottom")+
    
    #ylim(-3, 6) +
 ylab("Log2(FPKM+0.01)")+
    xlab("")+
    scale_fill_manual(values = c("skyblue","tomato","green3","purple"))
  #+labs(title="wilcox.test")
options(repr.plot.height=5, repr.plot.width=5)
NSUNS.TNM_M.draw

##############TNM_N的临床信息画图
NSUNS.fpkm.clinc.TNM_N<-dplyr::select(NSUNS.fpkm.clinc,c("TNM_N","NSUN1","NSUN2","NSUN3","NSUN4", "NSUN5","NSUN6","NSUN7"))%>%dplyr::filter(TNM_N %in% c("N0","N1","N2","N3"))%>% dplyr::mutate(TNM_N1= ifelse(TNM_N %in% c("N0"),"N0","N1-N3"))%>% pivot_longer(cols=-c("TNM_N","TNM_N1"), names_to = "Gene",values_to = "FPKM")
NSUNS.fpkm.clinc.TNM_N$TNM_N<-NSUNS.fpkm.clinc.TNM_N$TNM_N1
table(NSUNS.fpkm.clinc.TNM_N$TNM_N)
NSUNS.TNM_N.draw<- ggplot(data=NSUNS.fpkm.clinc.TNM_N,aes(x=Gene,y=log2(FPKM+0.01),fill=TNM_N))+ geom_boxplot()+
     theme_bw()+
 stat_summary(mapping=aes(x=Gene,y=log2(FPKM+0.01),color=TNM_N),fun="mean",geom="point",shape=21,size=3, fill="red",position=position_dodge(width=0.8))+
    stat_compare_means(aes(group=TNM_N),label="p.signif", method = "wilcox.test"
                       #label="p.signif",label.y = c(0,0.5)
    )+theme(text=element_text(size=14,face="plain",color="black"),
          axis.title=element_text(size=16,face="plain",color="black"),
          axis.text = element_text(size=12,face="plain",color="black",angle=45,hjust=0.9),
          legend.title = element_text(size=16,face="plain",color="black"),
          legend.text = element_text(size=12,face="plain",color="black"),
          legend.background = element_blank(),
          legend.position="bottom")+
    
    #ylim(-3, 6) +
 ylab("Log2(FPKM+0.01)")+
    xlab("")+
    scale_fill_manual(values = c("skyblue","tomato","green3","purple"))
  #+labs(title="wilcox.test")
options(repr.plot.height=5, repr.plot.width=5)
NSUNS.TNM_N.draw

##############Stage的临床信息画图
NSUNS.fpkm.clinc.Stage<-dplyr::select(NSUNS.fpkm.clinc,c("Stage","NSUN1","NSUN2","NSUN3","NSUN4", "NSUN5","NSUN6","NSUN7"))
NSUNS.fpkm.clinc.Stage[which(NSUNS.fpkm.clinc.Stage$Stage=="Stage IA"), ]$Stage='Stage I'
NSUNS.fpkm.clinc.Stage[which(NSUNS.fpkm.clinc.Stage$Stage=="Stage IB"), ]$Stage='Stage I'
NSUNS.fpkm.clinc.Stage[which(NSUNS.fpkm.clinc.Stage$Stage=="Stage IIA"), ]$Stage='Stage II'
NSUNS.fpkm.clinc.Stage[which(NSUNS.fpkm.clinc.Stage$Stage=="Stage IIB"), ]$Stage='Stage II'
NSUNS.fpkm.clinc.Stage[which(NSUNS.fpkm.clinc.Stage$Stage=="Stage IIIA"), ]$Stage='Stage III'
NSUNS.fpkm.clinc.Stage[which(NSUNS.fpkm.clinc.Stage$Stage=="Stage IIIB"), ]$Stage='Stage III'
NSUNS.fpkm.clinc.Stage[which(NSUNS.fpkm.clinc.Stage$Stage=="Stage IIIC"), ]$Stage='Stage III'
table(NSUNS.fpkm.clinc.Stage$Stage)

NSUNS.fpkm.clinc.Stage<-NSUNS.fpkm.clinc.Stage%>%dplyr::filter(Stage %in% c("Stage I","Stage II","Stage III","Stage IV"))%>% dplyr::mutate(Stage1= ifelse(Stage %in% c("Stage I","Stage II"),"I&II","III&IV"))%>% pivot_longer(cols=-c("Stage","Stage1"), names_to = "Gene",values_to = "FPKM")
NSUNS.fpkm.clinc.Stage$Stage<-NSUNS.fpkm.clinc.Stage$Stage1
table(NSUNS.fpkm.clinc.Stage$Stage)
NSUNS.Stage.draw<- ggplot(data=NSUNS.fpkm.clinc.Stage,aes(x=Gene,y=log2(FPKM+0.01),fill=Stage))+ geom_boxplot()+
     theme_bw()+
 stat_summary(mapping=aes(x=Gene,y=log2(FPKM+0.01),color=Stage),fun="mean",geom="point",shape=21,size=3, fill="red",position=position_dodge(width=0.8))+
    stat_compare_means(aes(group=Stage),label="p.signif", method = "wilcox.test"
                       #label="p.signif",label.y = c(0,0.5)
    )+theme(text=element_text(size=14,face="plain",color="black"),
          axis.title=element_text(size=16,face="plain",color="black"),
          axis.text = element_text(size=12,face="plain",color="black",angle=45,hjust=0.9),
          legend.title = element_text(size=16,face="plain",color="black"),
          legend.text = element_text(size=12,face="plain",color="black"),
          legend.background = element_blank(),
          legend.position="bottom")+
    
    #ylim(-3, 6) +
 ylab("Log2(FPKM+0.01)")+
    xlab("")+
    scale_fill_manual(values = c("skyblue","tomato","green3","purple"))
  #+labs(title="wilcox.test")
options(repr.plot.height=5, repr.plot.width=5)
NSUNS.Stage.draw

##############HER2.status的临床信息画图
NSUNS.fpkm.clinc.HER2.status<-dplyr::select(NSUNS.fpkm.clinc,c("HER2.status","NSUN1","NSUN2","NSUN3","NSUN4", "NSUN5","NSUN6","NSUN7"))%>%dplyr::filter(HER2.status %in% c("Negative","Positive"))%>% pivot_longer(cols=-c("HER2.status"), names_to = "Gene",values_to = "FPKM")
table(NSUNS.fpkm.clinc.HER2.status$HER2.status)
NSUNS.HER2.status.draw<- ggplot(data=NSUNS.fpkm.clinc.HER2.status,aes(x=Gene,y=log2(FPKM+0.01),fill=HER2.status))+ geom_boxplot()+
     theme_bw()+
 stat_summary(mapping=aes(x=Gene,y=log2(FPKM+0.01),color=HER2.status),fun="mean",geom="point",shape=21,size=3, fill="red",position=position_dodge(width=0.8))+
    stat_compare_means(aes(group=HER2.status),label="p.signif", method = "wilcox.test"
                       #label="p.signif",label.y = c(0,0.5)
    )+theme(text=element_text(size=14,face="plain",color="black"),
          axis.title=element_text(size=16,face="plain",color="black"),
          axis.text = element_text(size=12,face="plain",color="black",angle=45,hjust=0.9),
          legend.title = element_text(size=16,face="plain",color="black"),
          legend.text = element_text(size=12,face="plain",color="black"),
          legend.background = element_blank(),
          legend.position="bottom")+
    
    #ylim(-3, 6) +
 ylab("Log2(FPKM+0.01)")+
    xlab("")+
    scale_fill_manual(values = c("skyblue","tomato","green3","purple"))
  #+labs(title="wilcox.test")
options(repr.plot.height=5, repr.plot.width=5)
NSUNS.HER2.status.draw

##############ER.status的临床信息画图
NSUNS.fpkm.clinc.ER.status<-dplyr::select(NSUNS.fpkm.clinc,c("ER.status","NSUN1","NSUN2","NSUN3","NSUN4", "NSUN5","NSUN6","NSUN7"))%>%dplyr::filter(ER.status %in% c("Negative","Positive"))%>% pivot_longer(cols=-c("ER.status"), names_to = "Gene",values_to = "FPKM")
table(NSUNS.fpkm.clinc.ER.status$ER.status)
NSUNS.ER.status.draw<- ggplot(data=NSUNS.fpkm.clinc.ER.status,aes(x=Gene,y=log2(FPKM+0.01),fill=ER.status))+ geom_boxplot()+
     theme_bw()+
 stat_summary(mapping=aes(x=Gene,y=log2(FPKM+0.01),color=ER.status),fun="mean",geom="point",shape=21,size=3, fill="red",position=position_dodge(width=0.8))+
    stat_compare_means(aes(group=ER.status),label="p.signif", method = "wilcox.test"
                       #label="p.signif",label.y = c(0,0.5)
    )+theme(text=element_text(size=14,face="plain",color="black"),
          axis.title=element_text(size=16,face="plain",color="black"),
          axis.text = element_text(size=12,face="plain",color="black",angle=45,hjust=0.9),
          legend.title = element_text(size=16,face="plain",color="black"),
          legend.text = element_text(size=12,face="plain",color="black"),
          legend.background = element_blank(),
          legend.position="bottom")+
    
    #ylim(-3, 6) +
 ylab("Log2(FPKM+0.01)")+
    xlab("")+
    scale_fill_manual(values = c("skyblue","tomato","green3","purple"))
  #+labs(title="wilcox.test")
options(repr.plot.height=5, repr.plot.width=5)
NSUNS.ER.status.draw

############PR.status

NSUNS.fpkm.clinc.PR.status<-dplyr::select(NSUNS.fpkm.clinc,c("PR.status","NSUN1","NSUN2","NSUN3","NSUN4", "NSUN5","NSUN6","NSUN7"))%>%dplyr::filter(PR.status %in% c("Negative","Positive"))%>% pivot_longer(cols=-c("PR.status"), names_to = "Gene",values_to = "FPKM")
table(NSUNS.fpkm.clinc.PR.status$PR.status)
NSUNS.PR.status.draw<- ggplot(data=NSUNS.fpkm.clinc.PR.status,aes(x=Gene,y=log2(FPKM+0.01),fill=PR.status))+ geom_boxplot()+
     theme_bw()+
 stat_summary(mapping=aes(x=Gene,y=log2(FPKM+0.01),color=PR.status),fun="mean",geom="point",shape=21,size=3, fill="red",position=position_dodge(width=0.8))+
    stat_compare_means(aes(group=PR.status),label="p.signif", method = "wilcox.test"
                       #label="p.signif",label.y = c(0,0.5)
    )+theme(text=element_text(size=14,face="plain",color="black"),
          axis.title=element_text(size=16,face="plain",color="black"),
          axis.text = element_text(size=12,face="plain",color="black",angle=45,hjust=0.9),
          legend.title = element_text(size=16,face="plain",color="black"),
          legend.text = element_text(size=12,face="plain",color="black"),
          legend.background = element_blank(),
          legend.position="bottom")+
    
    #ylim(-3, 6) +
 ylab("Log2(FPKM+0.01)")+
    xlab("")+
    scale_fill_manual(values = c("skyblue","tomato","green3","purple"))
  #+labs(title="wilcox.test")
options(repr.plot.height=5, repr.plot.width=5)
NSUNS.PR.status.draw

