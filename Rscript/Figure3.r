####################自定义函数
corr<-function(x,y){
  cor_result=data.frame()
  for (i in colnames(x) ) {
    
    cor<-cor.test(as.numeric(x[  ,i]),as.numeric(y),method = "pearson")
    cor_result<-rbind(cor_result,cbind(i,round(cor$estimate[[1]],3),cor$p.value) )
  }
  colnames(cor_result)<-c("Immune_cell","coefficient","pvalue")
  cor_result$pvalue<-as.numeric(cor_result$pvalue)
  cor_result$pstar<-ifelse(cor_result$pvalue<0.05,
                           ifelse(cor_result$pvalue<0.01,ifelse(cor_result$pvalue<0.001,"***","**")
                                 ,"*" )
                           
                           ,"")
  cor_result
  
}
################数据处理
gene.TCGA<-TCGA.BRCA.fpkm[c("CTLA4","HAVCR2","IL2RA","LAG3","PDCD1","CD274","TIGIT"),]
gene.TCGA<-as.data.frame(t(gene.TCGA))
gene.TCGA$treat<-substr(rownames(gene.TCGA),14,15)
gene.TCGA<-dplyr::filter(gene.TCGA,treat=="01")
gene.TCGA$SampleID<-substr(rownames(gene.TCGA),1,15)

###################高低分组的免疫检查点基因表达分析以及免疫检查点基因表达与NSUNs基因表达的相关性分析
NSUN1.level<-dplyr::select(NSUNS.fpkm.clinc,c("OS.Status","OS.time","NSUN1","SampleID","treat"))%>%dplyr::filter(treat %in% c("01"))%>%dplyr::select(-c(treat)) %>%dplyr::mutate(NSUN1.level=ifelse(NSUN1>median(NSUN1),"High NSUN1","Low NSUN1"))

Immun.check.NSUN1<-merge(gene.TCGA,NSUN1.level,by="SampleID")%>% dplyr::select(-c("treat","SampleID","NSUN1","OS.Status","OS.time"))%>% pivot_longer(cols=-c("NSUN1.level"), names_to = "Gene",values_to = "Fpkm")

p.Immun.check.NSUN1<-ggplot(data=Immun.check.NSUN1,aes(x=Gene,y=log2(Fpkm+0.1),fill=NSUN1.level))+ geom_boxplot()+
     theme_bw()+
 stat_summary(mapping=aes(x=Gene,y=log2(Fpkm+0.1),color=NSUN1.level),fun="mean",geom="point",shape=21,size=3, fill="red",position=position_dodge(width=0.8))+
    stat_compare_means(aes(group=NSUN1.level),label="p.signif", method = "wilcox.test"
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
 ylab("Gene expression)")+
    xlab("TCGA")+
    scale_fill_manual(values = c("skyblue","tomato","green3"))
options(repr.plot.height=5, repr.plot.width=8)
p.Immun.check.NSUN1
NSUN1.Immun.check.cor<-merge(gene.TCGA,NSUN1.level,by="SampleID")
NSUN1.Immun.check.cor.coef<-corr(NSUN1.Immun.check.cor[,2:(ncol(NSUN1.Immun.check.cor)-5)],NSUN1.Immun.check.cor$NSUN1)
NSUN1.Immun.check.cor.coef$Risk.score="NSUN1"
NSUN1.Immun.check.cor.coef$coefficient=as.numeric(NSUN1.Immun.check.cor.coef$coefficient)

######################
NSUN2.level<-dplyr::select(NSUNS.fpkm.clinc,c("OS.Status","OS.time","NSUN2","SampleID","treat"))%>%dplyr::filter(treat %in% c("01"))%>%dplyr::select(-c(treat)) %>%dplyr::mutate(NSUN2.level=ifelse(NSUN2>median(NSUN2),"High NSUN2","Low NSUN2"))
NSUN2.infil<-merge(infil_CI,NSUN2.level,by="SampleID")%>% 
  pivot_longer(cols=matches("CIBERSORT"), names_to = "Celltype",values_to = "Fraction")
NSUN2.infil$Celltype<-gsub("_CIBERSORT", "",NSUN2.infil$Celltype)
NSUN2.infil$Celltype<-gsub("Myeloid.dendritic.cell.resting","Dendritic.cell.resting",NSUN2.infil$Celltype)
NSUN2.infil$Celltype<-gsub("Myeloid.dendritic.cell.activated","Dendritic.cell.activated",NSUN2.infil$Celltype)
NSUN2.infil.draw<- ggplot(data=NSUN2.infil,aes(x=Celltype,y=Fraction,fill=NSUN2.level))+ geom_boxplot()+
     theme_bw()+
 stat_summary(mapping=aes(x=Celltype,y=Fraction,color=NSUN2.level),fun="mean",geom="point",shape=21,size=3, fill="red",position=position_dodge(width=0.8))+
    stat_compare_means(aes(group=NSUN2.level),label="p.signif", method = "wilcox.test"
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
    scale_fill_manual(values = c("skyblue","tomato","green3"))
options(repr.plot.height=8, repr.plot.width=12)
NSUN2.infil.draw
NSUN2.cor<-merge(infil_CI,NSUN2.level,by="SampleID")
NSUN2.cor.coef<-corr(NSUN2.cor[,2:(ncol(NSUN2.cor)-4)],NSUN2.cor$NSUN2)
NSUN2.cor.coef$Risk.score="NSUN2"
NSUN2.cor.coef$coefficient=as.numeric(NSUN2.cor.coef$coefficient)
NSUN2.cor.coef$Immune_cell<-gsub("_CIBERSORT", "",NSUN2.cor.coef$Immune_cell)
NSUN2.cor.coef$Immune_cell<-gsub("Myeloid.dendritic.cell.resting","Dendritic.cell.resting",NSUN2.cor.coef$Immune_cell)
NSUN2.cor.coef$Immune_cell<-gsub("Myeloid.dendritic.cell.activated","Dendritic.cell.activated",NSUN2.cor.coef$Immune_cell)
################
NSUN3.level<-dplyr::select(NSUNS.fpkm.clinc,c("OS.Status","OS.time","NSUN3","SampleID","treat"))%>%dplyr::filter(treat %in% c("01"))%>%dplyr::select(-c(treat)) %>%dplyr::mutate(NSUN3.level=ifelse(NSUN3>median(NSUN3),"High NSUN3","Low NSUN3"))
NSUN3.infil<-merge(infil_CI,NSUN3.level,by="SampleID")%>% 
  pivot_longer(cols=matches("CIBERSORT"), names_to = "Celltype",values_to = "Fraction")
NSUN3.infil$Celltype<-gsub("_CIBERSORT", "",NSUN3.infil$Celltype)
NSUN3.infil$Celltype<-gsub("Myeloid.dendritic.cell.resting","Dendritic.cell.resting",NSUN3.infil$Celltype)
NSUN3.infil$Celltype<-gsub("Myeloid.dendritic.cell.activated","Dendritic.cell.activated",NSUN3.infil$Celltype)
NSUN3.infil.draw<- ggplot(data=NSUN3.infil,aes(x=Celltype,y=Fraction,fill=NSUN3.level))+ geom_boxplot()+
     theme_bw()+
 stat_summary(mapping=aes(x=Celltype,y=Fraction,color=NSUN3.level),fun="mean",geom="point",shape=21,size=3, fill="red",position=position_dodge(width=0.8))+
    stat_compare_means(aes(group=NSUN3.level),label="p.signif", method = "wilcox.test"
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
    scale_fill_manual(values = c("skyblue","tomato","green3"))
options(repr.plot.height=8, repr.plot.width=12)
NSUN3.infil.draw
NSUN3.cor<-merge(infil_CI,NSUN3.level,by="SampleID")
NSUN3.cor.coef<-corr(NSUN3.cor[,2:(ncol(NSUN3.cor)-4)],NSUN3.cor$NSUN3)
NSUN3.cor.coef$Risk.score="NSUN3"
NSUN3.cor.coef$coefficient=as.numeric(NSUN3.cor.coef$coefficient)
NSUN3.cor.coef$Immune_cell<-gsub("_CIBERSORT", "",NSUN3.cor.coef$Immune_cell)
NSUN3.cor.coef$Immune_cell<-gsub("Myeloid.dendritic.cell.resting","Dendritic.cell.resting",NSUN3.cor.coef$Immune_cell)
NSUN3.cor.coef$Immune_cell<-gsub("Myeloid.dendritic.cell.activated","Dendritic.cell.activated",NSUN3.cor.coef$Immune_cell)
############
NSUN4.level<-dplyr::select(NSUNS.fpkm.clinc,c("OS.Status","OS.time","NSUN4","SampleID","treat"))%>%dplyr::filter(treat %in% c("01"))%>%dplyr::select(-c(treat)) %>%dplyr::mutate(NSUN4.level=ifelse(NSUN4>median(NSUN4),"High NSUN4","Low NSUN4"))
NSUN4.infil<-merge(infil_CI,NSUN4.level,by="SampleID")%>% 
  pivot_longer(cols=matches("CIBERSORT"), names_to = "Celltype",values_to = "Fraction")
NSUN4.infil$Celltype<-gsub("_CIBERSORT", "",NSUN4.infil$Celltype)
NSUN4.infil$Celltype<-gsub("Myeloid.dendritic.cell.resting","Dendritic.cell.resting",NSUN4.infil$Celltype)
NSUN4.infil$Celltype<-gsub("Myeloid.dendritic.cell.activated","Dendritic.cell.activated",NSUN4.infil$Celltype)
NSUN4.infil.draw<- ggplot(data=NSUN4.infil,aes(x=Celltype,y=Fraction,fill=NSUN4.level))+ geom_boxplot()+
     theme_bw()+
 stat_summary(mapping=aes(x=Celltype,y=Fraction,color=NSUN4.level),fun="mean",geom="point",shape=21,size=3, fill="red",position=position_dodge(width=0.8))+
    stat_compare_means(aes(group=NSUN4.level),label="p.signif", method = "wilcox.test"
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
    scale_fill_manual(values = c("skyblue","tomato","green3"))
options(repr.plot.height=8, repr.plot.width=12)
NSUN4.infil.draw
NSUN4.cor<-merge(infil_CI,NSUN4.level,by="SampleID")
NSUN4.cor.coef<-corr(NSUN4.cor[,2:(ncol(NSUN4.cor)-4)],NSUN4.cor$NSUN4)
NSUN4.cor.coef$Risk.score="NSUN4"
NSUN4.cor.coef$coefficient=as.numeric(NSUN4.cor.coef$coefficient)
NSUN4.cor.coef$Immune_cell<-gsub("_CIBERSORT", "",NSUN4.cor.coef$Immune_cell)
NSUN4.cor.coef$Immune_cell<-gsub("Myeloid.dendritic.cell.resting","Dendritic.cell.resting",NSUN4.cor.coef$Immune_cell)
NSUN4.cor.coef$Immune_cell<-gsub("Myeloid.dendritic.cell.activated","Dendritic.cell.activated",NSUN4.cor.coef$Immune_cell)
#############################
NSUN5.level<-dplyr::select(NSUNS.fpkm.clinc,c("OS.Status","OS.time","NSUN5","SampleID","treat"))%>%dplyr::filter(treat %in% c("01"))%>%dplyr::select(-c(treat)) %>%dplyr::mutate(NSUN5.level=ifelse(NSUN5>median(NSUN5),"High NSUN5","Low NSUN5"))
NSUN5.infil<-merge(infil_CI,NSUN5.level,by="SampleID")%>% 
  pivot_longer(cols=matches("CIBERSORT"), names_to = "Celltype",values_to = "Fraction")
NSUN5.infil$Celltype<-gsub("_CIBERSORT", "",NSUN5.infil$Celltype)
NSUN5.infil$Celltype<-gsub("Myeloid.dendritic.cell.resting","Dendritic.cell.resting",NSUN5.infil$Celltype)
NSUN5.infil$Celltype<-gsub("Myeloid.dendritic.cell.activated","Dendritic.cell.activated",NSUN5.infil$Celltype)
NSUN5.infil.draw<- ggplot(data=NSUN5.infil,aes(x=Celltype,y=Fraction,fill=NSUN5.level))+ geom_boxplot()+
     theme_bw()+
 stat_summary(mapping=aes(x=Celltype,y=Fraction,color=NSUN5.level),fun="mean",geom="point",shape=21,size=3, fill="red",position=position_dodge(width=0.8))+
    stat_compare_means(aes(group=NSUN5.level),label="p.signif", method = "wilcox.test"
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
    scale_fill_manual(values = c("skyblue","tomato","green3"))
options(repr.plot.height=8, repr.plot.width=12)
NSUN5.infil.draw
NSUN5.cor<-merge(infil_CI,NSUN5.level,by="SampleID")
NSUN5.cor.coef<-corr(NSUN5.cor[,2:(ncol(NSUN5.cor)-4)],NSUN5.cor$NSUN5)
NSUN5.cor.coef$Risk.score="NSUN5"
NSUN5.cor.coef$coefficient=as.numeric(NSUN5.cor.coef$coefficient)
NSUN5.cor.coef$Immune_cell<-gsub("_CIBERSORT", "",NSUN5.cor.coef$Immune_cell)
NSUN5.cor.coef$Immune_cell<-gsub("Myeloid.dendritic.cell.resting","Dendritic.cell.resting",NSUN5.cor.coef$Immune_cell)
NSUN5.cor.coef$Immune_cell<-gsub("Myeloid.dendritic.cell.activated","Dendritic.cell.activated",NSUN5.cor.coef$Immune_cell)
#####################
NSUN6.level<-dplyr::select(NSUNS.fpkm.clinc,c("OS.Status","OS.time","NSUN6","SampleID","treat"))%>%dplyr::filter(treat %in% c("01"))%>%dplyr::select(-c(treat)) %>%dplyr::mutate(NSUN6.level=ifelse(NSUN6>median(NSUN6),"High NSUN6","Low NSUN6"))
NSUN6.infil<-merge(infil_CI,NSUN6.level,by="SampleID")%>% 
  pivot_longer(cols=matches("CIBERSORT"), names_to = "Celltype",values_to = "Fraction")
NSUN6.infil$Celltype<-gsub("_CIBERSORT", "",NSUN6.infil$Celltype)
NSUN6.infil$Celltype<-gsub("Myeloid.dendritic.cell.resting","Dendritic.cell.resting",NSUN6.infil$Celltype)
NSUN6.infil$Celltype<-gsub("Myeloid.dendritic.cell.activated","Dendritic.cell.activated",NSUN6.infil$Celltype)
NSUN6.infil.draw<- ggplot(data=NSUN6.infil,aes(x=Celltype,y=Fraction,fill=NSUN6.level))+ geom_boxplot()+
     theme_bw()+
 stat_summary(mapping=aes(x=Celltype,y=Fraction,color=NSUN6.level),fun="mean",geom="point",shape=21,size=3, fill="red",position=position_dodge(width=0.8))+
    stat_compare_means(aes(group=NSUN6.level),label="p.signif", method = "wilcox.test"
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
    scale_fill_manual(values = c("skyblue","tomato","green3"))
options(repr.plot.height=8, repr.plot.width=12)
NSUN6.infil.draw
NSUN6.cor<-merge(infil_CI,NSUN6.level,by="SampleID")
NSUN6.cor.coef<-corr(NSUN6.cor[,2:(ncol(NSUN6.cor)-4)],NSUN6.cor$NSUN6)
NSUN6.cor.coef$Risk.score="NSUN6"
NSUN6.cor.coef$coefficient=as.numeric(NSUN6.cor.coef$coefficient)
NSUN6.cor.coef$Immune_cell<-gsub("_CIBERSORT", "",NSUN6.cor.coef$Immune_cell)
NSUN6.cor.coef$Immune_cell<-gsub("Myeloid.dendritic.cell.resting","Dendritic.cell.resting",NSUN6.cor.coef$Immune_cell)
NSUN6.cor.coef$Immune_cell<-gsub("Myeloid.dendritic.cell.activated","Dendritic.cell.activated",NSUN6.cor.coef$Immune_cell)
#########################
NSUN7.level<-dplyr::select(NSUNS.fpkm.clinc,c("OS.Status","OS.time","NSUN7","SampleID","treat"))%>%dplyr::filter(treat %in% c("01"))%>%dplyr::select(-c(treat)) %>%dplyr::mutate(NSUN7.level=ifelse(NSUN7>median(NSUN7),"High NSUN7","Low NSUN7"))
NSUN7.infil<-merge(infil_CI,NSUN7.level,by="SampleID")%>% 
  pivot_longer(cols=matches("CIBERSORT"), names_to = "Celltype",values_to = "Fraction")
NSUN7.infil$Celltype<-gsub("_CIBERSORT", "",NSUN7.infil$Celltype)
NSUN7.infil$Celltype<-gsub("Myeloid.dendritic.cell.resting","Dendritic.cell.resting",NSUN7.infil$Celltype)
NSUN7.infil$Celltype<-gsub("Myeloid.dendritic.cell.activated","Dendritic.cell.activated",NSUN7.infil$Celltype)
NSUN7.infil.draw<- ggplot(data=NSUN7.infil,aes(x=Celltype,y=Fraction,fill=NSUN7.level))+ geom_boxplot()+
     theme_bw()+
 stat_summary(mapping=aes(x=Celltype,y=Fraction,color=NSUN7.level),fun="mean",geom="point",shape=21,size=3, fill="red",position=position_dodge(width=0.8))+
    stat_compare_means(aes(group=NSUN7.level),label="p.signif", method = "wilcox.test"
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
    scale_fill_manual(values = c("skyblue","tomato","green3"))
options(repr.plot.height=8, repr.plot.width=12)
NSUN7.infil.draw
NSUN7.cor<-merge(infil_CI,NSUN7.level,by="SampleID")
NSUN7.cor.coef<-corr(NSUN7.cor[,2:(ncol(NSUN7.cor)-4)],NSUN7.cor$NSUN7)
NSUN7.cor.coef$Risk.score="NSUN7"
NSUN7.cor.coef$coefficient=as.numeric(NSUN7.cor.coef$coefficient)
NSUN7.cor.coef$Immune_cell<-gsub("_CIBERSORT", "",NSUN7.cor.coef$Immune_cell)
NSUN7.cor.coef$Immune_cell<-gsub("Myeloid.dendritic.cell.resting","Dendritic.cell.resting",NSUN7.cor.coef$Immune_cell)
NSUN7.cor.coef$Immune_cell<-gsub("Myeloid.dendritic.cell.activated","Dendritic.cell.activated",NSUN7.cor.coef$Immune_cell)

##############  相关性热图的绘制
NSUNS.cor.coef<-rbind(NSUN1.cor.coef,NSUN2.cor.coef,NSUN3.cor.coef,NSUN4.cor.coef,NSUN5.cor.coef,NSUN6.cor.coef,NSUN7.cor.coef)
corheat<-ggplot(NSUNS.cor.coef,aes(Immune_cell,Risk.score))+geom_tile(aes(fill=coefficient) ,colour="white",size=1)+ scale_fill_gradient2(low = "blue", mid="white", high = "red")+
  geom_text(aes(label=pstar),col="black",size=5)+
   theme_bw()+ylab("")+xlab("")+
   #theme_minimal+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
        axis.text.x = element_text(angle=45,hjust=1),axis.ticks.x=element_blank(), panel.border = element_blank(),text=element_text(size=14,face="plain",color="black"))+
# axis.text.x = element_text(angle=45,hjust=1)  
#axis.text.x = element_blank(),axis.ticks.x=element_blank()
labs(fill=paste0("* p<0.05","\n\n","** p<0.01","\n\n","*** p<0.001","\n\n","Correlation"))
options(repr.plot.height=5, repr.plot.width=10)
corheat
