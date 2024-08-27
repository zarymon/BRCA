source("/mnt/nas_002/home/limma/02.project/04.BRCA/learning.r")
#####################################对差异表达基因进行单因素回归分析
###############
#OS.DEG.gene
   Cluster.DEG.OS_mutate.test<-data.frame()
   for (i in OS.DEG.gene){
   gene.level<-as.data.frame(t(gene[i,]))
   colnames(gene.level)=c("gene")
   gene.level$ID=substr(colnames(gene),1,12)
   gene.level<-na.omit(gene.level)
    if (median(gene.level$gene)==0){
        next}
   gene.level<- dplyr::mutate(gene.level,level=ifelse(gene >median(gene.level$gene),"high","low"))
   OS_data<-merge(gene.level,TCGA.BRCA$TCGA.clinc) %>% dplyr::select(c("level","OS.Status","OS.time"))
   #high_OS_data<-dplyr::filter(OS_data,level=="high")
   #low_OS_data<-dplyr::filter(OS_data,level=="low")
      # fit <- survfit(Surv(OS.time,OS.Status) ~ level,  data = OS_data) #2. 拟合生存曲线
       km<-survdiff(formula = Surv(OS.time,OS.Status) ~ level,  data = OS_data)
   
       p.value <- 1 - pchisq(km$chisq,length(km$n)-1)
       t<-coxph(formula = Surv(OS.time, OS.Status) ~ level, data = OS_data)
       t$coefficients
   
       Cluster.DEG.OS_mutate.test<-rbind(Cluster.DEG.OS_mutate.test,
                         cbind(name=i,
                               pvalue=p.value,coefficient=t$coefficients
                               ))
       }
   ####
Cluster.DEG.OS.gene<-dplyr::filter(Cluster.DEG.OS_mutate.test,pvalue<0.05) 

###########3
#############################################
###定义排序
##############################################建模分析
library("glmnet") 
library("survival")
TCGA.jianm.clinc<-TCGA.BRCA$TCGA.clinc %>% dplyr::select(c(ID,OS.Status,OS.time))
TCGA.jianm.gene<-gene
colnames(TCGA.jianm.gene)<-substr(colnames(TCGA.jianm.gene),1,12)
common_sampleL <- intersect(colnames(TCGA.jianm.gene),TCGA.jianm.clinc$ID)
rownames(TCGA.jianm.clinc)<-TCGA.jianm.clinc$ID
TCGA.jianm.clinc <-TCGA.jianm.clinc[common_sampleL,,drop=F]
TCGA.jianm.gene<-t(TCGA.jianm.gene[,common_sampleL,drop=F])

########################

inTrain<-sample(1:length(common_sampleL),round(length(common_sampleL)*3/8))
Train.clinc<-TCGA.jianm.clinc[inTrain,]
Test.clinc<-TCGA.jianm.clinc[-inTrain,]
Train.gene<-TCGA.jianm.gene[inTrain,]
Test.gene<-TCGA.jianm.gene[-inTrain,]
identical(rownames(Train.clinc),rownames(Train.gene))
identical(rownames(Test.clinc),rownames(Test.gene))
################
Metabric.jianm.clinc<-Metabric.clinical%>%dplyr::select(c(Patient_ID,OS.Status,OS.time))
Metabric.gene<-Metabric.BRCA.fpkm
Metabric.comm<-intersect(Metabric.jianm.clinc$Patient_ID,colnames(Metabric.gene))
rownames(Metabric.jianm.clinc)<-Metabric.jianm.clinc$Patient_ID
Metabric.gene<-t(Metabric.gene[,Metabric.comm,drop=F])
Metabric.jianm.clinc <-Metabric.jianm.clinc[Metabric.comm,,drop=F]   ######即便是取1列，也能保持矩阵的数据结构了，colnames 和rownames也能用
Metabric.jianm.clinc$ID=rownames(Metabric.jianm.clinc)
Metabric.jianm.clinc$OS.time<-as.numeric(Metabric.jianm.clinc$OS.time)
Metabric.jianm.clinc<- dplyr::select(Metabric.jianm.clinc,c(ID,OS.Status,OS.time))

identical(rownames(Metabric.gene),rownames(Metabric.jianm.clinc))

#################################
GSE20685<-readRDS("/mnt/nas_002/home/limma/02.project/04.BRCA/00.data/GSE20685/GSE20685.rds")
GSE20685.clinc<-GSE20685$clin %>% dplyr::select(c("event_death:ch1","follow_up_duration (years):ch1"))
#event_death:ch1 follow_up_duration (years):ch1
colnames(GSE20685.clinc)=c("OS.Status","OS.time")
GSE20685.clinc<-dplyr::filter(GSE20685.clinc,OS.Status %in% c("0","1"))

GSE20685.gene<-GSE20685$expr
rownames(GSE20685.gene)=GSE20685.gene$"Gene Symbol"
GSE20685.common_sampleL <- intersect(colnames(GSE20685.gene),rownames(GSE20685.clinc))
GSE20685.gene<-t(GSE20685.gene[,GSE20685.common_sampleL,drop=F])
GSE20685.clinc <-GSE20685.clinc[GSE20685.common_sampleL,,drop=F]   ######即便是取1列，也能保持矩阵的数据结构了，colnames 和rownames也能用
GSE20685.clinc$ID=rownames(GSE20685.clinc)
GSE20685.clinc$OS.time<-as.numeric(GSE20685.clinc$OS.time)
GSE20685.clinc<- dplyr::select(GSE20685.clinc,c(ID,OS.Status,OS.time))

identical(rownames(GSE20685.gene),rownames(GSE20685.clinc))
##################
###################################
set<-intersect(Cluster.DEG.OS.gene$name,colnames(Metabric.gene))
set<-intersect(set,colnames(GSE20685.gene))
#set<-Cluster.DEG.OS.gene$name
setwd("/mnt/nas_002/home/limma/02.project/04.BRCA/01.Analysis/00.jianmo")
Train.clinc$OS.time<-as.numeric(Train.clinc$OS.time)
Train.y<-data.matrix(Surv(time=as.double(Train.clinc$OS.time),event=as.double(Train.clinc$OS.Status)))
Train.x<-Train.gene[,set]
###################
Test.clinc$OS.time<-as.numeric(Test.clinc$OS.time)
Test.y<-data.matrix(Surv(time=as.double(Test.clinc$OS.time),event=as.double(Test.clinc$OS.Status)))
Test.x<-Test.gene[,set]
###########
TCGA.jianm.clinc$OS.time<-as.numeric(TCGA.jianm.clinc$OS.time)
TCGA.jianm.y<-data.matrix(Surv(time=as.double(TCGA.jianm.clinc$OS.time),event=as.double(TCGA.jianm.clinc$OS.Status)))
TCGA.jianm.x<-TCGA.jianm.gene[,set]
###################
GSE20685.clinc$OS.time<-as.numeric(GSE20685.clinc$OS.time)
GSE20685.clinc$OS.Status<-as.numeric(GSE20685.clinc$OS.Status)
GSE20685.y<-data.matrix(Surv(time=as.double(GSE20685.clinc$OS.time),event=as.double(GSE20685.clinc$OS.Status)))
GSE20685.x<-GSE20685.gene[,set]
#######################
Metabric.jianm.clinc$OS.time<-as.numeric(Metabric.jianm.clinc$OS.time)
Metabric.jianm.clinc$OS.Status<-as.numeric(Metabric.jianm.clinc$OS.Status)
Metabric.y<-data.matrix(Surv(time=as.double(Metabric.jianm.clinc$OS.time),event=as.double(Metabric.jianm.clinc$OS.Status)))
Metabric.x<-Metabric.gene[,set]

###################
seed=1
Algorithm="LASSO_COX"
cvfit<-lasso(TCGA.jianm.x,TCGA.jianm.y,seed,Algorithm,"TCGA_TCGA.jianm",TCGA.jianm.clinc,Test.x,Test.y,"TCGA_Test",Test.clinc,GSE20685.x,GSE20685.y,"GSE20685",GSE20685.clinc,Metabric.x,Metabric.y,"Metabric",Metabric.jianm.clinc)
#######################
seed=123
Algorithm="randomForestSRC"
cvfit<-randomForestSRC.ar(TCGA.jianm.x,TCGA.jianm.clinc,Algorithm,"TCGA_TCGA.jianm",Test.x,Test.clinc,"TCGA_Test",GSE20685.x,GSE20685.clinc,"GSE20685",Metabric.x,Metabric.jianm.clinc,"Metabric")

#########################
Algorithm="Coxboost"
cvfit<-Coxboost(TCGA.jianm.x,TCGA.jianm.clinc,Algorithm,"TCGA_TCGA.jianm",Test.x,Test.clinc,"TCGA_Test",GSE20685.x,GSE20685.clinc,"GSE20685",Metabric.x,Metabric.jianm.clinc,"Metabric")
#######################
Algorithm="gbm"
cvfit<-gbm.ar(TCGA.jianm.x,TCGA.jianm.clinc,Algorithm,"TCGA_TCGA.jianm",Test.x,Test.clinc,"TCGA_Test",GSE20685.x,GSE20685.clinc,"GSE20685",Metabric.x,Metabric.jianm.clinc,"Metabric")

#######################
seed=123
Algorithm="Ridge"
cvfit<-ridge.ar(TCGA.jianm.x,TCGA.jianm.y,seed,Algorithm,"TCGA_TCGA.jianm",TCGA.jianm.clinc,Test.x,Test.y,"TCGA_Test",Test.clinc,GSE20685.x,GSE20685.y,"GSE20685",GSE20685.clinc,Metabric.x,Metabric.y,"Metabric",Metabric.jianm.clinc)
#######################
seed=1
Algorithm="Enet0.2"

cvfit<-Enet.ar(TCGA.jianm.x,TCGA.jianm.y,seed,Algorithm,"TCGA_TCGA.jianm",TCGA.jianm.clinc,Test.x,Test.y,"TCGA_Test",Test.clinc,GSE20685.x,GSE20685.y,"GSE20685",GSE20685.clinc,Metabric.x,Metabric.y,"Metabric",Metabric.jianm.clinc,"0.2")
#######################
seed=1
Algorithm="Enet0.5"
cvfit<-Enet.ar(TCGA.jianm.x,TCGA.jianm.y,seed,Algorithm,"TCGA_TCGA.jianm",TCGA.jianm.clinc,Test.x,Test.y,"TCGA_Test",Test.clinc,GSE20685.x,GSE20685.y,"GSE20685",GSE20685.clinc,Metabric.x,Metabric.y,"Metabric",Metabric.jianm.clinc,"0.5")
#######################
seed=1
Algorithm="Enet0.8"
cvfit<-Enet.ar(TCGA.jianm.x,TCGA.jianm.y,seed,Algorithm,"TCGA_TCGA.jianm",TCGA.jianm.clinc,Test.x,Test.y,"TCGA_Test",Test.clinc,GSE20685.x,GSE20685.y,"GSE20685",GSE20685.clinc,Metabric.x,Metabric.y,"Metabric",Metabric.jianm.clinc,"0.8")


Algorithm="step.back"
cvfit<-step.back.ar(TCGA.jianm.x,TCGA.jianm.clinc,Algorithm,"TCGA_TCGA.jianm",Test.x,Test.clinc,"TCGA_Test",GSE20685.x,GSE20685.clinc,"GSE20685",Metabric.x,Metabric.jianm.clinc,"Metabric")
#####################
Algorithm="step.both"
cvfit<-step.both.ar(TCGA.jianm.x,TCGA.jianm.clinc,Algorithm,"TCGA_TCGA.jianm",Test.x,Test.clinc,"TCGA_Test",GSE20685.x,GSE20685.clinc,"GSE20685",Metabric.x,Metabric.jianm.clinc,"Metabric")


##################3
####################################4.LASSO_COX+
LASSO_COX<-read.table("/mnt/nas_002/home/limma/02.project/04.BRCA/01.Analysis/00.jianmo/LASSO_COX_selected_gene.xls",header=TRUE,sep="\t")
colnames(LASSO_COX)=c("names")
set<-LASSO_COX$names
source("/mnt/nas_002/home/limma/02.project/04.BRCA/learning.r")
setwd("/mnt/nas_002/home/limma/02.project/04.BRCA/01.Analysis/00.jianmo")
TCGA.jianm.clinc$OS.time<-as.numeric(TCGA.jianm.clinc$OS.time)
TCGA.jianm.y<-data.matrix(Surv(time=as.double(TCGA.jianm.clinc$OS.time),event=as.double(TCGA.jianm.clinc$OS.Status)))
TCGA.jianm.x<-TCGA.jianm.gene[,set]
###################
Test.clinc$OS.time<-as.numeric(Test.clinc$OS.time)
Test.y<-data.matrix(Surv(time=as.double(Test.clinc$OS.time),event=as.double(Test.clinc$OS.Status)))
Test.x<-Test.gene[,set]
###########
TCGA.jianm.clinc$OS.time<-as.numeric(TCGA.jianm.clinc$OS.time)
TCGA.jianm.y<-data.matrix(Surv(time=as.double(TCGA.jianm.clinc$OS.time),event=as.double(TCGA.jianm.clinc$OS.Status)))
TCGA.jianm.x<-TCGA.jianm.gene[,set]
#######################
GSE20685.clinc$OS.time<-as.numeric(GSE20685.clinc$OS.time)
GSE20685.clinc$OS.Status<-as.numeric(GSE20685.clinc$OS.Status)
GSE20685.y<-data.matrix(Surv(time=as.double(GSE20685.clinc$OS.time),event=as.double(GSE20685.clinc$OS.Status)))
GSE20685.x<-GSE20685.gene[,set]
#######################
Metabric.jianm.clinc$OS.time<-as.numeric(Metabric.jianm.clinc$OS.time)
Metabric.jianm.clinc$OS.Status<-as.numeric(Metabric.jianm.clinc$OS.Status)
Metabric.y<-data.matrix(Surv(time=as.double(Metabric.jianm.clinc$OS.time),event=as.double(Metabric.jianm.clinc$OS.Status)))
Metabric.x<-Metabric.gene[,set]
#######################################
##############3
seed=1
Algorithm="LASSOCOX_LASSO_COX"
cvfit<-lasso(TCGA.jianm.x,TCGA.jianm.y,seed,Algorithm,"TCGA_TCGA.jianm",TCGA.jianm.clinc,Test.x,Test.y,"TCGA_Test",Test.clinc,GSE20685.x,GSE20685.y,"GSE20685",GSE20685.clinc,Metabric.x,Metabric.y,"Metabric",Metabric.jianm.clinc)
#######################
seed=123
Algorithm="LASSOCOX_randomForestSRC"
cvfit<-randomForestSRC.ar(TCGA.jianm.x,TCGA.jianm.clinc,Algorithm,"TCGA_TCGA.jianm",Test.x,Test.clinc,"TCGA_Test",GSE20685.x,GSE20685.clinc,"GSE20685",Metabric.x,Metabric.jianm.clinc,"Metabric")

#########################
Algorithm="LASSOCOX_Coxboost"
cvfit<-Coxboost(TCGA.jianm.x,TCGA.jianm.clinc,Algorithm,"TCGA_TCGA.jianm",Test.x,Test.clinc,"TCGA_Test",GSE20685.x,GSE20685.clinc,"GSE20685",Metabric.x,Metabric.jianm.clinc,"Metabric")
#######################
Algorithm="LASSOCOX_gbm"
cvfit<-gbm.ar(TCGA.jianm.x,TCGA.jianm.clinc,Algorithm,"TCGA_TCGA.jianm",Test.x,Test.clinc,"TCGA_Test",GSE20685.x,GSE20685.clinc,"GSE20685",Metabric.x,Metabric.jianm.clinc,"Metabric")

#######################
seed=123
Algorithm="LASSOCOX_Ridge"
cvfit<-ridge.ar(TCGA.jianm.x,TCGA.jianm.y,seed,Algorithm,"TCGA_TCGA.jianm",TCGA.jianm.clinc,Test.x,Test.y,"TCGA_Test",Test.clinc,GSE20685.x,GSE20685.y,"GSE20685",GSE20685.clinc,Metabric.x,Metabric.y,"Metabric",Metabric.jianm.clinc)
#######################
seed=1
Algorithm="LASSOCOX_Enet0.2"

cvfit<-Enet.ar(TCGA.jianm.x,TCGA.jianm.y,seed,Algorithm,"TCGA_TCGA.jianm",TCGA.jianm.clinc,Test.x,Test.y,"TCGA_Test",Test.clinc,GSE20685.x,GSE20685.y,"GSE20685",GSE20685.clinc,Metabric.x,Metabric.y,"Metabric",Metabric.jianm.clinc,"0.2")
#######################
seed=1
Algorithm="LASSOCOX_Enet0.5"
cvfit<-Enet.ar(TCGA.jianm.x,TCGA.jianm.y,seed,Algorithm,"TCGA_TCGA.jianm",TCGA.jianm.clinc,Test.x,Test.y,"TCGA_Test",Test.clinc,GSE20685.x,GSE20685.y,"GSE20685",GSE20685.clinc,Metabric.x,Metabric.y,"Metabric",Metabric.jianm.clinc,"0.5")
#######################
seed=1
Algorithm="LASSOCOX_Enet0.8"
cvfit<-Enet.ar(TCGA.jianm.x,TCGA.jianm.y,seed,Algorithm,"TCGA_TCGA.jianm",TCGA.jianm.clinc,Test.x,Test.y,"TCGA_Test",Test.clinc,GSE20685.x,GSE20685.y,"GSE20685",GSE20685.clinc,Metabric.x,Metabric.y,"Metabric",Metabric.jianm.clinc,"0.8")

##########################

Algorithm="LASSOCOX_step.back"
cvfit<-step.back.ar(TCGA.jianm.x,TCGA.jianm.clinc,Algorithm,"TCGA_TCGA.jianm",Test.x,Test.clinc,"TCGA_Test",GSE20685.x,GSE20685.clinc,"GSE20685",Metabric.x,Metabric.jianm.clinc,"Metabric")
#####################
Algorithm="LASSOCOX_step.both"
cvfit<-step.both.ar(TCGA.jianm.x,TCGA.jianm.clinc,Algorithm,"TCGA_TCGA.jianm",Test.x,Test.clinc,"TCGA_Test",GSE20685.x,GSE20685.clinc,"GSE20685",Metabric.x,Metabric.jianm.clinc,"Metabric")
#########################

####################################4.randomForestSRC+
RSF<-read.table("/mnt/nas_002/home/limma/02.project/04.BRCA/01.Analysis/00.jianmo/randomForestSRC_selected_gene.xls",header=TRUE,sep="\t")
colnames(RSF)=c("names")
set<-RSF$names
source("/mnt/nas_002/home/limma/02.project/04.BRCA/learning.r")
setwd("/mnt/nas_002/home/limma/02.project/04.BRCA/01.Analysis/00.jianmo")
TCGA.jianm.clinc$OS.time<-as.numeric(TCGA.jianm.clinc$OS.time)
TCGA.jianm.y<-data.matrix(Surv(time=as.double(TCGA.jianm.clinc$OS.time),event=as.double(TCGA.jianm.clinc$OS.Status)))
TCGA.jianm.x<-TCGA.jianm.gene[,set]
###################
Test.clinc$OS.time<-as.numeric(Test.clinc$OS.time)
Test.y<-data.matrix(Surv(time=as.double(Test.clinc$OS.time),event=as.double(Test.clinc$OS.Status)))
Test.x<-Test.gene[,set]
#######################
GSE20685.clinc$OS.time<-as.numeric(GSE20685.clinc$OS.time)
GSE20685.clinc$OS.Status<-as.numeric(GSE20685.clinc$OS.Status)
GSE20685.y<-data.matrix(Surv(time=as.double(GSE20685.clinc$OS.time),event=as.double(GSE20685.clinc$OS.Status)))
GSE20685.x<-GSE20685.gene[,set]
#######################
Metabric.jianm.clinc$OS.time<-as.numeric(Metabric.jianm.clinc$OS.time)
Metabric.jianm.clinc$OS.Status<-as.numeric(Metabric.jianm.clinc$OS.Status)
Metabric.y<-data.matrix(Surv(time=as.double(Metabric.jianm.clinc$OS.time),event=as.double(Metabric.jianm.clinc$OS.Status)))
Metabric.x<-Metabric.gene[,set]
#######################################
##############3
seed=1
Algorithm="RSF_LASSO_COX"
cvfit<-lasso(TCGA.jianm.x,TCGA.jianm.y,seed,Algorithm,"TCGA_TCGA.jianm",TCGA.jianm.clinc,Test.x,Test.y,"TCGA_Test",Test.clinc,GSE20685.x,GSE20685.y,"GSE20685",GSE20685.clinc,Metabric.x,Metabric.y,"Metabric",Metabric.jianm.clinc)
#######################
seed=123
Algorithm="RSF_randomForestSRC"
cvfit<-randomForestSRC.ar(TCGA.jianm.x,TCGA.jianm.clinc,Algorithm,"TCGA_TCGA.jianm",Test.x,Test.clinc,"TCGA_Test",GSE20685.x,GSE20685.clinc,"GSE20685",Metabric.x,Metabric.jianm.clinc,"Metabric")

#########################
Algorithm="RSF_Coxboost"
cvfit<-Coxboost(TCGA.jianm.x,TCGA.jianm.clinc,Algorithm,"TCGA_TCGA.jianm",Test.x,Test.clinc,"TCGA_Test",GSE20685.x,GSE20685.clinc,"GSE20685",Metabric.x,Metabric.jianm.clinc,"Metabric")
#######################
Algorithm="RSF_gbm"
cvfit<-gbm.ar(TCGA.jianm.x,TCGA.jianm.clinc,Algorithm,"TCGA_TCGA.jianm",Test.x,Test.clinc,"TCGA_Test",GSE20685.x,GSE20685.clinc,"GSE20685",Metabric.x,Metabric.jianm.clinc,"Metabric")

#######################
seed=123
Algorithm="RSF_Ridge"
cvfit<-ridge.ar(TCGA.jianm.x,TCGA.jianm.y,seed,Algorithm,"TCGA_TCGA.jianm",TCGA.jianm.clinc,Test.x,Test.y,"TCGA_Test",Test.clinc,GSE20685.x,GSE20685.y,"GSE20685",GSE20685.clinc,Metabric.x,Metabric.y,"Metabric",Metabric.jianm.clinc)
#######################
seed=1
Algorithm="RSF_Enet0.2"

cvfit<-Enet.ar(TCGA.jianm.x,TCGA.jianm.y,seed,Algorithm,"TCGA_TCGA.jianm",TCGA.jianm.clinc,Test.x,Test.y,"TCGA_Test",Test.clinc,GSE20685.x,GSE20685.y,"GSE20685",GSE20685.clinc,Metabric.x,Metabric.y,"Metabric",Metabric.jianm.clinc,"0.2")
#######################
seed=1
Algorithm="RSF_Enet0.5"
cvfit<-Enet.ar(TCGA.jianm.x,TCGA.jianm.y,seed,Algorithm,"TCGA_TCGA.jianm",TCGA.jianm.clinc,Test.x,Test.y,"TCGA_Test",Test.clinc,GSE20685.x,GSE20685.y,"GSE20685",GSE20685.clinc,Metabric.x,Metabric.y,"Metabric",Metabric.jianm.clinc,"0.5")
#######################
seed=1
Algorithm="RSF_Enet0.8"
cvfit<-Enet.ar(TCGA.jianm.x,TCGA.jianm.y,seed,Algorithm,"TCGA_TCGA.jianm",TCGA.jianm.clinc,Test.x,Test.y,"TCGA_Test",Test.clinc,GSE20685.x,GSE20685.y,"GSE20685",GSE20685.clinc,Metabric.x,Metabric.y,"Metabric",Metabric.jianm.clinc,"0.8")

##########################

Algorithm="RSF_step.back"
cvfit<-step.back.ar(TCGA.jianm.x,TCGA.jianm.clinc,Algorithm,"TCGA_TCGA.jianm",Test.x,Test.clinc,"TCGA_Test",GSE20685.x,GSE20685.clinc,"GSE20685",Metabric.x,Metabric.jianm.clinc,"Metabric")
#####################
Algorithm="RSF_step.both"
cvfit<-step.both.ar(TCGA.jianm.x,TCGA.jianm.clinc,Algorithm,"TCGA_TCGA.jianm",Test.x,Test.clinc,"TCGA_Test",GSE20685.x,GSE20685.clinc,"GSE20685",Metabric.x,Metabric.jianm.clinc,"Metabric")


####################################4.Coxboost+
Coxboost<-read.table("/mnt/nas_002/home/limma/02.project/04.BRCA/01.Analysis/00.jianmo/Coxboost_selected_gene.xls",header=TRUE,sep="\t")
colnames(Coxboost)=c("names")
set<-Coxboost$names
source("/mnt/nas_002/home/limma/02.project/04.BRCA/learning.r")
setwd("/mnt/nas_002/home/limma/02.project/04.BRCA/01.Analysis/00.jianmo")
TCGA.jianm.clinc$OS.time<-as.numeric(TCGA.jianm.clinc$OS.time)
TCGA.jianm.y<-data.matrix(Surv(time=as.double(TCGA.jianm.clinc$OS.time),event=as.double(TCGA.jianm.clinc$OS.Status)))
TCGA.jianm.x<-TCGA.jianm.gene[,set]
###################
Test.clinc$OS.time<-as.numeric(Test.clinc$OS.time)
Test.y<-data.matrix(Surv(time=as.double(Test.clinc$OS.time),event=as.double(Test.clinc$OS.Status)))
Test.x<-Test.gene[,set]
#######################
GSE20685.clinc$OS.time<-as.numeric(GSE20685.clinc$OS.time)
GSE20685.clinc$OS.Status<-as.numeric(GSE20685.clinc$OS.Status)
GSE20685.y<-data.matrix(Surv(time=as.double(GSE20685.clinc$OS.time),event=as.double(GSE20685.clinc$OS.Status)))
GSE20685.x<-GSE20685.gene[,set]
#######################
Metabric.jianm.clinc$OS.time<-as.numeric(Metabric.jianm.clinc$OS.time)
Metabric.jianm.clinc$OS.Status<-as.numeric(Metabric.jianm.clinc$OS.Status)
Metabric.y<-data.matrix(Surv(time=as.double(Metabric.jianm.clinc$OS.time),event=as.double(Metabric.jianm.clinc$OS.Status)))
Metabric.x<-Metabric.gene[,set]
#######################################
##############3
seed=1
Algorithm="Coxboost_LASSO_COX"
cvfit<-lasso(TCGA.jianm.x,TCGA.jianm.y,seed,Algorithm,"TCGA_TCGA.jianm",TCGA.jianm.clinc,Test.x,Test.y,"TCGA_Test",Test.clinc,GSE20685.x,GSE20685.y,"GSE20685",GSE20685.clinc,Metabric.x,Metabric.y,"Metabric",Metabric.jianm.clinc)
#######################
seed=123
Algorithm="Coxboost_randomForestSRC"
cvfit<-randomForestSRC.ar(TCGA.jianm.x,TCGA.jianm.clinc,Algorithm,"TCGA_TCGA.jianm",Test.x,Test.clinc,"TCGA_Test",GSE20685.x,GSE20685.clinc,"GSE20685",Metabric.x,Metabric.jianm.clinc,"Metabric")

#########################
Algorithm="Coxboost_Coxboost"
cvfit<-Coxboost(TCGA.jianm.x,TCGA.jianm.clinc,Algorithm,"TCGA_TCGA.jianm",Test.x,Test.clinc,"TCGA_Test",GSE20685.x,GSE20685.clinc,"GSE20685",Metabric.x,Metabric.jianm.clinc,"Metabric")
#######################
Algorithm="Coxboost_gbm"
cvfit<-gbm.ar(TCGA.jianm.x,TCGA.jianm.clinc,Algorithm,"TCGA_TCGA.jianm",Test.x,Test.clinc,"TCGA_Test",GSE20685.x,GSE20685.clinc,"GSE20685",Metabric.x,Metabric.jianm.clinc,"Metabric")

#######################
seed=123
Algorithm="Coxboost_Ridge"
cvfit<-ridge.ar(TCGA.jianm.x,TCGA.jianm.y,seed,Algorithm,"TCGA_TCGA.jianm",TCGA.jianm.clinc,Test.x,Test.y,"TCGA_Test",Test.clinc,GSE20685.x,GSE20685.y,"GSE20685",GSE20685.clinc,Metabric.x,Metabric.y,"Metabric",Metabric.jianm.clinc)
#######################
seed=1
Algorithm="Coxboost_Enet0.2"

cvfit<-Enet.ar(TCGA.jianm.x,TCGA.jianm.y,seed,Algorithm,"TCGA_TCGA.jianm",TCGA.jianm.clinc,Test.x,Test.y,"TCGA_Test",Test.clinc,GSE20685.x,GSE20685.y,"GSE20685",GSE20685.clinc,Metabric.x,Metabric.y,"Metabric",Metabric.jianm.clinc,"0.2")
#######################
seed=1
Algorithm="Coxboost_Enet0.5"
cvfit<-Enet.ar(TCGA.jianm.x,TCGA.jianm.y,seed,Algorithm,"TCGA_TCGA.jianm",TCGA.jianm.clinc,Test.x,Test.y,"TCGA_Test",Test.clinc,GSE20685.x,GSE20685.y,"GSE20685",GSE20685.clinc,Metabric.x,Metabric.y,"Metabric",Metabric.jianm.clinc,"0.5")
#######################
seed=1
Algorithm="Coxboost_Enet0.8"
cvfit<-Enet.ar(TCGA.jianm.x,TCGA.jianm.y,seed,Algorithm,"TCGA_TCGA.jianm",TCGA.jianm.clinc,Test.x,Test.y,"TCGA_Test",Test.clinc,GSE20685.x,GSE20685.y,"GSE20685",GSE20685.clinc,Metabric.x,Metabric.y,"Metabric",Metabric.jianm.clinc,"0.8")

##########################

Algorithm="Coxboost_step.back"
cvfit<-step.back.ar(TCGA.jianm.x,TCGA.jianm.clinc,Algorithm,"TCGA_TCGA.jianm",Test.x,Test.clinc,"TCGA_Test",GSE20685.x,GSE20685.clinc,"GSE20685",Metabric.x,Metabric.jianm.clinc,"Metabric")
#####################
Algorithm="Coxboost_step.both"
cvfit<-step.both.ar(TCGA.jianm.x,TCGA.jianm.clinc,Algorithm,"TCGA_TCGA.jianm",Test.x,Test.clinc,"TCGA_Test",GSE20685.x,GSE20685.clinc,"GSE20685",Metabric.x,Metabric.jianm.clinc,"Metabric")

####################################Enet0.2+
Enet0.2<-read.table("/mnt/nas_002/home/limma/02.project/04.BRCA/01.Analysis/00.jianmo/Enet0.2_selected_gene.xls",header=TRUE,sep="\t")
colnames(Enet0.2)=c("names")
set<-Enet0.2$names
source("/mnt/nas_002/home/limma/02.project/04.BRCA/learning.r")
setwd("/mnt/nas_002/home/limma/02.project/04.BRCA/01.Analysis/00.jianmo")
TCGA.jianm.clinc$OS.time<-as.numeric(TCGA.jianm.clinc$OS.time)
TCGA.jianm.y<-data.matrix(Surv(time=as.double(TCGA.jianm.clinc$OS.time),event=as.double(TCGA.jianm.clinc$OS.Status)))
TCGA.jianm.x<-TCGA.jianm.gene[,set]
###################
Test.clinc$OS.time<-as.numeric(Test.clinc$OS.time)
Test.y<-data.matrix(Surv(time=as.double(Test.clinc$OS.time),event=as.double(Test.clinc$OS.Status)))
Test.x<-Test.gene[,set]
#######################
GSE20685.clinc$OS.time<-as.numeric(GSE20685.clinc$OS.time)
GSE20685.clinc$OS.Status<-as.numeric(GSE20685.clinc$OS.Status)
GSE20685.y<-data.matrix(Surv(time=as.double(GSE20685.clinc$OS.time),event=as.double(GSE20685.clinc$OS.Status)))
GSE20685.x<-GSE20685.gene[,set]
#######################
Metabric.jianm.clinc$OS.time<-as.numeric(Metabric.jianm.clinc$OS.time)
Metabric.jianm.clinc$OS.Status<-as.numeric(Metabric.jianm.clinc$OS.Status)
Metabric.y<-data.matrix(Surv(time=as.double(Metabric.jianm.clinc$OS.time),event=as.double(Metabric.jianm.clinc$OS.Status)))
Metabric.x<-Metabric.gene[,set]
#######################################
##############3
seed=1
Algorithm="Enet0.2_LASSO_COX"
cvfit<-lasso(TCGA.jianm.x,TCGA.jianm.y,seed,Algorithm,"TCGA_TCGA.jianm",TCGA.jianm.clinc,Test.x,Test.y,"TCGA_Test",Test.clinc,GSE20685.x,GSE20685.y,"GSE20685",GSE20685.clinc,Metabric.x,Metabric.y,"Metabric",Metabric.jianm.clinc)
#######################
seed=123
Algorithm="Enet0.2_randomForestSRC"
cvfit<-randomForestSRC.ar(TCGA.jianm.x,TCGA.jianm.clinc,Algorithm,"TCGA_TCGA.jianm",Test.x,Test.clinc,"TCGA_Test",GSE20685.x,GSE20685.clinc,"GSE20685",Metabric.x,Metabric.jianm.clinc,"Metabric")

#########################
Algorithm="Enet0.2_Coxboost"
cvfit<-Coxboost(TCGA.jianm.x,TCGA.jianm.clinc,Algorithm,"TCGA_TCGA.jianm",Test.x,Test.clinc,"TCGA_Test",GSE20685.x,GSE20685.clinc,"GSE20685",Metabric.x,Metabric.jianm.clinc,"Metabric")
#######################
Algorithm="Enet0.2_gbm"
cvfit<-gbm.ar(TCGA.jianm.x,TCGA.jianm.clinc,Algorithm,"TCGA_TCGA.jianm",Test.x,Test.clinc,"TCGA_Test",GSE20685.x,GSE20685.clinc,"GSE20685",Metabric.x,Metabric.jianm.clinc,"Metabric")

#######################
seed=123
Algorithm="Enet0.2_Ridge"
cvfit<-ridge.ar(TCGA.jianm.x,TCGA.jianm.y,seed,Algorithm,"TCGA_TCGA.jianm",TCGA.jianm.clinc,Test.x,Test.y,"TCGA_Test",Test.clinc,GSE20685.x,GSE20685.y,"GSE20685",GSE20685.clinc,Metabric.x,Metabric.y,"Metabric",Metabric.jianm.clinc)
#######################
seed=1
Algorithm="Enet0.2_Enet0.2"

cvfit<-Enet.ar(TCGA.jianm.x,TCGA.jianm.y,seed,Algorithm,"TCGA_TCGA.jianm",TCGA.jianm.clinc,Test.x,Test.y,"TCGA_Test",Test.clinc,GSE20685.x,GSE20685.y,"GSE20685",GSE20685.clinc,Metabric.x,Metabric.y,"Metabric",Metabric.jianm.clinc,"0.2")
#######################
seed=1
Algorithm="Enet0.2_Enet0.5"
cvfit<-Enet.ar(TCGA.jianm.x,TCGA.jianm.y,seed,Algorithm,"TCGA_TCGA.jianm",TCGA.jianm.clinc,Test.x,Test.y,"TCGA_Test",Test.clinc,GSE20685.x,GSE20685.y,"GSE20685",GSE20685.clinc,Metabric.x,Metabric.y,"Metabric",Metabric.jianm.clinc,"0.5")
#######################
seed=1
Algorithm="Enet0.2_Enet0.8"
cvfit<-Enet.ar(TCGA.jianm.x,TCGA.jianm.y,seed,Algorithm,"TCGA_TCGA.jianm",TCGA.jianm.clinc,Test.x,Test.y,"TCGA_Test",Test.clinc,GSE20685.x,GSE20685.y,"GSE20685",GSE20685.clinc,Metabric.x,Metabric.y,"Metabric",Metabric.jianm.clinc,"0.8")

##########################

Algorithm="Enet0.2_step.back"
cvfit<-step.back.ar(TCGA.jianm.x,TCGA.jianm.clinc,Algorithm,"TCGA_TCGA.jianm",Test.x,Test.clinc,"TCGA_Test",GSE20685.x,GSE20685.clinc,"GSE20685",Metabric.x,Metabric.jianm.clinc,"Metabric")
#####################
Algorithm="Enet0.2_step.both"
cvfit<-step.both.ar(TCGA.jianm.x,TCGA.jianm.clinc,Algorithm,"TCGA_TCGA.jianm",Test.x,Test.clinc,"TCGA_Test",GSE20685.x,GSE20685.clinc,"GSE20685",Metabric.x,Metabric.jianm.clinc,"Metabric")

####################################Enet0.5+
Enet0.5<-read.table("/mnt/nas_002/home/limma/02.project/04.BRCA/01.Analysis/00.jianmo/Enet0.5_selected_gene.xls",header=TRUE,sep="\t")
colnames(Enet0.5)=c("names")
set<-Enet0.5$names
source("/mnt/nas_002/home/limma/02.project/04.BRCA/learning.r")
setwd("/mnt/nas_002/home/limma/02.project/04.BRCA/01.Analysis/00.jianmo")
TCGA.jianm.clinc$OS.time<-as.numeric(TCGA.jianm.clinc$OS.time)
TCGA.jianm.y<-data.matrix(Surv(time=as.double(TCGA.jianm.clinc$OS.time),event=as.double(TCGA.jianm.clinc$OS.Status)))
TCGA.jianm.x<-TCGA.jianm.gene[,set]
###################
Test.clinc$OS.time<-as.numeric(Test.clinc$OS.time)
Test.y<-data.matrix(Surv(time=as.double(Test.clinc$OS.time),event=as.double(Test.clinc$OS.Status)))
Test.x<-Test.gene[,set]
#######################
GSE20685.clinc$OS.time<-as.numeric(GSE20685.clinc$OS.time)
GSE20685.clinc$OS.Status<-as.numeric(GSE20685.clinc$OS.Status)
GSE20685.y<-data.matrix(Surv(time=as.double(GSE20685.clinc$OS.time),event=as.double(GSE20685.clinc$OS.Status)))
GSE20685.x<-GSE20685.gene[,set]
#######################
Metabric.jianm.clinc$OS.time<-as.numeric(Metabric.jianm.clinc$OS.time)
Metabric.jianm.clinc$OS.Status<-as.numeric(Metabric.jianm.clinc$OS.Status)
Metabric.y<-data.matrix(Surv(time=as.double(Metabric.jianm.clinc$OS.time),event=as.double(Metabric.jianm.clinc$OS.Status)))
Metabric.x<-Metabric.gene[,set]
#######################################
##############3
seed=1
Algorithm="Enet0.5_LASSO_COX"
cvfit<-lasso(TCGA.jianm.x,TCGA.jianm.y,seed,Algorithm,"TCGA_TCGA.jianm",TCGA.jianm.clinc,Test.x,Test.y,"TCGA_Test",Test.clinc,GSE20685.x,GSE20685.y,"GSE20685",GSE20685.clinc,Metabric.x,Metabric.y,"Metabric",Metabric.jianm.clinc)
#######################
seed=123
Algorithm="Enet0.5_randomForestSRC"
cvfit<-randomForestSRC.ar(TCGA.jianm.x,TCGA.jianm.clinc,Algorithm,"TCGA_TCGA.jianm",Test.x,Test.clinc,"TCGA_Test",GSE20685.x,GSE20685.clinc,"GSE20685",Metabric.x,Metabric.jianm.clinc,"Metabric")

#########################
Algorithm="Enet0.5_Coxboost"
cvfit<-Coxboost(TCGA.jianm.x,TCGA.jianm.clinc,Algorithm,"TCGA_TCGA.jianm",Test.x,Test.clinc,"TCGA_Test",GSE20685.x,GSE20685.clinc,"GSE20685",Metabric.x,Metabric.jianm.clinc,"Metabric")
#######################
Algorithm="Enet0.5_gbm"
cvfit<-gbm.ar(TCGA.jianm.x,TCGA.jianm.clinc,Algorithm,"TCGA_TCGA.jianm",Test.x,Test.clinc,"TCGA_Test",GSE20685.x,GSE20685.clinc,"GSE20685",Metabric.x,Metabric.jianm.clinc,"Metabric")

#######################
seed=123
Algorithm="Enet0.5_Ridge"
cvfit<-ridge.ar(TCGA.jianm.x,TCGA.jianm.y,seed,Algorithm,"TCGA_TCGA.jianm",TCGA.jianm.clinc,Test.x,Test.y,"TCGA_Test",Test.clinc,GSE20685.x,GSE20685.y,"GSE20685",GSE20685.clinc,Metabric.x,Metabric.y,"Metabric",Metabric.jianm.clinc)
#######################
seed=1
Algorithm="Enet0.5_Enet0.2"

cvfit<-Enet.ar(TCGA.jianm.x,TCGA.jianm.y,seed,Algorithm,"TCGA_TCGA.jianm",TCGA.jianm.clinc,Test.x,Test.y,"TCGA_Test",Test.clinc,GSE20685.x,GSE20685.y,"GSE20685",GSE20685.clinc,Metabric.x,Metabric.y,"Metabric",Metabric.jianm.clinc,"0.2")
#######################
seed=1
Algorithm="Enet0.5_Enet0.5"
cvfit<-Enet.ar(TCGA.jianm.x,TCGA.jianm.y,seed,Algorithm,"TCGA_TCGA.jianm",TCGA.jianm.clinc,Test.x,Test.y,"TCGA_Test",Test.clinc,GSE20685.x,GSE20685.y,"GSE20685",GSE20685.clinc,Metabric.x,Metabric.y,"Metabric",Metabric.jianm.clinc,"0.5")
#######################
seed=1
Algorithm="Enet0.5_Enet0.8"
cvfit<-Enet.ar(TCGA.jianm.x,TCGA.jianm.y,seed,Algorithm,"TCGA_TCGA.jianm",TCGA.jianm.clinc,Test.x,Test.y,"TCGA_Test",Test.clinc,GSE20685.x,GSE20685.y,"GSE20685",GSE20685.clinc,Metabric.x,Metabric.y,"Metabric",Metabric.jianm.clinc,"0.8")

##########################

Algorithm="Enet0.5_step.back"
cvfit<-step.back.ar(TCGA.jianm.x,TCGA.jianm.clinc,Algorithm,"TCGA_TCGA.jianm",Test.x,Test.clinc,"TCGA_Test",GSE20685.x,GSE20685.clinc,"GSE20685",Metabric.x,Metabric.jianm.clinc,"Metabric")
#####################
Algorithm="Enet0.5_step.both"
cvfit<-step.both.ar(TCGA.jianm.x,TCGA.jianm.clinc,Algorithm,"TCGA_TCGA.jianm",Test.x,Test.clinc,"TCGA_Test",GSE20685.x,GSE20685.clinc,"GSE20685",Metabric.x,Metabric.jianm.clinc,"Metabric")

####################################Enet0.8+
Enet0.8<-read.table("/mnt/nas_002/home/limma/02.project/04.BRCA/01.Analysis/00.jianmo/Enet0.8_selected_gene.xls",header=TRUE,sep="\t")
colnames(Enet0.8)=c("names")
set<-Enet0.8$names
source("/mnt/nas_002/home/limma/02.project/04.BRCA/learning.r")
setwd("/mnt/nas_002/home/limma/02.project/04.BRCA/01.Analysis/00.jianmo")
TCGA.jianm.clinc$OS.time<-as.numeric(TCGA.jianm.clinc$OS.time)
TCGA.jianm.y<-data.matrix(Surv(time=as.double(TCGA.jianm.clinc$OS.time),event=as.double(TCGA.jianm.clinc$OS.Status)))
TCGA.jianm.x<-TCGA.jianm.gene[,set]
###################
Test.clinc$OS.time<-as.numeric(Test.clinc$OS.time)
Test.y<-data.matrix(Surv(time=as.double(Test.clinc$OS.time),event=as.double(Test.clinc$OS.Status)))
Test.x<-Test.gene[,set]
#######################
GSE20685.clinc$OS.time<-as.numeric(GSE20685.clinc$OS.time)
GSE20685.clinc$OS.Status<-as.numeric(GSE20685.clinc$OS.Status)
GSE20685.y<-data.matrix(Surv(time=as.double(GSE20685.clinc$OS.time),event=as.double(GSE20685.clinc$OS.Status)))
GSE20685.x<-GSE20685.gene[,set]
#######################
Metabric.jianm.clinc$OS.time<-as.numeric(Metabric.jianm.clinc$OS.time)
Metabric.jianm.clinc$OS.Status<-as.numeric(Metabric.jianm.clinc$OS.Status)
Metabric.y<-data.matrix(Surv(time=as.double(Metabric.jianm.clinc$OS.time),event=as.double(Metabric.jianm.clinc$OS.Status)))
Metabric.x<-Metabric.gene[,set]
#######################################
##############3
seed=1
Algorithm="Enet0.8_LASSO_COX"
cvfit<-lasso(TCGA.jianm.x,TCGA.jianm.y,seed,Algorithm,"TCGA_TCGA.jianm",TCGA.jianm.clinc,Test.x,Test.y,"TCGA_Test",Test.clinc,GSE20685.x,GSE20685.y,"GSE20685",GSE20685.clinc,Metabric.x,Metabric.y,"Metabric",Metabric.jianm.clinc)
#######################
seed=123
Algorithm="Enet0.8_randomForestSRC"
cvfit<-randomForestSRC.ar(TCGA.jianm.x,TCGA.jianm.clinc,Algorithm,"TCGA_TCGA.jianm",Test.x,Test.clinc,"TCGA_Test",GSE20685.x,GSE20685.clinc,"GSE20685",Metabric.x,Metabric.jianm.clinc,"Metabric")

#########################
Algorithm="Enet0.8_Coxboost"
cvfit<-Coxboost(TCGA.jianm.x,TCGA.jianm.clinc,Algorithm,"TCGA_TCGA.jianm",Test.x,Test.clinc,"TCGA_Test",GSE20685.x,GSE20685.clinc,"GSE20685",Metabric.x,Metabric.jianm.clinc,"Metabric")
#######################
Algorithm="Enet0.8_gbm"
cvfit<-gbm.ar(TCGA.jianm.x,TCGA.jianm.clinc,Algorithm,"TCGA_TCGA.jianm",Test.x,Test.clinc,"TCGA_Test",GSE20685.x,GSE20685.clinc,"GSE20685",Metabric.x,Metabric.jianm.clinc,"Metabric")

#######################
seed=123
Algorithm="Enet0.8_Ridge"
cvfit<-ridge.ar(TCGA.jianm.x,TCGA.jianm.y,seed,Algorithm,"TCGA_TCGA.jianm",TCGA.jianm.clinc,Test.x,Test.y,"TCGA_Test",Test.clinc,GSE20685.x,GSE20685.y,"GSE20685",GSE20685.clinc,Metabric.x,Metabric.y,"Metabric",Metabric.jianm.clinc)
#######################
seed=1
Algorithm="Enet0.8_Enet0.2"

cvfit<-Enet.ar(TCGA.jianm.x,TCGA.jianm.y,seed,Algorithm,"TCGA_TCGA.jianm",TCGA.jianm.clinc,Test.x,Test.y,"TCGA_Test",Test.clinc,GSE20685.x,GSE20685.y,"GSE20685",GSE20685.clinc,Metabric.x,Metabric.y,"Metabric",Metabric.jianm.clinc,"0.2")
#######################
seed=1
Algorithm="Enet0.8_Enet0.5"
cvfit<-Enet.ar(TCGA.jianm.x,TCGA.jianm.y,seed,Algorithm,"TCGA_TCGA.jianm",TCGA.jianm.clinc,Test.x,Test.y,"TCGA_Test",Test.clinc,GSE20685.x,GSE20685.y,"GSE20685",GSE20685.clinc,Metabric.x,Metabric.y,"Metabric",Metabric.jianm.clinc,"0.5")
#######################
seed=1
Algorithm="Enet0.8_Enet0.8"
cvfit<-Enet.ar(TCGA.jianm.x,TCGA.jianm.y,seed,Algorithm,"TCGA_TCGA.jianm",TCGA.jianm.clinc,Test.x,Test.y,"TCGA_Test",Test.clinc,GSE20685.x,GSE20685.y,"GSE20685",GSE20685.clinc,Metabric.x,Metabric.y,"Metabric",Metabric.jianm.clinc,"0.8")

##########################

Algorithm="Enet0.8_step.back"
cvfit<-step.back.ar(TCGA.jianm.x,TCGA.jianm.clinc,Algorithm,"TCGA_TCGA.jianm",Test.x,Test.clinc,"TCGA_Test",GSE20685.x,GSE20685.clinc,"GSE20685",Metabric.x,Metabric.jianm.clinc,"Metabric")
#####################
Algorithm="Enet0.8_step.both"
cvfit<-step.both.ar(TCGA.jianm.x,TCGA.jianm.clinc,Algorithm,"TCGA_TCGA.jianm",Test.x,Test.clinc,"TCGA_Test",GSE20685.x,GSE20685.clinc,"GSE20685",Metabric.x,Metabric.jianm.clinc,"Metabric")

####################################Ridge+
Ridge<-read.table("/mnt/nas_002/home/limma/02.project/04.BRCA/01.Analysis/00.jianmo/Ridge_selected_gene.xls",header=TRUE,sep="\t")
colnames(Ridge)=c("names")
set<-Ridge$names
source("/mnt/nas_002/home/limma/02.project/04.BRCA/learning.r")
setwd("/mnt/nas_002/home/limma/02.project/04.BRCA/01.Analysis/00.jianmo")
TCGA.jianm.clinc$OS.time<-as.numeric(TCGA.jianm.clinc$OS.time)
TCGA.jianm.y<-data.matrix(Surv(time=as.double(TCGA.jianm.clinc$OS.time),event=as.double(TCGA.jianm.clinc$OS.Status)))
TCGA.jianm.x<-TCGA.jianm.gene[,set]
###################
Test.clinc$OS.time<-as.numeric(Test.clinc$OS.time)
Test.y<-data.matrix(Surv(time=as.double(Test.clinc$OS.time),event=as.double(Test.clinc$OS.Status)))
Test.x<-Test.gene[,set]
#######################
GSE20685.clinc$OS.time<-as.numeric(GSE20685.clinc$OS.time)
GSE20685.clinc$OS.Status<-as.numeric(GSE20685.clinc$OS.Status)
GSE20685.y<-data.matrix(Surv(time=as.double(GSE20685.clinc$OS.time),event=as.double(GSE20685.clinc$OS.Status)))
GSE20685.x<-GSE20685.gene[,set]
#######################
Metabric.jianm.clinc$OS.time<-as.numeric(Metabric.jianm.clinc$OS.time)
Metabric.jianm.clinc$OS.Status<-as.numeric(Metabric.jianm.clinc$OS.Status)
Metabric.y<-data.matrix(Surv(time=as.double(Metabric.jianm.clinc$OS.time),event=as.double(Metabric.jianm.clinc$OS.Status)))
Metabric.x<-Metabric.gene[,set]
#######################################
##############3
seed=1
Algorithm="Ridge_LASSO_COX"
cvfit<-lasso(TCGA.jianm.x,TCGA.jianm.y,seed,Algorithm,"TCGA_TCGA.jianm",TCGA.jianm.clinc,Test.x,Test.y,"TCGA_Test",Test.clinc,GSE20685.x,GSE20685.y,"GSE20685",GSE20685.clinc,Metabric.x,Metabric.y,"Metabric",Metabric.jianm.clinc)
#######################
seed=123
Algorithm="Ridge_randomForestSRC"
cvfit<-randomForestSRC.ar(TCGA.jianm.x,TCGA.jianm.clinc,Algorithm,"TCGA_TCGA.jianm",Test.x,Test.clinc,"TCGA_Test",GSE20685.x,GSE20685.clinc,"GSE20685",Metabric.x,Metabric.jianm.clinc,"Metabric")

#########################
Algorithm="Ridge_Coxboost"
cvfit<-Coxboost(TCGA.jianm.x,TCGA.jianm.clinc,Algorithm,"TCGA_TCGA.jianm",Test.x,Test.clinc,"TCGA_Test",GSE20685.x,GSE20685.clinc,"GSE20685",Metabric.x,Metabric.jianm.clinc,"Metabric")
#######################
Algorithm="Ridge_gbm"
cvfit<-gbm.ar(TCGA.jianm.x,TCGA.jianm.clinc,Algorithm,"TCGA_TCGA.jianm",Test.x,Test.clinc,"TCGA_Test",GSE20685.x,GSE20685.clinc,"GSE20685",Metabric.x,Metabric.jianm.clinc,"Metabric")

#######################
seed=123
Algorithm="Ridge_Ridge"
cvfit<-ridge.ar(TCGA.jianm.x,TCGA.jianm.y,seed,Algorithm,"TCGA_TCGA.jianm",TCGA.jianm.clinc,Test.x,Test.y,"TCGA_Test",Test.clinc,GSE20685.x,GSE20685.y,"GSE20685",GSE20685.clinc,Metabric.x,Metabric.y,"Metabric",Metabric.jianm.clinc)
#######################
seed=1
Algorithm="Ridge_Enet0.2"

cvfit<-Enet.ar(TCGA.jianm.x,TCGA.jianm.y,seed,Algorithm,"TCGA_TCGA.jianm",TCGA.jianm.clinc,Test.x,Test.y,"TCGA_Test",Test.clinc,GSE20685.x,GSE20685.y,"GSE20685",GSE20685.clinc,Metabric.x,Metabric.y,"Metabric",Metabric.jianm.clinc,"0.2")
#######################
seed=1
Algorithm="Ridge_Enet0.5"
cvfit<-Enet.ar(TCGA.jianm.x,TCGA.jianm.y,seed,Algorithm,"TCGA_TCGA.jianm",TCGA.jianm.clinc,Test.x,Test.y,"TCGA_Test",Test.clinc,GSE20685.x,GSE20685.y,"GSE20685",GSE20685.clinc,Metabric.x,Metabric.y,"Metabric",Metabric.jianm.clinc,"0.5")
#######################
seed=1
Algorithm="Ridge_Enet0.8"
cvfit<-Enet.ar(TCGA.jianm.x,TCGA.jianm.y,seed,Algorithm,"TCGA_TCGA.jianm",TCGA.jianm.clinc,Test.x,Test.y,"TCGA_Test",Test.clinc,GSE20685.x,GSE20685.y,"GSE20685",GSE20685.clinc,Metabric.x,Metabric.y,"Metabric",Metabric.jianm.clinc,"0.8")

##########################

Algorithm="Ridge_step.back"
cvfit<-step.back.ar(TCGA.jianm.x,TCGA.jianm.clinc,Algorithm,"TCGA_TCGA.jianm",Test.x,Test.clinc,"TCGA_Test",GSE20685.x,GSE20685.clinc,"GSE20685",Metabric.x,Metabric.jianm.clinc,"Metabric")
#####################
Algorithm="Ridge_step.both"
cvfit<-step.both.ar(TCGA.jianm.x,TCGA.jianm.clinc,Algorithm,"TCGA_TCGA.jianm",Test.x,Test.clinc,"TCGA_Test",GSE20685.x,GSE20685.clinc,"GSE20685",Metabric.x,Metabric.jianm.clinc,"Metabric")
####################################Ridge+
step.back<-read.table("/mnt/nas_002/home/limma/02.project/04.BRCA/01.Analysis/00.jianmo/step.back_selected_gene.xls",header=TRUE,sep="\t")
colnames(step.back)=c("names")
set<-step.back$names
source("/mnt/nas_002/home/limma/02.project/04.BRCA/learning.r")
setwd("/mnt/nas_002/home/limma/02.project/04.BRCA/01.Analysis/00.jianmo")
TCGA.jianm.clinc$OS.time<-as.numeric(TCGA.jianm.clinc$OS.time)
TCGA.jianm.y<-data.matrix(Surv(time=as.double(TCGA.jianm.clinc$OS.time),event=as.double(TCGA.jianm.clinc$OS.Status)))
TCGA.jianm.x<-TCGA.jianm.gene[,set]
###################
Test.clinc$OS.time<-as.numeric(Test.clinc$OS.time)
Test.y<-data.matrix(Surv(time=as.double(Test.clinc$OS.time),event=as.double(Test.clinc$OS.Status)))
Test.x<-Test.gene[,set]
#######################
GSE20685.clinc$OS.time<-as.numeric(GSE20685.clinc$OS.time)
GSE20685.clinc$OS.Status<-as.numeric(GSE20685.clinc$OS.Status)
GSE20685.y<-data.matrix(Surv(time=as.double(GSE20685.clinc$OS.time),event=as.double(GSE20685.clinc$OS.Status)))
GSE20685.x<-GSE20685.gene[,set]
#######################
Metabric.jianm.clinc$OS.time<-as.numeric(Metabric.jianm.clinc$OS.time)
Metabric.jianm.clinc$OS.Status<-as.numeric(Metabric.jianm.clinc$OS.Status)
Metabric.y<-data.matrix(Surv(time=as.double(Metabric.jianm.clinc$OS.time),event=as.double(Metabric.jianm.clinc$OS.Status)))
Metabric.x<-Metabric.gene[,set]
#######################################
##############3
seed=1
Algorithm="step.back_LASSO_COX"
cvfit<-lasso(TCGA.jianm.x,TCGA.jianm.y,seed,Algorithm,"TCGA_TCGA.jianm",TCGA.jianm.clinc,Test.x,Test.y,"TCGA_Test",Test.clinc,GSE20685.x,GSE20685.y,"GSE20685",GSE20685.clinc,Metabric.x,Metabric.y,"Metabric",Metabric.jianm.clinc)
#######################
seed=123
Algorithm="step.back_randomForestSRC"
cvfit<-randomForestSRC.ar(TCGA.jianm.x,TCGA.jianm.clinc,Algorithm,"TCGA_TCGA.jianm",Test.x,Test.clinc,"TCGA_Test",GSE20685.x,GSE20685.clinc,"GSE20685",Metabric.x,Metabric.jianm.clinc,"Metabric")

#########################
Algorithm="step.back_Coxboost"
cvfit<-Coxboost(TCGA.jianm.x,TCGA.jianm.clinc,Algorithm,"TCGA_TCGA.jianm",Test.x,Test.clinc,"TCGA_Test",GSE20685.x,GSE20685.clinc,"GSE20685",Metabric.x,Metabric.jianm.clinc,"Metabric")
#######################
Algorithm="step.back_gbm"
cvfit<-gbm.ar(TCGA.jianm.x,TCGA.jianm.clinc,Algorithm,"TCGA_TCGA.jianm",Test.x,Test.clinc,"TCGA_Test",GSE20685.x,GSE20685.clinc,"GSE20685",Metabric.x,Metabric.jianm.clinc,"Metabric")

#######################
seed=123
Algorithm="step.back_Ridge"
cvfit<-ridge.ar(TCGA.jianm.x,TCGA.jianm.y,seed,Algorithm,"TCGA_TCGA.jianm",TCGA.jianm.clinc,Test.x,Test.y,"TCGA_Test",Test.clinc,GSE20685.x,GSE20685.y,"GSE20685",GSE20685.clinc,Metabric.x,Metabric.y,"Metabric",Metabric.jianm.clinc)
#######################
seed=1
Algorithm="step.back_Enet0.2"

cvfit<-Enet.ar(TCGA.jianm.x,TCGA.jianm.y,seed,Algorithm,"TCGA_TCGA.jianm",TCGA.jianm.clinc,Test.x,Test.y,"TCGA_Test",Test.clinc,GSE20685.x,GSE20685.y,"GSE20685",GSE20685.clinc,Metabric.x,Metabric.y,"Metabric",Metabric.jianm.clinc,"0.2")
#######################
seed=1
Algorithm="step.back_Enet0.5"
cvfit<-Enet.ar(TCGA.jianm.x,TCGA.jianm.y,seed,Algorithm,"TCGA_TCGA.jianm",TCGA.jianm.clinc,Test.x,Test.y,"TCGA_Test",Test.clinc,GSE20685.x,GSE20685.y,"GSE20685",GSE20685.clinc,Metabric.x,Metabric.y,"Metabric",Metabric.jianm.clinc,"0.5")
#######################
seed=1
Algorithm="step.back_Enet0.8"
cvfit<-Enet.ar(TCGA.jianm.x,TCGA.jianm.y,seed,Algorithm,"TCGA_TCGA.jianm",TCGA.jianm.clinc,Test.x,Test.y,"TCGA_Test",Test.clinc,GSE20685.x,GSE20685.y,"GSE20685",GSE20685.clinc,Metabric.x,Metabric.y,"Metabric",Metabric.jianm.clinc,"0.8")

##########################

Algorithm="step.back_step.back"
cvfit<-step.back.ar(TCGA.jianm.x,TCGA.jianm.clinc,Algorithm,"TCGA_TCGA.jianm",Test.x,Test.clinc,"TCGA_Test",GSE20685.x,GSE20685.clinc,"GSE20685",Metabric.x,Metabric.jianm.clinc,"Metabric")
#####################
Algorithm="step.back_step.both"
cvfit<-step.both.ar(TCGA.jianm.x,TCGA.jianm.clinc,Algorithm,"TCGA_TCGA.jianm",Test.x,Test.clinc,"TCGA_Test",GSE20685.x,GSE20685.clinc,"GSE20685",Metabric.x,Metabric.jianm.clinc,"Metabric")









