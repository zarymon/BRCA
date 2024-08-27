library(randomForest)
library(ImageGP)
library(ggplot2)
library(caret)
library("glmnet") 
library("survival")
library(MASS)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)
library(cowplot)
library(GSVA)
#library(msigdbr)
library(GSEABase)
library(limma)
library(pheatmap)
library(patchwork)
library(ggpubr)
library(ggsignif)
library(data.table)
library(survivalROC)
###########
lasso<-function(x,y,seed,Algorithm,name,TCGA.clinc,test.x,test.y,test.name,test.clinc,GSE20685.x,GSE20685.y,GSE20685.name,GSE20685.clinc,Metabric.x,Metabric.y,Metabric.name,Metabric.clinc){
  
  set.seed(seed)
  cvfit <- cv.glmnet(x, y, family = "cox",nfolds = 10)
  pdf(file = paste(Algorithm,".cvglmnet.pdf",sep=""), height = 5)
  plot(cvfit)
  dev.off()
  coef.min = coef(cvfit, s = "lambda.min") 
  ## lambda.min & lambda.1se 取一个
  cvfit$lambda.min
  active.min = which(coef.min@i != 0) ## 找出那些回归系数没有被惩罚为0的
  
  active<-as.numeric(coef.min)[active.min]
  lasso_geneids <-rownames(as.data.frame(which(coef.min[,1]!=0)))
  #lasso_geneids <-rownames(coef.min)[active.min]
  
  #lasso_geneids <- coef.min@Dimnames[[1]][coef.min@i+1] ## 提取基因名称
  lasso.coef<-coef.min@x
  names(lasso.coef)<-lasso_geneids
  lasso.coef<-as.data.frame(lasso.coef)
  saveRDS(cvfit,paste(Algorithm,".pre.rds",sep=""))
  write.table(lasso_geneids, paste(Algorithm,"_selected_gene.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
   write.table(lasso.coef, paste(Algorithm,"_lasso.coef.xls",sep=""),sep = "\t",quote = FALSE,row.names =TRUE)
  ###################
  fit <- glmnet(x, y, family = "cox", alpha = 1) # make sure alpha = 1
  plot(fit, xvar="lambda",label = F,las=1)
  abline(v = log(cvfit$lambda.min), lty = 3,
         lwd = 2,
         col = "black")
  pdf(file = paste(Algorithm,".glmnet.pdf",sep=""), height = 5)
  plot(fit, xvar="lambda",label = F,las=1)
  abline(v = log(cvfit$lambda.min), lty = 3,
         lwd = 2,
         col = "black")
  #plot(fit, xvar="dev",label = F)
  dev.off()
  ###############################################3.在train中进行predict
  ############################在train中进行predict
  #输出新数据的预测值，type参数允许选择预测的类型并提供预测值
  ###############################
  per<-predict(cvfit, newx = x,  s="lambda.min",type="link")
  c_index<-apply(per, 2, Cindex, y=y)
  per<-as.data.frame(per)
  colnames(per)<-c("predict")
  per$ID=rownames(per)
  Auc_data<-merge(TCGA.clinc,per,by="ID")
  write.table(Auc_data, paste(Algorithm,"Auc_data.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  #########################
  library(survivalROC)
  auc<-data.frame()
  for (i in 1:3){
  ROC= survivalROC(Stime=Auc_data$OS.time,##生存时间
                   status=Auc_data$OS.Status,## 终止事件    
                   marker = Auc_data$predict, ## marker value    
                   #other_markers = as.matrix(feng_36[,c("LAG3","IDO1")]), # 矫正的协变量
                   predict.time = i,## 预测时间截点
                   method="KM"
  )##span,NNE法的namda
   if (i==1){auc<-as.data.frame(ROC$AUC)}
   else{auc<-cbind(auc,as.data.frame(ROC$AUC))}
  }
  colnames(auc)=c("1year","2year","3year")
  #write.table(mutate(as.data.frame(ROC$AUC),data=name,Algorithm=Algorithm), paste(Algorithm,".ROC.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  #######################c_index/Auc
  #输出新数据的预测值，type参数允许选择预测的类型并提供预测值
  test.per<-predict(cvfit, newx = test.x,  s="lambda.min",type="link")
  test.c_index<-apply(test.per, 2, Cindex, y=test.y)
  test.per<-as.data.frame(test.per)
  colnames(test.per)<-c("predict")
  test.per$ID=rownames(test.per)
  test.Auc_data<-merge(test.clinc,test.per,by="ID")
  write.table(test.Auc_data, paste(Algorithm,"test.Auc_data.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  #########################
  test.auc<-data.frame()
  for (i in 1:3){
  test.ROC= survivalROC(Stime=test.Auc_data$OS.time,##生存时间
                        status=test.Auc_data$OS.Status,## 终止事件    
                        marker = test.Auc_data$predict, ## marker value    
                        #other_markers = as.matrix(feng_36[,c("LAG3","IDO1")]), # 矫正的协变量
                        predict.time = i,## 预测时间截点
                        method="KM"
  )##span,NNE法的namda
   if (i==1){test.auc<-as.data.frame(test.ROC$AUC)}
   else{test.auc<-cbind(test.auc,as.data.frame(test.ROC$AUC))}
  }
  colnames(test.auc)=c("test1year","test2year","test3year")
  #######################c_index/Auc
  #输出新数据的预测值，type参数允许选择预测的类型并提供预测值
  GSE20685.per<-predict(cvfit, newx = GSE20685.x,  s="lambda.min",type="link")
  GSE20685.c_index<-apply(GSE20685.per, 2, Cindex, y=GSE20685.y)
  GSE20685.per<-as.data.frame(GSE20685.per)
  colnames(GSE20685.per)<-c("predict")
  GSE20685.per$ID=rownames(GSE20685.per)
  GSE20685.Auc_data<-merge(GSE20685.clinc,GSE20685.per,by="ID")
  write.table(GSE20685.Auc_data, paste(Algorithm,"GSE20685.Auc_data.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  #########################
   GSE20685.auc<-data.frame()
  for (i in 1:3){
  GSE20685.ROC= survivalROC(Stime=GSE20685.Auc_data$OS.time,##生存时间
                        status=GSE20685.Auc_data$OS.Status,## 终止事件    
                        marker = GSE20685.Auc_data$predict, ## marker value    
                        #other_markers = as.matrix(feng_36[,c("LAG3","IDO1")]), # 矫正的协变量
                        predict.time = i,## 预测时间截点
                        method="KM"
  )##span,NNE法的namda
   if (i==1){GSE20685.auc<-as.data.frame(GSE20685.ROC$AUC)}
   else{GSE20685.auc<-cbind(GSE20685.auc,as.data.frame(GSE20685.ROC$AUC))}
  }
  colnames(GSE20685.auc)=c("GSE206851year","GSE206852year","GSE206853year")
  #######################c_index/Auc
  #输出新数据的预测值，type参数允许选择预测的类型并提供预测值
  Metabric.per<-predict(cvfit, newx = Metabric.x,  s="lambda.min",type="link")
  Metabric.c_index<-apply(Metabric.per, 2, Cindex, y=Metabric.y)
  Metabric.per<-as.data.frame(Metabric.per)
  colnames(Metabric.per)<-c("predict")
  Metabric.per$ID=rownames(Metabric.per)
  Metabric.Auc_data<-merge(Metabric.clinc,Metabric.per,by="ID")
  write.table(Metabric.Auc_data, paste(Algorithm,"Metabric.Auc_data.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  #########################
   Metabric.auc<-data.frame()
  for (i in 1:3){
  Metabric.ROC= survivalROC(Stime=Metabric.Auc_data$OS.time,##生存时间
                        status=Metabric.Auc_data$OS.Status,## 终止事件    
                        marker = Metabric.Auc_data$predict, ## marker value    
                        #other_markers = as.matrix(feng_36[,c("LAG3","IDO1")]), # 矫正的协变量
                        predict.time = i,## 预测时间截点
                        method="KM"
  )##span,NNE法的namda
   if (i==1){Metabric.auc<-as.data.frame(Metabric.ROC$AUC)}
   else{Metabric.auc<-cbind(Metabric.auc,as.data.frame(Metabric.ROC$AUC))}
  }
  colnames(Metabric.auc)=c("Metabric1year","Metabric2year","Metabric3year")
  ###################
  #Auc.out<-cbind(as.data.frame(ROC$AUC),as.data.frame(test.ROC$AUC),as.data.frame(GSE20685.auc))
   Auc.out<-cbind(auc,test.auc,GSE20685.auc,Metabric.auc)
  Auc.out$Algorithm=Algorithm
  write.table(Auc.out,paste(Algorithm,".Auc.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  c_index<-cbind(as.data.frame(c_index) ,as.data.frame(test.c_index),as.data.frame(GSE20685.c_index),as.data.frame(Metabric.c_index)) 
  colnames(c_index)=c(name,test.name,GSE20685.name,Metabric.name)
  c_index$Algorithm=Algorithm
  write.table(c_index,paste(Algorithm,".c_index.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  #write.table(test.c_index, paste(Algorithm,".c_index.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE,append=TRUE,col.names=FALSE)
  #write.table(c_index, paste(Algorithm,".c_index.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  #write.table(mutate(as.data.frame(test.ROC$AUC),data=test.name,Algorithm=Algorithm), paste(Algorithm,".ROC.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE,append=TRUE,col.names=FALSE)
}
ridge.ar<-function(x,y,seed,Algorithm,name,TCGA.clinc,test.x,test.y,test.name,test.clinc,GSE20685.x,GSE20685.y,GSE20685.name,GSE20685.clinc,Metabric.x,Metabric.y,Metabric.name,Metabric.clinc){
  set.seed(seed)
  cvfit <- cv.glmnet(x, y, family = "cox",nfolds = 10, alpha = 0)
  plot(cvfit)
  pdf(file = paste(Algorithm,".cvglmnet.pdf",sep=""), height = 5)
  plot(cvfit)
  dev.off()
  coef.min = coef(cvfit, s = "lambda.min") 
  ## lambda.min & lambda.1se 取一个
  cvfit$lambda.min
  active.min = which(abs(coef.min )>0.01 ) ## 找出那些回归系数没有被惩罚为0的
  
  #active<-as.numeric(coef.min)[active.min]
  lasso_geneids <-rownames(coef.min)[active.min]
  
  #lasso_geneids <- coef.min@Dimnames[[1]][coef.min@i+1] ## 提取基因名称
  saveRDS(cvfit,paste(Algorithm,".pre.rds",sep=""))
  write.table(lasso_geneids, paste(Algorithm,"_selected_gene.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  ###################
  fit <- glmnet(x, y, family = "cox", alpha = 0) # make sure alpha = 1
  plot(fit, xvar="lambda",label = F,las=1)
  abline(v = log(cvfit$lambda.min), lty = 3,
         lwd = 2,
         col = "black")
  pdf(file = paste(Algorithm,".glmnet.pdf",sep=""), height = 5)
  plot(fit, xvar="lambda",label = F,las=1)
  abline(v = log(cvfit$lambda.min), lty = 3,
         lwd = 2,
         col = "black")
  #plot(fit, xvar="dev",label = F)
  dev.off()
  ###############################################3.在train中进行predict
  ############################在train中进行predict
  #输出新数据的预测值，type参数允许选择预测的类型并提供预测值
  name="TCGA"
  ###############################
  per<-predict(cvfit, newx = x,  s="lambda.min",type="link")
  c_index<-apply(per, 2, Cindex, y=y)
  per<-as.data.frame(per)
  colnames(per)<-c("predict")
  per$ID=rownames(per)
  Auc_data<-merge(TCGA.clinc,per,by="ID")
 write.table(Auc_data, paste(Algorithm,"Auc_data.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  #########################
  library(survivalROC)
  auc<-data.frame()
  for (i in 1:3){
  ROC= survivalROC(Stime=Auc_data$OS.time,##生存时间
                   status=Auc_data$OS.Status,## 终止事件    
                   marker = Auc_data$predict, ## marker value    
                   #other_markers = as.matrix(feng_36[,c("LAG3","IDO1")]), # 矫正的协变量
                   predict.time = i,## 预测时间截点
                   method="KM"
  )##span,NNE法的namda
   if (i==1){auc<-as.data.frame(ROC$AUC)}
   else{auc<-cbind(auc,as.data.frame(ROC$AUC))}
  }
  colnames(auc)=c("1year","2year","3year")
 # write.table(mutate(as.data.frame(ROC$AUC),data=name,Algorithm=Algorithm), paste(Algorithm,".ROC.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  ###############################################4.在test中进行predict
  #######################c_index/Auc
  #输出新数据的预测值，type参数允许选择预测的类型并提供预测值
  test.per<-predict(cvfit, newx = test.x,  s="lambda.min",type="link")
  test.c_index<-apply(test.per, 2, Cindex, y=test.y)
  test.per<-as.data.frame(test.per)
  colnames(test.per)<-c("predict")
  test.per$ID=rownames(test.per)
  test.Auc_data<-merge(test.clinc,test.per,by="ID")
 write.table(test.Auc_data, paste(Algorithm,"test.Auc_data.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  #########################
 test.auc<-data.frame()
  for (i in 1:3){
  test.ROC= survivalROC(Stime=test.Auc_data$OS.time,##生存时间
                        status=test.Auc_data$OS.Status,## 终止事件    
                        marker = test.Auc_data$predict, ## marker value    
                        #other_markers = as.matrix(feng_36[,c("LAG3","IDO1")]), # 矫正的协变量
                        predict.time = i,## 预测时间截点
                        method="KM"
  )##span,NNE法的namda
   if (i==1){test.auc<-as.data.frame(test.ROC$AUC)}
   else{test.auc<-cbind(test.auc,as.data.frame(test.ROC$AUC))}
  }
  colnames(test.auc)=c("test1year","test2year","test3year")
  #######################c_index/Auc
  #输出新数据的预测值，type参数允许选择预测的类型并提供预测值
  GSE20685.per<-predict(cvfit, newx = GSE20685.x,  s="lambda.min",type="link")
  GSE20685.c_index<-apply(GSE20685.per, 2, Cindex, y=GSE20685.y)
  GSE20685.per<-as.data.frame(GSE20685.per)
  colnames(GSE20685.per)<-c("predict")
  GSE20685.per$ID=rownames(GSE20685.per)
  GSE20685.Auc_data<-merge(GSE20685.clinc,GSE20685.per,by="ID")
  write.table(GSE20685.Auc_data, paste(Algorithm,"GSE20685.Auc_data.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  #########################
  GSE20685.auc<-data.frame()
  for (i in 1:3){
  GSE20685.ROC= survivalROC(Stime=GSE20685.Auc_data$OS.time,##生存时间
                        status=GSE20685.Auc_data$OS.Status,## 终止事件    
                        marker = GSE20685.Auc_data$predict, ## marker value    
                        #other_markers = as.matrix(feng_36[,c("LAG3","IDO1")]), # 矫正的协变量
                        predict.time = i,## 预测时间截点
                        method="KM"
  )##span,NNE法的namda
   if (i==1){GSE20685.auc<-as.data.frame(GSE20685.ROC$AUC)}
   else{GSE20685.auc<-cbind(GSE20685.auc,as.data.frame(GSE20685.ROC$AUC))}
  }
  colnames(GSE20685.auc)=c("GSE206851year","GSE206852year","GSE206853year")
  #######################c_index/Auc
  #输出新数据的预测值，type参数允许选择预测的类型并提供预测值
  Metabric.per<-predict(cvfit, newx = Metabric.x,  s="lambda.min",type="link")
  Metabric.c_index<-apply(Metabric.per, 2, Cindex, y=Metabric.y)
  Metabric.per<-as.data.frame(Metabric.per)
  colnames(Metabric.per)<-c("predict")
  Metabric.per$ID=rownames(Metabric.per)
  Metabric.Auc_data<-merge(Metabric.clinc,Metabric.per,by="ID")
  write.table(Metabric.Auc_data, paste(Algorithm,"Metabric.Auc_data.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  #########################
  Metabric.auc<-data.frame()
  for (i in 1:3){
  Metabric.ROC= survivalROC(Stime=Metabric.Auc_data$OS.time,##生存时间
                        status=Metabric.Auc_data$OS.Status,## 终止事件    
                        marker = Metabric.Auc_data$predict, ## marker value    
                        #other_markers = as.matrix(feng_36[,c("LAG3","IDO1")]), # 矫正的协变量
                        predict.time = i,## 预测时间截点
                        method="KM"
  )##span,NNE法的namda
   if (i==1){Metabric.auc<-as.data.frame(Metabric.ROC$AUC)}
   else{Metabric.auc<-cbind(Metabric.auc,as.data.frame(Metabric.ROC$AUC))}
  }
  colnames(Metabric.auc)=c("Metabric1year","Metabric2year","Metabric3year")
  #################
   #Auc.out<-cbind(as.data.frame(ROC$AUC),as.data.frame(test.ROC$AUC),as.data.frame(GSE20685.auc))
  #colnames(Auc.out)=c(name,test.name,GSE20685.name)
  Auc.out<-cbind(auc,test.auc,GSE20685.auc,Metabric.auc)
  Auc.out$Algorithm=Algorithm
  write.table(Auc.out,paste(Algorithm,".Auc.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  c_index<-cbind(as.data.frame(c_index) ,as.data.frame(test.c_index),as.data.frame(GSE20685.c_index),as.data.frame(Metabric.c_index)) 
  colnames(c_index)=c(name,test.name,GSE20685.name,Metabric.name)
  c_index$Algorithm=Algorithm
  write.table(c_index,paste(Algorithm,".c_index.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
}
Enet.ar<-function(x,y,seed,Algorithm,name,TCGA.clinc,test.x,test.y,test.name,test.clinc,GSE20685.x,GSE20685.y,GSE20685.name,GSE20685.clinc,Metabric.x,Metabric.y,Metabric.name,Metabric.clinc,alpha){
  set.seed(seed)
  cvfit <- cv.glmnet(x, y, family = "cox",nfolds = 10, alpha = alpha)
  plot(cvfit)
  pdf(file = paste(Algorithm,".cvglmnet.pdf",sep=""), height = 5)
  plot(cvfit)
  dev.off()
  coef.min = coef(cvfit, s = "lambda.min") 
  ## lambda.min & lambda.1se 取一个
  cvfit$lambda.min
  active.min = which(abs(coef.min )>0.01 ) ## 找出那些回归系数没有被惩罚为0的
  
  #active<-as.numeric(coef.min)[active.min]
  lasso_geneids <-rownames(coef.min)[active.min]
  
  #lasso_geneids <- coef.min@Dimnames[[1]][coef.min@i+1] ## 提取基因名称
  saveRDS(cvfit,paste(Algorithm,".pre.rds",sep=""))
  write.table(lasso_geneids, paste(Algorithm,"_selected_gene.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  ###################
  fit <- glmnet(x, y, family = "cox", alpha = 0) # make sure alpha = 1
  plot(fit, xvar="lambda",label = F,las=1)
  abline(v = log(cvfit$lambda.min), lty = 3,
         lwd = 2,
         col = "black")
  pdf(file = paste(Algorithm,".glmnet.pdf",sep=""), height = 5)
  plot(fit, xvar="lambda",label = F,las=1)
  abline(v = log(cvfit$lambda.min), lty = 3,
         lwd = 2,
         col = "black")
  #plot(fit, xvar="dev",label = F)
  dev.off()
  ###############################################3.在train中进行predict
  ############################在train中进行predict
  #输出新数据的预测值，type参数允许选择预测的类型并提供预测值
  name="TCGA"
  ###############################
  per<-predict(cvfit, newx = x,  s="lambda.min",type="link")
  c_index<-apply(per, 2, Cindex, y=y)
  per<-as.data.frame(per)
  colnames(per)<-c("predict")
  per$ID=rownames(per)
  Auc_data<-merge(TCGA.clinc,per,by="ID")
 # write.table(Auc_data, paste(Algorithm,"Auc_data.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  #########################
  library(survivalROC)
   auc<-data.frame()
  for (i in 1:3){
  ROC= survivalROC(Stime=Auc_data$OS.time,##生存时间
                   status=Auc_data$OS.Status,## 终止事件    
                   marker = Auc_data$predict, ## marker value    
                   #other_markers = as.matrix(feng_36[,c("LAG3","IDO1")]), # 矫正的协变量
                   predict.time = i,## 预测时间截点
                   method="KM"
  )##span,NNE法的namda
   if (i==1){auc<-as.data.frame(ROC$AUC)}
   else{auc<-cbind(auc,as.data.frame(ROC$AUC))}
  }
  colnames(auc)=c("1year","2year","3year")
 # write.table(mutate(as.data.frame(ROC$AUC),data=name,Algorithm=Algorithm), paste(Algorithm,".ROC.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  ###############################################4.在test中进行predict
  #######################c_index/Auc
  #输出新数据的预测值，type参数允许选择预测的类型并提供预测值
  test.per<-predict(cvfit, newx = test.x,  s="lambda.min",type="link")
  test.c_index<-apply(test.per, 2, Cindex, y=test.y)
  test.per<-as.data.frame(test.per)
  colnames(test.per)<-c("predict")
  test.per$ID=rownames(test.per)
  test.Auc_data<-merge(test.clinc,test.per,by="ID")
 # write.table(test.Auc_data, paste(Algorithm,"test.Auc_data.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  #########################
   test.auc<-data.frame()
  for (i in 1:3){
  test.ROC= survivalROC(Stime=test.Auc_data$OS.time,##生存时间
                        status=test.Auc_data$OS.Status,## 终止事件    
                        marker = test.Auc_data$predict, ## marker value    
                        #other_markers = as.matrix(feng_36[,c("LAG3","IDO1")]), # 矫正的协变量
                        predict.time = i,## 预测时间截点
                        method="KM"
  )##span,NNE法的namda
   if (i==1){test.auc<-as.data.frame(test.ROC$AUC)}
   else{test.auc<-cbind(test.auc,as.data.frame(test.ROC$AUC))}
  }
  colnames(test.auc)=c("test1year","test2year","test3year")
  #######################c_index/Auc
  #输出新数据的预测值，type参数允许选择预测的类型并提供预测值
  GSE20685.per<-predict(cvfit, newx = GSE20685.x,  s="lambda.min",type="link")
  GSE20685.c_index<-apply(GSE20685.per, 2, Cindex, y=GSE20685.y)
  GSE20685.per<-as.data.frame(GSE20685.per)
  colnames(GSE20685.per)<-c("predict")
  GSE20685.per$ID=rownames(GSE20685.per)
  GSE20685.Auc_data<-merge(GSE20685.clinc,GSE20685.per,by="ID")
 #write.table(GSE20685.Auc_data, paste(Algorithm,"GSE20685.Auc_data.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  #########################
    GSE20685.auc<-data.frame()
  for (i in 1:3){
  GSE20685.ROC= survivalROC(Stime=GSE20685.Auc_data$OS.time,##生存时间
                        status=GSE20685.Auc_data$OS.Status,## 终止事件    
                        marker = GSE20685.Auc_data$predict, ## marker value    
                        #other_markers = as.matrix(feng_36[,c("LAG3","IDO1")]), # 矫正的协变量
                        predict.time = i,## 预测时间截点
                        method="KM"
  )##span,NNE法的namda
   if (i==1){GSE20685.auc<-as.data.frame(GSE20685.ROC$AUC)}
   else{GSE20685.auc<-cbind(GSE20685.auc,as.data.frame(GSE20685.ROC$AUC))}
  }
  colnames(GSE20685.auc)=c("GSE206851year","GSE206852year","GSE206853year")
  #######################c_index/Auc
  #输出新数据的预测值，type参数允许选择预测的类型并提供预测值
  Metabric.per<-predict(cvfit, newx = Metabric.x,  s="lambda.min",type="link")
  Metabric.c_index<-apply(Metabric.per, 2, Cindex, y=Metabric.y)
  Metabric.per<-as.data.frame(Metabric.per)
  colnames(Metabric.per)<-c("predict")
  Metabric.per$ID=rownames(Metabric.per)
  Metabric.Auc_data<-merge(Metabric.clinc,Metabric.per,by="ID")
 #write.table(Metabric.Auc_data, paste(Algorithm,"Metabric.Auc_data.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  #########################
    Metabric.auc<-data.frame()
  for (i in 1:3){
  Metabric.ROC= survivalROC(Stime=Metabric.Auc_data$OS.time,##生存时间
                        status=Metabric.Auc_data$OS.Status,## 终止事件    
                        marker = Metabric.Auc_data$predict, ## marker value    
                        #other_markers = as.matrix(feng_36[,c("LAG3","IDO1")]), # 矫正的协变量
                        predict.time = i,## 预测时间截点
                        method="KM"
  )##span,NNE法的namda
   if (i==1){Metabric.auc<-as.data.frame(Metabric.ROC$AUC)}
   else{Metabric.auc<-cbind(Metabric.auc,as.data.frame(Metabric.ROC$AUC))}
  }
  colnames(Metabric.auc)=c("Metabric1year","Metabric2year","Metabric3year")
  #################
   #Auc.out<-cbind(as.data.frame(ROC$AUC),as.data.frame(test.ROC$AUC),as.data.frame(GSE20685.auc))
  #colnames(Auc.out)=c(name,test.name,GSE20685.name)
    Auc.out<-cbind(auc,test.auc,GSE20685.auc,Metabric.auc)
  Auc.out$Algorithm=Algorithm
  write.table(Auc.out,paste(Algorithm,".Auc.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  c_index<-cbind(as.data.frame(c_index) ,as.data.frame(test.c_index),as.data.frame(GSE20685.c_index),as.data.frame(Metabric.c_index)) 
  colnames(c_index)=c(name,test.name,GSE20685.name,Metabric.name)
  c_index$Algorithm=Algorithm
  write.table(c_index,paste(Algorithm,".c_index.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
}
####################2.CoxBoost
Coxboost<-function(x,survival_data,Algorithm,name,test.x,test.survival_data,test.name,GSE20685.x,GSE20685.survival_data,GSE20685.name,Metabric.x,Metabric.survival_data,Metabric.name){library(CoxBoost)
  optim.res <- optimCoxBoostPenalty(time=survival_data$OS.time,status=survival_data$OS.Status,x=x,trace=TRUE,start.penalty=50)
  y<-Surv(time=as.double(survival_data$OS.time),event=as.double(survival_data$OS.Status))
  test.y<-Surv(time=as.double(test.survival_data$OS.time),event=as.double(test.survival_data$OS.Status))
  GSE20685.y<-Surv(time=as.double(GSE20685.survival_data$OS.time),event=as.double(GSE20685.survival_data$OS.Status))
  Metabric.y<-Surv(time=as.double(Metabric.survival_data$OS.time),event=as.double(Metabric.survival_data$OS.Status))
  CoxBoost <- CoxBoost(time=survival_data$OS.time,status=survival_data$OS.Status,x=x,
                       stepno=optim.res$cv.res$optimal.step,
                       penalty=optim.res$penalty)
  saveRDS(CoxBoost,"CoxBoost.rds")
  coef<-coef(CoxBoost,at.step=optim.res$cv.res$optimal.step )
  coef<-as.data.frame(coef)
  coef$ID<-rownames(coef)
  
  geneids<- coef[coef$coef!=0,,drop=F] %>% dplyr::select(ID,everything())
  #geneids
  write.table(geneids, paste(Algorithm,"_selected_gene.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  #########################3.构建模型
  per<-predict(CoxBoost,newdata=x,
               newtime=survival_data$OS.time,
               newstatus=survival_data$OS.Status,
               at.step=optim.res$cv.res$optimal.step,type="lp")
  per<-t(per)
  c_index<-apply(per, 2, Cindex, y=y)
  per<-as.data.frame(per)
  colnames(per)<-c("predict")
  per$ID=rownames(per)
  Auc_data<-cbind(survival_data,per)
  #write.table(Auc_data, paste(Algorithm,"Auc_data.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  #geneid<-CoxBoost$xnames[CoxBoost$coefficients[which.max(pre),]!=0] ####当at.step=0:optim.res$cv.res$optimal.step时使用
  #########################
  library(survivalROC)
  auc<-data.frame()
  for (i in 1:3){
  ROC= survivalROC(Stime=Auc_data$OS.time,##生存时间
                   status=Auc_data$OS.Status,## 终止事件    
                   marker = Auc_data$predict, ## marker value    
                   #other_markers = as.matrix(feng_36[,c("LAG3","IDO1")]), # 矫正的协变量
                   predict.time = i,## 预测时间截点
                   method="KM"
  )##span,NNE法的namda
   if (i==1){auc<-as.data.frame(ROC$AUC)}
   else{auc<-cbind(auc,as.data.frame(ROC$AUC))}
  }
  colnames(auc)=c("1year","2year","3year")
  ################################
  
  #######################c_index/Auc
  #输出新数据的预测值，type参数允许选择预测的类型并提供预测值
  test.per<-predict(CoxBoost,newdata=test.x,
                    newtime=test.survival_data$OS.time,
                    newstatus=test.survival_data$OS.Status,
                    at.step=optim.res$cv.res$optimal.step,type="lp")
  test.per<-t(test.per)
  test.c_index<-apply(test.per, 2, Cindex, y=test.y)
  test.per<-as.data.frame(test.per)
  colnames(test.per)<-c("predict")
  test.per$ID=rownames(test.per)
  test.Auc_data<-cbind(test.survival_data,test.per)
 # write.table(test.Auc_data, paste(Algorithm,"test.Auc_data.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  #########################
  test.auc<-data.frame()
  for (i in 1:3){
  test.ROC= survivalROC(Stime=test.Auc_data$OS.time,##生存时间
                        status=test.Auc_data$OS.Status,## 终止事件    
                        marker = test.Auc_data$predict, ## marker value    
                        #other_markers = as.matrix(feng_36[,c("LAG3","IDO1")]), # 矫正的协变量
                        predict.time = i,## 预测时间截点
                        method="KM"
  )##span,NNE法的namda
   if (i==1){test.auc<-as.data.frame(test.ROC$AUC)}
   else{test.auc<-cbind(test.auc,as.data.frame(test.ROC$AUC))}
  }
  colnames(test.auc)=c("test1year","test2year","test3year")
    #######################c_index/Auc
  #输出新数据的预测值，type参数允许选择预测的类型并提供预测值
  GSE20685.per<-predict(CoxBoost,newdata=GSE20685.x,
                    newtime=GSE20685.survival_data$OS.time,
                    newstatus=GSE20685.survival_data$OS.Status,
                    at.step=optim.res$cv.res$optimal.step,type="lp")
  GSE20685.per<-t(GSE20685.per)
  GSE20685.c_index<-apply(GSE20685.per, 2, Cindex, y=GSE20685.y)
  GSE20685.per<-as.data.frame(GSE20685.per)
  colnames(GSE20685.per)<-c("predict")
  GSE20685.per$ID=rownames(GSE20685.per)
  GSE20685.Auc_data<-cbind(GSE20685.survival_data,GSE20685.per)
 # write.table(GSE20685.Auc_data, paste(Algorithm,"GSE20685.Auc_data.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  #########################
  GSE20685.auc<-data.frame()
  for (i in 1:3){
  GSE20685.ROC= survivalROC(Stime=GSE20685.Auc_data$OS.time,##生存时间
                        status=GSE20685.Auc_data$OS.Status,## 终止事件    
                        marker = GSE20685.Auc_data$predict, ## marker value    
                        #other_markers = as.matrix(feng_36[,c("LAG3","IDO1")]), # 矫正的协变量
                        predict.time = i,## 预测时间截点
                        method="KM"
  )##span,NNE法的namda
   if (i==1){GSE20685.auc<-as.data.frame(GSE20685.ROC$AUC)}
   else{GSE20685.auc<-cbind(GSE20685.auc,as.data.frame(GSE20685.ROC$AUC))}
  }
  colnames(GSE20685.auc)=c("GSE206851year","GSE206852year","GSE206853year")
   #输出新数据的预测值，type参数允许选择预测的类型并提供预测值
  Metabric.per<-predict(CoxBoost,newdata=Metabric.x,
                    newtime=Metabric.survival_data$OS.time,
                    newstatus=Metabric.survival_data$OS.Status,
                    at.step=optim.res$cv.res$optimal.step,type="lp")
  Metabric.per<-t(Metabric.per)
  Metabric.c_index<-apply(Metabric.per, 2, Cindex, y=Metabric.y)
  Metabric.per<-as.data.frame(Metabric.per)
  colnames(Metabric.per)<-c("predict")
  Metabric.per$ID=rownames(Metabric.per)
  Metabric.Auc_data<-cbind(Metabric.survival_data,Metabric.per)
 # write.table(Metabric.Auc_data, paste(Algorithm,"Metabric.Auc_data.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  #########################
  Metabric.auc<-data.frame()
  for (i in 1:3){
  Metabric.ROC= survivalROC(Stime=Metabric.Auc_data$OS.time,##生存时间
                        status=Metabric.Auc_data$OS.Status,## 终止事件    
                        marker = Metabric.Auc_data$predict, ## marker value    
                        #other_markers = as.matrix(feng_36[,c("LAG3","IDO1")]), # 矫正的协变量
                        predict.time = i,## 预测时间截点
                        method="KM"
  )##span,NNE法的namda
   if (i==1){Metabric.auc<-as.data.frame(Metabric.ROC$AUC)}
   else{Metabric.auc<-cbind(Metabric.auc,as.data.frame(Metabric.ROC$AUC))}
  }
  colnames(Metabric.auc)=c("Metabric1year","Metabric2year","Metabric3year")
  
#Auc.out<-cbind(as.data.frame(ROC$AUC),as.data.frame(test.ROC$AUC),as.data.frame(GSE20685.auc))
 # colnames(Auc.out)=c(name,test.name,GSE20685.name)
  Auc.out<-cbind(auc,test.auc,GSE20685.auc,Metabric.auc)
  Auc.out$Algorithm=Algorithm
  write.table(Auc.out,paste(Algorithm,".Auc.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  c_index<-cbind(as.data.frame(c_index) ,as.data.frame(test.c_index),as.data.frame(GSE20685.c_index),as.data.frame(Metabric.c_index)) 
  colnames(c_index)=c(name,test.name,GSE20685.name,Metabric.name)
  c_index$Algorithm=Algorithm
  write.table(c_index,paste(Algorithm,".c_index.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
}
gbm.ar<-function(x,survival_data,Algorithm,name,test.x,test.survival_data,test.name,GSE20685.x,GSE20685.survival_data,GSE20685.name,Metabric.x,Metabric.survival_data,Metabric.name){
  set.seed(seed)
  library(gbm)
  y<-Surv(time=as.double(survival_data$OS.time),event=as.double(survival_data$OS.Status))
  test.y<-Surv(time=as.double(test.survival_data$OS.time),event=as.double(test.survival_data$OS.Status))
  GSE20685.y<-Surv(time=as.double(GSE20685.survival_data$OS.time),event=as.double(GSE20685.survival_data$OS.Status))
  Metabric.y<-Surv(time=as.double(Metabric.survival_data$OS.time),event=as.double(Metabric.survival_data$OS.Status))
  gbm.fit<- gbm.fit(x=x,y=y,distribution = 'coxph',
                    n.trees = 10000,
                    interaction.depth = 3,
                    n.minobsinnode = 10,
                    shrinkage = 0.001)
  #cv.folds = 10,n.cores = 6)
  saveRDS(gbm.fit,"gbm.fit.rds")
  # find index for number trees with minimum CV error
  best.iter <- gbm.perf(gbm.fit)
  #############################################3.构建模型
  per<-predict(gbm.fit, newdata = x, n.trees = best.iter)
  per<-t(as.matrix(t(per)))
  c_index<-apply(per, 2, Cindex, y=y)
  per<-as.data.frame(per)
  colnames(per)<-c("predict")
  per$ID=rownames(per)
  Auc_data<-cbind(survival_data,per)
  #write.table(Auc_data, paste(Algorithm,"Auc_data.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  #geneid<-CoxBoost$xnames[CoxBoost$coefficients[which.max(pre),]!=0
  #########################
  library(survivalROC)
  auc<-data.frame()
  for (i in 1:3){
  ROC= survivalROC(Stime=Auc_data$OS.time,##生存时间
                   status=Auc_data$OS.Status,## 终止事件    
                   marker = Auc_data$predict, ## marker value    
                   #other_markers = as.matrix(feng_36[,c("LAG3","IDO1")]), # 矫正的协变量
                   predict.time = i,## 预测时间截点
                   method="KM"
  )##span,NNE法的namda
   if (i==1){auc<-as.data.frame(ROC$AUC)}
   else{auc<-cbind(auc,as.data.frame(ROC$AUC))}
  }
  colnames(auc)=c("1year","2year","3year")
  

  #######################c_index/Auc
  #输出新数据的预测值，type参数允许选择预测的类型并提供预测值
  test.per<-predict(gbm.fit,newdata=test.x,n.trees = best.iter)
  test.per<-t(as.matrix(t(test.per)))
  test.c_index<-apply(test.per, 2, Cindex, y=test.y)
  test.per<-as.data.frame(test.per)
  colnames(test.per)<-c("predict")
  test.per$ID=rownames(test.per)
  test.Auc_data<-cbind(test.survival_data,test.per)
  #write.table(test.Auc_data, paste(Algorithm,"test.Auc_data.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  #########################
  test.auc<-data.frame()
  for (i in 1:3){
  test.ROC= survivalROC(Stime=test.Auc_data$OS.time,##生存时间
                        status=test.Auc_data$OS.Status,## 终止事件    
                        marker = test.Auc_data$predict, ## marker value    
                        #other_markers = as.matrix(feng_36[,c("LAG3","IDO1")]), # 矫正的协变量
                        predict.time = i,## 预测时间截点
                        method="KM"
  )##span,NNE法的namda
   if (i==1){test.auc<-as.data.frame(test.ROC$AUC)}
   else{test.auc<-cbind(test.auc,as.data.frame(test.ROC$AUC))}
  }
  colnames(test.auc)=c("test1year","test2year","test3year")
   #######################c_index/Auc
  #输出新数据的预测值，type参数允许选择预测的类型并提供预测值
  GSE20685.per<-predict(gbm.fit,newdata=GSE20685.x,n.trees = best.iter)
  GSE20685.per<-t(as.matrix(t(GSE20685.per)))
  GSE20685.c_index<-apply(GSE20685.per, 2, Cindex, y=GSE20685.y)
  GSE20685.per<-as.data.frame(GSE20685.per)
  colnames(GSE20685.per)<-c("predict")
  GSE20685.per$ID=rownames(GSE20685.per)
  GSE20685.Auc_data<-cbind(GSE20685.survival_data,GSE20685.per)
  #write.table(GSE20685.Auc_data, paste(Algorithm,"GSE20685.Auc_data.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  #########################
  GSE20685.auc<-data.frame()
  for (i in 1:3){
  GSE20685.ROC= survivalROC(Stime=GSE20685.Auc_data$OS.time,##生存时间
                        status=GSE20685.Auc_data$OS.Status,## 终止事件    
                        marker = GSE20685.Auc_data$predict, ## marker value    
                        #other_markers = as.matrix(feng_36[,c("LAG3","IDO1")]), # 矫正的协变量
                        predict.time = i,## 预测时间截点
                        method="KM"
  )##span,NNE法的namda
   if (i==1){GSE20685.auc<-as.data.frame(GSE20685.ROC$AUC)}
   else{GSE20685.auc<-cbind(GSE20685.auc,as.data.frame(GSE20685.ROC$AUC))}
  }
  colnames(GSE20685.auc)=c("GSE206851year","GSE206852year","GSE206853year")
   Auc.out<-cbind(auc,test.auc,GSE20685.auc)
    #######################c_index/Auc
  #输出新数据的预测值，type参数允许选择预测的类型并提供预测值
  Metabric.per<-predict(gbm.fit,newdata=Metabric.x,n.trees = best.iter)
  Metabric.per<-t(as.matrix(t(Metabric.per)))
  Metabric.c_index<-apply(Metabric.per, 2, Cindex, y=Metabric.y)
  Metabric.per<-as.data.frame(Metabric.per)
  colnames(Metabric.per)<-c("predict")
  Metabric.per$ID=rownames(Metabric.per)
  Metabric.Auc_data<-cbind(Metabric.survival_data,Metabric.per)
  #write.table(Metabric.Auc_data, paste(Algorithm,"Metabric.Auc_data.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  #########################
  Metabric.auc<-data.frame()
  for (i in 1:3){
  Metabric.ROC= survivalROC(Stime=Metabric.Auc_data$OS.time,##生存时间
                        status=Metabric.Auc_data$OS.Status,## 终止事件    
                        marker = Metabric.Auc_data$predict, ## marker value    
                        #other_markers = as.matrix(feng_36[,c("LAG3","IDO1")]), # 矫正的协变量
                        predict.time = i,## 预测时间截点
                        method="KM"
  )##span,NNE法的namda
   if (i==1){Metabric.auc<-as.data.frame(Metabric.ROC$AUC)}
   else{Metabric.auc<-cbind(Metabric.auc,as.data.frame(Metabric.ROC$AUC))}
  }
  colnames(Metabric.auc)=c("Metabric1year","Metabric2year","Metabric3year")
   Auc.out<-cbind(auc,test.auc,Metabric.auc)
  
  #Auc.out<-cbind(as.data.frame(ROC$AUC),as.data.frame(test.ROC$AUC),as.data.frame(GSE20685.auc))
  #colnames(Auc.out)=c(name,test.name,GSE20685.name)
  Auc.out<-cbind(auc,test.auc,GSE20685.auc,Metabric.auc)
  Auc.out$Algorithm=Algorithm
  write.table(Auc.out,paste(Algorithm,".Auc.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  c_index<-cbind(as.data.frame(c_index) ,as.data.frame(test.c_index),as.data.frame(GSE20685.c_index),as.data.frame(Metabric.c_index)) 
  colnames(c_index)=c(name,test.name,GSE20685.name,Metabric.name)
  c_index$Algorithm=Algorithm
  write.table(c_index,paste(Algorithm,".c_index.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
}
######################################randomForestSRC

randomForestSRC.ar<-function(x,survival_data,Algorithm,name,test.x,test.survival_data,test.name,GSE20685.x,GSE20685.survival_data,GSE20685.name,Metabric.x,Metabric.survival_data,Metabric.name){
y<-data.matrix(Surv(time=as.double(survival_data$OS.time),event=as.double(survival_data$OS.Status)))
merge<- cbind(survival_data,x) %>% dplyr::select(-c("ID"))
#h
#######################2.randomForestSRC
library(randomForestSRC)
##########PE (true OOB)      : 33.2984 
set.seed(seed)
rf.model <-rfsrc(Surv(OS.time, OS.Status) ~ . , data=merge, splitrule = "logrank",sampletye="swor",
                 # perf.type="none",
                 ntree = 1000, nodesize = 15,importance=TRUE,
                 #samplesize=517,
                 block.size=100, nsplit=10,save.memory = TRUE)
vs.pbc<-var.select(object=rf.model)#注意，一般用method="md"，但如果特征数量远大于样本数（10倍以上），要选用method="vh"，其他参数请看说明书按需设定。
geneids<-vs.pbc$topvars
saveRDS(rf.model,paste(Algorithm,".pre.rds",sep=""))
write.table(geneids, paste(Algorithm,"_selected_gene.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
#########################3.构建模型
pre<-predict(rf.model,newdata=merge,importance = TRUE)
per=pre$predicted
c_index<-apply(as.matrix(per), 2, Cindex, y=y)
per<-as.data.frame(per)
colnames(per)<-c("predict")
per$ID=rownames(per)
Auc_data<-cbind(survival_data,per)
#write.table(Auc_data, paste(Algorithm,"Auc_data.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
#geneid<-CoxBoost$xnames[CoxBoost$coefficients[which.max(pre),]!=0] ####当at.step=0:optim.res$cv.res$optimal.step时使用
#########################
library(survivalROC)
 auc<-data.frame()
  for (i in 1:3){
  ROC= survivalROC(Stime=Auc_data$OS.time,##生存时间
                   status=Auc_data$OS.Status,## 终止事件    
                   marker = Auc_data$predict, ## marker value    
                   #other_markers = as.matrix(feng_36[,c("LAG3","IDO1")]), # 矫正的协变量
                   predict.time = i,## 预测时间截点
                   method="KM"
  )##span,NNE法的namda
   if (i==1){auc<-as.data.frame(ROC$AUC)}
   else{auc<-cbind(auc,as.data.frame(ROC$AUC))}
  }
  colnames(auc)=c("1year","2year","3year")
################################
test.y<-data.matrix(Surv(time=as.double(test.survival_data$OS.time),event=as.double(test.survival_data$OS.Status)))
test.merge<- cbind(test.survival_data,test.x) %>% dplyr::select(-c("ID"))
test.merge$OS.Status<-as.factor(test.merge$OS.Status)
test.merge$OS.time<-as.numeric(test.merge$OS.time)

#######################c_index/Auc
#输出新数据的预测值，type参数允许选择预测的类型并提供预测值
test.pre<-predict(rf.model,newdata=test.merge,importance = TRUE)
test.per=test.pre$predicted
test.c_index<-apply(as.matrix(test.per), 2, Cindex, y=test.y)
test.per<-as.data.frame(test.per)
colnames(test.per)<-c("predict")
test.per$ID=rownames(test.per)
test.Auc_data<-cbind(test.survival_data,test.per)
#write.table(test.Auc_data, paste(Algorithm,"test.Auc_data.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
#########################
  test.auc<-data.frame()
  for (i in 1:3){
  test.ROC= survivalROC(Stime=test.Auc_data$OS.time,##生存时间
                        status=test.Auc_data$OS.Status,## 终止事件    
                        marker = test.Auc_data$predict, ## marker value    
                        #other_markers = as.matrix(feng_36[,c("LAG3","IDO1")]), # 矫正的协变量
                        predict.time = i,## 预测时间截点
                        method="KM"
  )##span,NNE法的namda
   if (i==1){test.auc<-as.data.frame(test.ROC$AUC)}
   else{test.auc<-cbind(test.auc,as.data.frame(test.ROC$AUC))}
  }
  colnames(test.auc)=c("test1year","test2year","test3year")
################################
GSE20685.y<-data.matrix(Surv(time=as.double(GSE20685.survival_data$OS.time),event=as.double(GSE20685.survival_data$OS.Status)))
GSE20685.merge<- cbind(GSE20685.survival_data,GSE20685.x) %>% dplyr::select(-c("ID"))
GSE20685.merge$OS.Status<-as.factor(GSE20685.merge$OS.Status)
GSE20685.merge$OS.time<-as.numeric(GSE20685.merge$OS.time)

#######################c_index/Auc
#输出新数据的预测值，type参数允许选择预测的类型并提供预测值
GSE20685.pre<-predict(rf.model,newdata=GSE20685.merge,importance = TRUE)
GSE20685.per=GSE20685.pre$predicted
GSE20685.c_index<-apply(as.matrix(GSE20685.per), 2, Cindex, y=GSE20685.y)
GSE20685.per<-as.data.frame(GSE20685.per)
colnames(GSE20685.per)<-c("predict")
GSE20685.per$ID=rownames(GSE20685.per)
GSE20685.Auc_data<-cbind(GSE20685.survival_data,GSE20685.per)
#write.table(GSE20685.Auc_data, paste(Algorithm,"GSE20685.Auc_data.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
#########################
  GSE20685.auc<-data.frame()
  for (i in 1:3){
  GSE20685.ROC= survivalROC(Stime=GSE20685.Auc_data$OS.time,##生存时间
                        status=GSE20685.Auc_data$OS.Status,## 终止事件    
                        marker = GSE20685.Auc_data$predict, ## marker value    
                        #other_markers = as.matrix(feng_36[,c("LAG3","IDO1")]), # 矫正的协变量
                        predict.time = i,## 预测时间截点
                        method="KM"
  )##span,NNE法的namda
   if (i==1){GSE20685.auc<-as.data.frame(GSE20685.ROC$AUC)}
   else{GSE20685.auc<-cbind(GSE20685.auc,as.data.frame(GSE20685.ROC$AUC))}
  }
  colnames(GSE20685.auc)=c("GSE206851year","GSE206852year","GSE206853year")
  ################################
Metabric.y<-data.matrix(Surv(time=as.double(Metabric.survival_data$OS.time),event=as.double(Metabric.survival_data$OS.Status)))
Metabric.merge<- cbind(Metabric.survival_data,Metabric.x) %>% dplyr::select(-c("ID"))
Metabric.merge$OS.Status<-as.factor(Metabric.merge$OS.Status)
Metabric.merge$OS.time<-as.numeric(Metabric.merge$OS.time)

#######################c_index/Auc
#输出新数据的预测值，type参数允许选择预测的类型并提供预测值
Metabric.pre<-predict(rf.model,newdata=Metabric.merge,importance = TRUE)
Metabric.per=Metabric.pre$predicted
Metabric.c_index<-apply(as.matrix(Metabric.per), 2, Cindex, y=Metabric.y)
Metabric.per<-as.data.frame(Metabric.per)
colnames(Metabric.per)<-c("predict")
Metabric.per$ID=rownames(Metabric.per)
Metabric.Auc_data<-cbind(Metabric.survival_data,Metabric.per)
#write.table(Metabric.Auc_data, paste(Algorithm,"Metabric.Auc_data.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
#########################
  Metabric.auc<-data.frame()
  for (i in 1:3){
  Metabric.ROC= survivalROC(Stime=Metabric.Auc_data$OS.time,##生存时间
                        status=Metabric.Auc_data$OS.Status,## 终止事件    
                        marker = Metabric.Auc_data$predict, ## marker value    
                        #other_markers = as.matrix(feng_36[,c("LAG3","IDO1")]), # 矫正的协变量
                        predict.time = i,## 预测时间截点
                        method="KM"
  )##span,NNE法的namda
   if (i==1){Metabric.auc<-as.data.frame(Metabric.ROC$AUC)}
   else{Metabric.auc<-cbind(Metabric.auc,as.data.frame(Metabric.ROC$AUC))}
  }
  colnames(Metabric.auc)=c("Metabric1year","Metabric2year","Metabric3year")

#################
  Auc.out<-cbind(auc,test.auc,GSE20685.auc,Metabric.auc)
  Auc.out$Algorithm=Algorithm
  write.table(Auc.out,paste(Algorithm,".Auc.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  c_index<-cbind(as.data.frame(c_index) ,as.data.frame(test.c_index),as.data.frame(GSE20685.c_index),as.data.frame(Metabric.c_index)) 
  colnames(c_index)=c(name,test.name,GSE20685.name,Metabric.name)
  c_index$Algorithm=Algorithm
  write.table(c_index,paste(Algorithm,".c_index.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
}
step.back.ar<-function(x,survival_data,Algorithm,name,test.x,test.survival_data,test.name,GSE20685.x,GSE20685.survival_data,GSE20685.name,Metabric.x,Metabric.survival_data,Metabric.name){
   y<-Surv(time=as.double(survival_data$OS.time),event=as.double(survival_data$OS.Status))
   test.y<-data.matrix(Surv(time=as.double(test.survival_data$OS.time),event=as.double(test.survival_data$OS.Status)))
   GSE20685.y<-data.matrix(Surv(time=as.double(GSE20685.survival_data$OS.time),event=as.double(GSE20685.survival_data$OS.Status)))
   Metabric.y<-Surv(time=as.double(Metabric.survival_data$OS.time),event=as.double(Metabric.survival_data$OS.Status))
  merge<- cbind(survival_data,x) %>% dplyr::select(-c("ID"))
  colnames(merge)<-gsub("-","_",colnames(merge))
  formula_string = paste0("Surv(OS.time, OS.Status)~",paste0(colnames(merge)[3:ncol(merge)],collapse = "+"))
  formula_expression= as.formula(formula_string)
  ####################################2.stepAIC
  fit <- coxph( formula_expression, data=merge )
  #summary(fit)
  step.fit.back <- stepAIC(fit,direction = "backward")

  #############################3.模型预测
  back.geneids<- as.data.frame(coefficients(step.fit.back))
  back.geneids$name<-rownames(back.geneids)

  colnames(back.geneids)=c("cofficients","name")
  back.geneids<-dplyr::select(back.geneids,c("name","cofficients"))
  write.table(back.geneids, paste(Algorithm,"_selected_gene.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  per<-predict(step.fit.back)
  per<-t(as.matrix(t(per)))
  c_index<-apply(per, 2, Cindex, y=y)
  per<-as.data.frame(per)
  colnames(per)<-c("predict")
  per$ID=rownames(per)
  Auc_data<-cbind(survival_data,per)
 # write.table(Auc_data, paste(Algorithm,"Auc_data.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  #geneid<-CoxBoost$xnames[CoxBoost$coefficients[which.max(pre),]!=0
  #########################
  auc<-data.frame()
  for (i in 1:3){
  ROC= survivalROC(Stime=Auc_data$OS.time,##生存时间
                   status=Auc_data$OS.Status,## 终止事件    
                   marker = Auc_data$predict, ## marker value    
                   #other_markers = as.matrix(feng_36[,c("LAG3","IDO1")]), # 矫正的协变量
                   predict.time = i,## 预测时间截点
                   method="KM"
  )##span,NNE法的namda
   if (i==1){auc<-as.data.frame(ROC$AUC)}
   else{auc<-cbind(auc,as.data.frame(ROC$AUC))}
  }
  colnames(auc)=c("1year","2year","3year")
  ###############################################4.在test中进行predict
  test.merge<- cbind(test.survival_data,test.x) %>% dplyr::select(-c("ID"))
  colnames(test.merge)<-gsub("-","_",colnames(test.merge))
  #######################c_index/Auc
  #输出新数据的预测值，type参数允许选择预测的类型并提供预测值
  test.per<-predict(step.fit.back,newdata=test.merge)
  test.per<-t(as.matrix(t(test.per)))
  test.c_index<-apply(test.per, 2, Cindex, y=test.y)
  test.per<-as.data.frame(test.per)
  colnames(test.per)<-c("predict")
  test.per$ID=rownames(test.per)
  test.Auc_data<-cbind(test.survival_data,test.per)
 # write.table(test.Auc_data, paste(Algorithm,"test.Auc_data.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  #########################
   test.auc<-data.frame()
  for (i in 1:3){
  test.ROC= survivalROC(Stime=test.Auc_data$OS.time,##生存时间
                        status=test.Auc_data$OS.Status,## 终止事件    
                        marker = test.Auc_data$predict, ## marker value    
                        #other_markers = as.matrix(feng_36[,c("LAG3","IDO1")]), # 矫正的协变量
                        predict.time = i,## 预测时间截点
                        method="KM"
  )##span,NNE法的namda
   if (i==1){test.auc<-as.data.frame(test.ROC$AUC)}
   else{test.auc<-cbind(test.auc,as.data.frame(test.ROC$AUC))}
  }
  colnames(test.auc)=c("test1year","test2year","test3year")
  ###############################################4.在GSE20685中进行predict
  GSE20685.merge<- cbind(GSE20685.survival_data,GSE20685.x) %>% dplyr::select(-c("ID"))
  colnames(GSE20685.merge)<-gsub("-","_",colnames(GSE20685.merge))
  #######################c_index/Auc
  #输出新数据的预测值，type参数允许选择预测的类型并提供预测值
  GSE20685.per<-predict(step.fit.back,newdata=GSE20685.merge)
  GSE20685.per<-t(as.matrix(t(GSE20685.per)))
  GSE20685.c_index<-apply(GSE20685.per, 2, Cindex, y=GSE20685.y)
  GSE20685.per<-as.data.frame(GSE20685.per)
  colnames(GSE20685.per)<-c("predict")
  GSE20685.per$ID=rownames(GSE20685.per)
  GSE20685.Auc_data<-cbind(GSE20685.survival_data,GSE20685.per)
 # write.table(GSE20685.Auc_data, paste(Algorithm,"GSE20685.Auc_data.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  #########################
  GSE20685.auc<-data.frame()
  for (i in 1:3){
  GSE20685.ROC= survivalROC(Stime=GSE20685.Auc_data$OS.time,##生存时间
                        status=GSE20685.Auc_data$OS.Status,## 终止事件    
                        marker = GSE20685.Auc_data$predict, ## marker value    
                        #other_markers = as.matrix(feng_36[,c("LAG3","IDO1")]), # 矫正的协变量
                        predict.time = i,## 预测时间截点
                        method="KM"
  )##span,NNE法的namda
   if (i==1){GSE20685.auc<-as.data.frame(GSE20685.ROC$AUC)}
   else{GSE20685.auc<-cbind(GSE20685.auc,as.data.frame(GSE20685.ROC$AUC))}
  }
  colnames(GSE20685.auc)=c("GSE206851year","GSE206852year","GSE206853year")
  
   ###############################################4.在Metabric中进行predict
  Metabric.merge<- cbind(Metabric.survival_data,Metabric.x) %>% dplyr::select(-c("ID"))
  colnames(Metabric.merge)<-gsub("-","_",colnames(Metabric.merge))
  #######################c_index/Auc
  #输出新数据的预测值，type参数允许选择预测的类型并提供预测值
  Metabric.per<-predict(step.fit.back,newdata=Metabric.merge)
  Metabric.per<-t(as.matrix(t(Metabric.per)))
  Metabric.c_index<-apply(Metabric.per, 2, Cindex, y=Metabric.y)
  Metabric.per<-as.data.frame(Metabric.per)
  colnames(Metabric.per)<-c("predict")
  Metabric.per$ID=rownames(Metabric.per)
  Metabric.Auc_data<-cbind(Metabric.survival_data,Metabric.per)
 # write.table(Metabric.Auc_data, paste(Algorithm,"Metabric.Auc_data.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  #########################
  Metabric.auc<-data.frame()
  for (i in 1:3){
  Metabric.ROC= survivalROC(Stime=Metabric.Auc_data$OS.time,##生存时间
                        status=Metabric.Auc_data$OS.Status,## 终止事件    
                        marker = Metabric.Auc_data$predict, ## marker value    
                        #other_markers = as.matrix(feng_36[,c("LAG3","IDO1")]), # 矫正的协变量
                        predict.time = i,## 预测时间截点
                        method="KM"
  )##span,NNE法的namda
   if (i==1){Metabric.auc<-as.data.frame(Metabric.ROC$AUC)}
   else{Metabric.auc<-cbind(Metabric.auc,as.data.frame(Metabric.ROC$AUC))}
  }
  colnames(Metabric.auc)=c("Metabric1year","Metabric2year","Metabric3year")
  
  ##############
   #Auc.out<-cbind(as.data.frame(ROC$AUC),as.data.frame(test.ROC$AUC),as.data.frame(GSE20685.auc))
  #colnames(Auc.out)=c(name,test.name,GSE20685.name)
 Auc.out<-cbind(auc,test.auc,GSE20685.auc,Metabric.auc)
  Auc.out$Algorithm=Algorithm
  write.table(Auc.out,paste(Algorithm,".Auc.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  c_index<-cbind(as.data.frame(c_index) ,as.data.frame(test.c_index),as.data.frame(GSE20685.c_index),as.data.frame(Metabric.c_index)) 
  colnames(c_index)=c(name,test.name,GSE20685.name,Metabric.name)
  c_index$Algorithm=Algorithm
  write.table(c_index,paste(Algorithm,".c_index.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
}
step.both.ar<-function(x,survival_data,Algorithm,name,test.x,test.survival_data,test.name,GSE20685.x,GSE20685.survival_data,GSE20685.name,Metabric.x,Metabric.survival_data,Metabric.name){
  y<-Surv(time=as.double(survival_data$OS.time),event=as.double(survival_data$OS.Status))
  test.y<-data.matrix(Surv(time=as.double(test.survival_data$OS.time),event=as.double(test.survival_data$OS.Status)))
   GSE20685.y<-data.matrix(Surv(time=as.double(GSE20685.survival_data$OS.time),event=as.double(GSE20685.survival_data$OS.Status)))
   Metabric.y<-Surv(time=as.double(Metabric.survival_data$OS.time),event=as.double(Metabric.survival_data$OS.Status))
  merge<- cbind(survival_data,x) %>% dplyr::select(-c("ID"))
  colnames(merge)<-gsub("-","_",colnames(merge))
  formula_string = paste0("Surv(OS.time, OS.Status)~",paste0(colnames(merge)[3:ncol(merge)],collapse = "+"))
  formula_expression= as.formula(formula_string)
  ####################################2.stepAIC
  fit <- coxph( formula_expression, data=merge )
  #summary(fit)
  step.fit.both <- stepAIC(fit,direction = "both")
  
  #############################3.模型预测
  both.geneids<- as.data.frame(coefficients(step.fit.both))
  both.geneids$name<-rownames(both.geneids)

  colnames(both.geneids)=c("cofficients","name")
  both.geneids<-dplyr::select(both.geneids,c("name","cofficients"))
  write.table(both.geneids, paste(Algorithm,"_selected_gene.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  per<-predict(step.fit.both)
  per<-t(as.matrix(t(per)))
  c_index<-apply(per, 2, Cindex, y=y)
  per<-as.data.frame(per)
  colnames(per)<-c("predict")
  per$ID=rownames(per)
  Auc_data<-cbind(survival_data,per)
  #write.table(Auc_data, paste(Algorithm,"Auc_data.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  #geneid<-CoxBoost$xnames[CoxBoost$coefficients[which.max(pre),]!=0
  #########################
  auc<-data.frame()
  for (i in 1:3){
  ROC= survivalROC(Stime=Auc_data$OS.time,##生存时间
                   status=Auc_data$OS.Status,## 终止事件    
                   marker = Auc_data$predict, ## marker value    
                   #other_markers = as.matrix(feng_36[,c("LAG3","IDO1")]), # 矫正的协变量
                   predict.time = i,## 预测时间截点
                   method="KM"
  )##span,NNE法的namda
   if (i==1){auc<-as.data.frame(ROC$AUC)}
   else{auc<-cbind(auc,as.data.frame(ROC$AUC))}
  }
  colnames(auc)=c("1year","2year","3year")
  ###############################################4.在test中进行predict
  test.merge<- cbind(test.survival_data,test.x) %>% dplyr::select(-c("ID"))
  colnames(test.merge)<-gsub("-","_",colnames(test.merge))
  #######################c_index/Auc
  #输出新数据的预测值，type参数允许选择预测的类型并提供预测值
  test.per<-predict(step.fit.both,newdata=test.merge)
  test.per<-t(as.matrix(t(test.per)))
  test.c_index<-apply(test.per, 2, Cindex, y=test.y)
  test.per<-as.data.frame(test.per)
  colnames(test.per)<-c("predict")
  test.per$ID=rownames(test.per)
  test.Auc_data<-cbind(test.survival_data,test.per)
 # write.table(test.Auc_data, paste(Algorithm,"test.Auc_data.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  #########################
  test.auc<-data.frame()
  for (i in 1:3){
  test.ROC= survivalROC(Stime=test.Auc_data$OS.time,##生存时间
                        status=test.Auc_data$OS.Status,## 终止事件    
                        marker = test.Auc_data$predict, ## marker value    
                        #other_markers = as.matrix(feng_36[,c("LAG3","IDO1")]), # 矫正的协变量
                        predict.time = i,## 预测时间截点
                        method="KM"
  )##span,NNE法的namda
   if (i==1){test.auc<-as.data.frame(test.ROC$AUC)}
   else{test.auc<-cbind(test.auc,as.data.frame(test.ROC$AUC))}
  }
  colnames(test.auc)=c("test1year","test2year","test3year")
###############################################4.在GSE20685中进行predict
  GSE20685.merge<- cbind(GSE20685.survival_data,GSE20685.x) %>% dplyr::select(-c("ID"))
  colnames(GSE20685.merge)<-gsub("-","_",colnames(GSE20685.merge))
  #######################c_index/Auc
  #输出新数据的预测值，type参数允许选择预测的类型并提供预测值
  GSE20685.per<-predict(step.fit.both,newdata=GSE20685.merge)
  GSE20685.per<-t(as.matrix(t(GSE20685.per)))
  GSE20685.c_index<-apply(GSE20685.per, 2, Cindex, y=GSE20685.y)
  GSE20685.per<-as.data.frame(GSE20685.per)
  colnames(GSE20685.per)<-c("predict")
  GSE20685.per$ID=rownames(GSE20685.per)
  GSE20685.Auc_data<-cbind(GSE20685.survival_data,GSE20685.per)
  #write.table(GSE20685.Auc_data, paste(Algorithm,"GSE20685.Auc_data.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  #########################
  GSE20685.auc<-data.frame()
  for (i in 1:3){
  GSE20685.ROC= survivalROC(Stime=GSE20685.Auc_data$OS.time,##生存时间
                        status=GSE20685.Auc_data$OS.Status,## 终止事件    
                        marker = GSE20685.Auc_data$predict, ## marker value    
                        #other_markers = as.matrix(feng_36[,c("LAG3","IDO1")]), # 矫正的协变量
                        predict.time = i,## 预测时间截点
                        method="KM"
  )##span,NNE法的namda
   if (i==1){GSE20685.auc<-as.data.frame(GSE20685.ROC$AUC)}
   else{GSE20685.auc<-cbind(GSE20685.auc,as.data.frame(GSE20685.ROC$AUC))}
  }
  colnames(GSE20685.auc)=c("GSE206851year","GSE206852year","GSE206853year")
  
  ###############################################4.在Metabric中进行predict
  Metabric.merge<- cbind(Metabric.survival_data,Metabric.x) %>% dplyr::select(-c("ID"))
  colnames(Metabric.merge)<-gsub("-","_",colnames(Metabric.merge))
  #######################c_index/Auc
  #输出新数据的预测值，type参数允许选择预测的类型并提供预测值
  Metabric.per<-predict(step.fit.both,newdata=Metabric.merge)
  Metabric.per<-t(as.matrix(t(Metabric.per)))
  Metabric.c_index<-apply(Metabric.per, 2, Cindex, y=Metabric.y)
  Metabric.per<-as.data.frame(Metabric.per)
  colnames(Metabric.per)<-c("predict")
  Metabric.per$ID=rownames(Metabric.per)
  Metabric.Auc_data<-cbind(Metabric.survival_data,Metabric.per)
  #write.table(Metabric.Auc_data, paste(Algorithm,"Metabric.Auc_data.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  #########################
  Metabric.auc<-data.frame()
  for (i in 1:3){
  Metabric.ROC= survivalROC(Stime=Metabric.Auc_data$OS.time,##生存时间
                        status=Metabric.Auc_data$OS.Status,## 终止事件    
                        marker = Metabric.Auc_data$predict, ## marker value    
                        #other_markers = as.matrix(feng_36[,c("LAG3","IDO1")]), # 矫正的协变量
                        predict.time = i,## 预测时间截点
                        method="KM"
  )##span,NNE法的namda
   if (i==1){Metabric.auc<-as.data.frame(Metabric.ROC$AUC)}
   else{Metabric.auc<-cbind(Metabric.auc,as.data.frame(Metabric.ROC$AUC))}
  }
  colnames(Metabric.auc)=c("Metabric1year","Metabric2year","Metabric3year")
 ##############
   #Auc.out<-cbind(as.data.frame(ROC$AUC),as.data.frame(test.ROC$AUC),as.data.frame(GSE20685.auc))
  #colnames(Auc.out)=c(name,test.name,GSE20685.name)
  Auc.out<-cbind(auc,test.auc,GSE20685.auc,Metabric.auc)
  Auc.out$Algorithm=Algorithm
  write.table(Auc.out,paste(Algorithm,".Auc.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
  c_index<-cbind(as.data.frame(c_index) ,as.data.frame(test.c_index),as.data.frame(GSE20685.c_index),as.data.frame(Metabric.c_index)) 
  colnames(c_index)=c(name,test.name,GSE20685.name,Metabric.name)
  c_index$Algorithm=Algorithm
  write.table(c_index,paste(Algorithm,".c_index.xls",sep=""),sep = "\t",quote = FALSE,row.names =FALSE)
}
