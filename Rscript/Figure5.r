#######################ROC曲线分析
####数据处理  TCGA.BRCA.fpkm来自Figure1
NSUN.FPKM<-TCGA.BRCA.fpkm[c("NOP2","NSUN2","NSUN3","NSUN4","NSUN5","NSUN6","NSUN7"),]
rownames(NSUN.FPKM)
rownames(NSUN.FPKM)<-c("NSUN1","NSUN2","NSUN3","NSUN4","NSUN5","NSUN6","NSUN7")
NSUNS.level<-as.data.frame(t(NSUN.FPKM))
rownames(NSUNS.level)<-colnames(NSUN.FPKM)
NSUNS.level$treat<-substr(rownames(NSUNS.level),14,15)
NSUNS.level<-dplyr::filter(NSUNS.level,treat %in% c("01","11"))
NSUNS.level$group="Primary tumor"
NSUNS.level[which(NSUNS.level$treat=="11"),]$group="Normal"
NSUNS.level<-NSUNS.level %>% pivot_longer(cols=-c("treat","group"), names_to = "Gene",values_to = "FPKM")
######ROC画图
library(pROC)
library(ggplot2)
NSUN1<-dplyr::filter(NSUNS.level,Gene=="NSUN1")
NSUN1.res<-roc(group~FPKM,data=NSUN1,aur=TRUE,
         ci=TRUE, # 显示95%CI
         #percent=TRUE, # 是否需要以百分比显示
         smooth=TRUE,# 是否平滑曲线
         levels=c('Normal','Primary tumor'),direction="<" #设置分组方向
         )
NSUN2<-dplyr::filter(NSUNS.level,Gene=="NSUN2")
NSUN2.res<-roc(group~FPKM,data=NSUN2,aur=TRUE,
         ci=TRUE, # 显示95%CI
         #percent=TRUE, # 是否需要以百分比显示
         smooth=TRUE,# 是否平滑曲线
         levels=c('Normal','Primary tumor'),direction="<" #设置分组方向
         )
NSUN3<-dplyr::filter(NSUNS.level,Gene=="NSUN3")
NSUN3.res<-roc(group~FPKM,data=NSUN3,aur=TRUE,
         ci=TRUE, # 显示95%CI
         #percent=TRUE, # 是否需要以百分比显示
         smooth=TRUE,# 是否平滑曲线
         levels=c('Normal','Primary tumor'),direction="<" #设置分组方向
         )
NSUN4<-dplyr::filter(NSUNS.level,Gene=="NSUN4")
NSUN4.res<-roc(group~FPKM,data=NSUN4,aur=TRUE,
         ci=TRUE, # 显示95%CI
         #percent=TRUE, # 是否需要以百分比显示
         smooth=TRUE,# 是否平滑曲线
         levels=c('Normal','Primary tumor'),direction="<" #设置分组方向
         )
NSUN5<-dplyr::filter(NSUNS.level,Gene=="NSUN5")
NSUN5.res<-roc(group~FPKM,data=NSUN5,aur=TRUE,
         ci=TRUE, # 显示95%CI
         #percent=TRUE, # 是否需要以百分比显示
         smooth=TRUE,# 是否平滑曲线
         levels=c('Normal','Primary tumor'),direction="<" #设置分组方向
         )

NSUN6<-dplyr::filter(NSUNS.level,Gene=="NSUN6")
NSUN6.res<-roc(group~FPKM,data=NSUN6,aur=TRUE,
         ci=TRUE, # 显示95%CI
         #percent=TRUE, # 是否需要以百分比显示
         smooth=TRUE,# 是否平滑曲线
         levels=c('Normal','Primary tumor'),direction="<" #设置分组方向
         )
NSUN7<-dplyr::filter(NSUNS.level,Gene=="NSUN7")
NSUN7.res<-roc(group~FPKM,data=NSUN7,aur=TRUE,
         ci=TRUE, # 显示95%CI
         #percent=TRUE, # 是否需要以百分比显示
         smooth=TRUE,# 是否平滑曲线
         levels=c('Normal','Primary tumor'),direction="<" #设置分组方向
         )

####
options(repr.plot.height=6, repr.plot.width=6)
#plot(NSUN1.res$specificities,NSUN1.res$sensitivities,col="red",xlim=c(1,0))
 plot(NSUN1.res$specificities,NSUN1.res$sensitivities,type="l",col="blue",lwd=2,xlim=c(0,1), ylim=c(0,1),   
     ylab="sensitivities")
#abline(0,1,col="gray",lty=2,,lwd=2)
lines(NSUN2.res$specificities,NSUN2.res$sensitivities, type="l",col="tomato",lwd=2,xlim=c(0,1), ylim=c(0,1))
lines(NSUN3.res$specificities,NSUN3.res$sensitivities, type="l",col="red",lwd=2,xlim=c(0,1), ylim=c(0,1))
lines(NSUN4.res$specificities,NSUN4.res$sensitivities, type="l",col="green3",lwd=2,xlim=c(0,1), ylim=c(0,1))
lines(NSUN5.res$specificities,NSUN5.res$sensitivities, type="l",col="cyan",lwd=2,xlim=c(0,1), ylim=c(0,1))
lines(NSUN6.res$specificities,NSUN6.res$sensitivities, type="l",col="magenta",lwd=2,xlim=c(0,1), ylim=c(0,1))
lines(NSUN7.res$specificities,NSUN7.res$sensitivities, type="l",col="yellow",lwd=2,xlim=c(0,1), ylim=c(0,1))
legend("bottomleft",c(paste("NSUN1:",round(NSUN1.res$auc,3),sep=" "),
                 paste("NSUN2:",round(NSUN2.res$auc,3),sep=" "),
                 paste("NSUN3:",round(NSUN3.res$auc,3),sep=" "),
                 paste("NSUN4:",round(NSUN4.res$auc,3),sep=" "),
                 paste("NSUN5:",round(NSUN5.res$auc,3),sep=" "),
                 paste("NSUN6:",round(NSUN6.res$auc,3),sep=" "),
                 paste("NSUN7:",round(NSUN7.res$auc,3),sep=" ")
                ),
                 x.intersp=1, y.intersp=0.8,
                 lty= 1 ,lwd= 2,col=c("blue","tomato","red","green3","cyan","magenta","yellow"),
                 bty = "n",# bty框的类型
                 seg.len=1,cex=0.8)


##########################生存分析
##############NSUN*.level来自Figure2
###################KM生存分析
library(survminer) 
library(survival) 
#2. 拟合生存曲线
NSUN7.fit <- survfit(Surv(OS.time,OS.Status) ~ NSUN7.level,  data = NSUN7.level) #2. 拟合生存曲线

#summary(fit) #查看生存分析结果
#3. 绘制基础曲线
NSUN1.fit <- survfit(Surv(OS.time,OS.Status) ~ NSUN1.level,  data = NSUN1.level) #2. 拟合生存曲线

#summary(fit) #查看生存分析结果
#3. 绘制基础曲线
NSUN1.sur3<-ggsurvplot(NSUN1.fit, # 创建的拟合对象
           data = NSUN1.level,  # 指定变量数据来源
           conf.int = TRUE, # 显示置信区间
           pval = TRUE, # 添加P值
           risk.table = TRUE, # 绘制累计风险曲线
           surv.median.line = "hv", # 添加中位生存时间线
           #add.all = TRUE, # 添加总患者生存曲线
            palette=c("#FF3030","#0000FF")) + xlab("Time (year)")
NSUN1.sur3
#ggsave("Cluster.KM.pdf",sur,width=6,height=6)
pdf(file = "/mnt/nas_002/home/limma/02.project/04.BRCA/01.Analysis/02.KM/NSUN1.KM.pdf",width =6,height = 6)
print(NSUN1.sur3)
dev.off()
NSUN2.fit <- survfit(Surv(OS.time,OS.Status) ~ NSUN2.level,  data = NSUN2.level) #2. 拟合生存曲线

#summary(fit) #查看生存分析结果
#3. 绘制基础曲线
NSUN2.sur3<-ggsurvplot(NSUN2.fit, # 创建的拟合对象
           data = NSUN2.level,  # 指定变量数据来源
           conf.int = TRUE, # 显示置信区间
           pval = TRUE, # 添加P值
           risk.table = TRUE, # 绘制累计风险曲线
           surv.median.line = "hv", # 添加中位生存时间线
           #add.all = TRUE, # 添加总患者生存曲线
            palette=c("#FF3030","#0000FF")) + xlab("Time (year)")
NSUN2.sur3
#ggsave("Cluster.KM.pdf",sur,width=6,height=6)
pdf(file = "/mnt/nas_002/home/limma/02.project/04.BRCA/01.Analysis/02.KM/NSUN2.KM.pdf",width =6,height = 6)
print(NSUN2.sur3)
dev.off()
NSUN3.fit <- survfit(Surv(OS.time,OS.Status) ~ NSUN3.level,  data = NSUN3.level) #2. 拟合生存曲线

NSUN3.sur3<-ggsurvplot(NSUN3.fit, # 创建的拟合对象
           data = NSUN3.level,  # 指定变量数据来源
           conf.int = TRUE, # 显示置信区间
           pval = TRUE, # 添加P值
           risk.table = TRUE, # 绘制累计风险曲线
           surv.median.line = "hv", # 添加中位生存时间线
           #add.all = TRUE, # 添加总患者生存曲线
            palette=c("#FF3030","#0000FF")) + xlab("Time (year)")
NSUN3.sur3
#ggsave("Cluster.KM.pdf",sur,width=6,height=6)
pdf(file = "/mnt/nas_002/home/limma/02.project/04.BRCA/01.Analysis/02.KM/NSUN3.KM.pdf",width =6,height = 6)
print(NSUN3.sur3)
dev.off()
NSUN4.fit <- survfit(Surv(OS.time,OS.Status) ~ NSUN4.level,  data = NSUN4.level) #2. 拟合生存曲线

NSUN4.sur3<-ggsurvplot(NSUN4.fit, # 创建的拟合对象
           data = NSUN4.level,  # 指定变量数据来源
           conf.int = TRUE, # 显示置信区间
           pval = TRUE, # 添加P值
           risk.table = TRUE, # 绘制累计风险曲线
           surv.median.line = "hv", # 添加中位生存时间线
           #add.all = TRUE, # 添加总患者生存曲线
            palette=c("#FF3030","#0000FF")) + xlab("Time (year)")
NSUN4.sur3
#ggsave("Cluster.KM.pdf",sur,width=6,height=6)
pdf(file = "/mnt/nas_002/home/limma/02.project/04.BRCA/01.Analysis/02.KM/NSUN4.KM.pdf",width =6,height = 6)
print(NSUN4.sur3)
dev.off()
NSUN5.fit <- survfit(Surv(OS.time,OS.Status) ~ NSUN5.level,  data = NSUN5.level) #2. 拟合生存曲线

NSUN5.sur3<-ggsurvplot(NSUN5.fit, # 创建的拟合对象
           data = NSUN5.level,  # 指定变量数据来源
           conf.int = TRUE, # 显示置信区间
           pval = TRUE, # 添加P值
           risk.table = TRUE, # 绘制累计风险曲线
           surv.median.line = "hv", # 添加中位生存时间线
           #add.all = TRUE, # 添加总患者生存曲线
            palette=c("#FF3030","#0000FF")) + xlab("Time (year)")
NSUN5.sur3
#ggsave("Cluster.KM.pdf",sur,width=6,height=6)
pdf(file = "/mnt/nas_002/home/limma/02.project/04.BRCA/01.Analysis/02.KM/NSUN5.KM.pdf",width =6,height = 6)
print(NSUN5.sur3)
dev.off()
NSUN6.fit <- survfit(Surv(OS.time,OS.Status) ~ NSUN6.level,  data = NSUN6.level) #2. 拟合生存曲线

#summary(fit) #查看生存分析结果
#3. 绘制基础曲线
NSUN6.sur3<-ggsurvplot(NSUN6.fit, # 创建的拟合对象
           data = NSUN6.level,  # 指定变量数据来源
           conf.int = TRUE, # 显示置信区间
           pval = TRUE, # 添加P值
           risk.table = TRUE, # 绘制累计风险曲线
           surv.median.line = "hv", # 添加中位生存时间线
           #add.all = TRUE, # 添加总患者生存曲线
            palette=c("#FF3030","#0000FF")) + xlab("Time (year)")
NSUN6.sur3
#ggsave("Cluster.KM.pdf",sur,width=6,height=6)
pdf(file = "/mnt/nas_002/home/limma/02.project/04.BRCA/01.Analysis/02.KM/NSUN6.KM.pdf",width =6,height = 6)
print(NSUN6.sur3)
dev.off()

NSUN7.sur3<-ggsurvplot(NSUN7.fit, # 创建的拟合对象
           data = NSUN7.level,  # 指定变量数据来源
           conf.int = TRUE, # 显示置信区间
           pval = TRUE, # 添加P值
           risk.table = TRUE, # 绘制累计风险曲线
           surv.median.line = "hv", # 添加中位生存时间线
           #add.all = TRUE, # 添加总患者生存曲线
            palette=c("#FF3030","#0000FF")) + xlab("Time (year)")
NSUN7.sur3
#ggsave("Cluster.KM.pdf",sur,width=6,height=6)
pdf(file = "/mnt/nas_002/home/limma/02.project/04.BRCA/01.Analysis/02.KM/NSUN7.KM.pdf",width =6,height = 6)
print(NSUN7.sur3)
dev.off()

############################## 单因素和多因素的分析,NSUNS.fpkm.clinc来自Figure1
Sing.mul.level<-NSUNS.fpkm.clinc%>%dplyr::filter(treat %in% c("01"))%>%dplyr::select(-c(treat))%>% dplyr::filter(Stage !="")
Sing.mul.level<-dplyr::select(Sing.mul.level,c("OS.Status","OS.time","Stage","age","TNM_T","TNM_M","TNM_N","NSUN1","NSUN2","NSUN3","NSUN4","NSUN5","NSUN6","NSUN7"))
clin.com<-Sing.mul.level
clin.com[which(clin.com$Stage=="Stage IA"), ]$Stage='Stage I'
clin.com[which(clin.com$Stage=="Stage IB"), ]$Stage='Stage I'
clin.com[which(clin.com$Stage=="Stage IIA"), ]$Stage='Stage II'
clin.com[which(clin.com$Stage=="Stage IIB"), ]$Stage='Stage II'
clin.com[which(clin.com$Stage=="Stage IIIA"), ]$Stage='Stage III'
clin.com[which(clin.com$Stage=="Stage IIIB"), ]$Stage='Stage III'
clin.com[which(clin.com$Stage=="Stage IIIC"), ]$Stage='Stage III'
clin.com[which(clin.com$TNM_T %in% c("T2","T3","T4")), ]$TNM_T='T2-T4'
clin.com[which(clin.com$TNM_N %in% c("N1","N2","N3")), ]$TNM_N='N1-N3'
#clin.com[which(clin.com$Stage %in% c("Stage I","Stage II")), ]$Stage='Stage I && II'
#clin.com[which(clin.com$Stage %in% c("Stage III","Stage IV")), ]$Stage='Stage III &&IV'

clin.com$Stage=as.character(clin.com$Stage)
#clin.com[which(clin.com$Stage %in% c("Stage I","Stage II")), ]$Stage='Stage I && II'
#clin.com[which(clin.com$Stage %in% c("Stage III","Stage IV")), ]$Stage='Stage III &&IV'
clin.com[which(clin.com$TNM_M %in% c("MX")),]$TNM_M<-NA
clin.com[which(clin.com$TNM_N %in% c("NX")),]$TNM_N<-NA
clin.com[which(clin.com$TNM_T %in% c("TX")),]$TNM_T<-NA
clin.com[which(clin.com$Stage %in% c("Stage X","")),]$Stage<-NA
clin.com$Stage <- factor(clin.com$Stage)
clin.com$TNM_T <- factor(clin.com$TNM_T)
clin.com$TNM_M <- factor(clin.com$TNM_M)
clin.com$TNM_N <- factor(clin.com$TNM_N)
clin.com<-droplevels (clin.com)

#######################################3
single_result<-data.frame()
for (a in colnames(clin.com[,4:ncol(clin.com)])){
  formula_string = paste0("Surv(OS.time, OS.Status)~",a)
  formula_expression= as.formula(formula_string)
    formula_expression
  cox <- coxph(formula_expression, data = clin.com)
  coxSummary = summary(cox)
  single_result=rbind(single_result,
                      cbind(Variable=rownames(coxSummary$conf.int),
                            HR=coxSummary$conf.int[,"exp(coef)"],
                            LowerCI=coxSummary$conf.int[,"lower .95"],
                            UpperCI=coxSummary$conf.int[,"upper .95"],
                            pvalue=coxSummary$coefficients[,"Pr(>|z|)"]))
}
single_result[,2:5] <- apply(single_result[,2:5],2,as.numeric)
############################通过单因素
single_result$LowerCI<- round(single_result$LowerCI,2)
single_result$UpperCI <- round(single_result$UpperCI,2)
single_result$HR <- round(single_result$HR,2)
single_result$pvalue <- round(single_result$pvalue,4)
single_result$HR_mean <- paste0(single_result$HR, "(", 
                                single_result$LowerCI, "~", single_result$UpperCI, ")")
single_result$pvalue <- ifelse(single_result$pvalue <0.001,"<0.001",single_result$pvalue)
single_result<-single_result %>% dplyr::select(Variable,HR_mean,pvalue,HR,everything())
#write.table(as.data.frame(single_result),file = "./singleCox.xls",col.names = TRUE,row.names=F,quote=FALSE,sep="\t")
########################多变因素
formula_string = paste0("Surv(OS.time, OS.Status)~",paste0((colnames(clin.com[,4:ncol(clin.com)])),collapse = "+"))

formula_expression= as.formula(formula_string)  
cox <- coxph(formula_expression, data = clin.com)
coxSummary = summary(cox)
coxSummary
multi_result<-data.frame()
multi_result=rbind(multi_result,cbind(Variable=rownames(coxSummary$conf.int),
                            HR=coxSummary$conf.int[,"exp(coef)"],
                            LowerCI=coxSummary$conf.int[,"lower .95"],
                            UpperCI=coxSummary$conf.int[,"upper .95"],
                            pvalue=coxSummary$coefficients[,"Pr(>|z|)"]))
multi_result[,2:5] <- apply(multi_result[,2:5],2,as.numeric)
############################
multi_result$LowerCI<- round(multi_result$LowerCI,2)
multi_result$UpperCI <- round(multi_result$UpperCI,2)
multi_result$HR <- round(multi_result$HR,2)
multi_result$pvalue <- round(multi_result$pvalue,4)
multi_result$HR_mean <- paste0(multi_result$HR, "(", 
                                multi_result$LowerCI, "~", multi_result$UpperCI, ")")
multi_result$pvalue <- ifelse(multi_result$pvalue <0.001,"<0.001",multi_result$pvalue)
multi_result<-multi_result %>% dplyr::select(Variable,HR_mean,pvalue,HR,everything())

##################################################单因素分析结果的可视化
single_result<-dplyr::filter(single_result,Variable !="TNM_TTX")%>%dplyr::filter(Variable !="TNM_MMX")%>%dplyr::filter(Variable !="TNM_NNX")%>%dplyr::filter(Variable !="StageStage X")
labeltext <- as.matrix(single_result[,c(1:3)])   
blacktext <- rep(F,nrow(labeltext))
blacktext[c(1,4,7,12,17,20,23)] <- TRUE
library(RColorBrewer)
Color <- brewer.pal(8, "Set2")                            
library(forestplot) 
singleCox <- 
  forestplot(single_result,labeltext,  # 图形文本部分 
             mean = HR,  # 图形 HR 部分 
             graph.pos=3, #箱线图所在的位置
             #箱线图
             lower = LowerCI, # 95%CI下限           
             upper = UpperCI, # 95%CI上限
             #箱线图中基准线的位置
             zero = 1.0, lwd.zero = 2, # 设置无效线的横坐标和线条宽度
             xlab = "<Hazard ratio>\n UnivariateCoxForest", # 设置x轴标题 
             xlog = TRUE,
             #xticks = c(0.125,0.25,0.5,1,2,4,8,16,32),  # 设置x轴刻度ticks
             xticks.digits = 3,
             clip = c(0.05,100),
             cex=0.9, lineheight = "auto",
             lwd.xaxis = 2, # 设置x轴的宽度
             boxsize = 0.5, # 设置Box的大小，保持大小一样
             lty.ci = 1, lwd.ci = 2, # 森林图置信区间的线条类型和宽度
             ci.vertices = TRUE,  # 森林图置信区间两端添加小竖线，默认FALSE
             ci.vertices.height = 0.2, # 设置森林图置信区间两端小竖线的高度，默认0.1
             align = "l",  # 文字对齐方式，"l"、"r"和"c"分别为左对齐，右对齐和居中对齐
             is.summary= blacktext,# 向量长度等于图表行数，TRUE为行加粗，且该行下添加一条直线，但在未设置颜色时不显示。
             col=fpColors(box=Color[1],lines = 'black',zero = '#7AC5CD'),# 设置坐标轴标题大小
             txt_gp = fpTxtGp(ticks = gpar(cex = 1), xlab = gpar(cex = 1.1)))
singleCox
##################################多因素分析结果的可视化
multi_result<-dplyr::filter(as.data.frame(multi_result),Variable !="TNM_TTX")%>%dplyr::filter(Variable !="TNM_MMX")%>%dplyr::filter(Variable !="TNM_NNX")%>%dplyr::filter(Variable !="StageStage X")
labeltext <- as.matrix(multi_result[,c(1:3)])   
blacktext <- rep(F,nrow(labeltext))
blacktext[c(1,4,7,12,17,20,23)] <- TRUE
library(RColorBrewer)
Color <- brewer.pal(8, "Set2")
library(forestplot) 
multi.Cox <- 
  forestplot(multi_result,labeltext,  # 图形文本部分 
             mean = HR,  # 图形 HR 部分 
             graph.pos=3, #箱线图所在的位置
             #箱线图
             lower = LowerCI, # 95%CI下限           
             upper = UpperCI, # 95%CI上限
             #箱线图中基准线的位置
             zero = 1.0, lwd.zero = 2, # 设置无效线的横坐标和线条宽度
             xlab = "<Hazard ratio>\n MultivariateCoxForest", # 设置x轴标题 
             xlog = TRUE,
             #xticks = c(0.125,0.25,0.5,1,2,4,8,16,32),  # 设置x轴刻度ticks
             xticks.digits = 3,
             clip = c(0.05,100),
             cex=0.9, lineheight = "auto",
             lwd.xaxis = 2, # 设置x轴的宽度
             boxsize = 0.5, # 设置Box的大小，保持大小一样
             lty.ci = 1, lwd.ci = 2, # 森林图置信区间的线条类型和宽度
             ci.vertices = TRUE,  # 森林图置信区间两端添加小竖线，默认FALSE
             ci.vertices.height = 0.2, # 设置森林图置信区间两端小竖线的高度，默认0.1
             align = "l",  # 文字对齐方式，"l"、"r"和"c"分别为左对齐，右对齐和居中对齐
             is.summary= blacktext,# 向量长度等于图表行数，TRUE为行加粗，且该行下添加一条直线，但在未设置颜色时不显示。
             col=fpColors(box=Color[1],lines = 'black',zero = '#7AC5CD'),# 设置坐标轴标题大小
             txt_gp = fpTxtGp(ticks = gpar(cex = 1), xlab = gpar(cex = 1.1)))
multi.Cox
pdf("/mnt/nas_002/home/limma/02.project/04.BRCA/01.Analysis/04.Fig/Fig5j.multi.Cox.forest.pdf")
print(multi.Cox)
dev.off()
