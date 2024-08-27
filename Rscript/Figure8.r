##############################
ROC<- read.table("/mnt/nas_002/home/limma/02.project/04.BRCA/01.Analysis/00.jianmo/All.Auc.xls",sep = "\t",header=TRUE)
rownames(ROC)=ROC$Algorithm
colnames(ROC)=c("TCGA","GSE20685","Metabric","Ave","Algorithm")
ROC<-dplyr::select(ROC,-Algorithm)
################3
annotation_col = data.frame(Evaluate = factor(colnames(ROC)))
rownames(annotation_col)=colnames(ROC)
ann_colors = list(Evaluate = c(TCGA= "red", GSE20685 = "red", Metabric= "red",Ave="Purple"))
#ann_colors = list(Evaluate = c(R_TCGA = "red", Test1.TCGA= "tomato",Test2.GSE39582="Orchid",Test3.GSE17537="Purple"))
###########

bk1 <- c(seq(1,0.85,by=-0.01))
bk2 <- c(seq(0.85,0.65,by=-0.01))
bk3 <- c(seq(0.65,0.4,by=-0.01))
bk<-c(bk1,bk2,bk3)
#clour = c(colorRampPalette(colors = c("white","DarkCyan"))(length(bk1)),colorRampPalette(colors = c("white","white"))(length(bk2)),colorRampPalette(colors = c("white","tomato"))(length(bk2)),colorRampPalette(colors = c("white","white"))(length(bk4)))
clour = c(colorRampPalette(colors = c("white","white"))(length(bk3)),colorRampPalette(colors = c("white","tomato"))(length(bk2)),colorRampPalette(colors = c("white","DarkCyan"))(length(bk1)))
#clour = c(colorRampPalette(colors = c("white","white"))(length(bk2)),colorRampPalette(colors = c("white","tomato"))(length(bk1)),colorRampPalette(colors = c("white","white"))(length(bk4)),colorRampPalette(colors = c("white","DarkCyan"))(length(bk3)))
library(pheatmap)
options(repr.plot.height=15, repr.plot.width=10)
p<-pheatmap(ROC,border=FALSE,cluster_rows = F,cluster_cols = F,cellwidth = 30,
           # cellheight= 12,show_rownames = T,
            display_numbers = TRUE, #设定在每个热图格子中显示相应的数值
            number_format = "%.2f", #设定数值的显示格式
            number_color = "Black", #设置热图格子中数值的颜色
            gaps_col = c(3),
            annotation_col = annotation_col,annotation_colors = ann_colors,color=clour) 
            #color =colorRampPalette(c("blue","white","red"))(100))

p
pdf("/mnt/nas_002/home/limma/02.project/04.BRCA/01.Analysis/00.jianmo/Figure/ROC.index.heatmap.pdf",height = 15,width=10)
print(p)
dev.off()
#################################

coef <-read.table("/mnt/nas_002/home/limma/02.project/04.BRCA/01.Analysis/00.jianmo/LASSO_COX_lasso.coef.xls",sep = "\t",header=TRUE)
colnames(coef)=c("coef")
coef$geneid=rownames(coef)
library(ggplot2) 
p = ggplot(coef,aes(coef,geneid))
p= p+ geom_point(aes(size=abs(coef),color=coef))
p= p+ geom_vline(xintercept = 0)+
 #scale_color_gradient(low="yellow",high = "red")+
  geom_segment(aes(x=0, y=geneid, xend=coef, yend=geneid,color=coef))+scale_color_gradient(low="blue",high = "red")+
  theme(text=element_text(size=16,face="plain",color="black"),
        axis.title=element_text(size=16,face="plain",color="black"),
        axis.text = element_text(size=10,face="plain",color="black"),
        legend.title = element_text(size=16,face="plain",color="black"),
        legend.text = element_text(size=10,face="plain",color="black"),
        legend.background = element_blank(),
        legend.position="bottom")
p = p+labs(size="coef",x="coefficient",y="Gene",title="")+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
options(repr.plot.height=5, repr.plot.width=5)
p
pdf("/mnt/nas_002/home/limma/02.project/04.BRCA/01.Analysis/00.jianmo/Figure/Coefficen.pdf")
print(p)
dev.off()

###########################1.多因素风险分析
clin.com<-TCGA.BRCA$TCGA.clinc%>%dplyr::select(c("ID","gender","Stage","age","TNM_T","TNM_M","TNM_N"))
clin.com<-merge(Risk.score,clin.com,by="ID")%>% dplyr::select(-c("level","Rank","gender","TNM_M"))
clin.com[which(clin.com$TNM_T %in% c("TX","Ti")),]$TNM_T<-NA
clin.com[which(clin.com$TNM_N %in% c("NX")),]$TNM_N<-NA
#clin.com[which(clin.com$TNM_M %in% c("MX")),]$TNM_M<-NA
clin.com[which(clin.com$Stage=="Stage IA"), ]$Stage='Stage I'
clin.com[which(clin.com$Stage=="Stage IB"), ]$Stage='Stage I'
clin.com[which(clin.com$Stage=="Stage IIA"), ]$Stage='Stage II'
clin.com[which(clin.com$Stage=="Stage IIB"), ]$Stage='Stage II'
clin.com[which(clin.com$Stage=="Stage IIIA"), ]$Stage='Stage III'
clin.com[which(clin.com$Stage=="Stage IIIB"), ]$Stage='Stage III'
clin.com[which(clin.com$Stage=="Stage IIIC"), ]$Stage='Stage III'
clin.com[which(clin.com$Stage %in% c("Stage X","")),]$Stage<-NA
clin.com<-droplevels(clin.com)
table(clin.com$Stage)
########################多变因素
formula_string = paste0("Surv(OS.time, OS.Status)~",paste0((colnames(clin.com[,4:ncol(clin.com)])),collapse = "+"))

formula_expression= as.formula(formula_string)  
cox <- coxph(formula_expression, data = clin.com)
coxSummary = summary(cox)
coxSummary
multi_result=cbind(Variable=rownames(coxSummary$conf.int),
                            HR=coxSummary$conf.int[,"exp(coef)"],
                            LowerCI=coxSummary$conf.int[,"lower .95"],
                            UpperCI=coxSummary$conf.int[,"upper .95"],
                            pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
multi_result[,2:5] <- apply(multi_result[,2:5],2,as.numeric)
####
multi_result<-as.data.frame(multi_result)
multi_result$LowerCI<- round(as.numeric(multi_result$LowerCI),2)
multi_result$UpperCI <- round(as.numeric(multi_result$UpperCI),2)
multi_result$HR <- round(as.numeric(multi_result$HR),2)
multi_result$pvalue <- round(as.numeric(multi_result$pvalue),4)
multi_result$HR_mean <- paste0(multi_result$HR, "(", 
                                multi_result$LowerCI, "~", multi_result$UpperCI, ")")
#multi_result$pvalue <- ifelse(multi_result$pvalue <0.001,"<0.001",multi_result$pvalue)
multi_result<-multi_result %>% dplyr::select(Variable,HR_mean,pvalue,HR,everything())
write.table(as.data.frame(multi_result),file = "/mnt/nas_002/home/limma/02.project/04.BRCA/01.Analysis/00.jianmo/Multi.Cox.xls",col.names = TRUE,row.names=F,quote=FALSE,sep="\t")
multi_result

labeltext <- as.matrix(multi_result[,c(1:3)])   
#blacktext <- rep(F,nrow(labeltext))
#blacktext[c(1,4,7,12,17,20,23)] <- TRUE
library(RColorBrewer)
#Color <- brewer.pal(8, "Set2")
Color <- c("skyblue")   
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
             is.summary= FALSE,# 向量长度等于图表行数，TRUE为行加粗，且该行下添加一条直线，但在未设置颜色时不显示。
             col=fpColors(box='#cc340c',lines = "black",zero = 'grey'),# 设置坐标轴标题大小
             txt_gp = fpTxtGp(ticks = gpar(cex = 1), xlab = gpar(cex = 1.1)))
multi.Cox
pdf("/mnt/nas_002/home/limma/02.project/04.BRCA/01.Analysis/00.jianmo/Fig5j.multi.Cox.forest.pdf")
print(multi.Cox)
dev.off()



###############################2.KM生存分析

library(survminer)
library(survival)
library(cowplot)
# 根据risk score 分成 高/低 风险2组样本，以中位数分割
Risk.score<-read.table("/mnt/nas_002/home/limma/02.project/04.BRCA/01.Analysis/00.jianmo/LASSO_COXTCGA.Auc_data.xls",header=TRUE,sep="\t")
Risk.score<- dplyr::mutate(Risk.score,level=ifelse(predict > median(predict), "high", "low")) %>% dplyr::arrange(predict)
Rank = 1:length(Risk.score$level)
Risk.score <- data.frame(Risk.score,Rank)
fit <- survfit(Surv(OS.time,OS.Status==1) ~ level, data =  Risk.score) #2. 拟合生存曲线
TCGA.sur3<-ggsurvplot(fit, # 创建的拟合对象
                 data =Risk.score,  # 指定变量数据来源
                 conf.int = TRUE, # 显示置信区间
                 pval = TRUE, # 添加P值
                 risk.table = TRUE, # 绘制累计风险曲线
                 surv.median.line = "hv", # 添加中位生存时间线
                 #add.all = TRUE, # 添加总患者生存曲线
                 palette=c("#156077","#f46f20")) + xlab("Time (year)")
TCGA.sur3
pdf("/mnt/nas_002/home/limma/02.project/04.BRCA/01.Analysis/00.jianmo/Figure/TCGA.sur3.Risk.Survplot.pdf")
print(TCGA.sur3)
dev.off()

# 根据risk score 分成 高/低 风险2组样本，以中位数分割
GSE20685.Risk.score<-read.table("/mnt/nas_002/home/limma/02.project/04.BRCA/01.Analysis/00.jianmo/LASSO_COXGSE20685.Auc_data.xls",header=TRUE,sep="\t")
GSE20685.Risk.score<- dplyr::mutate(GSE20685.Risk.score,level=ifelse(predict > median(predict), "high", "low")) %>% dplyr::arrange(predict)
Rank = 1:length(GSE20685.Risk.score$level)
GSE20685.Risk.score <- data.frame(GSE20685.Risk.score,Rank)
fit <- survfit(Surv(OS.time,OS.Status==1) ~ level, data =  GSE20685.Risk.score) #2. 拟合生存曲线
GSE20685.sur3<-ggsurvplot(fit, # 创建的拟合对象
                 data =GSE20685.Risk.score,  # 指定变量数据来源
                 conf.int = TRUE, # 显示置信区间
                 pval = TRUE, # 添加P值
                 risk.table = TRUE, # 绘制累计风险曲线
                 surv.median.line = "hv", # 添加中位生存时间线
                 #add.all = TRUE, # 添加总患者生存曲线
                 palette=c("#156077","#f46f20")) + xlab("Time (year)")
GSE20685.sur3
pdf("/mnt/nas_002/home/limma/02.project/04.BRCA/01.Analysis/00.jianmo/Figure//GSE20685.sur3.Risk.Survplot.pdf")
print(GSE20685.sur3)
dev.off()

Metabric.Risk.score<-read.table("/mnt/nas_002/home/limma/02.project/04.BRCA/01.Analysis/00.jianmo/LASSO_COXMetabric.Auc_data.xls",header=TRUE,sep="\t")
Metabric.Risk.score<- dplyr::mutate(Metabric.Risk.score,level=ifelse(predict > median(predict), "high", "low")) %>% dplyr::arrange(predict)
Rank = 1:length(Metabric.Risk.score$level)
Metabric.Risk.score <- data.frame(Metabric.Risk.score,Rank)
fit <- survfit(Surv(OS.time,OS.Status==1) ~ level, data =  Metabric.Risk.score) #2. 拟合生存曲线
Metabric.sur3<-ggsurvplot(fit, # 创建的拟合对象
                 data =Metabric.Risk.score,  # 指定变量数据来源
                 conf.int = TRUE, # 显示置信区间
                 pval = TRUE, # 添加P值
                 risk.table = TRUE, # 绘制累计风险曲线
                 surv.median.line = "hv", # 添加中位生存时间线
                 #add.all = TRUE, # 添加总患者生存曲线
                 palette=c("#156077","#f46f20")) + xlab("Time (year)")
Metabric.sur3
pdf("/mnt/nas_002/home/limma/02.project/04.BRCA/01.Analysis/00.jianmo/Figure//Metabric.sur3.Risk.Survplot.pdf")
print(Metabric.sur3)
dev.off()
######################################3.ROC 曲线分析
options(repr.plot.height=8, repr.plot.width=8)
library(survivalROC)
cutoff <- 1
ROC1= survivalROC(Stime=Risk.score$OS.time,##生存时间
                     status=Risk.score$OS.Status,## 终止事件    
                     marker = Risk.score$predict, ## marker value    
                     #other_markers = as.matrix(feng_36[,c("LAG3","IDO1")]), # 矫正的协变量
                     predict.time = 1,## 预测时间截点
                     method="KM"
                     )##span,NNE法的namda
#再加一个3年的线
ROC2=survivalROC(Stime=Risk.score$OS.time,##生存时间
                     status=Risk.score$OS.Status,## 终止事件    
                     marker = Risk.score$predict, ## marker value    
                     predict.time = 2,## 预测时间截点
                     method="KM")
#再加一个5年的线
ROC3=survivalROC(Stime=Risk.score$OS.time,##生存时间
                     status=Risk.score$OS.Status,## 终止事件    
                     marker = Risk.score$predict, ## marker value    
                     predict.time = 3,## 预测时间截点
                     method="KM")
plotROC.1<- function() {plot(ROC1$FP, ROC1$TP, 
     type="l",col="red",lwd=2,xlim=c(0,1), ylim=c(0,1),   
    # xlab=paste( "FP", "\n", "AUC = ",round(ROC$AUC,3)), 
     ylab="TP")
abline(0,1)
aucText=paste0("1 years"," (AUC=",sprintf("%.3f",ROC1$AUC),")")
lines(ROC2$FP, ROC2$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="cyan4",lwd = 2)
aucText2=paste0("2 years"," (AUC=",sprintf("%.3f",ROC2$AUC),")") #这个后面添加legend用

lines(ROC3$FP, ROC3$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="blue",lwd = 2)
aucText3=paste0("3 years"," (AUC=",sprintf("%.3f",ROC3$AUC),")") #这个后面添加legend用
#添加legend
legend("bottomright", c(aucText,aucText2,aucText3),
       lwd=2,bty="n",col=c("red","cyan4","blue"),cex=1.2)
                        }
plotROC.1()
pdf("/mnt/nas_002/home/limma/02.project/04.BRCA/01.Analysis/00.jianmo/Figure/TCGA.ROC.pdf")
plotROC.1()
dev.off()

######################################ROC 曲线分析

GSE20685.ROC1= survivalROC(Stime=GSE20685.Risk.score$OS.time,##生存时间
                     status=GSE20685.Risk.score$OS.Status,## 终止事件    
                     marker = GSE20685.Risk.score$predict, ## marker value    
                     #other_markers = as.matrix(feng_36[,c("LAG3","IDO1")]), # 矫正的协变量
                     predict.time = 1,## 预测时间截点
                     method="KM"
                     )##span,NNE法的namda
#再加一个3年的线
GSE20685.ROC2=survivalROC(Stime=GSE20685.Risk.score$OS.time,##生存时间
                     status=GSE20685.Risk.score$OS.Status,## 终止事件    
                     marker = GSE20685.Risk.score$predict, ## marker value    
                     predict.time = 2,## 预测时间截点
                     method="KM")
#再加一个5年的线
GSE20685.ROC3=survivalROC(Stime=GSE20685.Risk.score$OS.time,##生存时间
                     status=GSE20685.Risk.score$OS.Status,## 终止事件    
                     marker = GSE20685.Risk.score$predict, ## marker value    
                     predict.time = 3,## 预测时间截点
                     method="KM")
plotROC.1<- function() {plot(GSE20685.ROC1$FP, GSE20685.ROC1$TP, 
     type="l",col="red",lwd=2,xlim=c(0,1), ylim=c(0,1),   
    # xlab=paste( "FP", "\n", "AUC = ",round(GSE20685.ROC$AUC,3)), 
     ylab="TP")
abline(0,1)
aucText=paste0("1 years"," (AUC=",sprintf("%.3f",GSE20685.ROC1$AUC),")")

lines(GSE20685.ROC2$FP, GSE20685.ROC2$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="cyan4",lwd = 2)
aucText2=paste0("2 years"," (AUC=",sprintf("%.3f",GSE20685.ROC2$AUC),")") #这个后面添加legend用

lines(GSE20685.ROC3$FP, GSE20685.ROC3$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="blue",lwd = 2)
aucText3=paste0("3 years"," (AUC=",sprintf("%.3f",GSE20685.ROC3$AUC),")") #这个后面添加legend用
#添加legend
legend("bottomright", c(aucText,aucText2,aucText3),
       lwd=2,bty="n",col=c("red","cyan4","blue"),cex=1.2)
                        }
plotROC.1()
pdf("//mnt/nas_002/home/limma/02.project/04.BRCA/01.Analysis/00.jianmo/Figure/GSE20685.ROC.pdf")
plotROC.1()
dev.off()

######################################ROC 曲线分析

Metabric.ROC1= survivalROC(Stime=Metabric.Risk.score$OS.time,##生存时间
                     status=Metabric.Risk.score$OS.Status,## 终止事件    
                     marker = Metabric.Risk.score$predict, ## marker value    
                     #other_markers = as.matrix(feng_36[,c("LAG3","IDO1")]), # 矫正的协变量
                     predict.time = 1,## 预测时间截点
                     method="KM"
                     )##span,NNE法的namda
#再加一个3年的线
Metabric.ROC2=survivalROC(Stime=Metabric.Risk.score$OS.time,##生存时间
                     status=Metabric.Risk.score$OS.Status,## 终止事件    
                     marker = Metabric.Risk.score$predict, ## marker value    
                     predict.time = 2,## 预测时间截点
                     method="KM")
#再加一个5年的线
Metabric.ROC3=survivalROC(Stime=Metabric.Risk.score$OS.time,##生存时间
                     status=Metabric.Risk.score$OS.Status,## 终止事件    
                     marker = Metabric.Risk.score$predict, ## marker value    
                     predict.time = 3,## 预测时间截点
                     method="KM")
plotROC.1<- function() {plot(Metabric.ROC1$FP, Metabric.ROC1$TP, 
     type="l",col="red",lwd=2,xlim=c(0,1), ylim=c(0,1),   
    # xlab=paste( "FP", "\n", "AUC = ",round(Metabric.ROC$AUC,3)), 
     ylab="TP")
abline(0,1)
aucText=paste0("1 years"," (AUC=",sprintf("%.3f",Metabric.ROC1$AUC),")")

lines(Metabric.ROC2$FP, Metabric.ROC2$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="cyan4",lwd = 2)
aucText2=paste0("2 years"," (AUC=",sprintf("%.3f",Metabric.ROC2$AUC),")") #这个后面添加legend用

lines(Metabric.ROC3$FP, Metabric.ROC3$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="blue",lwd = 2)
aucText3=paste0("3 years"," (AUC=",sprintf("%.3f",Metabric.ROC3$AUC),")") #这个后面添加legend用
#添加legend
legend("bottomright", c(aucText,aucText2,aucText3),
       lwd=2,bty="n",col=c("red","cyan4","blue"),cex=1.2)
                        }
plotROC.1()
pdf("//mnt/nas_002/home/limma/02.project/04.BRCA/01.Analysis/00.jianmo/Figure/Metabric.ROC.pdf")
plotROC.1()
dev.off()

################################ 免疫微环境分析
my_comparisons <- list(c("high","low"))
estimate<-read.table("/mnt/nas_002/home/limma/02.project/04.BRCA/00.data/00.TCGA/breast_cancer_RNAseqV2.txt",header = TRUE, sep="\t")
estimate$ID=substr(estimate$ID,1,15)
Risk.score.new<-Risk.score
estimate$ID=substr(estimate$ID,1,12)
estimate<-merge(estimate,Risk.score.new,by="ID")%>% dplyr::select(-c("OS.Status","OS.time","predict","Rank"))

estimate1<-estimate%>% pivot_longer(cols=-c("ID","level"), names_to = "type",values_to = "Score")
esti_p<- ggplot(data=estimate1,aes(x=type,y=Score,fill=level))+ geom_boxplot()+
     theme_bw()+
# stat_summary(mapping=aes(x=type,y=Score,color=group),fun="mean",geom="point",shape=21,size=3, fill="grey",position=position_dodge(width=0.8))+
    stat_compare_means(aes(group=level),label="p.signif", method = "wilcox.test"
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
    scale_fill_manual(values = c("#156077","#f46f20"))
  #+labs(title="wilcox.test")
options(repr.plot.height=6, repr.plot.width=6)
esti_p
ggsave("//mnt/nas_002/home/limma/02.project/04.BRCA/01.Analysis/00.jianmo/Figure/ESTI.pdf",esti_p,height=6,width=6)

##############################
####################免疫评分
infil<-read.table("/mnt/nas_002/home/limma/02.project/02.COAD/infiltration_estimation_for_tcga.csv",header = TRUE, sep=",",row.names = 1)
infil<-dplyr::select(infil,ends_with("CIBERSORT"))
infil$treat<-substr(row.names(infil),14,15)
infil$ID<-substr(row.names(infil),1,12)
infil.group<-Risk.score.new

Tide_infil<-merge(infil,infil.group,by="ID")%>% 
  pivot_longer(cols=matches("CIBERSORT"), names_to = "Celltype",values_to = "Fraction")
Tide_infil$Celltype<-gsub("_CIBERSORT", "",Tide_infil$Celltype)



tme_p<- ggplot(data=Tide_infil,aes(x=Celltype,y=Fraction,fill=level))+ geom_boxplot()+
     theme_bw()+
 #stat_summary(mapping=aes(x=Celltype,y=Fraction,color=level),fun="mean",geom="point",shape=21,size=3, fill="grey",position=position_dodge(width=0.8))+
    stat_compare_means(aes(group=level),label="p.signif", method = "wilcox.test"
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
    scale_fill_manual(values = c("#156077","#f46f20"))
options(repr.plot.height=8, repr.plot.width=12)
tme_p

ggsave("//mnt/nas_002/home/limma/02.project/04.BRCA/01.Analysis/00.jianmo/Figure/TME.pdf",tme_p,height=8,width=12)

#################免疫检查点分析
check.gene<-gene[c("CTLA4","HAVCR2","IL2RA","LAG3","PDCD1","TIGIT"),]
check.gene<-as.data.frame(t(check.gene))
check.gene$ID<-substr(rownames(check.gene),1,12)
check.gene<-merge(check.gene,Risk.score.new,by="ID")%>%dplyr::select(-c("OS.Status","OS.time","predict","Rank"))
check.gene<-check.gene %>% pivot_longer(cols=-c("ID","level"), names_to = "Gene",values_to = "Fpkm")
head(check.gene)

p.check.gene<-ggplot(data=check.gene,aes(x=Gene,y=log2(Fpkm+0.1),fill=level))+ geom_boxplot()+theme_bw()+stat_compare_means(aes(group=level),label="p.signif", method = "wilcox.test" )

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
  scale_fill_manual(values = c("#156077","#f46f20"))
#+labs(title="wilcox.test")

options(repr.plot.height=5, repr.plot.width=8)
p.check.gene
ggsave("//mnt/nas_002/home/limma/02.project/04.BRCA/01.Analysis/00.jianmo/Figure/check.gene.pdf",p.check.gene,height=5,width=8)

####################hub基因分析
hub.gene.data<-gene[c("NFKBIA","PTEN","TPI1","GAPDH","NSUN2"),]
hub.gene.data<-as.data.frame(t(hub.gene.data))
hub.gene.data$ID<-substr(rownames(hub.gene.data),1,12)
hub.gene.data<-merge(hub.gene.data,Risk.score.new,by="ID")%>%dplyr::select(-c("OS.Status","OS.time","Rank"))
hub.gene<-hub.gene.data %>%dplyr::select(-c( "predict"))%>%pivot_longer(cols=-c("ID","level"), names_to = "Gene",values_to = "Fpkm")
head(hub.gene)

p.hub.gene<-ggplot(data=hub.gene,aes(x=Gene,y=log2(Fpkm+0.1),fill=level))+ geom_boxplot()+theme_bw()+stat_compare_means(aes(group=level),label="p.signif", method = "wilcox.test" )

p.hub.gene<-p.hub.gene+ theme_bw()+
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
  scale_fill_manual(values = c("#156077","#f46f20"))
#+labs(title="wilcox.test")

options(repr.plot.height=5, repr.plot.width=8)
p.hub.gene
ggsave("//mnt/nas_002/home/limma/02.project/04.BRCA/01.Analysis/00.jianmo/Figure/hub.gene.pdf",p.hub.gene,height=5,width=8)



########################hub基因与风险评分的相关性分析
or<-cor.test(as.numeric(hub.gene.data$predict),as.numeric(hub.gene.data$NFKBIA),method = "pearson")
cor
NFKBIA.data<-hub.gene.data[,c("NFKBIA","predict")]
NFKBIA.data$cor=paste0("NFKBIA",signif(cor$p.value,3),sep="")
colnames(NFKBIA.data)<-c("gene","predict","cor")

cor<-cor.test(as.numeric(hub.gene.data$predict),as.numeric(hub.gene.data$PTEN),method = "pearson")
cor
PTEN.data<-hub.gene.data[,c("PTEN","predict")]
PTEN.data$cor=paste0("PTEN",signif(cor$p.value,3),sep="")
colnames(PTEN.data)<-c("gene","predict","cor")
cor<-cor.test(as.numeric(hub.gene.data$predict),as.numeric(hub.gene.data$TPI1),method = "pearson")
cor
TPI1.data<-hub.gene.data[,c("TPI1","predict")]
TPI1.data$cor=paste0("TPI1",signif(cor$p.value,3),sep="")
colnames(TPI1.data)<-c("gene","predict","cor")
cor<-cor.test(as.numeric(hub.gene.data$predict),as.numeric(hub.gene.data$GAPDH),method = "pearson")
cor
GAPDH.data<-hub.gene.data[,c("GAPDH","predict")]
GAPDH.data$cor=paste0("GAPDH",signif(cor$p.value,3),sep="")
colnames(GAPDH.data)<-c("gene","predict","cor")

cor<-cor.test(as.numeric(hub.gene.data$predict),as.numeric(hub.gene.data$NSUN2),method = "pearson")
cor
NSUN2.data<-hub.gene.data[,c("NSUN2","predict")]
NSUN2.data$cor=paste0("NSUN2",signif(cor$p.value,3),sep="")
colnames(NSUN2.data)<-c("gene","predict","cor")
cor.data<-rbind(NFKBIA.data,PTEN.data,TPI1.data,GAPDH.data,NSUN2.data)
head(cor.data)

TCGA.Risk.score.cor<-ggplot(cor.data,aes(cor.data$predict,cor.data$gene))+geom_point(size=3,colour = 'tomato3')+geom_smooth(method="lm")+#geom_vline(xintercept = 40)+
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
    )+labs(x ="Risk.score" ,y="expression")+facet_wrap(.~cor,ncol=5,scales = "free_y")
options(repr.plot.height=4, repr.plot.width=16)
TCGA.Risk.score.cor
