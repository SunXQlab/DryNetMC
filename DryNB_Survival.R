source("http://bioc.ism.ac.jp/biocLite.R")
biocLite("AnnotationDbi")
install.packages("AnnotationDbi")
install.packages("bit")
biocLite("org.Hs.eg.db")
library(bit)
library(AnnotationDbi)
library(org.Hs.eg.db)
###### Survival package
install.packages("lattice")
library(lattice)
install.packages("survival")
library(survival)

######## plot KM curves
install.packages("reshape2")
library(reshape2)
install.packages("data.table")
library(data.table)
install.packages("zoo")
library("zoo")
install.packages("survminer")
library("survminer")
install.packages("survival")
library("survival")

#### ROC package
# install.packages("gplots")
# library(gplots)
# install.packages("ROCR")
# library(ROCR)
# install.packages("dplyr")
# library("dplyr")
# install.packages("lubridate")
# library(lubridate)
# install.packages("robustbase")
# install.packages("caret")
# library(caret)
install.packages("pROC")
library(pROC)



install.packages("glmnet")
library(glmnet)  # manually load this package by downloading and installing from Toll. 


setwd("Path/Survival analysis")

#### TCGA TargetedTherapy
Gene_LGG=read.csv("Data/TCGA_LGG_GeneExpression.csv")
rownames(Gene_LGG)=Gene_LGG[,1]
colnames(Gene_LGG)=gsub(".", "-", colnames(Gene_LGG), fixed = TRUE)
Gene_LGG=Gene_LGG[,-1]

Survival_LGG=read.csv("Data/TCGA_LGG_Survival data.csv")
Survival_LGG=Survival_LGG[,c(1,2,4,93)]
rownames(Survival_LGG)=Survival_LGG[,1]
colnames(Survival_LGG)=c('SampleID','OS_event','OS_time','Target_Therapy')

#####
Gene_GBM=read.csv("Data/TCGA_GBM_GeneExpression.csv")
rownames(Gene_GBM)=Gene_GBM[,1]
colnames(Gene_GBM)=gsub(".", "-", colnames(Gene_GBM), fixed = TRUE)
Gene_GBM=Gene_GBM[,-1]

Survival_GBM=read.csv("Data/TCGA_GBM_Survival data.csv")
Survival_GBM=Survival_GBM[,c(1,33,32,109)]
rownames(Survival_GBM)=Survival_GBM[,1]
colnames(Survival_GBM)=c('SampleID','OS_event','OS_time','Target_Therapy')

Gene_TCGA=cbind(Gene_GBM[intersect(rownames(Gene_GBM),rownames(Gene_LGG)),],Gene_LGG[intersect(rownames(Gene_GBM),rownames(Gene_LGG)),])
colnames(Gene_TCGA)=gsub(".", "-", colnames(Gene_TCGA), fixed = TRUE)
Survival_TCGA=rbind(Survival_GBM,Survival_LGG)

Survival_TCGA=Survival_TCGA[Survival_TCGA[,4]=='YES',]
dim(Survival_TCGA)


DNG=read.csv("Data/Quantifing each node.csv")  # Differential network genes
Score=as.matrix(DNG[,c(1,9)])  # Importance Score
rownames(Score)=Score[,1]
Score=Score[,-1]

Gene=names(Score)

#####################
#####################
### 100 computations
#####################
#####################

Time=100
sample_ratio=(5:8)*0.1
# AUC_Training_LASSO=matrix(0,length(sample_ratio),Time)
# AUC_Test_LASSO=matrix(0,length(sample_ratio),Time)
# # perf_Training_LASSO=data.frame(Time)
# # perf_Test_LASSO=data.frame(Time)

AUC_Training_wLASSO=matrix(0,length(sample_ratio),Time)
AUC_Test_wLASSO=matrix(0,length(sample_ratio),Time)

sr_ind=0  # index of sample ratio

for (p in sample_ratio)
{
  sr_ind=sr_ind+1
  
  for (i in 1:Time)
  {
    #########  Training set - TCGA dataset 1
    ind <- sample(2, nrow(Survival_TCGA), replace = TRUE, prob = c(p, 1-p))
    
    Survival_Gene=Gene_TCGA[names(Score),which(ind==1)]
    
    # Survival_Gene_w=Gene_TCGA[,which(ind==1)]
    
    colnames(Survival_Gene)=gsub(".", "-", colnames(Survival_Gene), fixed = TRUE)
    
    
    Survival=Survival_TCGA[intersect(rownames(Survival_TCGA),colnames(Survival_Gene)),]
    Survival=as.matrix(Survival)
    Survival_Gene=Survival_Gene[,intersect(rownames(Survival),colnames(Survival_Gene))]
    
    status=as.numeric(Survival[,2])
    time=as.numeric(Survival[,3])
    Genes=rownames(na.omit(Survival_Gene))
    Gene_marker=matrix(as.double(as.matrix(Survival_Gene[Genes,])),dim(Survival_Gene))
    
    time=time[!is.na(status)]
    Gene_marker=Gene_marker[,!is.na(status)]
    status=status[!is.na(status)]
    
    status[status==1&time>365*3]=0  # 3-year survival status
    
    S=as.numeric(Score)
    Gene_marker_w=Gene_marker#*(S-min(S))/max(S)  # Importance score-weighted expression
    
    
    ##Dynamic Network-based LASSO
    S=as.numeric(Score)
    Gene_marker_w=Gene_marker  # Importance score-weighted expression
    cvfit_wLASSO<-cv.glmnet(t(Gene_marker_w),as.double(status),family="binomial",alpha=1,penalty.factor=S/max(S))
  
    Prediction = predict(cvfit_wLASSO,newx = t(Gene_marker_w),s="lambda.min",type="response")
    pred <- prediction(Prediction,status)
    # perf_Training_wLASSO[i] <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
    roc0=roc(status,as.numeric(Prediction))
    AUC_Training_wLASSO[sr_ind,i]=as.numeric(auc(status,as.numeric(Prediction)))
    
    ######### Validation set - TCGA dataset 2
    Survival_Gene=Gene_TCGA[names(Score),]
    Survival_Gene=Survival_Gene[,which(ind==2)]
    
    colnames(Survival_Gene)=gsub(".", "-", colnames(Survival_Gene), fixed = TRUE)
    
    
    Survival=Survival_TCGA[intersect(rownames(Survival_TCGA),colnames(Survival_Gene)),]
    Survival=as.matrix(Survival)
    Survival_Gene=Survival_Gene[,intersect(rownames(Survival),colnames(Survival_Gene))]
    
    status=as.numeric(Survival[,2])
    time=as.numeric(Survival[,3])
    Genes=rownames(na.omit(Survival_Gene))
    Gene_marker=matrix(as.double(as.matrix(Survival_Gene[Genes,])),dim(Survival_Gene))
    
    time=time[!is.na(status)]
    Gene_marker=Gene_marker[,!is.na(status)]
    status=status[!is.na(status)]
    
    status[status==1&time>365*3]=0  # 3-year survival status
    
   
    ### Validation for Dynamic Network-based LASSO
    S=as.numeric(Score)
    Gene_marker_w=Gene_marker#*(S-min(S))/max(S)  # Importance score-weighted expression
    
    Prediction = predict(cvfit_wLASSO,newx = t(Gene_marker_w),s="lambda.min",type="response")
    pred <- prediction(Prediction,status)
    # perf_Test_wLASSO[i] <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
    roc1=roc(status,as.numeric(Prediction))
    AUC_Test_wLASSO[sr_ind,i]=as.numeric(auc(status,as.numeric(Prediction)))
    
  }
}


# AUC_Training_wLASSO

AUC_Test_wLASSO

mean(AUC_Test_wLASSO)

#### Barplots of AUCs_Test dataset of Dynamic Network-based LASSO
dev.new()
boxplot(t(AUC_Training_wLASSO),ylim=c(0,1))


dev.new()
boxplot(t(AUC_Test_wLASSO),ylim=c(0,1))


############### Differential expression profiles of 7 genes in the Normal and Tumor tissues

setwd("Path/Survival analysis") 

Gene_LGG=read.csv("Data/TCGA_LGG_GeneExpression.csv")
rownames(Gene_LGG)=Gene_LGG[,1]
colnames(Gene_LGG)=gsub(".", "-", colnames(Gene_LGG), fixed = TRUE)
Gene_LGG=Gene_LGG[,-1]
# 
Survival_LGG=read.csv("Data/TCGA_LGG_Survival data.csv")
Survival_LGG=Survival_LGG[,c(1,2,4,87)]
rownames(Survival_LGG)=Survival_LGG[,1]
colnames(Survival_LGG)=c('SampleID','OS_event','OS_time','Target_Therapy')

#####
Gene_GBM=read.csv("Data/TCGA_GBM_GeneExpression.csv")
rownames(Gene_GBM)=Gene_GBM[,1]
colnames(Gene_GBM)=gsub(".", "-", colnames(Gene_GBM), fixed = TRUE)
Gene_GBM=Gene_GBM[,-1]

Survival_GBM=read.csv("Data/TCGA_GBM_Survival data.csv")
Survival_GBM=Survival_GBM[,c(1,33,32,106)]
rownames(Survival_GBM)=Survival_GBM[,1]
colnames(Survival_GBM)=c('SampleID','OS_event','OS_time','Sample Type')


Gene_TCGA=Gene_GBM
colnames(Gene_TCGA)=gsub(".", "-", colnames(Gene_TCGA), fixed = TRUE)
Survival_TCGA=Survival_GBM[intersect(colnames(Gene_GBM),rownames(Survival_GBM)),]

Gene_Normal=Gene_TCGA[,Survival_TCGA[,4]=='Solid Tissue Normal']
Gene_Tumor=Gene_TCGA[,Survival_TCGA[,4]!='Solid Tissue Normal']
dim(Gene_Normal)
dim(Gene_Tumor)

Gene_TCGA=cbind(Gene_LGG,Gene_GBM)
Gene_Tumor=Gene_TCGA[,setdiff(colnames(Gene_TCGA),colnames(Gene_Normal))]

MarkerGenes=c("KIF2C", "CCNA2", "NDC80", "KIF11", "KIF23", "ANLN", "CENPM")

# ID=ensemble2symbol[intersect(rownames(ensemble2symbol),GeneName),] 

Marker_Normal=matrix(as.numeric(as.matrix(Gene_Normal[MarkerGenes,])),dim(Gene_Normal[MarkerGenes,]))
Marker_Tumor= matrix(as.numeric(as.matrix(Gene_Tumor[MarkerGenes,])),dim(Gene_Tumor[MarkerGenes,]))

pvalue=matrix(0,1,length(MarkerGenes))
for (i in 1:length(MarkerGenes))
{
  Wtest=wilcox.test(Marker_Normal[i,],Marker_Tumor[i,])
  pvalue[i]=Wtest$p.value
}

Data=as.data.frame(cbind(Marker_Normal,Marker_Tumor))
Data=t(na.omit(Data))

type=as.factor(c(rep('Normal',5),rep('Tumor',697)))
x=cbind(Data,(type))

# install.packages("vioplot")
library(vioplot)

D=x[x[,8]==1,1:7]
E=x[x[,8]==2,1:7]
dev.new()
vioplot(D[,1],D[,2],D[,3],D[,4],D[,5],D[,6],D[,7],at=seq(1,19,3), names=paste(as.character(MarkerGenes), 'Normal',sep='_'), xlim=c(-0.5,16),col="DarkKhaki", border="gray51", lty=1, lwd=1, rectCol="gray", 
        colMed="white", pchMed=19)
vioplot(E[,1],E[,2],E[,3],E[,4],E[,5],E[,6],E[,7],at=seq(2,20,3), add=T, names=paste(as.character(MarkerGenes), 'Tumor',sep='_'), xlim=c(-0.5,16),col="DeepSkyBlue3", border="gray51", lty=1, lwd=1, rectCol="gray", 
        colMed="white", pchMed=19)

##################################################################################################################
##################################################################################################################

#################  K-M analysis

##################  ensemble ID  to gene symbol   #################################
ensembl2gene <- toTable(org.Hs.egENSEMBL2EG)
gene2symbol <- toTable(org.Hs.egSYMBOL)
ensemble2symbol <- merge(ensembl2gene, gene2symbol, by = 'gene_id')[2:3]
ensemble2symbol=as.matrix(ensemble2symbol)
rownames(ensemble2symbol)=ensemble2symbol[,2]

###TCGA data
setwd("Path/Survival analysis")

Gene_LGG=read.csv("Data/TCGA_LGG_GeneExpression.csv")
rownames(Gene_LGG)=Gene_LGG[,1]
colnames(Gene_LGG)=gsub(".", "-", colnames(Gene_LGG), fixed = TRUE)
Gene_LGG=Gene_LGG[,-1]

Survival_LGG=read.csv("Data/TCGA_LGG_Survival data.csv")
Survival_LGG=Survival_LGG[,c(1,2,4)]
rownames(Survival_LGG)=Survival_LGG[,1]
colnames(Survival_LGG)=c('SampleID','OS_event','OS_time')

#####
Gene_GBM=read.csv("Data/TCGA_GBM_GeneExpression.csv")
rownames(Gene_GBM)=Gene_GBM[,1]
colnames(Gene_GBM)=gsub(".", "-", colnames(Gene_GBM), fixed = TRUE)
Gene_GBM=Gene_GBM[,-1]

Survival_GBM=read.csv("Data/TCGA_GBM_Survival data.csv")
Survival_GBM=Survival_GBM[,c(1,33,32)]
rownames(Survival_GBM)=Survival_GBM[,1]
colnames(Survival_GBM)=c('SampleID','OS_event','OS_time')

Gene_TCGA=cbind(Gene_GBM[intersect(rownames(Gene_GBM),rownames(Gene_LGG)),],Gene_LGG[intersect(rownames(Gene_GBM),rownames(Gene_LGG)),])
colnames(Gene_TCGA)=gsub(".", "-", colnames(Gene_TCGA), fixed = TRUE)
Survival_TCGA=rbind(Survival_GBM,Survival_LGG)


######################## Select genes for signature   #####################################

MNB=c("KIF2C", "CCNA2", "NDC80", "KIF11", "KIF23", "ANLN", "CENPM")  # Cell Cycle-related gene signature

Survival_Gene=Gene_TCGA[t(MNB),]
colnames(Survival_Gene)=gsub(".", "-", colnames(Survival_Gene), fixed = TRUE)


Survival=Survival_TCGA[intersect(rownames(Survival_TCGA),colnames(Survival_Gene)),]
Survival=as.matrix(Survival)
Survival_Gene=Survival_Gene[,intersect(rownames(Survival),colnames(Survival_Gene))]

### Training dataset
status=as.numeric(Survival[,2])
time=as.numeric(Survival[,3])
Genes=rownames(na.omit(Survival_Gene))
Gene_marker=matrix(as.double(as.matrix(Survival_Gene[Genes,])),dim(Survival_Gene))



time=time[!is.na(status)]
Gene_marker=Gene_marker[,!is.na(status)]
status=status[!is.na(status)]

### COX Model
surv=Surv(as.double(time),as.double(status))
fit <- coxph(surv ~ t(Gene_marker), method="efron")
A=coef(fit)


# 
# surv=Surv(as.double(time),as.double(status))
# 
# cvfit<-cv.glmnet(t(Gene_marker),surv,family="cox",alpha=1)
# plot(cvfit)
# A=coef(cvfit, s = "lambda.min")
# A=coef(cvfit, s="lambda.1se")
# which(A!=0)
# sum(A!=0)
# MNB[which(A!=0)]
# 
# fit <- coxph(Surv(time, status) ~ t(Gene_marker[which(A!=0),]), method="breslow")

########  ROC for Training set
predicted=matrix(0,1,dim(Gene_marker)[2])
for (i in 1:dim(Gene_marker)[2])
{
  predicted[i]=sum(A*Gene_marker[,i])
}

predicted0=as.numeric(predicted)
predicted0=predicted0[!is.na(status)]
status0=as.matrix(na.omit(status))
pred <- prediction(predicted0,status0)
perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
performance(pred,"auc") # shows calculated AUC for model
dev.new()
plot(perf,colorize=FALSE, col="DarkCyan") # plot ROC curve
lines(c(0,1),c(0,1),col = "gray", lty = 4 )

roc0=roc((status0),predicted0,smooth=TRUE,plot=TRUE)
roc0
AUC=auc(roc0,predicted0)
AUC
## optimal combination
opt <- which.max(rowSums(cbind(roc0$sensitivities,roc0$specificities)))
## optimal cut-off point 
sort(predicted0,F)[opt]

##########  K-M survival curves for MNB in TCGA datasets
groups=matrix(0,1,length(status))
groups[predicted<=sort(predicted,F)[opt]]=1
groups[predicted>sort(predicted,F)[opt]]=2
# groups[predicted<=median(predicted)]=1
# groups[predicted>median(predicted)]=2
groups=t(groups)
groups=as.numeric(groups)


fit<- survfit(Surv(time/365, status) ~ groups, data = as.data.frame(Survival_Gene))
# Drawing survival curves
dev.new()
ggsurvplot(fit, legend = c(0.2, 0.2))




# Add risk table
# and change risk table y text colors by strata
ggsurvplot(fit, pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, risk.table.y.text.col = TRUE)




