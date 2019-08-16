

library(glmnet)

setwd("Path/Validation_Network")
AA=read.csv('ODE_Coefficients.csv',fill= T,header = F)




AUC_DryNetMC=matrix(0,1,100)
AUC_PCCNet=matrix(0,1,100)
AUC_ODELASSO=matrix(0,1,100)
AUC_GENIE3=matrix(0,1,100)
AUC_ODELASSOP=matrix(0,1,100)
AUC_GRENITS=matrix(0,1,100)

Prec_DryNetMC=matrix(0,1,100)
Prec_PCCNet=matrix(0,1,100)
Prec_ODELASSO=matrix(0,1,100)
Prec_GENIE3=matrix(0,1,100)
Prec_ODELASSOP=matrix(0,1,100)

Rec_DryNetMC=matrix(0,1,100)
Rec_PCCNet=matrix(0,1,100)
Rec_ODELASSO=matrix(0,1,100)
Rec_GENIE3=matrix(0,1,100)
Rec_ODELASSOP=matrix(0,1,100)

F1_DryNetMC=matrix(0,1,100)
F1_PCCNet=matrix(0,1,100)
F1_ODELASSO=matrix(0,1,100)
F1_GENIE3=matrix(0,1,100)
F1_ODELASSOP=matrix(0,1,100)

for (Net_ind in 1:100)
{
  cat(Net_ind)
  
  h=0.1
  time=seq(0,48,h)
  BB=AA
  
  x=matrix(0,dim(BB)[1],length(time)+1)
  CC=matrix(as.numeric(unlist(AA)),dim(AA))*0.2
  C=as.numeric(matrix(1,5,1)*(-0.1))
  x[,1]=runif(dim(BB)[1])
  # randomly generate networks
  element1=sample(1:5, 1)
  element2=sample(1:5, 1)
  CC[element1,element2]=runif(1, -0.1, 0.1)
  diag(CC)=0
  
  ## Calculate AUC of different methods
  ind=1
  for (t in time)
  {
    x[,ind+1]=x[,ind]+h*(CC%*%x[,ind]+C*x[,ind])
    
    ind=ind+1
  }
  
  # 
  TimeP=c(0,6,12,24,48)
  x_P=x[,TimeP*10+1]
  
  write.table(x_P,"Simulated_Gene_Expression_Matrix.txt")
  
  #####  hermite interpolation 
  # library(pracma)
  
  x_con=matrix(0,dim(x)[1],length(0:48))
  for (i in 1:dim(x)[1])
  {
    x_con[i,]=pchip(c(0,6,12,24,48),x_P[i,],c(0:48))  # Piecewise cubic hermite interpolation 
  }
  
  
  ######### PPI information
  PPI=sign(rand(5,5)>0.6)+sign(abs(CC))
  PPI=sign(PPI)
  diag(PPI)=0
  ##### Correlation
  Netsize=5
  PCC= matrix(NA,nrow=Netsize,ncol=Netsize)  # Partial Pearson correlation coefficient 
  PCC_p= matrix(NA,nrow=Netsize,ncol=Netsize)  # p value 
  
  for (i in 1:Netsize)
  {
    for (j in 1:Netsize)
    {
      
      B=cor.test(as.numeric(x_con[i,]),as.numeric(x_con[j,]),method="pearson")
      PCC[i,j]=B$estimate
      PCC_p[i,j]=B$p.value
    }
  }
  
  
  CorrNet=matrix(as.numeric(PCC_p<0.05),5,5)   #(abs(PCC)>0.5&
  E=(PCC_p<0.05)*PPI  #abs(PCC)>0.75&
  E=PPI
  ############## ODE network
  y=t(diff(t(x_con)))
  x=x_con[,1:48]
  
  A = matrix(0,nrow=dim(x)[1],ncol=dim(x)[1]+1)  # 
  
  
  for (i in 1:dim(x)[1])
  {
    if (sum(E[i,]==0)!=61)
    {
      cvfit=cv.glmnet(t(x),y[i,],family="gaussian",alpha = 1,exclude=which(E[i,]==0))  # range of lambda: 10^(-5) to 10^(-1)
      Coef = coef(cvfit)
      Coef_min = coef(cvfit,s="lambda.min")
      A[i,]=as.numeric(Coef_min) #AA
      
    }
    
  }
  
  
  B=A[,-c(1)]
  
  ####### Compare B with sign(abs(CC))
  # library(ROCR)
  # library(pROC)
  pred <- prediction(as.vector(abs(B)),as.vector(sign(abs(CC))))
  perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
  performance(pred,"auc") # shows calculated AUC for model
  # dev.new()
  # plot(perf,colorize=FALSE, col="red",lwd=2) # plot ROC curve
  # lines(c(0,1),c(0,1),col = "gray", lty = 4 )
  
  roc1=roc(as.vector(sign(abs(CC))),as.vector(abs(B)))
  # roc
  AUC1=roc1$auc
  AUC_DryNetMC[Net_ind]=AUC1
  
  opt <- which.max(rowSums(cbind(roc1$sensitivities,roc1$specificities)))
  
  # ACC=unlist(performance(pred,"acc")@y.values)[opt]
  # ACC_DryNetMC[Net_ind]=ACC
  
  Prec=unlist(performance(pred,"prec")@y.values)[opt]
  Prec_DryNetMC[Net_ind]=Prec
  
  Rec=unlist(performance(pred,"rec")@y.values)[opt]
  Rec_DryNetMC[Net_ind]=Rec
  
  F1=unlist(performance(pred,"f")@y.values)[opt]
  F1_DryNetMC[Net_ind]=F1
  

  
  ####### Compare CorrNet=PCC*(PCC_p<0.05) with sign(abs(CC))
  CorrNet=PCC*(PCC_p<0.05)
  # CorrNet=1/(1e-3+PCC_p)
  diag(CorrNet)=0
  
  pred <- prediction(as.vector(abs(CorrNet)),as.vector(sign(abs(CC))))
  perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
  performance(pred,"auc") # shows calculated AUC for model
  # par(new=TRUE)
  # plot(perf,colorize=FALSE, col="blue",lwd=2) # plot ROC curve
  # lines(c(0,1),c(0,1),col = "gray", lty = 4 )
  
  roc2=roc(as.vector(sign(abs(CC))),as.vector(abs(CorrNet)))
  # roc1
  AUC2=roc2$auc
  AUC_PCCNet[Net_ind]=AUC2
  opt <- which.max(rowSums(cbind(roc2$sensitivities,roc2$specificities)))
  
  # ACC=unlist(performance(pred,"acc")@y.values)[opt]
  # ACC_PCCNet[Net_ind]=ACC
  
  Prec=unlist(performance(pred,"prec")@y.values)[opt]
  Prec_PCCNet[Net_ind]=Prec
  
  Rec=unlist(performance(pred,"rec")@y.values)[opt]
  Rec_PCCNet[Net_ind]=Rec
  
  F1=unlist(performance(pred,"f")@y.values)[opt]
  F1_PCCNet[Net_ind]=F1
  
  ## ODE_LASSO
  y=t(diff(t(x_P)))
  x=x_P[,1:4]
  
  A = matrix(0,nrow=dim(x)[1],ncol=dim(x)[1]+1)  # 
  
  
  for (i in 1:dim(x)[1])
  {

      cvfit=cv.glmnet(t(x),y[i,],family="gaussian",alpha = 1)  # range of lambda: 10^(-5) to 10^(-1)
      Coef = coef(cvfit)
      Coef_min = coef(cvfit,s="lambda.min")
      A[i,]=as.numeric(Coef_min) #AA
  }
  
  B=A[,-c(1)]
  
  ####### Compare B with sign(abs(CC))
  pred <- prediction(as.vector(abs(B)),as.vector(sign(abs(CC))))
  perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
  performance(pred,"auc") # shows calculated AUC for model
  # dev.new()
  # plot(perf,colorize=FALSE, col="red",lwd=2) # plot ROC curve
  # lines(c(0,1),c(0,1),col = "gray", lty = 4 )
  
  roc3=roc(as.vector(sign(abs(CC))),as.vector(abs(B)))
  # roc
  AUC3=roc3$auc
  AUC_ODELASSO[Net_ind]=AUC3
  
  opt <- which.max(rowSums(cbind(roc3$sensitivities,roc3$specificities)))
  # 
  # ACC=unlist(performance(pred,"acc")@y.values)[opt]
  # ACC_ODELASSO[Net_ind]=ACC
  
  Prec=unlist(performance(pred,"prec")@y.values)[opt]
  Prec_ODELASSO[Net_ind]=Prec
  
  Rec=unlist(performance(pred,"rec")@y.values)[opt]
  Rec_ODELASSO[Net_ind]=Rec
  
  F1=unlist(performance(pred,"f")@y.values)[opt]
  F1_ODELASSO[Net_ind]=F1
  
  
  ## ODE_LASSO_Pri
  y=t(diff(t(x_P)))
  x=x_P[,1:4]
  
  A = matrix(0,nrow=dim(x)[1],ncol=dim(x)[1]+1)  # 
  
  
  for (i in 1:dim(x)[1])
  {
    if (sum(E[i,]==0)!=61)
    {
      cvfit=cv.glmnet(t(x),y[i,],family="gaussian",alpha = 1,exclude=which(E[i,]==0))  # range of lambda: 10^(-5) to 10^(-1)
      Coef = coef(cvfit)
      Coef_min = coef(cvfit,s="lambda.min")
      A[i,]=as.numeric(Coef_min) #AA
      
    }
    
  }
  
  B=A[,-c(1)]
  
  ####### Compare B with sign(abs(CC))
  pred <- prediction(as.vector(abs(B)),as.vector(sign(abs(CC))))
  perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
  performance(pred,"auc") # shows calculated AUC for model
  # dev.new()
  # plot(perf,colorize=FALSE, col="red",lwd=2) # plot ROC curve
  # lines(c(0,1),c(0,1),col = "gray", lty = 4 )
  
  roc3=roc(as.vector(sign(abs(CC))),as.vector(abs(B)))
  # roc
  AUC3=roc3$auc
  AUC_ODELASSOP[Net_ind]=AUC3
  
  opt <- which.max(rowSums(cbind(roc3$sensitivities,roc3$specificities)))
  
  # ACC=unlist(performance(pred,"acc")@y.values)[opt]
  # ACC_ODELASSOP[Net_ind]=ACC
  
  Prec=unlist(performance(pred,"prec")@y.values)[opt]
  Prec_ODELASSOP[Net_ind]=Prec
  
  Rec=unlist(performance(pred,"rec")@y.values)[opt]
  Rec_ODELASSOP[Net_ind]=Rec
  
  F1=unlist(performance(pred,"f")@y.values)[opt]
  F1_ODELASSOP[Net_ind]=F1
  
  
  ## GENIE3
  # install.packages("randomForest")
  # source("GENIE3.R")
  x_P=read.expr.matrix("Simulated_Gene_Expression_Matrix.txt",form="rows.are.genes")
  B=GENIE3(as.data.frame(x_P))
  B=t(B)
  ####### Compare B with sign(abs(CC))
  pred <- prediction(as.vector(abs(B)),as.vector(sign(abs(CC))))
  perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
  performance(pred,"auc") # shows calculated AUC for model
  # dev.new()
  # plot(perf,colorize=FALSE, col="red",lwd=2) # plot ROC curve
  # lines(c(0,1),c(0,1),col = "gray", lty = 4 )
  
  roc4=roc(as.vector(sign(abs(CC))),as.vector(abs(B)))
  # roc
  AUC4=roc4$auc
  AUC_GENIE3[Net_ind]=AUC4
  
  opt <- which.max(rowSums(cbind(roc4$sensitivities,roc4$specificities)))
  
  # ACC=unlist(performance(pred,"acc")@y.values)[opt]
  # ACC_GENIE3[Net_ind]=ACC
  
  Prec=unlist(performance(pred,"prec")@y.values)[opt]
  Prec_GENIE3[Net_ind]=Prec
  
  Rec=unlist(performance(pred,"rec")@y.values)[opt]
  Rec_GENIE3[Net_ind]=Rec
  
  F1=unlist(performance(pred,"f")@y.values)[opt]
  F1_GENIE3[Net_ind]=F1
  
  
  ### GRENITS
  output.folder <- "Path/Validation_Network"
  LinearNet(output.folder, as.data.frame(x_P) )
  analyse.output(output.folder)
  dir(output.folder)
  prob.file <- paste(output.folder, "/NetworkProbability_Matrix.txt", sep = "")
  prob.mat <- read.table(prob.file)
  print(prob.mat)
  
  M=as.matrix(prob.mat)
  
  pred <- prediction(as.vector(M),as.vector(sign(abs(CC))))
  perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
  performance(pred,"auc") # shows calculated AUC for model
  # par(new=TRUE)
  # plot(perf,colorize=FALSE, col="darkgoldenrod",lwd=2) # plot ROC curve
  # lines(c(0,1),c(0,1),col = "gray", lty = 4 )
  # 
  roc1=roc(as.vector(sign(abs(CC))),as.vector(M))
  AUC6=roc1$auc

  AUC_GRENITS[Net_ind]=AUC6
  
}


Wtest=wilcox.test(AUC_DryNetMC,AUC_PCCNet,alternative = 'greater',paired = T)  #"two.sided" (default)  #,paired = T
pvalue=Wtest$p.value
pvalue

Wtest=wilcox.test(AUC_DryNetMC,AUC_ODELASSO,alternative = 'greater',paired = T)  #"two.sided" (default)  #,paired = T
pvalue=Wtest$p.value
pvalue

Wtest=wilcox.test(AUC_DryNetMC,AUC_ODELASSOP,alternative = 'greater',paired = T)  #"two.sided" (default)  #,paired = T
pvalue=Wtest$p.value
pvalue

Wtest=wilcox.test(AUC_DryNetMC,AUC_GENIE3,alternative = 'greater',paired = T)  #"two.sided" (default)  #,paired = T
pvalue=Wtest$p.value
pvalue

Wtest=wilcox.test(AUC_DryNetMC,AUC_GRENITS,alternative = 'greater',paired = T)  #"two.sided" (default)  #,paired = T
pvalue=Wtest$p.value
pvalue

dataset <- data.frame(value = c(AUC_PCCNet,AUC_GRENITS, AUC_GENIE3,AUC_ODELASSO,AUC_ODELASSOP,AUC_DryNetMC), group = factor(rep(c("PCCNet", "GRENITS", "GENIE3", "OdeLasso","OdeLassoP","DryNetMC"), times = c(length(AUC_PCCNet), length(AUC_GRENITS), length(AUC_GENIE3),length(AUC_ODELASSO),length(AUC_ODELASSOP), length(AUC_DryNetMC)))))

dev.new()
boxplot( value ~ group,  notch = F, dataset, border = c( "red", "darkgoldenrod", "purple","#009E73","green", "blue"),cex = 1,cex.axis=1,pars = list(boxwex = 0.5, staplewex = 0.5, outwex = 0.5),ylim = c(0, 1))  #,col.axis = "#009E73"



dataset <- data.frame(value = c(Prec_PCCNet,Prec_GENIE3,Prec_ODELASSO,Prec_ODELASSOP,Prec_DryNetMC), group = factor(rep(c("PCCNet", "GENIE3", "OdeLasso","OdeLassoP","DryNetMC"), times = c(length(AUC_PCCNet), length(AUC_GENIE3),length(AUC_ODELASSO),length(AUC_ODELASSOP), length(AUC_DryNetMC)))))

dev.new()
boxplot( value ~ group,  notch = F, dataset, border = c( "red","purple","#009E73","green", "blue"),cex = 1,cex.axis=1,pars = list(boxwex = 0.5, staplewex = 0.5, outwex = 0.5))  #,col.axis = "#009E73"

dataset <- data.frame(value = c(Rec_PCCNet,Rec_GENIE3,Rec_ODELASSO,Rec_ODELASSOP,Rec_DryNetMC), group = factor(rep(c("PCCNet", "GENIE3", "OdeLasso","OdeLassoP","DryNetMC"), times = c(length(AUC_PCCNet), length(AUC_GENIE3),length(AUC_ODELASSO),length(AUC_ODELASSOP), length(AUC_DryNetMC)))))

dev.new()
boxplot( value ~ group,  notch = F, dataset, border = c( "red","purple","#009E73","green", "blue"),cex = 1,cex.axis=1,pars = list(boxwex = 0.5, staplewex = 0.5, outwex = 0.5))  #,col.axis = "#009E73"

dataset <- data.frame(value = c(F1_PCCNet,F1_GENIE3,F1_ODELASSO,F1_ODELASSOP,F1_DryNetMC), group = factor(rep(c("PCCNet", "GENIE3", "OdeLasso","OdeLassoP","DryNetMC"), times = c(length(AUC_PCCNet), length(AUC_GENIE3),length(AUC_ODELASSO),length(AUC_ODELASSOP), length(AUC_DryNetMC)))))

dev.new()
boxplot( value ~ group,  notch = F, dataset, border = c( "red","purple","#009E73","green", "blue"),cex = 1,cex.axis=1,pars = list(boxwex = 0.5, staplewex = 0.5, outwex = 0.5))  #,col.axis = "#009E73"


# dev.new()
# # library(ggpubr)
# # library(digest)
# p<-ggboxplot(dataset, "group", "value",
#              color = "Samples", palette =c("#009E73","#00AFBB", "#FC4E07","#E7B800"),  #, "#E7B800"
#              linetype = 1,
#              title = FALSE, xlab = FALSE,
#              font.label = list(size = 50, face = "plain"),
#              order = c("PCCNet", "ODE_LASSO","GENIE3", "DryNetMC"),
#              # fill = "Samples", palette =c("#FC4E07", "#00AFBB"),
#              add = "jitter")
# 
# p
# 
# compar<-list(c("Normal","BRCA1_Mutated"))
# 
# p+stat_compare_means(comparisons = compar, method = "wilcox.test",  method.args = list(alternative = "less"),aes(label = paste0("p = ", ..p.format..)))+stat_compare_means(label.y = 1, show.legend = FALSE)   #
# 
# 
