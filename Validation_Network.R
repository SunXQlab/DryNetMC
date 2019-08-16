
library(glmnet)

library(pracma)

#### ROC package
# install.packages("gplots")
library(gplots)
# install.packages("ROCR")
library(ROCR)
# install.packages("dplyr")
library("dplyr")
# install.packages("lubridate")
library(lubridate)
# install.packages("robustbase")
# install.packages("caret")
library(caret)
library(pROC)

# install.packages("timeROC")
# install.packages("prodlim")
# library(prodlim)
# library(quantreg)
# install.packages("quantreg")
# install.packages("polspline")
# library("polspline")
# library(timeROC)

setwd("Path/Validation_Network")
A=read.csv('ODE_Coefficients.csv',fill= T,header = F)


h=0.1
time=seq(0,48,h)
BB=A

  x=matrix(0,dim(BB)[1],length(time)+1)
  CC=matrix(as.numeric(unlist(A)),dim(A))*0.2
  C=as.numeric(matrix(1,5,1)*(-0.1))
  ind=1
  x[,1]=as.vector(rand(1,dim(BB)[1]))
  
  for (t in time)
  {
    x[,ind+1]=x[,ind]+h*(CC%*%x[,ind]+C*x[,ind])

    ind=ind+1
  }

  # 
  TimeP=c(0,6,12,24,48)
  x_P=x[,TimeP*10+1]

 
  
  #####  hermite interpolation 
  # library(pracma)
  
  x_con=matrix(0,dim(x)[1],length(0:48))
  for (i in 1:dim(x)[1])
  {
    x_con[i,]=pchip(c(0,6,12,24,48),x_P[i,],c(0:48))  # Piecewise cubic hermite interpolation 
  }
  
  dev.new()
  plot(1:482,x[1,], col="DodgerBlue",ylim=c(-1,1),pch=16, cex=0.2,lty=1,lwd=2,xaxt="n",yaxt="n",ann=F) 
  lines(1:482,x[2,],col = "brown" ,lty=1,lwd=2)
  lines(1:482,x[3,],col = "DarkGoldenrod4" ,lty=1,lwd=2)
  lines(1:482,x[4,],col = "green" ,lty=1,lwd=2)
  lines(1:482,x[5,],col = "DarkSlateGray" ,lty=1,lwd=2)
  axis(1,TimeP*10+1, c(0,6,12,24,48))
  axis(2,seq(-1,1,0.5), seq(-1,1,0.5))
  # legend("bottomleft", inset=.05,title="Original data", c("A","B","C","D","E"), lty=c(1,1,1,1,1), lwd=c(2,2,2,2,2), col=c("DodgerBlue", "brown","DarkGoldenrod4","green","DarkSlateGray"))
  
  par(new=TRUE)
  plot(TimeP*10+1,x[1,TimeP*10+1], col="DodgerBlue",ylim=c(-1,1),pch=16, cex=1.5,xaxt="n",yaxt="n",ann="n") 
  par(new=TRUE)
  plot(TimeP*10+1,x[2,TimeP*10+1], col="brown",ylim=c(-1,1),pch=16, cex=1.5,xaxt="n",yaxt="n",ann="n") 
  par(new=TRUE)
  plot(TimeP*10+1,x[3,TimeP*10+1], col="DarkGoldenrod4",ylim=c(-1,1),pch=16, cex=1.5,xaxt="n",yaxt="n",ann="n") 
  par(new=TRUE)
  plot(TimeP*10+1,x[4,TimeP*10+1], col="green",ylim=c(-1,1),pch=16, cex=1.5,xaxt="n",yaxt="n",ann="n") 
  par(new=TRUE)
  plot(TimeP*10+1,x[5,TimeP*10+1], col="DarkSlateGray",ylim=c(-1,1),pch=16, cex=1.5,xaxt="n",yaxt="n",ann="n") 
  
  # legend("bottom", inset=.05, title="Sampling", c("A","B","C","D","E"), pch=c(16,16,16,16,16), cex=c(1,1,1,1,1), col=c("DodgerBlue", "brown","DarkGoldenrod4","green","DarkSlateGray"))
  
par(new=TRUE)
plot(1:49,x_con[1,], col="DodgerBlue",ylim=c(-1,1),pch=16, cex=0.2,lty=3,lwd=2,xaxt="n",yaxt="n",ann="n") 
lines(1:49,x_con[2,],col = "brown" ,lty=3,lwd=2)
lines(1:49,x_con[3,],col = "DarkGoldenrod4" ,lty=3,lwd=2)
lines(1:49,x_con[4,],col = "green" ,lty=3,lwd=2)
lines(1:49,x_con[5,],col = "DarkSlateGray" ,lty=3,lwd=2)
# legend("bottomright", inset=.05, title="Interpolation", c("A","B","C","D","E"), lty=c(3,3,3,3,3), col=c("DodgerBlue", "brown","DarkGoldenrod4","green","DarkSlateGray"))


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

pred <- prediction(as.vector(abs(B)),as.vector(sign(abs(CC))))
perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
performance(pred,"auc") # shows calculated AUC for model
dev.new()
plot(perf,colorize=FALSE, col="red",lwd=2) # plot ROC curve
lines(c(0,1),c(0,1),col = "gray", lty = 4 )

roc=roc(as.vector(sign(abs(CC))),as.vector(abs(B)))
roc
AUC=roc$auc
AUC

opt <- which.max(rowSums(cbind(roc$sensitivities,roc$specificities)))
## optimal cut-off point 
sort(as.vector(abs(B)),F)[opt]


####### Compare CorrNet=PCC*(PCC_p<0.05) with sign(abs(CC))
CorrNet=PCC*(PCC_p<0.05)
# CorrNet=1/(1e-3+PCC_p)
diag(CorrNet)=0

pred <- prediction(as.vector(abs(CorrNet)),as.vector(sign(abs(CC))))
perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
performance(pred,"auc") # shows calculated AUC for model
par(new=TRUE)
plot(perf,colorize=FALSE, col="blue",lwd=2) # plot ROC curve
lines(c(0,1),c(0,1),col = "gray", lty = 4 )

roc1=roc(as.vector(sign(abs(CC))),as.vector(abs(CorrNet)))
roc1
AUC=roc1$auc
AUC

### Compare with GRENITS (Bayesian method)

library(GRENITS)
# help(GRENITS)

# ### An Example 
# data(Athaliana_ODE)
# output.folder <- paste(tempdir(), "/Example_LinearNet", sep="")
# LinearNet(output.folder, Athaliana_ODE)
# analyse.output(output.folder)
# dir(output.folder)
# prob.file <- paste(output.folder, "/NetworkProbability_Matrix.txt", sep = "")
# prob.mat <- read.table(prob.file)
# print(prob.mat)


### Performance on ODE-simulated data
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
par(new=TRUE)
plot(perf,colorize=FALSE, col="blue",lwd=2) # plot ROC curve
lines(c(0,1),c(0,1),col = "gray", lty = 4 )

roc1=roc(as.vector(sign(abs(CC))),as.vector(M))
roc1
AUC=roc1$auc
AUC
