

setwd("F:/胶质瘤诱导分化组学数据分析/Code-cluster_based/Results") 
install.packages("iterators")
library(glmnet)  # manually load this package by downloading and installing from Toll. 
library(foreach)
library(Matrix)
############ 


# DEG_R_cluster_con=read.csv('DEG_R_cluster_con.csv',fill= T,header = F)
# rowname=DEG_R_cluster_con[,1]
# DEG_R_cluster_con=DEG_R_cluster_con[,-1]
# DEG_R_cluster_con=matrix(as.numeric(as.matrix(DEG_R_cluster_con)),dim(DEG_R_cluster_con)[1],dim(DEG_R_cluster_con)[2])
# rownames(DEG_R_cluster_con)=rowname[-1]
# DEG_R_cluster_con=na.omit(DEG_R_cluster_con)


install.packages("pracma")
library(pracma)

DEG_R_cluster_con=matrix(0,dim(TC_gene_2_Normalized)[1],length(0:48))
for (i in 1:dim(TC_gene_2_Normalized)[1])
{
  DEG_R_cluster_con[i,]=pchip(c(0,6,12,24,48),TC_gene_2_Normalized[i,],c(0:48))  # Piecewise cubic hermite interpolation 
}
rownames(DEG_R_cluster_con)=rownames(TC_gene_2_Normalized)




Node_gene=DEG_R_cluster_con[NetGene,]
# Node_gene=Node_gene[,-1]
Node_gene=na.omit(Node_gene)
Netsize=dim(Node_gene)[1]  #75

############
Net_Gene=read.csv('Net_Gene.csv',fill= T,header = F)
Net_Gene=matrix(as.matrix(Net_Gene),dim(Net_Gene)[1],1)
x=DEG_R_cluster_con[Net_Gene,]
y=t(diff(t(x)))
x=x[,1:48]


### Select genes in the network

edge_R=read.table("edge matrix for Resistant genes.txt")
EM=edge_R
e=EM
e=e[2:62,2:62]

dim(e)

# A=rownames(edge_R[rowSums(abs(edge_R))!=0, ])
# B=colnames(edge_R[,colSums(abs(edge_R))!=0])
# Net_Gene=union(A,B)
# 
# edge_R=edge_R[Net_Gene,Net_Gene]
# dim(edge_R)  #51
# 
# x=x[Net_Gene,]
# dim(x)
# 
# EM=edge_R
# e=EM


A = matrix(0,nrow=dim(x)[1],ncol=dim(x)[1]+1)  # 
Lambda=matrix(NA,nrow=1,ncol=dim(x)[1])
Error_CV_mean=matrix(NA,nrow=1,ncol=dim(x)[1])

for (i in 1:dim(x)[1])
{

  if (sum(e[i,]==0)!=61)
  {
    x_del=x
    cvfit=cv.glmnet(t(x_del),y[i,],family="gaussian",alpha = 1,exclude=which(e[i,]==0))  # range of lambda: 10^(-5) to 10^(-1)
    Coef = coef(cvfit)
    Coef_min = coef(cvfit,s="lambda.min")
    #AA[AA!=0]=as.numeric(Coef_min)
    A[i,]=as.numeric(Coef_min) #AA
    Lambda[i]=cvfit$lambda.min
    Error_CV_mean[i] = min(cvfit$cvm)
  }
  
}

plot(cvfit)
title("Fitting",line=2.5)

plot(1:length(Lambda),log10(Lambda),type="b",main="Lambda",sub='',xlim=c(1,65), ylim=c(0, 0.05),xlab="Lambda NO.",ylab='Lambda Values')
# title("Lambda",line=2.5)
# bins=seq(min(log10(Lambda)),max(log10(Lambda)),by=0.1)
hist(log10(Lambda),freq = T, col=rainbow(4,alpha=0.5)[3],breaks = 30,xlim=c(-4.5,-2.5),ylim=c(0, 8))
# lines(density(log10(Lambda)), col="red")

hist(Error_CV_mean,freq = T,col=rainbow(4,alpha=0.5)[2], breaks = 30,xlim=c(0,0.0030),ylim=c(0, 15))

setwd("F:/胶质瘤诱导分化组学数据分析/Code-cluster_based/Results") 
write.table(A, file="Edge Coefficient_Resistant Network.xls",col.names = NA,sep = "\t")

A=read.table("Edge Coefficient_Resistant Network.xls",header = T,sep="\t",fileEncoding = "UTF-8",row.names=NULL)


max(max(A))
min(min(A))
sum(sum(A[,-1]!=0))   # 764
B=A[,-1]
B=as.matrix(B)
dim(B)


## AIC to select significant interactions

Theta=NULL
BIC=NULL
RSS=NULL


timestart<-Sys.time()
A = matrix(0,nrow=dim(x)[1],ncol=dim(x)[1]+1)  # 
variable=1:dim(x)[1]
for (i in variable)#)
{
  if (sum(e[i,]==0)!=61)
  {
    x_del=x
    cvfit=cv.glmnet(t(x_del),y[i,],family="gaussian",alpha = 1,exclude=which(e[i,]==0))  # range of lambda: 10^(-5) to 10^(-1): lambda=10^seq(-1,-5,-0.1),
    Coef = coef(cvfit)
    Coef_min = coef(cvfit,s="lambda.min")
    #AA[AA!=0]=as.numeric(Coef_min)
    A[i,]=as.numeric(Coef_min) #AA
  }
  
  
}

ind=1

for (theta in 10^seq(-5,0,0.5))  # seq(.)里的(to - from)/by 
{
  
  # Y_predict=matrix(NA,nrow=dim(x)[1],ncol=dim(x)[2])
  Y_predict=y
  
  for (i in variable)
  {    
    D=A[,-1]
    D[abs(D)<=theta]=0
    x_new=x
    x_new[which(D[i,]==0),]=0
    
    if (sum(e[i,]==0)!=61)
    {
      x_del=x
      cvfit=cv.glmnet(t(x_del),y[i,],family="gaussian",alpha = 1,exclude=which(e[i,]==0))  # range of lambda: 10^(-5) to 10^(-1): lambda=10^seq(-1,-5,-0.1),
      
      Y_predict[i,]=predict(cvfit,t(x_new),s="lambda.min")  #for predict of glmnet: type="response",exact=TRUE)[,1]
      # Y_predict=D%*%x+A[,1]
    }
    
  }
  
  
  p=sum(D!=0)
  E=y-Y_predict   #Y_predict=D%*%x+A[,1]   #(t(t(x)%*%t(A[,-1]))+A[,1])
  N=dim(E)[1]*dim(E)[2]
  R=sum(E^2)  #/(dim(E)[1]*dim(E)[2])
  RSS[ind]=R
  # BIC[ind]=log(1/N*R+1+log(2*pi))+2*p/N
  BIC[ind]=N*log(R/N)+p*log(N)
  Theta[ind]=theta
  ind=ind+1
}

MSE=(RSS/N)^0.5

# par(mfrow=c(2,2))
plot(log(Theta,10),MSE,type="b",main="MSE",sub='',xlim=c(-5,0), ylim=c(min(MSE), max(MSE)),xlab="theta",ylab='MSE',pch=21,cex=1.5,bg="blue")
plot(log(Theta,10),BIC,type="b",main="BIC",sub='',xlim=c(-5,0), ylim=c(min(BIC), max(BIC)),xlab="theta",ylab='BIC',pch=21,col="gray",cex=2,bg="blue",lty=4,lwd=2)


# R画图设置参见:  https://blog.csdn.net/zhyoulun/article/details/46430807  

timeend<-Sys.time()
runningtime<-timeend-timestart
print(runningtime) 


B[abs(B)<=0.01]=0
sum(sum(B!=0))   # 100
rownames(B)=Net_Gene
colnames(B)=Net_Gene

write.table(B, file="Significant Coefficient_Resistant Network.xls",quote=F, row.names=F,col.names = F,sep = "\t")

### judge if ODE network is a subnetwork of PPI 
for (i in 1:dim(x)[1])
{
  print(i)
  print(all(which(B[i,]!=0)%in%which(e[i,]!=0)))
}

### Select genes in the network
# a=rownames(edge_R[rowSums(abs(B))!=0, ])
# b=colnames(edge_R[,colSums(abs(B))!=0])
# DirectedNet_Gene=union(a,b)
# 
# B_Net=B[DirectedNet_Gene,DirectedNet_Gene]
# dim(B_Net)  # 56,56
# 
# write.table(B_Net, file="Coefficient_GeneNetwork_Resistant.txt",col.names = NA,sep = "\t") 


#### save gene expression data in the initial network
Node_gene_Directed=DEG_R_cluster_con[DirectedNet_Gene,]
write.table(Node_gene_Directed, file="Node_R_gene_Directed.xls",col.names = NA,sep = "\t") 
Node_gene_Directed["AURKB",]
### Reform txt data for cytoscape 
BB=B#_Net

C=matrix(NA,nrow=sum(sum(BB!=0)),ncol=3)
row=1
for (i in 1:dim(BB)[1])
{
  for (j in 1:dim(BB)[1])
  {
    if (BB[i,j]!=0)
    {
      C[row,1]=as.matrix(rownames(BB)[i])
      C[row,3]=as.matrix(colnames(BB)[j])
      C[row,2]=BB[i,j]
      row=row+1
    }
  }
}

write.table(C, file="Cytoscape_Net_Resistant_new.txt",quote=F,row.names=F,col.names = F,sep = "\t")





