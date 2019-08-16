


install.packages("iterators")
library(glmnet)  # manually load this package by downloading and installing from Toll. 

setwd("Path/Results") 
# DEG_S_cluster_con=read.csv("DEG_S_cluster_con.csv",fill= T,header = F)
# rowname=DEG_S_cluster_con[,1]
# DEG_S_cluster_con=DEG_S_cluster_con[,-1]
# DEG_S_cluster_con=matrix(as.numeric(as.matrix(DEG_S_cluster_con)),dim(DEG_S_cluster_con)[1],dim(DEG_S_cluster_con)[2])
# rownames(DEG_S_cluster_con)=rowname

install.packages("pracma")
library(pracma)

DEG_S_cluster_con=matrix(0,dim(TC_gene_1_Normalized)[1],length(0:48))
for (i in 1:dim(TC_gene_1_Normalized)[1])
{
  DEG_S_cluster_con[i,]=pchip(c(0,6,12,24,48),TC_gene_1_Normalized[i,],c(0:48))  # Piecewise cubic hermite interpolation 
}
rownames(DEG_S_cluster_con)=rownames(TC_gene_1_Normalized)



PPI_pair=read.csv("PPI.csv",header = FALSE)
NetGene12=PPI_pair[,1:2]
NetGene=union(NetGene12[,1],NetGene12[,2])
Netsize=length(NetGene)  #75 


Net_Gene=read.csv('Net_Gene.csv',fill= T,header = F)
Net_Gene=matrix(as.matrix(Net_Gene),dim(Net_Gene)[1],1)
x=DEG_S_cluster_con[Net_Gene,]
y=t(diff(t(x)))
x=x[,1:48]


edge_S=read.table("edge matrix for Sensitive genes.txt")
EM=edge_S
e=EM


dim(e)
e=e[2:62,2:62]

A = matrix(0,nrow=dim(x)[1],ncol=dim(x)[1]+1)  # 
Lambda=matrix(NA,nrow=1,ncol=dim(x)[1])
Error_CV_mean=matrix(NA,nrow=1,ncol=dim(x)[1])

for (i in 1:dim(x)[1])
{
  #penalty.factor=matrix(1,1,48)
  #penalty.factor[43]=Inf
  #exclude = which(e[i,]==0)  # set the link that are not in PPI to be 0.
  #exclude=matrix(0,1,44)
  #AA=matrix(1000,1,45)
  #AA[which(e[i,]==0)]=0
  if (sum(e[i,]==0)!=61)
  {
  x_del=x
  # x_del[which(e[i,]==0),] = rnorm(48,10,1)*1e20  # reset x to exclude variables that are not in PPI. #x[-which(e[i,]==0),]
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
bins=seq(min(log10(Lambda)),max(log10(Lambda)),by=0.1)
hist(log10(Lambda),freq = T, col=rainbow(4,alpha=0.5)[3],breaks = 30,xlim=c(-5,-2),ylim=c(0, 15))
# lines(density(log10(Lambda)), col="red")

hist(Error_CV_mean,freq = T,col=rainbow(4,alpha=1)[2], breaks = 30,xlim=c(0,0.0035),ylim=c(0, 25))


write.table(A, file="Edge Coefficient_Sensitive Network.xls",col.names = NA,sep = "\t")

A=read.table("Coefficient_Sensitive Network.xls",header = T,sep="\t",fileEncoding = "UTF-8",row.names=NULL)


max(max(A[,-1]))
min(min(A[,-1]))
hist(A[,-1],freq = F, breaks = 50)
sum(sum(A[,-1]!=0))   # 764
B=A[,-c(1,2)]
B=matrix(as.numeric(B),dim(B)[1],dim(B)[1])
dim(B)

## BIC to select significant interactions
Theta=NULL
BIC=NULL
RSS=NULL


timestart<-Sys.time()
A = matrix(0,nrow=dim(x)[1],ncol=dim(x)[1]+1)  # 
variable=1:dim(x)[1]
for (i in variable)#)
{
  if (sum(e[i,]==0)!=dim(e)[1])
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

for (theta in 10^seq(-5,0,0.1))  # seq(.)里的(to - from)/by 
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
plot(log(Theta,10),MSE,type="b",main="MSE",sub='',xlim=c(-5,0), ylim=c(min(MSE), max(MSE)),xlab="theta",ylab='MSE')
plot(log(Theta,10),BIC,type="b",main="BIC",sub='',xlim=c(-5,0), ylim=c(min(BIC), max(BIC)),xlab="theta",ylab='BIC',pch=21,col="gray",cex=2,bg="blue",lty=4,lwd=2)

timeend<-Sys.time()
runningtime<-timeend-timestart
print(runningtime) 



##########################
B[abs(B)<=0.01]=0
sum(sum(B!=0))   # 205
rownames(B)=t(Net_Gene)
colnames(B)=t(Net_Gene)

write.table(B, file="Significant Coefficient_Sensitive Network.xls",quote=F, row.names=F,col.names = F,sep = "\t")

### judge if ODE network is a subnetwork of PPI 
for (i in 1:dim(x)[1])
{
  print(i)
  print(all(which(B[i,]!=0)%in%which(e[i,]!=0)))
}

### Select genes in the network
a=rownames(edge_S[rowSums(abs(B))!=0, ])
b=colnames(edge_S[,colSums(abs(B))!=0])
DirectedNet_Gene=union(a,b)

B_Net=B[DirectedNet_Gene,DirectedNet_Gene]
dim(B_Net)  # 56,56

write.table(B_Net, file="Coefficient_GeneNetwork_Sensitive.txt",col.names = NA,sep = "\t") 


#### save gene expression data in the initial network
Node_gene_Directed=DEG_cluster_con[DirectedNet_Gene,]
write.table(Node_gene_Directed, file="Node_S_gene_Directed.xls",col.names = NA,sep = "\t") 
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
      C[row,1]=as.matrix(colnames(BB)[j])
      C[row,3]=as.matrix(rownames(BB)[i])
      C[row,2]=BB[i,j]
      row=row+1
    }
  }
}

write.table(C, file="Cytoscape_Net_Sensitive_new.txt",quote=F,row.names=F,col.names = F,sep = "\t")




