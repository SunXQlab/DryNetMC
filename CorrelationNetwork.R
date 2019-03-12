

##### read DEGs 
setwd("Path") 

install.packages("ggm")
library(ggm)

source("pcor.R")

if(!require("ppcor")){
  install.packages("ppcor", repos='http://cran.us.r-project.org')
  library(ppcor)
}

DEG_S=read.table("Path/Results/Drug Sensitive Genes.txt",fill= FALSE)
DEG_R=read.table("Path/Results/Drug Resistant Genes.txt",fill= FALSE) 

names_SR=union(rownames(DEG_S),rownames(DEG_R))


Gene_RPKM=read.csv("Path/GBM_gene_RPKM.csv",fill= T,header = F)
Gene_RPKM=na.omit(Gene_RPKM)
Gene_S=Gene_RPKM[,2:6]  ## Sensitive cell line
Gene_R=Gene_RPKM[,7:11] ## Resistant cell line
Gene_S=matrix(as.numeric(as.matrix(Gene_S)),dim(Gene_S)[1],dim(Gene_S)[2])
Gene_R=matrix(as.numeric(as.matrix(Gene_R)),dim(Gene_R)[1],dim(Gene_R)[2])
rownames(Gene_S)=Gene_RPKM[,1]
rownames(Gene_R)=Gene_RPKM[,1]



Node_gene_1=Gene_S[names_SR,]
Node_gene_2=Gene_R[names_SR,]

write.table(Node_gene_1,file="Path/Results/PPI Node",sep = "\t", col.names = NA)


##### Initial PPI network
PPI_pair=read.csv("Path/Results/PPI.csv",header = FALSE)
NetGene12=PPI_pair[,1:2]
NetGene=union(NetGene12[,1],NetGene12[,2])
Netsize=length(NetGene)  #75 

PPI=matrix(0,nrow=Netsize,ncol=Netsize)
for (i in 1:Netsize)
{
  for (j in 1:Netsize)
  {
    for (k in 1:dim(NetGene12)[1])
    {
      #if (is.na(prod(which(as.matrix(PPI_pair)[k,]==c(names_SR[i],names_SR[j]))))==0)
      if (sum(c(NetGene[i],NetGene[j])==as.matrix(NetGene12)[k,])==2)
        PPI[i,j]=PPI[i,j]+1
      else
        PPI[i,j]=PPI[i,j]+0  
      
    }
    
  }  
}


for (i in 1:Netsize)
{
  for (j in 1:Netsize)
  {
    PPI[i,j]=max(PPI[i,j],PPI[j,i])
  }
}

rownames(PPI)=NetGene
colnames(PPI)=NetGene
sum(PPI)  #=  edges

setwd("Path/Results")  
write.table(PPI, file="PPI matrix.xls",sep='\t') 


PC_S= matrix(NA,nrow=Netsize,ncol=Netsize)  # Partial Pearson correlation coefficient 
PC_p_S= matrix(NA,nrow=Netsize,ncol=Netsize)  # p value 

PC_R= matrix(NA,nrow=Netsize,ncol=Netsize)  # Partial Pearson correlation coefficient 
PC_p_R= matrix(NA,nrow=Netsize,ncol=Netsize)  # p value 



##   Correlation for Sensitive genes

for (i in 1:Netsize)
{
  for (j in 1:Netsize)
  {

    B=cor.test(as.numeric(Node_gene_1[i,]),as.numeric(Node_gene_1[j,]),method="pearson")
    PC_S[i,j]=B$estimate
    PC_p_S[i,j]=B$p.value
  }
}


##  Correlation for Resistant genes
for (i in 1:Netsize)
{
  for (j in 1:Netsize)
  {
    
    B=cor.test(as.numeric(Node_gene_2[i,]),as.numeric(Node_gene_2[j,]),method="pearson")
    PC_R[i,j]=B$estimate
    PC_p_R[i,j]=B$p.value
  }
}




#####  Save PPC 
setwd("Path/Results")  
write.table(PC_S, file="Pearson correlation for sensitive genes.txt") 
write.table(PC_R, file="Pearson correlation for Resistant genes.txt") 
write.table(PC_p_S, file="p value of Pearson correlation for sensitive genes.txt") 
write.table(PC_p_R, file="p value of Pearson correlation for Resistant genes.txt") 

edge_S=matrix(NA,nrow=Netsize,ncol=Netsize)
edge_R=matrix(NA,nrow=Netsize,ncol=Netsize)
rownames(edge_S)=NetGene; colnames(edge_S)=NetGene; 
rownames(edge_R)=NetGene; colnames(edge_R)=NetGene; 


for (i in 1:Netsize){
  for (j in 1:Netsize){
    if (PC_S[i,j]>0.75 & PC_p_S[i,j]<0.05 & PPI[i,j]==1)
      edge_S[i,j]=1
    else if (PC_S[i,j]<-0.75 & PC_p_S[i,j]<0.05 & PPI[i,j]==1)
      edge_S[i,j]=-1
    else
      edge_S[i,j]=0
  }
}


for (i in 1:Netsize){
  for (j in 1:Netsize){
    if (PC_R[i,j]>0.75 & PC_p_R[i,j]<0.05 & PPI[i,j]==1)
      edge_R[i,j]=1
    else if (PC_R[i,j]<-0.75 & PC_p_R[i,j]<0.05 & PPI[i,j]==1)
      edge_R[i,j]=-1
    else
      edge_R[i,j]=0
  }
}


sum(edge_S!=0)   # 620
sum(edge_R!=0)   # 206

### Select genes in the network
A=rownames(edge_S[rowSums(abs(edge_S))!=0, ])
B=colnames(edge_S[,colSums(abs(edge_S))!=0])
Net_Gene=union(A,B)

edge_S=edge_S[Net_Gene,Net_Gene]
edge_R=edge_R[Net_Gene,Net_Gene]
dim(edge_S)  #61
dim(edge_R)  #61



write.table(edge_S, file="edge matrix for Sensitive genes.txt",col.names = NA,sep = "\t") 
write.table(edge_R, file="edge matrix for Resistant genes.txt",col.names = NA,sep = "\t") 


#### save gene expression data in the initial network
Node_gene_S=Gene_S[Net_Gene,]
Node_gene_R=Gene_R[Net_Gene,]

setwd("Path/Results") 
write.table(Node_gene_S, file="Node_S.xls",col.names = NA,sep = "\t") 
write.table(Node_gene_R, file="Node_R.xls",col.names = NA,sep = "\t") 

