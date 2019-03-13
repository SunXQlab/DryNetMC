
# Single-cell simulation and data analysis



DNG=read.csv("Path/Survival analysis/Data/Quantifing each node.csv")  # Differential network genes
Score=as.matrix(DNG[,c(1,9)])  # Importance Score
Nodes=Score[,1]
setwd("Path/Results")  

E_D=read.table("Path/Results/Differential Network Coefficient.txt")
E_D=as.matrix(E_D)
BB=sign(E_D[Nodes,Nodes])

# install.packages("pracma")  # Manually install
library('pracma')

num=10000
h=0.1
time=seq(0,48,h)
sc_GeneExp=matrix(0,dim(BB)[1],num)
for (ii in 1:num)
{
  x=matrix(0,dim(BB)[1],length(time)+1)
  # CC=randn(42,42)*BB*0.01
  CC=E_D[Nodes,Nodes]*(1+randn(42,42)*0.001)
  C=as.vector(rand(42,1)*(-0.01))
  noise=as.vector(randn(42,1)*0.01)
  # f <- function(t,xx) 
  # {
  #   dx = (CC%*%xx)+C*xx+noise  # *(1-xx/10)*xx
  #   return(dx)
  # }
  # a=0
  # b=48
  # y0=as.vector(rand(1,dim(BB)[1]))
  # n=100
  # sol <- rk4sys(f, a, b, y0, n)
  # sc_GeneExp[,ii]=sol$y[n+1,]
  ind=1
  x[,1]=as.vector(rand(1,dim(BB)[1]))
  
  for (t in time)
  {
    x[,ind+1]=x[,ind]+h*(CC%*%x[,ind]+C*x[,ind])+noise*sqrt(h)
    x[x<0]=0
    x[x>100]=100
    ind=ind+1
  }
  sc_GeneExp[,ii]=x[,ind]
  
  
  
}

# sol$y[n,]
# sol$y[n-1,]
sc_GeneExp[,ii]
sc_GeneExp[,ii-1]
hist(sc_GeneExp[10,])

Gene=sc_GeneExp
###########################

### DM

# # source("http://bioconductor.org/biocLite.R")
# 
# # biocLite("destiny")
# # biocLite("rgl")
# # install.packages("stringi")
# 
# library(stringi)
# library(destiny)
# library(rgl)
# library(Biobase)

# dm=DiffusionMap(t(Gene),n_eigs = 2)
# dm=eigenvectors(dm)
# dev.new()
# plot(dm,col = "#00FF66FF", type = "p", pch = 16,cex=0.5)

dm=DiffusionMap(t(Gene))
dpt <- DPT(dm)
dev.new()
plot(dpt, col_by = 'branch')

dev.new()
plot(dpt, root = 3, paths_to = c(1,2), col_by = 'branch')


dev.new()
plot(dpt, col_by = 'branch', divide = 2, dcs = c(-1,-3,2), pch = 20)


#####################
#########tSNE
# install.packages("Rtsne")
library(Rtsne)

# whole genome data
Gene=sc_GeneExp

tsne = Rtsne(t(Gene), dims = 2, perplexity=220, theta=0,verbose=TRUE, max_iter = 50)

dev.new()
plot(tsne$Y, col = "#00FF66FF", type = "p", pch = 16,cex=1) #col = "#00FF66FF",

####  Heatmap

# install.packages("gplots") #下载gplots程序包
library(gtools) #加载gplots程序
library(gplots) #加载gplots程序

x=sc_GeneExp
rownames(x)=Nodes

x_Normalized=(x-apply(x,1,min))/(apply(x,1,max)-apply(x,1,min))
x=x_Normalized

dev.new()
heatmap.2(x,col=cm.colors(255),distfun =function(x) dist(x,method ='euclidean'),hclustfun =function(x) hclust(x,method ='ward'), dendrogram = "row", scale="row", Rowv=T,Colv=F,key=T,keysize=1.5,trace="none",cexCol=0.5,cexRow=0.35,na.color=par("bg"))

# dev.new()
# heatmap.2(t(x),col=cm.colors(255),distfun =function(x) 1-cor(t(x)),hclustfun =function(x) hclust(t(x),method ='ward'), dendrogram = "row", scale="row", Rowv=T,Colv=F,key=T,keysize=1.5,trace="none",cexCol=0.5,cexRow=0.5,na.color=par("bg"))
# 
# 
dev.new()
heatmap.2(x,col=bluered(255),distfun = dist,hclustfun =function(x) hclust(x,method ='ward'),scale="row", Rowv=T,Colv=NA,dendrogram="row",key=T,keysize=1.5,trace="none",cexCol=0.5,cexRow=0.45,na.color=par("bg"))

# 
# hv1=heatmap.2(t(x),col=redblue(100),scale="none",Rowv=NA,Colv=T,dendrogram="col",key=T,keysize=1.5,trace="none",cexCol=0.5,cexRow=1)
# hv2=heatmap.2(x,distfun = dist,hclustfun = hclust,scale="row",na.rm=TRUE,col=redblue(100),Rowv=NA,Colv=T,dendrogram="col",key=T,keysize=1.5,trace="none",cexCol=0.5,cexRow=1)
# 
# dev.off

library(pheatmap)
dev.off()

dev.new()
pheatmap(x,cluster_row=T, cluster_cols=T, clustering_distance_rows='euclidean',clustering_method = "ward", cellwidth = 0.35, cellheight = 5,color = colorRampPalette(c("CornflowerBlue", "white", "firebrick3"))(50), fontsize=9, fontsize_row=6,labRow=NA) #自定义颜色


