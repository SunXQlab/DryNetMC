


########### Differential Network #############

setwd("F:/胶质瘤诱导分化组学数据分析/Code-cluster_based/Results") 
Net_Gene=read.csv('Net_Gene.csv',fill= T,header = F)
Edge_S=read.table("F:/胶质瘤诱导分化组学数据分析/Code/ElasticNet/Significant Coefficient_Sensitive Network.xls")
Edge_R=read.table("F:/胶质瘤诱导分化组学数据分析/Code/ElasticNet/Significant Coefficient_Resistant Network.xls")


E_S=E_S[,-c(1,2)]
E_R=E_R[,-c(1,2)]
# rownames(E_S)=t(Net_Gene)
rownaems(E_R)=t(Net_Gene)
# colnames(E_S)=t(Net_Gene)
colnames(E_R)=t(Net_Gene)

E_S_matrix=matrix(as.numeric(as.matrix(E_S)),dim(E_S)[1],dim(E_S)[2])
rownames(E_S_matrix)=rownames(E_S)
colnames(E_S_matrix)=colnames(E_S)



E_R_matrix=matrix(as.numeric(as.matrix(E_R)),dim(E_R)[1],dim(E_R)[2])
rownames(E_R_matrix)=rownames(E_R)
colnames(E_R_matrix)=colnames(E_R)

E_D=matrix(0,dim(E_R_matrix)[1],dim(E_R_matrix)[2])
rownames(E_D)=t(Net_Gene)
colnames(E_D)=t(Net_Gene)




for (i in 1:dim(E_R_matrix)[1])
{
  for (j in 1:dim(E_R_matrix)[2])
  {
    if (E_R_matrix[i,j]!=0 & E_S_matrix[i,j]==0 & abs(E_R_matrix[i,j])>=0.01)
    {E_D[i,j]=E_R_matrix[i,j]}
  }
}

dim(E_D)

setwd("F:/胶质瘤诱导分化组学数据分析/Code-cluster_based/Results")  
write.table(E_D, file="Differential Network Coefficient.txt",col.names = NA,sep = "\t")

### Reform txt data for cytoscape 
E_D=read.table("Differential Network Coefficient.txt")
E_D=as.matrix(E_D)
BB=E_D

C=matrix(NA,nrow=sum(sum(BB!=0)),ncol=3)
row=1
for (i in 1:dim(BB)[1])
{
  for (j in 1:dim(BB)[1])
  {
    if (BB[i,j]!=0)
    {
      C[row,1]=rownames(BB)[i]
      C[row,3]=colnames(BB)[j]
      C[row,2]=BB[i,j]
      row=row+1
    }
  }
}

write.table(C, file="Cytoscape_Differential Network.txt",quote=F,row.names=F,col.names = F,sep = "\t")


########################### Network topology analysis ###############################

library(igraph) # Load the igraph package
setwd("F:/胶质瘤诱导分化组学数据分析/Code-cluster_based/Results")  
C=read.table("Cytoscape_Differential Network.txt")
Net_Diff=as.data.frame(cbind(as.character(C[,1]),as.character(C[,3]),C[,2]))  # Input data for R analysis, data.frame

gene=NULL
ind=1
for (i in variable)#)
{
  if (sum(E_D[i,]==0)!=61 | (sum(E_D[,i]==0)!=61)) 
  {
    gene[ind]=rownames(E_D)[i]
    ind=ind+1
  }
}


net=graph_from_data_frame(Net_Diff, directed = T,vertices = NULL)  # 
triad_census(net) # for directed networks

V(net)  #
E(net)



#Diameter
diam <- get_diameter(net, directed=T)

# Degree
deg <- degree(net, mode="all")
# plot(net, vertex.size=deg*3)
hist(deg, breaks=1:max(deg), main="Histogram of node degree")
deg.dist <- degree_distribution(net, cumulative=T, mode="all")
plot( x=0:max(deg), y=1-deg.dist, pch=19, cex=1.2, col="orange",
      xlab="Degree", ylab="Cumulative Frequency")

## Centrality & centralization

degree(net, mode="in")
degree(net, mode="out")
centr_degree(net, mode="in", normalized=T)


eigen_centrality(net, directed=T, weights=NA)
centr_eigen(net, directed=T, normalized=T)

## Hubs and Authorities
net=graph_from_data_frame(Net_Diff, directed = T, vertices = NULL)  # network toplogy without direction

hs <- hub_score(net, weights=T)$vector
which(hs==max(hs)) # PLK1

as <- authority_score(net, weights=NA)$vector
which(as==max(as)) # PLK1

par(mfrow=c(1,2))
plot(net, vertex.size=hs*50, main="Hubs")
plot(net, vertex.size=as*30, main="Authorities")

hs_sort=sort(hs, decreasing = T)
as_sort=sort(as, decreasing = T)

# cliques
net.sym <- as.undirected(net, mode= "each",
                         edge.attr.comb=list(weight="sum", "ignore"))

cliques(net.sym) # list of cliques
sapply(cliques(net.sym), length) # clique sizes
largest_cliques(net.sym) # cliques with max number of nodes
vcol <- rep("grey80", vcount(net.sym))
vcol[unlist(largest_cliques(net.sym))] <- "gold"
plot(as.undirected(net.sym), vertex.label=V(net.sym)$name, vertex.color=vcol)

## Community detection

cfg <- cluster_fast_greedy(as.undirected(net))
plot(cfg, as.undirected(net))
cfg$membership #label the group ID of each node
# V(net)$community <- cfg$membership
# colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)
# plot(net, vertex.color=colrs[V(net)$community])

ceb <- cluster_edge_betweenness(net)
dendPlot(ceb, mode="hclust")
plot(ceb, net)

## neighbor

neighbors(net, "ANLN", mode="all")


############## Adaptation dynamics

setwd("F:/胶质瘤诱导分化组学数据分析/Code-cluster_based/Results") 

Node_S=read.table("Node_S.xls") 
Node_R=read.table("Node_R.xls") 


rownames(Node_S)=Node_S[,1]
Node_S=Node_S[-1,-1]

rownames(Node_R)=Node_R[,1]
Node_R=Node_R[-1,-1]

Node_S=matrix(as.numeric(as.matrix(unlist(Node_S))),dim(Node_S)[1],dim(Node_S)[2])
Node_R=matrix(as.numeric(as.matrix(Node_R)),dim(Node_R)[1],dim(Node_R)[2])

## Sensitive AI
AI_S=NULL
AI_S=abs((Node_S[V(net_S),5]-Node_S[V(net_S),1])/Node_S[V(net_S),1])
AI_S=as.numeric(AI_S)
names(AI_S)=t(rownames(as.matrix(V(net_S))))

## Resistant AI
AI_R=NULL
AI_R=abs((Node_R[V(net_R),5]-Node_R[V(net_R),1])/Node_R[V(net_R),1])
AI_R=as.numeric(AI_R)
names(AI_R)=t(rownames(as.matrix(V(net_R))))

hist(AI_S)
hist(AI_R)
wilcox.test(AI_S,AI_R,alternative="greater", paired=F)  # p-value = 1.029e-15

install.packages("vioplot")
library(vioplot)
vioplot(AI_S,AI_R,names=c("Sensitive Network","Resistant Network"),col="gold",ylab="Relative Response")
title("Dynamic changes in gene expressions")

## Differntial AI
AI=NULL
AI=t(abs((Node_S[V(net),5]-Node_S[V(net),1])/(Node_R[V(net),5]-Node_R[V(net),1])))*(Node_R[V(net),1]/Node_S[V(net),1])
AI=as.numeric(AI)
names(AI)=t(rownames(as.matrix(V(net))))

AI_sort=sort(AI,decreasing=T)


############### Entropy of each node
net=graph_from_data_frame(Net_Diff, directed = T, vertices = NULL) 

S=matrix(0,1,length(V(net)))

for (i in 1:length(V(net)))
{
  neigh_out <- neighbors(net, V(net)[i], mode="out")
  neigh_in <- neighbors(net, V(net)[i], mode="in")
  aij_out=E_D[names(V(net))[i],names(neigh_out)]
  aij_in=E_D[names(neigh_in),names(V(net))[i]]
  # aij=c(aij_out,aij_in)
  aij=aij_out
  S[i]=-sum(abs(aij)/sum(abs(aij))*log(abs(aij)/sum(abs(aij))))#/log(deg[V(net)[i]])
}
colnames(S)=names(V(net))
S_sort=sort(S,decreasing=T)

#### Entropy of node in sensitive network
library(igraph) # Load the igraph package
setwd("F:/胶质瘤诱导分化组学数据分析/Code-cluster_based/Results") 
Net_Gene=read.csv('Net_Gene.csv',fill= T,header = F)
E_S=read.table("Coefficient_GeneNetwork_Sensitive.txt",header = T,sep="\t",fileEncoding = "UTF-8",row.names=1)
E_R=read.table("Significant Coefficient_Resistant Network.xls",header = F,sep="\t",fileEncoding = "UTF-8")


# rownames(E_S)=t(Net_Gene)
rownames(E_R)=t(Net_Gene)
# colnames(E_S)=t(Net_Gene)
colnames(E_R)=t(Net_Gene)

E_S_matrix=matrix(as.numeric(as.matrix(E_S)),dim(E_S)[1],dim(E_S)[2])
rownames(E_S_matrix)=rownames(E_S)
colnames(E_S_matrix)=colnames(E_S)



E_R_matrix=matrix(as.numeric(as.matrix(E_R)),dim(E_R)[1],dim(E_R)[2])
rownames(E_R_matrix)=rownames(E_R)
colnames(E_R_matrix)=colnames(E_R)



setwd("F:/胶质瘤诱导分化组学数据分析/Code-cluster_based/Results")  
C_S=read.table("Cytoscape_Net_Sensitive_New.txt")
Net_S=as.data.frame(cbind(as.character(C_S[,1]),as.character(C_S[,3]),C_S[,2]))  # Input data for R analysis, data.frame
net_S=graph_from_data_frame(Net_S, directed = T,vertices = NULL) 
C_S=as.matrix(C_S,dim(C_S)[1],dim(C_S)[2])
S_S=matrix(0,1,dim(E_S_matrix))
deg <- degree(net_S, mode="out")

for (i in 1:dim(E_S_matrix)[1])
{
  neigh_out <- colnames(E_S_matrix)[E_S_matrix[i,]!=0]
  neigh_in <- rownames(E_S_matrix)[E_S_matrix[,i]!=0]
  aij_out=E_S_matrix[i,neigh_out]
  aij_in=E_S_matrix[neigh_in,i]
  aij=as.numeric(c(aij_out,aij_in))
  if (sum(abs(aij))!=0) {
    S_S[i]=-sum(abs(aij)/sum(abs(aij))*log2(abs(aij)/sum(abs(aij))))#/log(sum(abs(aij)!=0))
  }
}
names(S_S)=rownames(E_S_matrix)

setwd("F:/胶质瘤诱导分化组学数据分析/Code-cluster_based/Results")  
C_R=read.table("Cytoscape_Net_Resistant_new.txt")
Net_R=as.data.frame(cbind(as.character(C_R[,1]),as.character(C_R[,3]),C_R[,2]))  # Input data for R analysis, data.frame
net_R=graph_from_data_frame(Net_R, directed = T,vertices = NULL) 
C_R=as.matrix(C_R,dim(C_R)[1],dim(C_R)[2])
S_R=matrix(0,1,length(V(net_R)))
id=matrix(0,1,length(V(net_R)))
deg <- degree(net_R, mode="out")
ind=1
for (i in 1:dim(E_R_matrix)[1])
{
  # neigh_out <- neighbors(net_R, V(net_R)[i], mode="out")
  # neigh_in <- neighbors(net_R, V(net_R)[i], mode="in")
  neigh_out <- colnames(E_R_matrix)[E_R_matrix[i,]!=0]
  neigh_in <- rownames(E_R_matrix)[E_R_matrix[,i]!=0]
  aij_out=E_R_matrix[i,neigh_out]
  aij_in=E_R_matrix[neigh_in,i]
  aij=as.numeric(c(aij_out,aij_in))
  if (sum(abs(aij))!=0) {
    S_R[ind]=-sum(abs(aij)/sum(abs(aij))*log2(abs(aij)/sum(abs(aij))))#/log(sum(abs(aij)!=0))
    ind=ind+1
    id[ind]=rownames(E_R_matrix)[i]
  }
  
}

names(S_R)=t(id[-1])

S_delta=S_R[V(net)]-S_S[V(net)]
S_delta_sort=sort(S_delta,decreasing=T)

wilcox.test(S_S,S_R, alternative="less", paired=F)  # p-value = 1.419e-12

vioplot(t(S_S),t(S_R),names=c("Sensitive Network","Resistant Network"),col="green")
title("Entropy")



################### Rank of nodes according to their attributes
hs_rank=rank(hs)
S_rank=rank(S)
S_delta_rank=rank(S_delta)
AI_rank=rank(AI)

Node_rank=hs_rank+S_delta_rank+AI_rank

ranksum_rank=43-rank(Node_rank)

Node_attribute=NULL
Node_attribute=rbind(hs,S_delta,AI,hs_rank,S_delta_rank,AI_rank,Node_rank,ranksum_rank)

write.table(t(Node_attribute), file="Quantifing each node.txt",quote=F,row.names=T,col.names = T,sep = "\t")

DiffNet_Node_S=Node_S[V(net),]
DiffNet_Node_R=Node_R[V(net),]

DiffNet_Node=cbind(DiffNet_Node_S,DiffNet_Node_R)
write.table(DiffNet_Node, file="DiffNet_Node.txt",quote=F,row.names=F,col.names = T,sep = "\t")

cor.test(hs, S_delta,method="spearman")
cor.test(hs_rank, AI_rank)
cor.test(AI_rank, S_delta_rank)

####### Calculating hs
AM=as_adjacency_matrix(net)

AM=E_D[V(net),V(net)]
AM=abs(AM)
dim(AM)
EG1=eigen(AM*t(AM))
which(EG1$values==max(EG1$values))
hs=EG1$vectors[,1]

EG2=eigen(t(AM)*AM)
which(EG2$values==max(EG2$values))
as=EG2$vectors[,1]


