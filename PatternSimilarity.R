

install.packages("vioplot")
library(vioplot)

#######  Select DEGs from RNA seq Data
setwd("F:/Files from office/胶质瘤诱导分化组学数据分析/Data")  
Gene_RPKM=read.csv("F:/Files from office/胶质瘤诱导分化组学数据分析/Data/GBM_gene_RPKM.csv")
Gene_RPKM=na.omit(Gene_RPKM)
Gene_S=Gene_RPKM[,2:6]  ## Sensitive cell line
Gene_R=Gene_RPKM[,7:11] ## Resistant cell line
Gene_U=Gene_RPKM[,12:16] ## test cell line
Gene_S=matrix(as.numeric(as.matrix(Gene_S)),dim(Gene_S)[1],dim(Gene_S)[2])
Gene_R=matrix(as.numeric(as.matrix(Gene_R)),dim(Gene_R)[1],dim(Gene_R)[2])
Gene_U=matrix(as.numeric(as.matrix(Gene_U)),dim(Gene_U)[1],dim(Gene_U)[2])

rownames(Gene_S)=Gene_RPKM[,1]
rownames(Gene_R)=Gene_RPKM[,1]
rownames(Gene_U)=Gene_RPKM[,1]

###### DEG of sensitive cells
timestart<-Sys.time()
DEG_S=c()
RNames_DEG_S=c()
for (i in 1:dim(Gene_S)[1]) 
{ for (j in 1:5)
{ for (k in 1:max(j-1,1))
{if ((max(Gene_S[i,])>=10) & (Gene_S[i,k]>0) & (Gene_S[i,j]>0) & (Gene_S[i,j]/Gene_S[i,k]>5 || Gene_S[i,j]/Gene_S[i,k]<0.2)){
  DEG_S=rbind(DEG_S,Gene_S[i,])
  RNames_DEG_S=rbind(RNames_DEG_S,rownames(Gene_S)[i])
}
}
}
}
rownames(DEG_S)=RNames_DEG_S
DEG_S=unique(DEG_S)
# DEG_S=DEG_S[-1,]
timeend<-Sys.time()
runningtime<-timeend-timestart
print(runningtime)


summary(DEG_S)
dim(DEG_S)  #163,5

plot(c(0,6,12,24,48),DEG_S["STC1",],type="b",main="DEG_S",sub='',xlim=c(0,50), ylim=c(min(DEG_S["STC1",]), max(DEG_S["STC1",])),xlab="Time (Hours)",ylab='Gene Expression')


####  DEG of resistance 

timestart<-Sys.time()
DEG_R=c()
RNames_DEG_R=c()
for (i in 1:dim(Gene_R)[1]) 
{ for (j in 1:5)
{ for (k in 1:max(j-1,1))
{if ((max(Gene_R[i,])>=10) & (Gene_R[i,k]>0) & (Gene_R[i,j]>0) & (Gene_R[i,j]/Gene_R[i,k]>5 || Gene_R[i,j]/Gene_R[i,k]<0.2)){
  DEG_R=rbind(DEG_R,Gene_R[i,])
  RNames_DEG_R=rbind(RNames_DEG_R,rownames(Gene_R)[i])
}
}
}
}
rownames(DEG_R)=RNames_DEG_R
DEG_R=unique(DEG_R)
# DEG_S=DEG_S[-1,]
timeend<-Sys.time()
runningtime<-timeend-timestart
print(runningtime)


summary(DEG_R)
dim(DEG_R)  #101,5

# plot(c(0,6,12,24,48),DEG_R["STC1",],type="b",main="DEG_R",sub='',xlim=c(0,50), ylim=c(min(DEG_R["STC1",]), max(DEG_R["STC1",])),xlab="Time (Hours)",ylab='Gene Expression')

########## 



#####  Save DEGs 
# setwd("F:/胶质瘤诱导分化组学数据分析/Code-cluster_based/Results")  
# write.table(DEG_S, file="Drug Sensitive Genes.txt") 
# write.table(DEG_R, file="Drug Resistant Genes.txt") 
# 
# 
# DEG_S=read.table("F:/胶质瘤诱导分化组学数据分析/Code-cluster_based/Results/Drug Sensitive Genes.txt")
# DEG_R=read.table("F:/胶质瘤诱导分化组学数据分析/Code-cluster_based/Results/Drug Resistant Genes.txt")



names_SR=union(rownames(DEG_S),rownames(DEG_R))

setwd("F:/Files from office/胶质瘤诱导分化组学数据分析/Code-cluster_based/Results") 
Gene_Annotation=read.csv("Gene Annotation_ranked.csv",sep = ",") 
GeneName=Gene_Annotation[,1]

# TC_gene_1=Gene_S[names_SR,]  # Temporally changing genes of sensitivity
# TC_gene_2=Gene_R[names_SR,]  # Temporally changing genes of resistance
# 
# TC_gene_3=Gene_U[names_SR,]  # Temporally changing genes of test cell line

TC_gene_1=Gene_S[GeneName,]  # Temporally changing genes of sensitivity
TC_gene_2=Gene_R[GeneName,]  # Temporally changing genes of resistance
TC_gene_3=Gene_U[GeneName,]  # Temporally changing genes of test cell line

 setwd("F:/胶质瘤诱导分化组学数据分析/Code-cluster_based/Results")  
 write.table(TC_gene_1, file="TC_gene_1.txt",quote=F,row.names=F,col.names = F,sep = "\t") 
 write.table(TC_gene_2, file="TC_gene_2.txt",quote=F,row.names=F,col.names = F,sep = "\t") 
 write.table(TC_gene_3, file="TC_gene_3.txt",quote=F,row.names=F,col.names = F,sep = "\t") 


### Euclid similarity measurement
D_1=sqrt(rowSums((TC_gene_3-TC_gene_1)^2))
D_2=sqrt(rowSums((TC_gene_3-TC_gene_2)^2))

wilcox.test(D_1,D_2,alternative="less", paired=T)  # Wilcoxon signed rank test with continuity correction;  p-value =0.01625
vioplot(D_1,D_2,names=c("Distance to Sensitive cells","Distance to Resistant cells"),col=c("green","brown"))
title("Temporal Pattern Similarity_Euclid Distance")

temp <- locator(1) # 在图表上，你喜欢的地方点击一下，文字就出来了
text(temp,"p-value = 0.01625")

### Manhattan Distance
D_1=rowSums(abs(TC_gene_3-TC_gene_1)^1)
D_2=rowSums(abs(TC_gene_3-TC_gene_2)^1)

wilcox.test(D_1,D_2,alternative="less", paired=T)  # Wilcoxon signed rank test with continuity correction;  p-value = 0.02675
vioplot(D_1,D_2,names=c("Distance to Sensitive cells","Distance to Resistant cells"),col=c("green","brown"))
title("Temporal Pattern Similarity_Manhattan Distance")

temp <- locator(1) # 在图表上，你喜欢的地方点击一下，文字就出来了
text(temp,"p-value = 0.02675")

### Slope Distance
Time=t(matrix(c(6,6,12,24),dim(TC_gene_1)[2]-1,dim(TC_gene_1)[1]))
Slope_1=t(diff(t(TC_gene_1)))/Time
Slope_2=t(diff(t(TC_gene_2)))/Time
Slope_3=t(diff(t(TC_gene_3)))/Time
  
D_1=sqrt(rowSums((Slope_3-Slope_1)^2))
D_2=sqrt(rowSums((Slope_3-Slope_2)^2))

wilcox.test(D_1,D_2,alternative="less", paired=T)  # Wilcoxon signed rank test with continuity correction;  p-value =0.01625
vioplot(D_1,D_2,names=c("Distance to Sensitive cells","Distance to Resistant cells"),col=c("green","brown"))
title("Temporal Pattern Similarity")

### DTW, dynamic time wrapping 
install.packages("dtw")
library(dtw)

for (i in 1:42)
{
  D1[i]=dtw(TCG3[i,],TCG1[i,])$distance
  D2[i]=dtw(TCG3[i,],TCG2[i,])$distance
 
}

# Or ## calculate DTW distance in matlab "PatternSimilarity" and read the results from there saved
# D1=read.table("F:/胶质瘤诱导分化组学数据分析/Code-cluster_based/Results/DTWDistance1.txt")
# D2=read.table("F:/胶质瘤诱导分化组学数据分析/Code-cluster_based/Results/DTWDistance2.txt")
# D1=as.numeric(D1)
# D2=as.numeric(D2)

wilcox.test(D1,D2,alternative="less", paired=T)  # Wilcoxon signed rank test with continuity correction;  p-value = 0.01523
vioplot(D1,D2,names=c("DTW Distance to Sensitive cells","DTW Distance to Resistant cells"),col=c("green","brown"))
title("Temporal Pattern Similarity_DTW Distance")

temp <- locator(1)
text(temp,"p-value = 0.01523")


##########  Heatmap of gene expression in the differential network

DNG=read.csv("F:/Files from office/胶质瘤诱导分化组学数据分析/Code-cluster_based/Survival analysis/Data/Quantifing each node.csv")  # Differential network genes
Score=as.matrix(DNG[,c(1,9)])  # Importance Score
Nodes=Score[,1]

x=cbind(TC_gene_1,TC_gene_2,TC_gene_3)
rownames(x)=Nodes
colnames(x)=c('S0h','S6h','S12h','S24h','S48h','R0h','R6h','R12h','R24h','R48h','U0h','U6h','U12h','U24h','U48h')

x_Normalized=(x-apply(x,1,min))/(apply(x,1,max)-apply(x,1,min)+1e-3)
x=x_Normalized
# x[is.nan(x)]==0

library(pheatmap)
dev.off()

dev.new()
pheatmap(x,cluster_row=T, cluster_cols=F, clustering_distance_rows='euclidean',clustering_method = "ward", cellwidth = 10, cellheight = 8,color = colorRampPalette(c("CornflowerBlue", "white", "firebrick3"))(50), fontsize=9, fontsize_row=6,labRow=NA) #自定义颜色

dev.new()
pheatmap(t(x),cluster_row=F, cluster_cols=T, clustering_distance_rows='euclidean',clustering_method = "ward", cellwidth = 10, cellheight = 8,color = colorRampPalette(c("CornflowerBlue", "white", "firebrick3"))(50), fontsize=9, fontsize_row=6,labRow=NA) #自定义颜色


