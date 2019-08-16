########## Existing DEG detection methods
###DEseq2
source ("http://bioconductor.org/biocLite.R")
biocLite("DESeq2") 
install.packages("backports")
library(DESeq2) # 加载DESeq2包

setwd("F:/Files from office/胶质瘤诱导分化组学数据分析/")
Gene_RPKM=read.csv("Data/GBM_RNA-seq_countTab_v2.csv")
Gene_RPKM=na.omit(Gene_RPKM)
Gene_S=Gene_RPKM[,2:6]  ## Sensitive cell line
Gene_R=Gene_RPKM[,7:11] ## Resistant cell line
Gene_S=matrix(as.numeric(as.matrix(Gene_S)),dim(Gene_S)[1],dim(Gene_S)[2])
Gene_R=matrix(as.numeric(as.matrix(Gene_R)),dim(Gene_R)[1],dim(Gene_R)[2])
rownames(Gene_S)=Gene_RPKM[,1]
rownames(Gene_R)=Gene_RPKM[,1]
colnames(Gene_S)=colnames(Gene_RPKM)[2:6] 
colnames(Gene_R)=colnames(Gene_RPKM)[7:11]

countData=cbind(Gene_S,Gene_R)

condition <- factor(c("S","S","S","S","S","R","R","R","R","R")) # 定义condition
colData <- data.frame(row.names=colnames(countData), condition) # 样品信息矩阵

dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), design= ~ condition ) # 构建dds矩阵
head(dds) # 查看dds矩阵的前6行

dds <- DESeq(dds) # 对dds进行Normalize
resultsNames(dds) # 查看结果的名称
res <- results(dds) # 使用results()函数获取结果，并赋值给res
head (res, n=5) # 查看res矩阵的前5行

mcols(res,use.names= TRUE) # 查看res矩阵每一列的含义
summary(res) # 对res矩阵进行总结

table(res$padj<0.05) # 取padj小于0.05的数据
res <- res[order(res$padj),] # 按照padj的大小将res重新排列
diff_gene_deseq2 <- subset(res,padj < 0.05 & (log2FoldChange >1 | log2FoldChange < -1)) # 获取padj小于0.05，表达倍数取以2为对数后绝对值大于1的差异表达基因，赋值给diff_gene_deseq2
head (diff_gene_deseq2, n=5) # 查看diff_gene_deseq2矩阵的前5行


DEG_FC_order=diff_gene_deseq2[order(abs(diff_gene_deseq2[,2]),decreasing=T),]
head(DEG_FC_order)

Gene_S=Gene_S[intersect(rownames(Gene_S),rownames(DEG_FC_order)),]
Gene_R=Gene_R[intersect(rownames(Gene_R),rownames(DEG_FC_order)),]

# Top 5: OCIAD2, FLG, SPP1, SOX3, NEFL
DEG=rownames(DEG_FC_order)[1:5]




# install.packages("vioplot")
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

## 

GeneName=c('KIF2C', 'CCNA2', 'NDC80', 'KIF11', 'KIF23')

# GeneName=c('OCIAD2', 'FLG', 'SPP1', 'SOX3', 'NEFL')

TC_gene_1=Gene_S[GeneName,]  # Temporally changing genes of sensitivity
TC_gene_2=Gene_R[GeneName,]  # Temporally changing genes of resistance
TC_gene_3=Gene_U[GeneName,]  # Temporally changing genes of test cell line

# setwd("F:/Files from office/胶质瘤诱导分化组学数据分析/Code-cluster_based/Results")
# write.table(TC_gene_1, file="TC_gene_1.txt",quote=F,row.names=F,col.names = F,sep = "\t")
# write.table(TC_gene_2, file="TC_gene_2.txt",quote=F,row.names=F,col.names = F,sep = "\t")
# write.table(TC_gene_3, file="TC_gene_3.txt",quote=F,row.names=F,col.names = F,sep = "\t")


### DTW, dynamic time wrapping 
# install.packages("dtw")
# install.packages("proxy")
library(proxy)
library(dtw)

D_1=matrix(0,5,1)
D_2=matrix(0,5,1)
for (i in 1:5)
{
  DTW=dtw(TC_gene_1[i,],TC_gene_3[i,])
  D_1[i,]=DTW$distance
  DTW=dtw(TC_gene_2[i,],TC_gene_3[i,])
  D_2[i,]=DTW$distance
}


wilcox.test(D_1,D_2,alternative="greater", paired=T)  # Wilcoxon signed rank test with continuity correction;  p-value = 0.01523
# vioplot(D1,D2,names=c("DTW Distance to Sensitive cells","DTW Distance to Resistant cells"),col=c("green","brown"))
# title("Temporal Pattern Similarity_DTW Distance")

dataset <- data.frame(value = c(D_2,D_1), group = factor(rep(c("Distance to resistant cells","Distance to sensitive cells"), times = c(length(D_1), length(D_2)))))
dev.new()
boxplot( t(value) ~ t(group),  notch = F, dataset, outline = FALSE, border = c( "#009E73","purple"),cex = 1,cex.axis=1,pars = list(boxwex = 0.5, staplewex = 0.5, outwex = 0.5))  #,col.axis = "#009E73"
temp <- locator(1) # 在图表上，你喜欢的地方点击一下，文字就出来了
text(temp,"p-value = 0.03125")
title("Temporal Pattern Similarity_DTW Distance")


D0=mean(D_1-D_2)
## Random 5 genes for 100 times

D=matrix(0,1,1000)
D1_mean=matrix(0,1,1000)
D2_mean=matrix(0,1,1000)

ind=1
while (ind <= 1000)
{
  GeneName=rownames(DEG_FC_order)[sample(1:dim(as.matrix(DEG_FC_order))[1],5)]
  if (length(intersect(GeneName,rownames(Gene_S)))==5)
  {
  TC_gene_1=Gene_S[GeneName,]  # Temporally changing genes of sensitivity
  TC_gene_2=Gene_R[GeneName,]  # Temporally changing genes of resistance
  TC_gene_3=Gene_U[GeneName,]  # Temporally changing genes of test cell line
  
  D_1=matrix(0,5,1)
  D_2=matrix(0,5,1)
  for (j in 1:5)
  {
    DTW=dtw(TC_gene_1[j,],TC_gene_3[j,])
    D_1[j,]=DTW$distance
    DTW=dtw(TC_gene_2[j,],TC_gene_3[j,])
    D_2[j,]=DTW$distance
  }
  D1_mean[ind]=mean(D_1)
  D2_mean[ind]=mean(D_2)
  D[ind]=mean(D_1-D_2)
  ind=ind+1
  }
}
# vioplot 
# install.packages("vioplot")
library(vioplot)
dev.new()
vioplot(t(D1_mean),t(D2_mean),names=c("Mean distance to sensitive cells","Mean distance to resistant cells"),col=c("#009E73","purple"))

dataset <- data.frame(value = c(D1_mean,D2_mean), group = factor(rep(c("Mean distance to sensitive cells","Mean distance to resistant cells"), times = c(length(D1_mean), length(D2_mean)))))
dev.new()
boxplot( t(value) ~ t(group),  notch = F, dataset, outline = FALSE, border = c( "#009E73","purple"),cex = 1,cex.axis=1,pars = list(boxwex = 0.5, staplewex = 0.5, outwex = 0.5))  #,col.axis = "#009E73"

wilcox.test(D1_mean,D2_mean,alternative="two.sided", paired=T)  # Wilcoxon signed rank test with continuity correction;  p-value = 0.01523

### plot Hist 
dev.new()
hist(D,breaks=500, prob=FALSE, col="#7F7FFF", xlab="D_S-D_R", main="",xaxs = "i", yaxs ="i",xlim = c(-1000, 1000))  #,xlim = c(0.4, 0.9),ylim = c(0, 10)
# lines(density(D), col="darkblue", lwd=2)
par(new=TRUE)
lines(rep(D0,200),0:199,col="#FF0000", lwd=3)

p.value<-length(D[D<D0])/1000  # 1e-3

p.value
temp <- locator(1) # 在图表上，你喜欢的地方点击一下，文字就出来了
text(temp,"p-value = 0.025")

### Correlation Network-derived 5 genes

#Performs Gene Sets Net Correlation Analysis (GSNCA) test to detect differentially coexpressed gene sets.
# https://rdrr.io/bioc/GSAR/man/GSNCAtest.html

source("https://bioconductor.org/biocLite.R")
biocLite("GSAR")

library("GSAR")
library(MASS)

Gene_selected=intersect(rownames(Gene_S),rownames(DEG_FC_order[1:1000,]))

gp <- rbind(t(Gene_S[Gene_selected,]),t(Gene_R[Gene_selected,]))
gp=matrix(as.numeric(gp),dim(gp))

for (i in 1:dim(gp)[2])
{
  if ( sd(gp[1:5,i])<1e-3 |  sd(gp[6:10,i])<1e-3)
  {
    gp[,i]=NA
  }
}

gp=t(na.omit(t(gp)))
dim(gp)

# dataset <- aperm(gp, c(2,1))
# ## first 20 samples belong to group 1
# ## second 20 samples belong to group 2
# pvalue <- GSNCAtest(object=dataset, group=c(rep(1,5),rep(2,5)),check.sd=TRUE, min.sd=1e-3, max.skip=100) 


objt=aperm(gp, c(2,1))
## calculate the weight matrix
Wmat <- as.matrix(dist(objt, method = "euclidean", diag = TRUE, upper = TRUE, p = 2))
## create a weighted undirectional graph from the weight matrix
gr <- graph_from_adjacency_matrix(Wmat, weighted = TRUE, mode = "undirected")
## find the minimum spanning tree
MST <- mst(gr)
# MST <- findMST2((objt))
# MST <- findMST2.PPI(gr)
## ranking according to HDP (according to the high directed preorder (HDP) traversal of the tree)
HDP.ranks <- HDP.ranking(MST)

G1=Gene_selected[which(HDP.ranks==c(1))]
G2=Gene_selected[which(HDP.ranks==c(2))]
G3=Gene_selected[which(HDP.ranks==c(3))]
G4=Gene_selected[which(HDP.ranks==c(4))]
G5=Gene_selected[which(HDP.ranks==c(5))]

Gene_GSNCA_5=c(G1,G2,G3,G4,G5)
##  "C5orf58"       "COL4A3"        "MUC16"         "DKFZp686O1327" "GAL"  

GeneName=Gene_GSNCA_5
TC_gene_1=Gene_S[GeneName,]  # Temporally changing genes of sensitivity
TC_gene_2=Gene_R[GeneName,]  # Temporally changing genes of resistance
TC_gene_3=Gene_U[GeneName,]  # Temporally changing genes of test cell line

### DTW, dynamic time wrapping 
D_1=matrix(0,5,1)
D_2=matrix(0,5,1)
for (i in 1:5)
{
  DTW=dtw(TC_gene_1[i,],TC_gene_3[i,])
  D_1[i,]=DTW$distance
  DTW=dtw(TC_gene_2[i,],TC_gene_3[i,])
  D_2[i,]=DTW$distance
}


wilcox.test(D_1,D_2,alternative="two.sided", paired=T)  # Wilcoxon signed rank test with continuity correction;  p-value = 0.01523
# vioplot(D1,D2,names=c("DTW Distance to Sensitive cells","DTW Distance to Resistant cells"),col=c("green","brown"))
# title("Temporal Pattern Similarity_DTW Distance")

dataset <- data.frame(value = c(D_2,D_1), group = factor(rep(c("Distance to resistant cells","Distance to sensitive cells"), times = c(length(D_1), length(D_2)))))
dev.new()
boxplot( t(value) ~ t(group),  notch = F, dataset, outline = FALSE, border = c( "#009E73","purple"),cex = 1,cex.axis=1,pars = list(boxwex = 0.5, staplewex = 0.5, outwex = 0.5))  #,col.axis = "#009E73"
temp <- locator(1) # 在图表上，你喜欢的地方点击一下，文字就出来了
text(temp,"p-value = 1")
title("Temporal Pattern Similarity_DTW Distance")
