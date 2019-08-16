# DryNetMC
Differential Regulatory Network-based Modeling and Characterization (DryNetMC) for Cancer Drug Resistance 

--Modeling and characterization of the dynamic gene regulatory networks underlying cancer drug resistance based on time course RNA-seq data



--- Work protocol---


The main codes were implemented in R (R version 3.5.1). The pipeline can be run in the following order:

1.	DEG.R

This piece of code performs analyzing conventional differentially expressed genes (DEGs) and temporally changing genes (TCGs) defined in our study, as well as normalization and heat map visualization. The following R packages are required to be installed: “plotrix”, "gplots", and “pheatmap”.  

2.	Clustering_visualization.m

This piece of code performs clustering and visualizing TCGs. 

3.	CorrelationNetwork.R

This piece of code performs constructing initial PPI networks and correlation networks. The following R packages should be installed: "ggm" and "ppcor".

4.	ODENet_Sensitive.R and ODENet_Resistant.R

This piece of code performs dynamic network reconstruction for sensitive cells or resistant cells. The output is a gene interaction matrix which is reformed for Cytoscape visualization. The following R packages should be installed:  "iterators" and " pracma".  


5.	Differential_Network.R 

This piece of code performs differential network analysis, visualization, and characterization (quantifying topological hubs, adaptation dynamics and local network entropy). The following R packages are required to be installed: “igraph” and "vioplot".

6.	FB_Motif_Detection.R

This piece of code detects various types of 2-node and 3-node feedback motifs in the sensitive network and resistant network. 

7.	PatternSimilarity.R

This piece of code performs calculation and comparison of the distance between the tested cell line and the sensitive or resistant cell line. The following R packages are required to be installed: “dtw” and "vioplot".  ‘PatternSimilarity.m’ can be alternatively used.

8.	D2NB_Survival.R

This piece of code performs D2NB identification, differential expression analysis of the identified genes and K-M analysis. The following R packages are required to be installed: "AnnotationDbi", "bit", “org.Hs.eg.db", "lattice", "survival", "reshape2", "data.table", "zoo", "survminer", "survival", "glmnet" and "pROC".  

9.	Validation_Network.R

This piece of code performs validation of the dynamic network reconstruction method based on a synthetic dataset in comparison to the correlation network method. The following R packages are required to be installed: “pracma”, “glmnet” and “ROCR”. 


The case study gene expression data are saved in “GBM_gene_RPKM.csv”. The file path in the codes for loading this dataset should be customized to the working folder.  


