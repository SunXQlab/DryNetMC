%%%%%%%%%%%%%% Patterm Similarity
cd('F:\胶质瘤诱导分化组学数据分析\Code-cluster_based\Results')
TCG1=textread('TC_gene_1.txt');
TCG2=textread('TC_gene_2.txt');
TCG3=textread('TC_gene_3.txt');

D1=zeros(1,42);
D2=zeros(1,42);
cd('F:\胶质瘤诱导分化组学数据分析\Code-cluster_based\ElasticNet')
for  i=1:42
    [Dist,D,k,w]=dtw(TCG3(i,:),TCG1(i,:));
    D1(i)=Dist;
    [Dist,D,k,w]=dtw(TCG3(i,:),TCG2(i,:));
    D2(i)=Dist;
end

save DTWDistance1.txt -ascii D1
save DTWDistance2.txt -ascii D2

[p,h,stats]=ranksum(D1,D2)  %,'alpha',0.01,'tail','left');

[p,h,stats] = signrank(D1,D2) %,'alpha',0.01,'tail','left','method','exact')

[h,p,ks2stat] = kstest2(D1,D2);

hist(D1)
figure,hist(D2)

figure,
boxplot([D1;D2]');% ,'plotstyle','compact')