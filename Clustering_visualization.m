cd('F:\胶质瘤诱导分化组学数据分析\Code-cluster_based\Results')
TCG1=xlsread('TC_gene_1_Normalized.xls');
TCG2=xlsread('TC_gene_2_Normalized.xls');

%%% Interpolation
TCG1=pchip([0 6 12 24 48],TCG1,0:1:48);
%%% Interpolation
TCG2=pchip([0 6 12 24 48],TCG2,0:1:48);

%%%%%%%%%%%%%%%  Sensitive genes
corrDist = pdist(TCG1,'euclidean');
clusterTree = linkage(corrDist,'average');

clusters = cluster(clusterTree,'maxclust',9);

close all
figure
for c = 1:9
    subplot(3,3,c);
    plot(0:1:48,TCG1((clusters == c),:)','Color',[0.8 0.8 0.8],'MarkerSize',2);
    hold on,
    plot([0:1:48],mean(TCG1((clusters == c),:)',2),'k-','LineWidth',3);
    plot([0 6 12 24 48],mean(TCG1((clusters == c),[0 6 12 24 48]+1)',2),'b .','MarkerSize',30);
    axis tight
end
suptitle('Hierarchical Clustering of Sensitive TCG Profiles');


%%%%%%%%%%%%%%%%%% Resistant genes
corrDist = pdist(TCG2,'euclidean');
clusterTree = linkage(corrDist,'average');

clusters = cluster(clusterTree,'maxclust',9);

% close all
figure
for c = 1:9
    subplot(3,3,c);
    plot(0:1:48,TCG2((clusters == c),:)','Color',[0.8 0.8 0.8],'MarkerSize',2);
    hold on,
    plot([0:1:48],mean(TCG2((clusters == c),:)',2),'k-','LineWidth',3);
    plot([0 6 12 24 48],mean(TCG2((clusters == c),[0 6 12 24 48]+1)',2),'b .','MarkerSize',30);
    axis tight
end
suptitle('Hierarchical Clustering of Resistant TCG Profiles');


%%%% K means clustering 
rng('default');

TCG1_new=TCG1(sum(TCG1,2)~=0,:);
TCG2_new=TCG2(sum(TCG2,2)~=0,:);



[cidx, ctrs] = kmeans(TCG1_new,9,'dist','sqeuclidean','rep',5,'disp','final');
figure
for c = 1:9
    subplot(3,3,c);
    plot(0:1:48,TCG1_new((cidx == c),:)','Color',[0.8 0.8 0.8],'MarkerSize',2);
    hold on,
    plot([0:1:48],mean(TCG1_new((cidx == c),:)',2),'k-','LineWidth',3);
    plot([0 6 12 24 48],mean(TCG1_new((cidx == c),[0 6 12 24 48]+1)',2),'b .','MarkerSize',30);
    axis tight
end
suptitle('K-Means Clustering of Sensitive TCGProfiles');


[cidx, ctrs] = kmeans(TCG2_new,9,'dist','sqeuclidean','rep',5,'disp','final');
figure
for c = 1:9
    subplot(3,3,c);
    plot(0:1:48,TCG2_new((cidx == c),:)','Color',[0.8 0.8 0.8],'MarkerSize',2);
    hold on,
    plot([0:1:48],mean(TCG2_new((cidx == c),:)',2),'k-','LineWidth',3);
    plot([0 6 12 24 48],mean(TCG2_new((cidx == c),[0 6 12 24 48]+1)',2),'b .','MarkerSize',30);
    axis tight
end
suptitle('K-Means Clustering of Resistant TCGProfiles');