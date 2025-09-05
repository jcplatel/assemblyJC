
%% Clustering

[NCell,NRace] = size(Race);
% NClini=10;
if ~exist ('NClini','var')
    NClini=0;
end
%Race=gpuArray(Race);
[IDX2,sCl,M,S,NClini] = kmeansopttest(Race,kmean_iter,'var',NClini);%increase to 50,  100 ?    
% M = CovarM(Race);
% IDX2 = kmedoids(M,NCl);
NCl = max(IDX2);
[~,x2] = sort(IDX2);%cluster de SCE
MSort = M(x2,x2);
disp(['nClusters: '  num2str(NCl)])


%% Save Clusters
 %save([namefull,'Clusters.mat'],'IDX2')   % liste de tous les SCE et dans quel cluster ils sont
    
%% Remove cluster non-statistically significant: basically do some sce permutation 
% and some clustering with the same number of cluster and if the best silhouette
% of one cluster is higher than real data then remove this cluster.
% since they are sorted remove the worst one 

% kmeans_surrogate=500; %or more...
sClrnd = zeros(1,kmeans_surrogate);
cRace = parallel.pool.Constant(Race);
parfor i = 1:kmeans_surrogate  
    sClrnd(i) = kmeansoptrnd(cRace.Value,10,NCl); 
end


%NClOK = sum(sCl>max(sClrnd));
NClOK =sum(sCl>prctile(sClrnd,95)); %use max, use 99%  ?
sClOK = sCl(1:NClOK)';
disp(['nClustersOK: ' num2str(NClOK)])


%%%new add
%Race clusters
% R = cell(0);        % which sce for each cluster
% CellScore = zeros(NCell,NClOK);  % number of time a cell participate in a sce of a given cluster
% CellScoreN = zeros(NCell,NClOK);
% for i = 1:NClOK
%     R{i} = find(IDX2==i);
%     CellScore(:,i) = sum(Race(:,R{i}),2);
%     CellScoreN(:,i) = CellScore(:,i)/length(R{i});
% end
% %Assign cells to cluster with which it most likely spikes in term of
% %percentage
% 
% [~,CellCl] = max(CellScoreN,[],2); %début orthogonalisation
% %Remove cells with less than 2 spikes in a given cluster   ###should
% %export this raw CellCl  
% CellCl(max(CellScore,[],2)<2) = 0;
% [X1,x1] = sort(CellCl);
% 
% %write list of cells in each cluster
% 
assemblyraw= cell(0);
% k = 0;
% for i = 1:NCl
%     k = k+1;
%     assemblyraw{k} = transpose(find(CellCl==i));
% end

%%%%%%%

RaceOK = Race(:,IDX2<=NClOK);
NRaceOK = size(RaceOK,2);
% disp(['nSCEOK: ' num2str(NRaceOK)])    


if NClOK>1
     RACE_Orthojc
     %%call rasterplot here
else
    NCl=NClOK;
     assemblyortho= cell(0);
     assemblystat= cell(0);
   % assemblyraw=[]
end

%% recalcul silhouette cluster finaux
[~,x2] = sort(IDX2);%cluster de SCE
a = find (IDX2(x2)>NCl,1, 'first')-1;
x2=x2(1:a);
MSort_ok = M(x2,x2);
IDX2_ok = IDX2(x2);
s = silh(MSort_ok, IDX2_ok);

sClOK = zeros(1,NCl);
for i = 1:NCl
    sClOK(i) = median(s(IDX2_ok==i));%original
end
mean_sClOK=mean(sClOK);
%%

if savefig==1 &&  NClOK>1
    figure('visible','off');
    % figure
    subplot(1,2,1)
    imagesc(MSort)
    colormap jet
    axis image
    xlabel('sorted SCE #')     %was RACE
    ylabel('sorted SCE #')     %was RACE
    
    subplot(1,2,2)
    % imagesc(Race(x1,x2),[-1 1.2])
    imagesc(Race(x1,RList),[-1 1.2])
    axis image
    xlabel('sorted SCE #')     %was RACE
    ylabel('sorted Cell #')
    
    exportgraphics(gcf,strcat(namefull ,'clusters.png'),'Resolution',300)
    close gcf
end

if NCl ==0
    NClOK=0;
     assemblyortho= cell(0);
     assemblystat= cell(0);
end
disp(['NCl: ' num2str(NCl)])  
% ncluster(sce_n_cells_threshold)=NCl;

%write all data in excel

% exportdata
% save([namefull,'results.mat'])  

