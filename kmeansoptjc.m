function [IDXs,sCl,M,S] = kmeansopt(E,N,type)

%E p parameters (cells) by N Events
%N number of trials per cluster number
%p p-value for cluster validation
%IDXs  index cluster chaque sce for best cluster value
%sCl: median silhouette each cluster of best cluster value
%M: CovarM(E);
%S  : silhouette moyenne pour chaque nombre de cluster de 2 à 20 


Ne = size(E,2);

%% Covariance matrix
if strcmp(type,'var')
    M = CovarM(E);
end

%% k-means loop
%
rng("default")
% parfor k = 1:N*18
%     NCl = floor((k-1)/N) + 2;
% %     IDX = kmeans(E',NCl)'; %Normal K-means on distance metric
%     IDX = kmeans(M,NCl); %Kmeans on distance of covariance metric
%     s = silh(M,IDX);
%     IDX0(k,:) = IDX;
%     S(k) = mean(s);
% end

stream = RandStream('mlfg6331_64');  
% Random number stream 
options = statset('UseParallel',1,'UseSubstreams',1,  'Streams',stream);
tic
for k = 2:20
    %NCl = floor((k-1)/N) + 2;
%     IDX = kmeans(E',NCl)'; %Normal K-means on distance metric
    IDX = kmeans(M,k,'Options',options,'Replicates',1000,'MaxIter',1000); %Kmeans on distance of covariance metric %'Options',statset('UseParallel',1)
    s = silh(M,IDX);
    IDX0(k,:) = IDX;
    S(k) = mean(s);
end
toc

%Best clustering for each cluster number
% IDX1 = zeros(18,Ne);
% for i = 2:20
%     [~,idx] =max(S(i));
%     %tmp = IDX0((i-1)*N+(1:N),:);
%     [~,idx] = max(S((i-1)*N+(1:N)));
%     IDX1(i,:) = tmp(idx,:);
% end

%% keep best silhouette
 [~,NCl] = max(S);
% NCl = floor((ClOK-1)/N) + 2;
 IDX = IDX0(NCl,:);
 s = silh(M,IDX);
 sCl = zeros(1,NCl);
 for i = 1:NCl
     sCl(i) = median(s(IDX==i));
 end

[sCl,xCl] = sort(sCl,'descend');
IDXs =IDX;
% IDXs = zeros(1,Ne);
% for i = 1:NCl
%     IDXs = IDXs + (IDX == xCl(i))*i;
% end
