function [IDXs,sCl,M,S] = kmeansopt(E,N,type)

%E p parameters (cells) by N Events
%N number of trials per cluster number
%p p-value for cluster validation

Ne = size(E,2);

%% Covariance matrix
if strcmp(type,'var')
    M = CovarM(E);
end
max_cluster=30;
%% k-means loop
% rng("default")
parfor k = 1:N*max_cluster
    NCl = floor((k-1)/N) + 2;
    % IDX = kmeans(E',NCl)'; %Normal K-means on distance metric
    % IDX = kmeans(M,NCl,"MaxIter",300,'OnlinePhase','on');%,'distance','cityblock');    % Kmeans on distance of covariance metric
    IDX = kmeans(M,NCl,"MaxIter",300);%,'distance','cityblock');    % Kmeans on distance of covariance metric
    s = silh(M,IDX);
    IDX0(k,:) = IDX;
    S(k) = mean(s);
end
% stream = RandStream('mlfg6331_64');  % Random number stream
% options = statset('UseParallel',1,'UseSubstreams',1, 'Streams',stream);
% % toc
% tic
% for k = 2:20
%     [IDX,~,sumD]  = kmeans(M,k,'Replicates',100,'Options',options,"MaxIter",300,'OnlinePhase','on','display','final'); 
%     sumDk{k}=sumD;
%     s=silhouette(M,IDX);
%     IDX0(k,:) = IDX;
%      S(k) = mean(s);
% end
% toc
% % 
% stream = RandStream('mlfg6331_64');  
% % Random number stream 
% options = statset('UseParallel',1,'UseSubstreams',1,  'Streams',stream);
% for k = 2:20
%     %NCl = floor((k-1)/N) + 2;
% %     IDX = kmeans(E',NCl)'; %Normal K-means on distance metric
%     IDX = kmeans(M,k,'Options',options,'Replicates',1000); %Kmeans on distance of covariance metric %'Options',statset('UseParallel',1)
%     s = silh(M,IDX);
%     IDX0(k,:) = IDX;
%     S(k) = mean(s);
% end


%Best clustering for each cluster number

% IDX1 = zeros(18,Ne);  %removed 2023-11-19
% for i = 1:18
%     tmp = IDX0((i-1)*N+(1:N),:);
%     [~,idx] = max(S((i-1)*N+(1:N))); % maybe 95% should be better
%     IDX1(i,:) = tmp(idx,:);
% end

%% keep best silhouette  %%%seem redundant...
for n=1:max_cluster
    med_S(n)=median(S(((n-1)*N+1):n*N),"omitmissing");
    max_S(n)=max(S(((n-1)*N+1):n*N));
    per_S(n)=prctile(S(((n-1)*N+1):n*N),95);
end
% % [~,ClOK] = max(med_S);
% best_K=max(max_S(3:max_cluster));
% [~,ClOK] = max();
%put k<4 at 0 to find max starting at k=4
S(1:(2*N)+1)=0;
% [~,ClOK] = max(S((2*N)+1:end)); % maybe 95% should be better  
[~,ClOK] = max(S);
% test = prctile(S,95); 
NCl = floor((ClOK-1)/N) + 2;
% NCl=ClOK-1;
IDX = IDX0(ClOK,:);
s = silh(M,IDX);
sCl = zeros(1,NCl);
for i = 1:NCl
    sCl(i) = median(s(IDX==i));
end

%sort RACE/silhouette of best cluster
[sCl,xCl] = sort(sCl,'descend');
IDXs = zeros(1,Ne);
for i = 1:NCl
    IDXs = IDXs + (IDX == xCl(i))*i;
end

