function [IDXs,sCl,M,S,NClini] = kmeansopttest(E,N,type,NClini) % on clusterize les SCE sur la base des cellules qui y participe

%E p parameters (cells) by N Events
%N number of trials per cluster number
%p p-value for cluster validation

Ne = size(E,2);

%% Covariance matrix
if strcmp(type,'var')
    M = CovarM(E);
end

%% k-means loop
% rng("default")
%loop to find best number of cluster
% M=gpuArray(M);
% if NClini==0   %calculate if does not exist
%     for NCl=5:18
%         parfor k = 1:100
%             % NCl = floor((k-1)/N) + 2;
%             % IDX = kmeans(E',NCl)'; %Normal K-means on distance metric
%             IDX = kmeans(M,NCl,'distance','correlation');    % Kmeans on distance of covariance metric
%             s = silh(M,IDX);
%             S(k) = mean(s);%original
%         end
%         SS(NCl)=mean(S);
%     end
%     [~ ,NClini]=max(SS);
% end

% S=[];
% s=[];
% IDX0=[];
S=zeros(N,1);
IDX0=zeros(N,Ne);
NCl=NClini;
cM = parallel.pool.Constant(M);
parfor k = 1:N
    Mloc = cM.Value; 
    % NCl = floor((k-1)/N) + 2;
    % IDX = kmeans(E',NCl)'; %Normal K-means on distance metric
    IDX = kmeans(Mloc,NCl,"MaxIter",300,'OnlinePhase','on');%,'distance','cityblock');    % Kmeans on distance of covariance metric
    % s = silh(M,IDX);
    IDX0(k,:) = IDX;
    % S(k) = mean(s);%original
    S(k) = mean(silh(Mloc,IDX));%original
end


[~,ClOK] = max(S); 
% test = prctile(S,95); 

IDX = IDX0(ClOK,:);
s = silh(M,IDX);
sCl = zeros(1,NCl);
for i = 1:NCl
    sCl(i) = median(s(IDX==i));%original
    % sCl(i) = mean(s(IDX==i));
end

%sort RACE/silhouette of best cluster
[sCl,xCl] = sort(sCl,'descend');
IDXs = zeros(1,Ne);
for i = 1:NCl
    IDXs = IDXs + (IDX == xCl(i))*i;
end

