function sCl = kmeansoptrnd(E,N,NCl)   %Race,10,NCl

[Np,Ne] = size(E);

%% Randomization
ERnd = zeros(Np,Ne);

for i = 1:Ne
    ERnd(:,i) = E(randperm(Np),i);    %random perm of sce for each cell
end

%% Covariance matrix
M = zeros(Ne,Ne);
parfor i = 1:Ne
    for j = 1:Ne
        M(i,j) = covnorm(ERnd(:,i),ERnd(:,j),0);
    end
end
M(isnan(M)) = 0;


%% k-means

for k = 1:N
    % IDX = kmeans(M,NCl,"Replicates",100); %modified on 2023-10-22
    IDX = kmeans(M,NCl);
    s = silh(M,IDX);
    IDX0(k,:) = IDX;
    S(k) = median(s);
end

% parfor k = 1:N
%     IDX = kmeans(M,NCl,),'Options',statset('UseParallel',1),'Replicates',N);
%     s = silh(M,IDX);
%     IDX0(k,:) = IDX;
%     S(k) = median(s);
% end

% keep best silhouette of permuted cell/SCE matrix
[~,ClOK] = max(S);      %max of silhouette median to find best clustering among the repeat
IDX = IDX0(ClOK,:);         %choose and write best clustering of sce
s = silh(M,IDX);                %calculate silhouette of this best clustering with a value by SCE
sCl = zeros(1,NCl);
for i = 1:NCl
    sCl(i) = median(s(IDX==i));         %calculate silhouette for each cluster
end
sCl = max(sCl);         