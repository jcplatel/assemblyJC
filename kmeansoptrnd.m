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
% stream = RandStream('mlfg6331_64');  % Random number stream
% options = statset('UseParallel',1,'UseSubstreams',1, 'Streams',stream);

for k = 1:N
    IDX = kmeans(M,NCl); %modified on 2023-10-22
    % IDX = kmeans(M,NCl,options,"MaxIter",300,'OnlinePhase','on','display','final');  %'OnlinePhase','on', 'TolFun', 1e-4=default
    s = silh(M,IDX);
    IDX0(k,:) = IDX;
    % S(k) = median(s);
    S(k) = mean(s);
end


% keep best silhouette of permuted cell/SCE matrix
[~,ClOK] = max(S);      %max of silhouette median to find best clustering among the repeat
IDX = IDX0(ClOK,:);         %choose and write best clustering of sce
s = silh(M,IDX);                %calculate silhouette of this best clustering with a value by SCE
sCl = zeros(1,NCl);
for i = 1:NCl
    % sCl(i) = median(s(IDX==i));         %calculate silhouette for each cluster
     sCl(i) = mean(s(IDX==i));
end
sCl = max(sCl);         