function sCl = kmeansoptrndjc(E,N,NCl)

[Np,Ne] = size(E);

%% Randomization
ERnd = zeros(Np,Ne);
for i = 1:Ne
    ERnd(:,i) = E(randperm(Np),i);
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
 % rng default
% parfor k = 1:N
%     IDX = kmeans(M,NCl);
%     s = silh(M,IDX);
%     IDX0(k,:) = IDX;
%     S(k) = median(s);
% end

stream = RandStream('mlfg6331_64');  
options = statset('UseParallel',1,'UseSubstreams',1,  'Streams',stream);


% IDX = kmeans(M,NCl,'Options',options,'Replicates',1000);
IDX = kmeans(M,NCl,options,"MaxIter",300,'OnlinePhase','on','Replicates',100);  %'OnlinePhase','on', 'TolFun', 1e-4=default
   
s= silhouette(M,IDX);
% s = silh(M,IDX);
% IDX0(k,:) = IDX;
% S = s;

% keep best silhouette
% [~,ClOK] = max(S);
% IDX = IDX0(ClOK,:);
%s = silh(M,IDX);
sCl = zeros(1,NCl);
for i = 1:NCl
    sCl(i) = median(s(IDX==i));
    % sCl(i) = mean(s(IDX==i));
end
sCl = max(sCl);
