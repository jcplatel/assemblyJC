%% Helper: silhouette avec matrice de distances
function s = silhouette_from_dist(distMat, labels)
n = size(distMat,1);
s = zeros(n,1);
uniqueLabs = unique(labels);
for i = 1:n
    lab = labels(i);
    inIdx = find(labels == lab);
    outIdx = find(labels ~= lab);
    if numel(inIdx) == 1
        s(i) = 0; continue;
    end
    a = mean(distMat(i, inIdx(inIdx~=i))); % intra
    b = inf;
    for other = uniqueLabs(:)'
        if other == lab, continue; end
        otherIdx = find(labels == other);
        b = min(b, mean(distMat(i, otherIdx)));
    end
    s(i) = (b - a) / max(a,b);
end
end