%% Helper: compute modularity for weighted adjacency matrix A and labels
function Q = compute_modularity(A, labels)
% A: symmetric adjacency (weights)
% labels: integer labels (n x 1)
A = (A + A')/2;
n = size(A,1);
m = sum(A(:))/2;
if m == 0
    Q = 0;
    return;
end
% degree
k = sum(A,2);
Q = 0;
unique_lbl = unique(labels(:));
for i = 1:numel(unique_lbl)
    members = find(labels == unique_lbl(i));
    subA = A(members,members);
    kin = sum(sum(subA));
    % expected weight = sum degrees product / (2m)
    deg_sum = sum(k(members));
    Q = Q + (kin/(2*m) - (deg_sum/(2*m))^2);
end
end