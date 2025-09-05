% Helper: Benjamini-Hochberg FDR (returns adj p-values)
function [h, crit_p, adj_p, sorted_p] = fdr_bh(pvals, q)
% Simple implementation of Benjamini-Hochberg that returns adjusted p-values
% pvals: vector of p-values
% q: desired FDR level (used here only for returned h)
p = pvals(:);
[m, ~] = size(p);
[sorted_p, sort_ids] = sort(p);
adj_sorted = zeros(size(p));
for i = 1:m
    adj_sorted(i) = sorted_p(i) * m / i;
end
% enforce monotonicity
for i = m-1:-1:1
    adj_sorted(i) = min(adj_sorted(i), adj_sorted(i+1));
end
% cap at 1
adj_sorted = min(adj_sorted, 1);
% unsort
adj_p = zeros(size(p));
adj_p(sort_ids) = adj_sorted;
h = adj_p <= q;
crit_p = max(sorted_p(h)) 
if any(h) 
else 
    NaN; %#ok<STCON>
end