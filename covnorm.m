function C=covnorm(A,B,win)

N=length(A);
% C=xcov(A,B,win)/std(A)/std(B)/N;
A_centered = A - mean(A);
B_centered = B - mean(B);
C = (A_centered' * B_centered) / (std(A) * std(B) * N);