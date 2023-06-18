function M = CovarM(E)

Ne = size(E,2);
M = zeros(Ne,Ne);

parfor i = 1:Ne  %was parfor
    for j = 1:Ne
        M(i,j) = covnorm(E(:,i),E(:,j),0);
    end
end
M(isnan(M)) = 0;