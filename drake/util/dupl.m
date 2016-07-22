function D = dupl(n)
%DUPL Returns the duplication matrix D_n such that vec(S) = D_n*vech(S).
%This is useful for taking derivatives with respect to symmetric matrices.

a = tril(ones(n));
i = find(a);
a(i) = 1:length(i);
a = a + tril(a,-1)';
j = vec(a);

D = sparse(1:(n*n), j, ones(n*n,1), n*n, n*(n+1)/2);

end

