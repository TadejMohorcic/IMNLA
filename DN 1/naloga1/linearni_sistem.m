function M = linearni_sistem(a,n,N)
% function M = linearni_sistem(a,n,N) poišče matriko M, ki predstavlja
% razpršen sistem Mx=b za reševanje robnega problema neke parcialne
% diferencialne enačbe na domeni [-a,a]*[-a,a]. Izhodni podatek je matrika
% M. Vhodni podatki pa so polovica širine intervala a, število delilnih
% točk n in kaznovalna matrika N.
%
% Tadej Mohorčič, 2023

h = 2*a/(n-1);
vec = reshape(N,1,[]);
col = find(vec);

A = -4*eye(n) + diag(repelem(1,n-1),1) + diag(repelem(1,n-1),-1);
Acell = repmat({A},1,n);
B = reshape(N,[],1);
M = sparse(blkdiag(Acell{:}) + h^2*diag(B));

l = size(M,1);

for i = (n+1):l
    M(i,i-n) = 1;
    M(i-n,i) = 1;
end

end