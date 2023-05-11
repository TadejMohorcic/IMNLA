function [L,l] = konvergenca(a,n,K,b,tol)
% function [L,l] = konvergenca(a,n,K,b,tol) računa vpliv parametra k pri
% kaznovalni metodi reševanja neke parcialne diferencialne enačbe na
% območju [-a,a]^2. Vhodni podatki so širina intervala 2*a, število
% delilnih točk n, vrstica parametrov k K, desna stran sistema b in
% toleranca tol. Izhodni podatki so vrstica norm L in število korakov l,
% preden je razlika zadnjih dveh norm manjša kot tol.
%
% Tadej Mohorčič, 2023

razlika = inf;
l = 2;
L = zeros(1,length(K));

N = kaznovalna_metoda(a,n,K(1));
M = linearni_sistem(a,n,N);
res = M\b;
L(1) = norm(res);

while (razlika > tol) && (l < length(K))
    N = kaznovalna_metoda(a,n,K(l));
    M = linearni_sistem(a,n,N);
    res = M\b;
    L(l) = norm(res);
    razlika = abs(L(l) - L(l-1));
    l = l+1;
end

end