function M = kaznovalna_metoda(a,n,k)
% function M = kaznovalna_metoda(a,n,k) poišče matriko M, v kateri so na
% mestih (i,j) vrednosti enake k, če je točka (x_i,y_j) od izhodišča
% oddaljena za manj kot 1/10, in 0 sicer. Izhodni podatek je matrika M.
% Vhodni podatki so polovica širine intervala a, število delilnih točk n in
% vrednost k.
%
% Tadej Mohorčič, 2023

x = linspace(-a,a,n);
y = linspace(-a,a,n);
M = zeros(n);

for i = 1:n
    for j = 1:n
        r = x(i)^2 + y(j)^2;
        if r < 1/10
            M(i,j) = k;
        else
            M(i,j) = 0;
    end
end

end