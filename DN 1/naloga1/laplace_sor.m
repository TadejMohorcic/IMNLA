function [U,k] = laplace_sor(V,K,h,tol,maxIter,w)
% function [U,k] = laplace_sor(V,K,h,tol,maxIter,w) poišče rešitev sistema
% Ax = b z metodo SOR (Successive over-relaxation), kjer je A Laplaceova
% matrika v dveh dimenzijah. Izhodna podatka sta matrika U z vrednostmi
% u(x_i,y_j) rešitve u Laplaceove enačbe in število korakov k, da je metoda
% skonvergirala. Vhodni podatki so matrika V robnih pogojev za območje
% [-a,a]*[-a,a], kaznovalna matrika K nekega območja pod [-a,a]*[-a,a],
% razdalja med delilnimi točkami h, toleranca tol, največje dovoljeno
% število korakov maxIter in parameter w.
%
% Tadej Mohorčič, 2023

n = size(V,1);
razlika = inf;
k = 0;
U = V;

while (razlika > tol) && (k < maxIter)
    for i = 2:n-1
        for j = 2:n-1
            U(i,j) = w*(h^2 - U(i-1,j) - U(i,j-1) - U(i+1,j) - U(i,j+1))/(h^2*K(i,j) - 4) + (1-w)*U(i,j);
        end
    end
    razlika = max(max(abs(V-U)));
    V = U;
    k = k+1;

end