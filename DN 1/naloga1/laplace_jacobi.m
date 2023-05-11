function [U,k] = laplace_jacobi(V,K,h,tol,maxIter)
% function [U,k] = laplace_jacobi(V,K,h,tol,maxIter) poišče rešitev sistema
% Ax = b z Jacobijevo iteracijo, kjer je A Laplaceova matrika v dveh
% dimenzijah. Izhodna podatka sta matrika U z vrednostmi u(x_i,y_j) rešitve
% Laplaceove enačbe in število korakov k, da je metoda skonvergirala.
% Vhodni podatki so matrika V robnih pogojev za območje [-a,a]*[-a,a],
% kaznovalna matrika K nekega območja pod [-a,a]*[-a,a], razdalja med
% delilnimi točkami h, toleranca tol in največje dovoljeno število korakov
% maxIter.
%
% Tadej Mohorčič, 2023

n = size(V,1);
razlika = inf;
k = 0;
U = V;

while (razlika > tol) && (k < maxIter)
    for i = 2:n-1
        for j = 2:n-1
            U(i,j) = (h^2 - V(i-1,j) - V(i,j-1) - V(i+1,j) - V(i,j+1))/(h^2*K(i,j) - 4);
        end
    end
    razlika = max(max(abs(V-U)));
    V = U;
    k = k+1;

end