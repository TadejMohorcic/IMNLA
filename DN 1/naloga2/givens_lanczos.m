function [x,res,c,s,af,bt,f,g,h] = givens_lanczos(A,b,x0,tol,k)
% function [x,res] = givens_lanczos(A,b,x0,tol,k) išče rešitev sistema Ax =
% b s pomočjo direktnega Lanczosovega algoritma z Givensovimi rotacijami.
% Vhodni podatki so matrika A, rešitev b, začetni približek x0, toleranca
% in število korakov. Funkcija vrne x bodisi po k korakih bodisi ko
% norma ostanka pade pod toleranco, in pa vektor rešitev res.
%
% Tadej Mohorčič, 2023

y = A\b;

af = zeros(1,k);
bt = zeros(1,k);
c = zeros(1,k);
s = zeros(1,k);
v = zeros(size(A,1),k);
p = zeros(size(A,1),k);

x = zeros(1,k);
res = zeros(1,k);

% če je k = 1 izvedemo kodo po svoje
if k == 1
    r0 = b-A*x0;
    v(:,1) = r0/norm(r0);
    
    z = A*v(:,1);
    af(1) = v(:,1)'*z;
    z = z - af(1)*v(:,1);
    bt(1) = norm(z);
    
    r = sqrt(af(1)^2 + bt(1)^2);
    c(1) = af(1)/r;
    s(1) = -bt(1)/r;
    
    f = af(1);
    p(:,1) = v(:,1)/f;
    ksi = norm(r0);
    x = x0 + p(:,1)*ksi;
    res(1) = -s(1)*norm(r0);
    return

% tudi če je k = 2 se izvede po svoje
elseif k == 2
    %- k = 1 ---------------
    r0 = b-A*x0;
    v(:,1) = r0/norm(r0);
    
    z = A*v(:,1);
    af(1) = v(:,1)'*z;
    z = z - af(1)*v(:,1);
    bt(1) = norm(z);

    r = sqrt(af(1)^2 + bt(1)^2);
    c(1) = af(1)/r;
    s(1) = -bt(1)/r;

    f = c(1)*af(1) - s(1)*bt(1);

    p(:,1) = v(:,1)/f;
    ksi = c(1)*norm(r0);
    x = x0 + p(:,1)*ksi;
    res(1) = -s(1)*norm(r0);

    v(:,2) = z/bt(1);

    %- k = 2 ---------------
    z = A*v(:,2)-bt(1)*v(:,1);
    af(2) = v(:,2)'*z;
    z = z - af(2)*v(:,2);
    bt(2) = norm(z);

    r = sqrt((s(1)*bt(1) + c(1)*af(2))^2 + bt(2)^2);
    c(2) = (s(1)*bt(1) + c(1)*af(2))/r;
    s(2) = -bt(2)/r;


    f = s(1)*bt(1) + c(1)*af(2);
    g = c(1)*bt(1) - s(1)*af(2);

    p(:,2) = (v(:,2) - g*p(:,1))/f;
    ksi = s(1)*norm(r0);

    x = x + p(:,2)*ksi;
    res(2) = abs(res(1)/(s(1)*bt(1) + c(1)*af(2)))*bt(2);
    return

else
    %- k = 1 ---------------
    r0 = b-A*x0;
    v(:,1) = r0/norm(r0);
    
    z = A*v(:,1);
    af(1) = v(:,1)'*z;
    z = z - af(1)*v(:,1);
    bt(1) = norm(z);

    r = sqrt(af(1)^2 + bt(1)^2);
    c(1) = af(1)/r;
    s(1) = -bt(1)/r;

    f = c(1)*af(1) - s(1)*bt(1);

    p(:,1) = v(:,1)/f;
    ksi = c(1)*norm(r0);
    x = x0 + p(:,1)*ksi;
    res(1) = -s(1)*norm(r0);

    v(:,2) = z/bt(1);

    %- k = 2 ---------------
    z = A*v(:,2)-bt(1)*v(:,1);
    af(2) = v(:,2)'*z;
    z = z - af(2)*v(:,2);
    bt(2) = norm(z);

    r = sqrt((s(1)*bt(1) + c(1)*af(2))^2 + bt(2)^2);
    c(2) = (s(1)*bt(1) + c(1)*af(2))/r;
    s(2) = -bt(2)/r;

    f = s(1)*c(2)*bt(1) + c(1)*c(2)*af(2) - s(2)*bt(2);
    g = c(1)*bt(1) - s(1)*af(2);

    p(:,2) = (v(:,2) - g*p(:,1))/f;
    ksi = (s(1)*c(2)*ksi)/c(1);
    x = x + p(:,2)*ksi;
    res(2) = abs(res(1)/(s(1)*bt(1) + c(1)*af(2)))*bt(2);
    y = s(1)*norm(r0);

    v(:,3) = z/bt(2);
    %- k > 2 ---------------
    for i = 3:k
        z = A*v(:,i) - bt(i-1)*v(:,i-1);
        af(i) = v(:,i)'*z;
        z = z - af(i)*v(:,i);
        bt(i) = norm(z);
        
        r = sqrt((s(i-1)*c(i-2)*bt(i-1) + c(i-1)*af(i))^2 + bt(i)^2);
        c(i) = (s(i-1)*c(i-2)*bt(i-1) + c(i-1)*af(i))/r;
        s(i) = -bt(i)/r;

        f = s(i-1)*c(i-2)*bt(i-1) + c(i-1)*af(i);
        g = c(i-1)*c(i-2)*bt(i-1) - s(i-1)*af(i);
        h = -s(i-2)*bt(i-1);

        p(:,i) = (v(:,i) - g*p(:,i-1) - h*p(:,i-2))/f;
        ksi = (s(i-1)*ksi)/c(i-1);
        y = s(i-1)*y;
        res(i) = abs(y/f)*bt(i);
        if res(i) < tol || i == k
            x = x + p(:,i)*ksi;
            return
        end

        f = c(i)*s(i-1)*c(i-2)*bt(i-1) + c(i-1)*c(i)*af(i) - s(i)*bt(i);

        p(:,i) = (v(:,i) - g*p(:,i-1) - h*p(:,i-2))/f;
        ksi = c(i)*ksi;
        x = x + p(:,i)*ksi;

        v(:,i+1) = z/bt(i);
    end

end