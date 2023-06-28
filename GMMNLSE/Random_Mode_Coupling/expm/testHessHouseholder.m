function [A,tau,P] = testHessHouseholder(A)
%--------------------------------------------------------------------------
% Syntax:       H = myHessHouseholder(mat);
%               [H Q] = myHessHouseholder(mat);
%
% Inputs:       mat is an arbitrary N x N matrix
%
% Outputs:      H is an upper Hessenberg matrix that is similar (in the
%               linear algebra sense) to mat. In particular, this means
%               that H has the same eigenvalues as input mat.
%
%               Q is the similarity transformation such that
%               mat = Q * H * Q';
%
% Description:  This function uses the Householder reduction algorithm to 
%               compute an upper Hessenberg matrix H that is similar
%               (in the linear algebra sense) to mat with respect to
%               similarity transformation Q. That is, mat = Q * H * Q'.
%
%               NOTE: In particular, H has same eigenvalues as input mat.
%
% Author:       Brian Moore
%               brimoor@umich.edu
%
% Date:         September 3, 2012
%--------------------------------------------------------------------------

% Check input matrix size
[m,n] = size(A);
if m ~= n
    error('Input matrix must be square');
end
if ishermitian(A)
    A_ishermitian = true;
else
    A_ishermitian = false;
end

tol = eps^2*max(abs(real(A)),abs(imag(A)));

tau = ones(m,1);
U = cell(1,n-2);

% Perform Householder transformation to upper Hessenberg form
for r = 2:(n-1)
    u = A(r:n,r-1);
    V = norm(u)^2;

    if V>tol
        if A(r,r-1) == 0
            delta = V;
            A(r,r-1) = -sqrt(V);
            u(1) = sqrt(V);
            tau(r) = -A(r,r-1);
        else
            root = sqrt(abs(A(r,r-1))^2*V);
            delta = V + root;
            ratio = V/root;
            tau(r) = -ratio*conj(A(r,r-1));
            u(1) = (ratio+1)*A(r,r-1);
            A(r,r-1) = -V*A(r,r-1)/root;
        end
    end
    
    v = u/delta;
    if A_ishermitian
        A(1:r-2,r:n) = 0;
        A(r-1,r) = conj(A(r,r-1));
        A(r:n,r:n) = H*A(r:n,r:n)*H;
    else
        H = eye(n-r+1) - u*v';
        A(1:r-1,r:n) = A(1:r-1,r:n)*H;
        A(r:n,r:n) = H*A(r:n,r:n)*H;
    end
    A((r+1):n,r-1) = 0;
end
tau(n) = conj(A(n,n-1));
%{
for i = 2:n
    bb = abs(tau(i));
    b(i-1) = bb;
    if bb == 0
        bb = 1;
        tau(i) = 1;
    end
    tau(i) = tau(i)*tau(i-1)/bb;
    
    for j = (i+1):n
        A(i,j) = A(i,j)*tau(i);
    end
    for j = 1:(i-1)
        A(j,i) = A(j,i)*conj(tau(i));
    end
end

T = diag(tau);
P = eye(n);
for j = (n-2):-1:1
    P((j+1):n,:) = P((j+1):n,:) - U{j} * (U{j}' * P((j+1):n,:));
end
P = P*T';

A = T*A*T';
%}