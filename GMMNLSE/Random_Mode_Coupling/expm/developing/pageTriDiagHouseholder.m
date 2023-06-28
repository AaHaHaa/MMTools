function [T,Q] = pageTriDiagHouseholder(M)
%PAGETRIDIAGHOUSEHOLDER It finds the tridiagonal matrix of M after 
%pagewised Householder transformation
%
% Inputs:       M is a multidimensional array with symmetric (if real) or
%               Hermitian (if complex) NxN matrix of each page i, (:,:,i).
%
% Outputs:      In each page, T is a tridiagonal matrix that is similar to  
%               the corresponding page of M; that is to say, T has the same
%               symmetry as M. In particular, this means that T has the 
%               same eigenvalues as M.
%               
%               Q is the similarity transformation such that M = Q*T*Q' of 
%               each page
%
% Description:  This function uses the Householder reduction algorithm to 
%               compute a tridiagonal matrix T that is similar to M with 
%               respect to similarity transformation Q. That is, M=Q*T*Q'.
%
% Author:       Yi-Hao Chen
%               yc2368@cornell.edu
%
% Date:         5/14/2023

if isequal(class(M),'gpuArray')
    use_gpu = true;
else
    use_gpu = false;
end

symTol = 1e-8;

% Check input matrix size
[m,n,numM] = size(M);
if m ~= n
    error('Input matrix must be square');
end
% Make sure input matrix has any symmetry structure
cc_M = pagectranspose(M); % complex-conjugate transpose
check_Hermitian = abs(M - cc_M);
check_skewHermitian = abs(M + cc_M);
if max(check_Hermitian(:)) < symTol
    symmetry = 1; % symmetry (if real) or Hermitian (if complex)
elseif max(check_skewHermitian(:)) < symTol
    symmetry = -1; % anti-symmetry (if real) or skew-Hermitian (if complex)
else
    error('Input matrix does not have any symmetry (symmetry and anti-symmetry for real matrices, or Hermitian and skew-Hermitian for complex matrices) to the working precision.');
end

% Initialize the tridiagonal matrices
T = M;
if use_gpu
    Q = eye(m,'gpuArray');
else
    Q = eye(m);
end

for i = 1:m-2
    % Property of Householder matrix P:
    % For any two vector x and y with the following conditions:
    %   1. they have the same length, |x|=|y|
    %   2. x'*y = y'*x,
    % the Householder matrix P = I-2*v*v', where v = (x-y)/|x-y| satisfies
    %    P*x = y.
    off_diagonal = -T(i+1,i,:).*sqrt(sum(abs(T(i+1:end,i,:)).^2))./abs(T(i+1,i,:)); % negative sign here is to avoid "loss of significance" with floating-point operations
    v = T(i+1:end,i,:);
    v(1,1,:) = v(1,1,:) - off_diagonal;
    v = v./sqrt(sum(abs(v).^2,1));
    if use_gpu
        P = eye(m-i,'gpuArray') - 2*pagefun(@mtimes, v,pagefun(@ctranspose,v));
        Pm = repmat(eye(m,'gpuArray'),1,1,numM); Pm(i+1:end,i+1:end,:) = P;
        Q = pagefun(@mtimes, Q,Pm);
    else
        P = eye(m-i) - 2*pagemtimes(v,pagectranspose(v));
        Pm = repmat(eye(m),1,1,numM); Pm(i+1:end,i+1:end,:) = P;
        Q = pagemtimes(Q,Pm);
    end
    
    T(i+1,i,:) = off_diagonal;
    T(i,i+1,:) = symmetry*conj(T(i+1,i,:));
    T(i+2:end,i,:) = 0;
    T(i,i+2:end,:) = 0;
    if use_gpu
        T(i+1:end,i+1:end,:) = pagefun(@mtimes, pagefun(@mtimes, P,T(i+1:end,i+1:end,:)),pagefun(@ctranspose,P));
    else
        T(i+1:end,i+1:end,:) = pagemtimes(pagemtimes(P,T(i+1:end,i+1:end,:)),pagectranspose(P));
    end
end

end