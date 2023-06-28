function [V,e] = pageTriDiagEig(T,Qtridiag)
%PAGETRIDIAGEIG  It applies the pagewised eigendecomposition of a 
%tridiagonal matrix with QR transformation
%
% Inputs:       T is a multidimensional array with (anti-)symmetric (if 
%               real) or (skew-)Hermitian (if complex) NxN matrix of each 
%               page i, (:,:,i).
%
% Outputs:      V is a multidimensional array with eigenvector-matrix in
%               each page
%
%               e is a multidimensioal array with eigenvalues in each page
%
% Description:  This function uses the QR decomposition for 
%               eigendecomposition.
%               
%               M = Qtridiag *     T      * Qtridiag'
%                 = Qtridiag * Q * e * Q' * Qtridiag', Q = Q1'*Q2'*... where Q1, Q2,... are Q matrices of each QR transform
%                 =            V * e * V'
%
% Author:       Yi-Hao Chen
%               yc2368@cornell.edu
%
% Date:         5/15/2023

if isequal(class(T),'gpuArray')
    use_gpu = true;
else
    use_gpu = false;
end

qrTol = 1e-10;

% Check input matrix size
[m,numM] = size(T,[1,3]);

% Initialize the tridiagonal matrices
R = T;
if use_gpu
    V = repmat(eye(m,'gpuArray'),1,1,numM);
    e = zeros(1,m,numM,'gpuArray');
else
    V = repmat(eye(m),1,1,numM);
    e = zeros(1,m,numM);
end

while m > 1
    if use_gpu
        Qm = repmat(eye(m,'gpuArray'),1,1,numM);
    else
        Qm = repmat(eye(m),1,1,numM);
    end
    
    % Apply Wilkinson shift: use the eigenvalue of the bottom-right 2x2 
    % matrix, which is closer to R(m,m,:), for the QR shift
    meanR = (R(m-1,m-1,:) + R(m,m,:))/2;
    detR = R(m-1,m-1,:).*R(m,m,:) - R(m-1,m,:).*R(m,m-1,:);
    tmpR = sqrt(meanR.^2-detR);
    lambda = cat(1,meanR+tmpR, meanR-tmpR); % the eigenvalues of a 2x2 matrix can be computed analytically
    [~,idx] = min(abs(lambda-R(m,m,:)),[],1);
    qr_shift = lambda(idx + permute((0:numM-1)*2,[1,3,2]));

    R = R - qr_shift.*eye(m);
    % QR transformation with Given's matrices
    for i = 1:m-1
        % Make zeros with Given's matrices
        r = sqrt(sum(abs(R([i,i+1],i,:)).^2,1));
        c = R(i  ,i,:)./r;
        s = R(i+1,i,:)./r;
        G = cat(1,cat(2,  conj(c), conj(s)),...
                  cat(2,      -s ,      c )); % complex version of Given's matrix (from the patent US8473539B1)
        if use_gpu
            Qm(:,[i,i+1],:) = pagefun(@mtimes, Qm(:,[i,i+1],:),pagefun(@ctranspose,G));

            R([i,i+1],:,:) = pagefun(@mtimes, G,R([i,i+1],:,:));
        else
            Qm(:,[i,i+1],:) = pagemtimes(Qm(:,[i,i+1],:),pagectranspose(G));

            R([i,i+1],:,:) = pagemtimes(G,R([i,i+1],:,:));
        end
    end
    if use_gpu
        R = pagefun(@mtimes,R,Qm) + qr_shift.*eye(m,'gpuArray');
    else
        R = pagemtimes(R,Qm) + qr_shift.*eye(m);
    end
    
    if use_gpu
        V(:,1:m,:) = pagefun(@mtimes,V(:,1:m,:),Qm);
    else
        V(:,1:m,:) = pagemtimes(V(:,1:m,:),Qm);
    end
    
    if all(abs(R(m,m,:)) > abs(R(m,m-1,:))/qrTol) && all(abs(R(m-1,m-1,:)) > abs(R(m,m-1,:))/qrTol)
        e(1,m,:) = R(m,m,:);
        m = m-1;
        R = R(1:m,1:m,:);
    end
end
e(1,1,:) = R(1,1,:);
if use_gpu
    V(:,[1,2],:) = pagefun(@mtimes,V(:,[1,2],:),Qm);
else
    V(:,[1,2],:) = pagemtimes(V(:,[1,2],:),Qm);
end

% Include the previous similarity transformation for obtaining the 
% tridiagonal matrix M
if use_gpu
    V = pagefun(@mtimes,Qtridiag,V);
else
    V = pagemtimes(Qtridiag,V);
end

end