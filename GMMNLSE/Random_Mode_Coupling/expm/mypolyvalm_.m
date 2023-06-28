function [Y,nmult,c,Xpwr] = mypolyvalm_(c,X,ascending_c,optim_c,Xpwr,use_gpu)
%MYPOLYVALM_ Evaluate polynomial with matrix argument with CPU or GPU and
% can also compute for "multiple" matrices "at once". while matrices are 
% stored as pages of a multidimensional array.
%
% mypolyvalm_(c,X) efficiently calculates a polynomial with scalar coefficients
% c(1), c(2), ... and square-matrix argument X. For a polynomial of order p
% (i.e., numel(c)==p+1) the number of matrix multiplies is approximately
% 2*(sqrt(p+1)-1).
%
% syntax:
%   Y = mypolyvalm_(c,X);
%   Y = mypolyvalm_(c,X,ascending_c,optim_c,Xpwr);
%   [Y,nmult,c,Xpwr] = mypolyvalm_(c,...);
%
% inputs:
%
%   c: size-[m,n] matrix.
%   The polynomial coefficients are c(:). The shape of c affects the
%   computation algorithm; see Notes.
%
% optional inputs:
%
%   X: a multidimensional array with each page a square matrix, default = [].
%   polynomial argument. (Set X = [] to just get outputs other than Y.)
%   Size: (n,n,...)
%
%   ascending_c: true or false, default = false.
%   ascending_c = true if the c(:) coefficients are in ascending monomial
%   order, false if in descending order. (The default, false, is compatible
%   with polyvalm.)
%
%   optim_c: true or false, default = true.
%   If true, c will be automatically reshaped (with the last column possibly
%   zero-padded) to minimize the number of matrix multiplies.
%
%   Xpwr: a multimensional array with the 1st index representing i-th power
%   of X. default = [].
%   It stores some or all precomputed X powers. Xpwr(i,:) must be either []
%   (uninitialized) or X^(i-1).
%
% outputs:
%
%   Y: a multidimensioanl array with each page a square matrix, size-matched to X.
%   The polynomial order is p = numel(c)-1.
%
%   For a matrix YY in each page and its corresponding input matrix XX,
%   If ~ascending_c, YY = c(1)*XX^p + c(2)*XX^(p-1) + ... + c(p)*XX + c(p+1)*I.
%   If ascending_c, YY = c(1)*I + c(2)*XX^2 + ... + c(p+1)*XX^p.
%   (I is an identity matrix.)
%
%   nmult: integer.
%   number of matrix multiplies (i.e., BLAS3 matrix-matrix multiplies).
%
%   c: matrix.
%   same as input, but in ascending order and, if optim_c, possibly reshaped
%   and zero-padded. (CAUTION: Set ascending_c = true when reusing c.)
%
%   Xpwr: a multimensional array with the 1st index representing i-th power
%   of X. default = [].
%   It stores some or all precomputed X powers. Xpwr(i,:) must be either []
%   (uninitialized) or X^(i-1).
%   "Same as input, but possibly length-extended and with more entries
%   initialized. (See Notes.)"
%
% Notes:
%
%   Y is calculated by first putting c in ascending order (i.e., if ~ascend_c,
%   then c(:) = c(end:-1:1)), and then evaluating the following nested sum:
%     Y = sum|k=1:n (sum|j=1:m X^(j-1)*c(j,k))*(X^m)^(k-1)
%   The inner sum is evaluated by pre-computing the powers X, X^2, X^3, ...
%   X^(m-1). If n>1 then the power X^m is also precomputed and used to
%   evaluate the outer sum by Horner's method.
%
%   If the polynomial is to be evaluated with multiple X arguments but the
%   same c coefficients, then the optimum c matrix can be precomputed:
%     [~,nmult,c] = mypolyvalm_(c,[],ascending_c,true);
%   This can then be applied to specific X values:
%     Y = mypolyvalm_(c,X,true,false);
%   Also, if the coefficient values change but the number of coefficients does
%   not, then the non-zero c coefficients can be changed and used with optim_c
%   = false.
%
%   If the polynomial is to be evaluated with a single X argument but with
%   multiple c coefficient arrays, then the Xpwr buffer can be recycled to
%   avoid redundant power calculations:
%     [Y,~,~,Xpwr] = mypolyvalm_(c,X,ascending_c,optim_c,Xpwr);
%   (The algorithm for optimizing c does not consider the multiplies saved by
%   providing the Xpwr argument. With all X powers precomputed and c formatted
%   as a column vector, the algorithm would not do any matrix multiplies.)
%
% Examples:
%
%   Determine an efficient procedure for evaluating an order-12 polynomial.
%     >> p = 12; c = ones(1,p+1); [~,nmult,c] = mypolyvalm_(c)
%     nmult =
%          5
%     c =
%          1     1     1     1     1
%          1     1     1     1     0
%          1     1     1     1     0
%   The polynomial algorithm has the following form (using syntax for scalar
%   X):
%     Y = [1,X,X^2]*c*[1;X^3;X^6;X^9;X^12]
%     = c(1,1)+X*c(2,1)+X^2*c(3,1) +
%      (c(1,2)+X*c(2,2)+X^2*c(3,2) +
%      (c(1,3)+X*c(2,3)+X^2*c(3,3) +
%      (c(1,4)+X*c(2,4)+X^2*c(3,4) +
%       c(1,5)*X^3)*X^3)*X^3)*X^3
%   5 matrix multiplies are required for the following terms:
%     X^2, X^3, (c(1,4)+...)*X^3, (c(1,3)+...)*X^3, (c(1,2)+...)*X^3
%
%   Compare runtime and results for polyvalm and mypolyvalm_:
%     >> p = 12; c = rand(1,p+1); X = rand(1000);
%     >> tic, Y = polyvalm(c,X); toc
%     Elapsed time is 0.668133 seconds.
%     >> tic, Y_ = mypolyvalm_(c,X); toc
%     Elapsed time is 0.509375 seconds.
%     >> max(abs(Y_(:)-Y(:)))/max(abs(Y(:)))
%     ans =
%          1.235038181253709e-15
%   Compare runtime with X powers and c shape stored for efficient reuse.
%     >> tic, [Y_,nmult,c_,Xpwr] = mypolyvalm_(c,X,false,true); toc
%     Elapsed time is 0.456861 seconds.
%     >> c_(1:p+1) = rand(p+1,1); % The shape of c_ is unchanged.
%     >> tic, Y_ = polyvalm_(c_,X,true,false,Xpwr); toc
%     Elapsed time is 0.319544 seconds.
%     >> c = c_(p+1:-1:1); tic, Y = polyvalm(c,X); toc
%     Elapsed time is 0.537951 seconds.
%     >> max(abs(Y_(:)-Y(:)))/max(abs(Y(:)))
%     ans =
%          1.235038181253709e-15
%
%
% Author: Ken Johnson  kjinnovation@earthlink.net  3/6/2014
%
%
% Acknowledgment of Federal support and disclaimer:
%
% This material is based upon work supported by the Department of Energy under
% Award Number DE-SC0011372.
%
% This report was prepared as an account of work sponsored by an agency of the
% United States Government.  Neither the United States Government nor any
% agency thereof, nor any of their employees, makes any warranty, express or
% implied, or assumes any legal liability or responsibility for the accuracy,
% completeness, or usefulness of any information, apparatus, product, or
% process disclosed, or represents that its use would not infringe privately
% owned rights.  Reference herein to any specific commercial product, process,
% or service by trade name, trademark, manufacturer, or otherwise does not
% necessarily constitute or imply its endorsement, recommendation, or favoring
% by the United States Government or any agency thereof.  The views and
% opinions of authors expressed herein do not necessarily state or reflect
% those of the United States Government or any agency thereof.
%
%
% =========================================================================
% Author: Yi-Hao Chen
% email: jollycyclopes@gmail.com
% date: 6/6/2018
%
% Acknowledgement:
%   This code is modified from "polyvalm_" by Ken Johnson.
%

if nargin<6
    use_gpu = false;
    if nargin<5
        Xpwr = [];
        if nargin<4
            optim_c = true;
            if nargin<3
                ascending_c = false;
                if nargin<2
                    X = [];
                end
            end
        end
    end
end

sc = size(c);

if use_gpu
    if ~isequal(class(c),'gpuArray') && length(sc)>2
        c = gpuArray(c);
    else
        if isequal(class(c),'gpuArray')
            c = gather(c);
        end
    end
    if ~isequal(class(X),'gpuArray')
        X = gpuArray(X);
    end
else
    if isequal(class(c),'gpuArray')
        c = gather(c);
    end
    if isequal(class(X),'gpuArray')
        X = gather(X);
    end
end

sX = size(X);
sz = sX(1);
nX = prod(sX(3:end)); % the number of targeted matrices
lX = length(sX);
if nX>1 && prod(sc(3:end))>1 && sc(2)>1 && optim_c
    error('"c" needs to be in a column for each matrix calculation, that is, "c" is (poly_order,1,...), under optimization.');
end
if sc(1)==1 && sc(2)~=1
    c = permute(c,[2 1 3:lX]);
    [sc(2),sc(1)] = deal(sc(1),sc(2));
end
if ~ascending_c % change it to ascending power series
    if sc(2)>1
        if use_gpu && length(sc)>2
            c = pagefun(@rot90,c,2);
        else
            c = rot90(c,2);
        end
    else
        c = flipud(c);
    end
end

p = find(c(:,:,1),1,'last')-1; % polynomial order or []

% p = 0,1
if isempty(p)
    % c is all-zero or empty.
    if use_gpu
        Y = zeros(sX,'gpuArray');
    else
        Y = zeros(sX);
    end
    nmult = 0;
    Xpwr = [];
    return
elseif p==0
    % Y = c(1)*I.
    Y = repmat(diag(repmat(c(1),1,sz)),[1,1,sX(3:end)]);
    if use_gpu
        Y = gpuArray(Y);
    end
    nmult = 0;
    Xpwr = repmat(eye(sz,sz),[1,1,sX(3:end)]);
    return
end

if optim_c
    % Modify c by reshaping and zero-padding to minimize the number of matrix
    % multiplies.
    len = p+1;
    c = reshape(c(1:len,:),[len,1,sc(3:end)]); % column vector, truncate trailing zeros
    % c will be partitioned into length-m columns, with the last column
    % zero-padded.
    m = 1:len; % all possible m values
    n = ceil(len./m); % number of columns for each m
    % Count the number of maxtix mult's for evaulating
    %   Y = sum|k=1:n (sum|j=1:m X^(j-1)*c(j,k))*(X^m)^(k-1)
    % For powers X^2, X^3, ... X^m: count = m-1
    % For outer sum (Horner's method): count = n-1
    % Two special cases: (1) If c is a column vector, then the number of
    % mult's is reduced by 1 because X^m is not needed. (2) If c has more than
    % 1 column and the last column is all zeros except for the first element,
    % then the number of mult's is reduced by 1 because the first step of
    % Horner's method only involves a scalar multiply.
    nmult = (m-1)+(n-1)-(n==1 | (n>1 & mod(len,m)==1));
    % Pick m to minimize nmult.
    m = find(nmult==min(nmult),1);
    n = n(m);
    nmult = nmult(m);
    c = cat(1,c,zeros([m*n-len,1,sc(3:end)]));
    c = reshape(c,[m,n,sc(3:end)]);
else
    m = sc(1);
    n = sc(2);
    nmult =(m-1)+(n-1);
    if n==1 || (n>1 && all(c(2:end,end,:)==0))
        nmult = nmult-1;
    end
end

if isempty(X)
    Y = [];
    return
end
% Compute needed X powers, X^j = Xpwr{j+1}, j = 1:m. (Unused Xpwr entries are
% left empty.)
sXpwr = size(Xpwr);
get_elements = arrayfun(@(x) 1:x,sX,'UniformOutput',false);
if isempty(Xpwr)
    if use_gpu
        Xpwr = zeros([m+1-(n==1),sX],'gpuArray');
    else
        Xpwr = zeros([m+1-(n==1),sX]);
    end
    % Xpwr(1,...) is a diagonal matrix
    if nX>1
        get_ndim = arrayfun(@(x) (0:x-1)',sX,'UniformOutput',false);
        diag_idx = (1:sz+1:sz^2)';
        previous_ndim = cumprod(sX);
        for i = 3:length(get_elements)
            next_ndim = shiftdim((get_ndim{i})*previous_ndim(i-1),-(i-2));
            diag_idx = bsxfun(@plus,diag_idx,next_ndim);
        end
    else
        diag_idx = (1:sz+1:sz^2)';
    end
    Xpwr(1,diag_idx(:)) = 1;
    
    for j = 2:m+1-(n==1)
        if use_gpu
            Xpwr(j,get_elements{:}) = pagefun(@mtimes, X,reshape(Xpwr(j-1,get_elements{:}),sX));
        else
            for i = 1:nX
                Xpwr(j,:,:,i) = shiftdim(X(:,:,i)*squeeze(Xpwr(j-1,:,:,i)),-1);
            end
        end
    end
elseif sXpwr(1) ~= m+1-(n==1)
    if sXpwr(1) > m+1-(n==1)
        Xpwr = Xpwr(1:(m+1-(n==1)),get_elements{:});
    else
        Xpwr = cat(1,Xpwr,zeros([m+1-(n==1)-sXpwr(1),sX]));
        for j = sXpwr(1)+1:m+1-(n==1)
            if use_gpu
                Xpwr(j,get_elements{:}) = pagefun(@mtimes, X,reshape(Xpwr(j-1,get_elements{:}),sX));
            else
                for i = 1:nX
                    Xpwr(j,:,:,i) = shiftdim(X(:,:,i)*squeeze(Xpwr(j-1,:,:,i)),-1);
                end
            end
        end
    end
end
X_ = reshape(Xpwr(end,get_elements{:}),sX);
if n==1 || m==1
    Xpwr_ = Xpwr;
else
    Xpwr_ = Xpwr(1:end-1,get_elements{:});
end

% Compute the sum
if n==1
    m_ = find(c(:,1,1),1,'last');
    if isempty(m_)
        if use_gpu
            Y = zeros(sX,'gpuArray');
        else
            Y = zeros(sX);
        end
    elseif m_==1 % Y = I*c
        if nX>1
            Y = cell2mat(arrayfun(@(cx) diag(repmat(cx,1,sz)), c,'UniformOutput',false));
        else
            Y = diag(repmat(c(1),1,sz));
        end
    else
        c_ = reshape(c(:,1,:),[m,1,1,sc(3:end)]);
        Y = reshape(sum(c_.*Xpwr_),sX);
    end
else % n>1
    % k=n
    m_ = find(c(:,n,1),1,'last');
    if isempty(m_) % Y{k} = [] (alias for all-zeros)
        if use_gpu
            Y = zeros(sX,'gpuArray');
        else
            Y = zeros(sX);
        end
    elseif m_==1
        if length(sc) == 2
            Y = c(1,n).*X_;
        else
            Y = bsxfun(@times, c(1,n,get_elements{3:end}),X_);
        end
    else
        c_ = reshape(c(:,n,:),[m,1,1,sc(3:end)]);
        Y = reshape(sum(bsxfun(@times, c_,Xpwr_)),sX);
        if use_gpu
            Y = pagefun(@mtimes, X_,Y);
        else
            try
                Y = mmx_mult(X_,Y);
            catch
                for i = 1:nX
                    Y(:,:,i) = X_(:,:,i)*Y(:,:,i);
                end
            end
        end
    end
    % k<n
    for k = n-1:-1:1
        m_ = find(c(:,k,1),1,'last');
        if m_==1
            Y = reshape(c(1,k,:),[1,1,sc(3:end)]).*repmat(eye(sz,sz),1,1,sX(3:end)) + Y; % Y = c(1,k,1)*eye(sz,sz) + Y;
        else
            c_ = reshape(c(:,k,:),[m,1,1,sc(3:end)]);
            Y = reshape(sum(bsxfun(@times, c_,Xpwr_)),sX) + Y;
        end
        if k>1
            if use_gpu
                Y = pagefun(@mtimes, X_,Y);
            else
                try
                    Y = mmx_mult(X_,Y);
                catch
                    for i = 1:nX
                        Y(:,:,i) = X_(:,:,i)*Y(:,:,i);
                    end
                end
            end
        end
    end
end