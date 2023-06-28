function [expA,m,s] = myexpm_(A,m,s,precise,minusI,use_gpu)
%myEXPM_ Matrix exponential algorithm with improved numerical precision which
%can use CPU or GPU for computations and can also compute for "multiple"
%matrices "at once". while matrices are stored as pages of a multidimensional array.
%
% myexpm_(A) calculates the matrix exponential of square matrices on "each 
% page" of a multidimensional array A, that is, A has the size
% (n,n,?,?,...), using a scaling and squaring algorithm with a Pade  
% approximation. An enhancement over the standard expm implementation  
% preserves numerical precision for arbitrarily small scaling factors.
%
% For CPU computing, it relies on "mmx" by Yuval on MATLAB File Exchange
% for a faster matrix multiplication. Without it, this code use "for-loop"
% for matrix multiplications.
% https://www.mathworks.com/matlabcentral/fileexchange/37515-mmx-multithreaded-matrix-operations-on-n-d-matrices
%
% Syntax:
%   expA = myexpm_(A);
%   expA = myexpm_(A,m,s,precise);
%   expA = myexpm_(A,m,s,precise,minusI);
%   expA = myexpm_(A,m,s,precise,minusI,use_gpu);
%   [expA,m,s] = myexpm_(A);
%
% inputs:
%
%   A: a multidimensional array with each page a square matrix (real or
%   complex). Size: (n,n,...)
%
% optional inputs:
%
%   m: positive integer or [], default = [].
%   Order of Pade approximation. If m is unspecified or [] it will be
%   determined automatically. (See Note.)
%
%   s: non-negative integer or [], default = [].
%   Number of squaring steps. (The scale factor is 2^-s.) If s is unspecified
%   or [] it will be determined automatically. (See Note.)
%
%   precise: true or false; default = true.
%   If false, the standard algorithm (no precision enhancement) is used.
%
%   minusI: true or false; default = false.
%   If true, the expA result is exp(A)-I (i.e., with identity matrix
%   subtracted). This is analogous to MATLAB's expm1 function. (The 'm1' means
%   'minus 1' - expm1 is unrelated to matrix exponentiation.)
%
%   use_gpu: true or false; default = false.
%   If true, use GPU for computations, otherwise use CPU.
%
% Note:
%   The determination criteria for m and s are not the same as expm; set m and
%   s manually to match expm. If only one of m or s is specified, the other
%   will be determined to achieve a relative approximation error at the
%   numerical precision limit (eps). If both are unspecified, they will be
%   determined to achieve the target approximation error at minimum cost in
%   terms of the number of matrix multiply operations.
%
% outputs:
%
%   expA: exponential of A, or exp(A)-I if minusI is true.
%
%   m, s: same as inputs (or auto-generated values)
%
% Test cases (using MATLAB R2017a):
%
% Test case 1
% Adapted from Eq 1.2 in Awad Al-Mohy and Nicholas Higham, A New Scaling and
% Squaring Algorithm for the Matrix Exponential, SIMAX 31, 970-989, 2009.
%
%   >> b = 1e8; A = [1,b;0,-1]; expA = [exp(1),b*sinh(1);0,exp(-1)];
%   >> myexpm(A)-expA
%   ans =
%                      0                   0
%                      0                   0
%   >> myexpm_(A,13,25,false)-expA
%   ans =
%      1.0e-03 *
%     -0.000000026279423  -0.353783369064331
%                      0   0.000000003573086
%   >> myexpm_(A,13,25)-expA
%   ans =
%      1.0e-07 *
%                      0   0.447034835815430
%                      0                   0
%   >> myexpm_(A)-expA
%   ans =
%      1.0e-15 *
%     -0.444089209850063                   0
%                      0                   0
%   >> myexpm_(A,[],[],true,false,true)-expA % use GPU
%   ans =
%      1.0e-07 *
%                      0   0.298023223876953
%                      0                   0
%
%   >> [~,m,s] = myexpm_(A)
%   m =
%       14
%   s =
%       53
%
% Test case 2
% Adapted from "A Balancing Act for the Matrix Exponential"
% Posted by Cleve Moler, July 23, 2012
% http://blogs.mathworks.com/cleve/2012/07/23/a-balancing-act-for-the-matrix-exponential/
%   >> a = 2e10; b = 4e8/6; c = 200/3; d = 3; e = 1e-8;
%   >> A = [0 e 0; -(a+b) -d a; c 0 -c]; expA = expmdemo3(A);
%   >> num2str(expm(A)-expA)
%   ans =
%    ' 1.5266e-14  2.5849e-23   7.494e-15'
%    '-8.2888e-07 -2.5865e-15 -7.7486e-07'
%    ' 1.4766e-14  2.4195e-23  7.0499e-15'
%   >> num2str(myexpm_(A)-expA)
%   ans =
%    '8.2712e-15  1.5716e-23   4.996e-15'
%    '6.7148e-07  2.2482e-15  7.1246e-07'
%    '7.7716e-15  1.4269e-23  4.6074e-15'
%   >> [~,m,s] = myexpm_(A)
%   m =
%       15
%   s =
%       69
%
%
% Author: Ken Johnson  kjinnovation@earthlink.net
% Originally posted 2/28/2014; last code revision 3/6/2014.
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
% =========================================================================
% Author: Yi-Hao Chen
% email: jollycyclopes@gmail.com
% date: 6/6/2018
%
% Acknowledgement:
%   This code is modified from "expm_" by Ken Johnson.
%
% Examples:
%   >> A = rand(10,10,3,6,2)+1i*rand(10,10,3,6,2);
%   >> expA = myexpm_(A);
%   >> expA = myexpm_(A,[],[],true,false,true); % use GPU with "precise"
%   turned on
%
%   Result:
%       size(expA) = [10,10,3,6,2]
%       expA(:,:,1,1,1) = expm(A(:,:,1,1,1))
%

sA = size(A);
sz = sA(1);
nA = prod(sA(3:end)); % the number of targeted matrices
slogA = [1 1 sA(3:end)];

if nargin<6
    use_gpu = false;
    if nargin<5
        minusI = false;
        if nargin<4
            precise = true;
            if nargin<3
                s = [];
                if nargin<2
                    m = [];
                end
            end
        end
    end
end

if use_gpu && ~isequal(class(A),'gpuArray')
    A = gpuArray(A);
elseif ~use_gpu && isequal(class(A),'gpuArray')
    A = gather(A);
end
input_A_is_single = false;
if use_gpu
    Atype = classUnderlying(A);
else
    Atype = class(A);
end
if isequal(Atype,'single')
    input_A_is_single = true;
    A = double(A);
end

if isempty(s)
    if ~isempty(m)
        % m is specified; determine s to achieve the target approximation
        % error.
        sm = size(m);
        if ~isequal(sm(3:end),sA(3:end)) && ~isscalar(m)
            error('"m" should be of the same size, starting with the 3rd dimension, as "A" or just a scalar.');
        end
        if isscalar(m)
            m = m*ones(slogA);
        end
        % Use the Frobenius norm (root-sum-square) because it is symmetric
        % under matrix transposition (ensuring that expm_(A) is similarly
        % symmetric under transposition).
        %
        % log2_normA = log2(norm(A,'fro'));
        normA = sum(reshape(abs(A),[sz^2,1,sA(3:end)]).^2);
        % MATLAB "log2,sqrt" are slow with GPU!
        if use_gpu
            log2_normA = log2(sqrt(gather(normA)));
        else
            log2_normA = log2(normA);
        end
        if use_gpu && nA>1
            log2_normA = gpuArray(log2_normA);
        end
        s = ceil(calc_s(m,log2_normA));
    else
        % m and s are both unspecified; determine both to minimize cost (n
        % terms of matrix muliplies) while achieving the target approximation
        % error.
        %
        % log2_normA = log2(norm(A,'fro'));
        normA = sum(reshape(abs(A),[sz^2,1,sA(3:end)]).^2);
        % MATLAB "log2,sqrt" are slow with GPU!
        if use_gpu
            log2_normA = log2(sqrt(gather(normA)));
        else
            log2_normA = log2(normA);
        end
        if use_gpu
            if nA>1
                log2_normA = gpuArray(log2_normA);
            end
            i_ = true(slogA,'gpuArray');
        else
            i_ = true(slogA);
        end
        mj_ = 1;
        sj_ = calc_s(mj_,log2_normA); % fractional, to facilitate optimization
        cost = calc_cost(mj_,sj_(:));
        m = ones(slogA);
        if use_gpu
            s = zeros(slogA,'gpuArray');
        else
            s = zeros(slogA);
        end
        while mj_<100
            mj = mj_+1;
            sj = calc_s(mj,log2_normA(i_));
            cost_ = calc_cost(mj,sj(:));
            i = cost_<cost;
            if (nA==1 && ~i) || ~all(i)
                ii_ = i_;
                i_(i_) = i; % remained
                ii_(ii_) = ~i; % finished
                m(ii_) = mj_;
                s(ii_) = sj_(~i);
            end
            if (nA==1 && ~i) || ~any(i)
                break
            end
            
            mj_ = mj;
            sj_ = sj(i);
            if any(~i)
                cost = cost_(i);
            else
                cost = cost_;
            end
        end
        s = ceil(s);
    end
elseif isempty(m)
    % s is specified; determine m to achieve the target approximation error.
    ss = size(s);
    if ~isequal(ss(3:end),sA(3:end)) && ~isscalar(s)
        error('"s" should be of the same size, starting with the 3rd dimension, as "A" or just a scalar.');
    end
    %
    % log2_normA = log2(norm(A,'fro'));
    normA = sum(reshape(abs(A),[sz^2,1,sA(3:end)]).^2);
    % MATLAB "log2" and "sqrt" are slow with GPU!
    if use_gpu
        log2_normA = log2(sqrt(gather(normA)));
    else
        log2_normA = log2(normA);
    end
    if use_gpu && nA>1
        log2_normA = gpuArray(log2_normA);
    end
    log2_eps = log2(eps);
    mj = 1;
    m = ones(slogA);
    i_ = true(slogA);
    while mj<100
        if isscalar(s)
            i = calc_log2_err(mj+1,s,log2_normA(i_))>log2_eps;
        else
            i = calc_log2_err(mj+1,s(i_),log2_normA(i_))>log2_eps;
        end
        if (nA==1 && ~i) || ~all(i)
            ii_ = i_;
            i_(i_) = i; % remained
            ii_(ii_) = ~i; % finished
            m(ii_) = mj;
        end
        if (nA==1 && ~i) || ~any(i)
            break
        end
            
        mj = mj+1;
    end
else
    if isscalar(m)
        m = m*ones(slogA);
    end
end
A = bsxfun(@rdivide, A,2.^s);

% Create an identity matrix
% I = repmat(eye(sz,sz),[1 1 sA(3:end]) is slow if there are lots of
% matrices.
if nA>1
    if use_gpu
        I = zeros(sA,'gpuArray');
    else
        I = zeros(sA);
    end
    get_ndim = arrayfun(@(x) (0:x-1)',sA,'UniformOutput',false);
    diag_idx = (1:sz+1:sz^2)';
    previous_ndim = cumprod(sA);
    for i = 3:length(get_ndim)
        next_ndim = shiftdim((get_ndim{i})*previous_ndim(i-1),-(i-2));
        diag_idx = bsxfun(@plus,diag_idx,next_ndim);
    end
    I(diag_idx(:)) = 1;
else
    if use_gpu
        I = eye(sA,'gpuArray');
    else
        I = eye(sA);
    end
end

% Pade approximation: expA = p(A,m)/p(-A,m). Compute the polynomial
% coefficients b(0), b(1) ... b(m) of p(A,m) from the recursion relations:
%   p(A,0) = I
%   p(A,1) = I+A/2
%   p(A,m) = p(A,m-1)+p(A,m-2)*A^2/(4*(2*m-1)*(2*m-3))
mu = unique(m,'sorted')';
max_m = mu(end);
b1_ = 1; % p(A,0) coef's [b(0)] ('1' connotes 1-based indexing)
b1 = [1,1/2]; % p(A,1) coef's [b(0),b(1)]
for j = 2:max_m
    % [p(A,j-1),p(A,j)] coef's:
    [b1_,b1] = deal(b1,[b1,0]+[0,0,b1_]/(4*(2*j-1)*(2*j-3)));
end
% Calculate even- and odd-order polynomials:
%   p_odd = (b(1)*I+b(3)*A^2+...)*A
%   p_even = b(0)+b(2)*A^2+b(4)*A^4+...
if use_gpu
    AA = pagefun(@mtimes, A,A);
else
    AA = zeros(sA);
    try % use "mmx" by Yuval on MATLAB File Exchange
        AA = mmx_mult(A,A);
    catch
        for k = 1:nA
            AA(:,:,k) = A(:,:,k)*A(:,:,k);
        end
    end
end
if use_gpu
    p_even = zeros(sA,'gpuArray');
    p_odd = zeros(sA,'gpuArray');
else
    p_even = zeros(sA);
    p_odd = zeros(sA);
end
for k = mu
    mi = (m==k);
    AAk = AA(:,:,mi);
    [pe,~,~,AApwr] = mypolyvalm_(b1(1:2:k+1),AAk,true,true,[],use_gpu);
    po = mypolyvalm_(b1(2:2:k+1),AAk,true,true,AApwr,use_gpu);
    if use_gpu
        Ak = A(:,:,mi);
        po = pagefun(@mtimes, po,Ak);
        po = reshape(po,sz,sz,[]);
    else
        try % use "mmx" by Yuval on MATLAB File Exchange
            Ak = A(:,:,mi);
            po = mmx_mult(po,Ak);
            po = reshape(po,sz,sz,[]);
        catch
            j = 1;
            mk = find(mi)';
            for h = mk
                po(:,:,j) = po(:,:,j)*A(:,:,h);
                j = j+1;
            end
        end
    end

    p_even(:,:,mi) = reshape(pe,sz,sz,[]);
    p_odd(:,:,mi) = po;
end
if ~precise
    % Standard implementation of Pade approximation.
    if use_gpu
        expA = pagefun(@mrdivide, p_even+p_odd, p_even-p_odd);
    else
        try % use "mmx" by Yuval on MATLAB File Exchange
            expA = multslash(p_even-p_odd,p_even+p_odd);
        catch
            expA = zeros(sA);
            for k = 1:nA
                expA(:,:,k) = (p_even(:,:,k)+p_odd(:,:,k))/(p_even(:,:,k)-p_odd(:,:,k));
            end
        end
    end
    % Undo scaling by repeated squaring
    max_s = max(s(:));
    if use_gpu
        for cs = 1:max_s
            ri = find(s>=cs);
            expA(:,:,ri) = pagefun(@mtimes, expA(:,:,ri),expA(:,:,ri));
        end
    else
        try % use "mmx" by Yuval on MATLAB File Exchange
            for cs = 1:max_s
                ri = find(s>=cs);
                expA(:,:,ri) = mmx_mult(expA(:,:,ri),expA(:,:,ri));
            end
        catch
            for k = 1:nA
                expA(:,:,k) = expA(:,:,k)^(2^s(k));
            end
        end
    end
    if minusI
        expA = expA-I;
    end
else
    % Modified algorithm - separate out dominant I term from expA to avoid
    % precision loss.
    if use_gpu % exp(A)-I
        expA = 2*pagefun(@mrdivide, p_odd, p_even-p_odd);
    else
        try % use "mmx" by Yuval on MATLAB File Exchange
            expA = 2*multslash(p_even-p_odd,p_odd);
        catch
            expA = zeros(sA);
            for k = 1:nA
                expA(:,:,k) = 2*p_odd(:,:,k)/(p_even(:,:,k)-p_odd(:,:,k));
            end
        end
    end
    % Undo scaling by repeated squaring
    max_s = max(s(:));
    if use_gpu % (exp(A)^2-I) = (exp(A)-I)^2+2*(exp(A)-I)
        for cs = 1:max_s
            ri = find(s>=cs);
            expA(:,:,ri) = pagefun(@mtimes, expA(:,:,ri),expA(:,:,ri))+2*expA(:,:,ri);
        end
    else
        try % use "mmx" by Yuval on MATLAB File Exchange
            for cs = 1:max_s
                ri = find(s>=cs);
                expA(:,:,ri) = mmx_mult(expA(:,:,ri),expA(:,:,ri))+2*expA(:,:,ri);
            end
        catch
            for k = 1:nA
                if isscalar(s)
                    s = s*ones(slogA);
                end
                for j = 1:s(k)
                    expA(:,:,k) = expA(:,:,k)*expA(:,:,k)+2*expA(:,:,k);
                end
            end
        end
    end
    if ~minusI
        expA = expA+I;
    end
end

if input_A_is_single
    expA = single(expA);
end

function s = calc_s(m,log2_normA)
% Calculate the scaling power, s, for the scale-and-square algorithm. (s is
% returned as a fractional value to facilitate cost optimization.)
%
% Outline of methodology:
% Denote by r(x) the order-m Pade rational approximation to exp(x).
% For small x, the approximation error is estimated as
%   r(x)-exp(x) = c*x^(2*m+1)
% where
%   c = (-1)^(m+1)*((m!)^2/((2*m)!*(2*m+1)!))
% (from Eq 1.3 in Nicholas Higham, The Scaling and Squaring Method for the
% Matrix Exponential Revisited, SIAM REVIEW  Vol. 51 (2009), pp. 747-764.)
%
% The error in the n-th power of exp(x) is estimated as
%   r(x)^n-exp(x)^n = n*exp(x)^(n-1)*(r(x)-exp(x))
%   = n*exp(x)^(n-1)*c*x^(2*m+1)
% Set x = A*2^-s and n = 2^s:
%   r(A*2^-s)^(2^s)-exp(A) = exp(A-A*2^-s)*c*A^(2*m+1)*2^(-2*m*s)
% The relative error in the Pade scale-and-square approximant to exp(A) is
%   r(A*2^-s)^(2^s)*exp(-A)-I = exp(-A*2^-s)*c*A^(2*m+1)*2^(-2*m*s)
% The leading factor exp(-A*2^-s) on the right is close to I, and the
% remainder of the expression will be bounded by eps,
%   abs(c)*norm(A^(2*m+1))*2^(-2*m*s) <= eps
% Using the relation norm(A^(2*m+1)) <= norm(A)^(2*m+1), it suffices to
% require
%   abs(c)*norm(A)^(2*m+1)*2^(-2*m*s) <= eps
% s is defined by making the above relation an equality:
%   s = log2(abs(c)*norm(A)^(2*m+1)/eps)/(2*m)
if isscalar(m)
    log2_abs_c = -2*sum(log2(m+1:2*m))-log2(2*m+1);
else
    log2_abs_c = arrayfun(@(m) -2*sum(log2(m+1:2*m))-log2(2*m+1),m);
end
s = max(0,(log2_abs_c-log2(eps)+(2*m+1).*log2_normA)./(2*m));

end

function log2_err = calc_log2_err(m,s,log2_normA)
% Calculate an approximate bound on the log2 relative error of the Pade
% approximation to exp(A):
%   err = abs(c)*norm(A)^(2*m+1)*2^(-2*m*s)
% with c defined in calc_s.
log2_abs_c = arrayfun(@(m) -2*sum(log2(m+1:2*m))-log2(2*m+1),m);
log2_err = log2_abs_c+(2*m+1).*log2_normA-2*m.*s;

end

function cost = calc_cost(m,s)
% Estimate the cost, in terms of matrix multiplies, for calculataing the
% exponential matrix. There are two calls to polyvalm_ (for p_even and
% p_odd), each for a polynomial of order m/2 (approximately). The cost of
% each call is approximately 2*(sqrt(m/2+1)-1). There are two additional
% multiplies for AA = A*A, and p_odd = (...)*A. The squaring operation
% entails s multiplications.
cost = 4*(sqrt(m/2+1)-1)+2+s;

end

end