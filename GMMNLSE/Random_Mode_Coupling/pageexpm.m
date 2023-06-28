function [E,jPwr] = pageexpm(D,options)
%PAGEEXPM Paged matrix exponential
%
%   For a square matrix argument D, pageexpm(D) is equivalent to expm(D)
%   (except that pageexpm uses a computational algorithm different from
%   MATLAB's expm function).
%
%   For a general array D (any number of dimensions) with
%   size(D,1)==size(D,2), the result E = pageexpm(D) is size-matched to D
%   with
%     E(:,:,j) = expm(D(:,:,j))
%   pageexpm is implemented using paged matrix operations (pagemtimes,
%   pagemldivide) for efficient operation on multicore and multiprocessor
%   computers:
%     https://www.mathworks.com/discovery/matlab-multicore.html
%
%   pageexpm requires MATLAB R2022a or later.
%
% Syntax:
%   E = pageexpm(D);
%   E = pageexpm(D,options);
%
% Inputs:
%
%   D: numeric array with size(D,1)==size(D,2)
%
%   options (optional, default = struct): struct with the following fields,
%   all optional:
%
%     options.minusI: true or false (default = false). If true, the result
%     will have an identity matrix subtracted off, equivalent to
%       E = pageexpm(D)-eye(size(D,1))
%     but without precision loss in the E diagonal elements E(j,j,:). This
%     is the matrix analog of function expm1. Use this option to avoid
%     precision loss with small arguments (e.g. for finite-difference
%     differentiation or numerical integration by parts).
%
%     options.RelTol: relative error tolerance, positive real scalar
%     (default = eps). For a matrix result E the error's norm is of order
%     RelTol*norm(E,'fro').
%
%     options.PadeOrder: Padé approximation order, positive odd integer
%     (default = 7).
%
% Outputs:
%
%   E: matrix exponential of D, numeric array, size-matched to D.
%
%   jPwr: scaling power (number of squaring steps used in scale-and-square
%   algorithm).

% Version 16-Mar-2022
% Copyright (c) 2022, Kenneth C. Johnson.
% KJ Innovation

% Documentation reference:
%   "Matrix Exponential Computational Algorithm" (15-Mar-2022)
%   https://vixra.org/abs/2203.0072
% Note: Equation numbers("Eq ...") in the code comments are from the
% above-cited reference.

if ~isfloat(D)
    error('pageexpm:inputType','Input must be single or double array.');
end

size_ = size(D);
size1 = size_(1);
if size1~=size_(2)
     error('pageexpm:inputMustBeSquare', ...
         ['Input D must be square in dimensions 1 and 2 ' ...
         '(size(D,1)==size(D,2).']);
end
size_(1:2) = [];

minusI = false;
RelTol = eps;
n = 7; % PadeOrder
if nargin>=2
    fieldcount = 0;
    if isfield(options,'minusI')
        fieldcount = fieldcount+1;
        minusI = options.minusI;
        if ~(islogical(minusI) && isscalar(minusI))
            error('pageexpm:invalidOption', ...
                'options.minusI must be true or false (logical scalar).');
        end
    end
    if isfield(options,'RelTol')
        fieldcount = fieldcount+1;
        RelTol = options.RelTol;
        if ~(isfloat(RelTol) && isreal(RelTol) ...
                && isscalar(RelTol) && RelTol>0)
            error('pageexpm:invalidOption', ...
                'options.RelTol must be a positive real scalar.');
        end
    end
    if isfield(options,'PadeOrder')
        fieldcount = fieldcount+1;
        n = options.PadeOrder;
        if ~(isnumeric(n) && isreal(n) && isscalar(n) ...
                && n==fix(n) && n>0 && mod(n,2)==1)
            error('pageexpm:invalidOption', ...
                'options.PadeOrder must be a positive odd integer.');
        end
    end
    if fieldcount<length(fieldnames(options))
        name = setdiff(fieldnames(options), ...
            {'minusI','RelTol','PadeOrder'});
            warning('pageexpm:unrecognizedOption', ...
                ['Unrecognized options field: ''' name{1} '''']);
    end
end

if isempty(D)
    E = D;
    jPwr = 0;
    return
end

if size1==1
    if minusI
        E = expm1(D);
    else
        E = exp(D);
    end
    jPwr = 0;
    return
end

% Get the P polynomial coefficients c, separated into even-order cEven and
% odd-order cOdd, for Padé order n. (Eqs 10, 11, 53)
[cEven,cOdd] = PadeCoef(n); % local function

[N,M] = size(cEven); % size(cOdd) is the same.
% Pre-compute the even powers of D in the P polynomial: DD{j} = (D^2)^j
% (for use in local function polyval_; see commnet header in polyval_).
if M>1
    DD{N} = []; % allocation
    DD{1} = pagemtimes(D,D);
    for j = 2:N
        DD{j} = pagemtimes(DD{1},DD{j-1}); % (D^2)*j
    end
elseif N>1
    % Case M==1: DD{N} is not needed.
    DD{N-1} = []; % allocation (omit DD{N})
    DD{1} = pagemtimes(D,D);
    for j = 2:N-1
        DD{j} = pagemtimes(DD{1},DD{j-1});
    end
else % M==1, N==1
    % n==1. D^2 powers are not needed for polynomial evaluation (Eqs 55),
    % but include D^2 in DD for Eqs 50 and 46 (jPwr initialization).
    DD = {pagemtimes(D,D)};
end
if any(isnan(DD{end}) | isinf(DD{end}))
    error('pageexpm:Overflow', ...
        'Numeric overflow from large D and/or large options.PadeOrder.')
end
% normD and normDD are used in Eqs 50 and 46.
normDD = sqrt(sum(real(DD{1}).^2+imag(DD{1}).^2,[1,2])); ...
    % = norm(DD{1},'fro') if ismatrix(DD{1})
if all(normDD(:)==0)
    % D^2 is zero (although D might be nonzero).
    E = eye(size1)+D;
    jPwr = 0;
    return
end
normD = sqrt(normDD);
% normD and normDD are size-[1,1,...]
normDpwr = []; % place-holder for leading factor in Eq 50

% Apply scale-and-square algorithm with jPwr = scaling power (i.e., number
% of squaring steps).
jPwr = calc_jPwr; % nested function
h_ = 2^(-jPwr-1); % 'h_' avoids conflict with 'h' in nested functions.
% Calculate PEven, POdd (Eqs 57).
if n>1
    % Incorporate h_^(2*j) factor in DD{j}.
    hhDD = DD;
    hh_ = h_^2;
    hhDD{1} = hh_*DD{1}; % (h_*D)^2
    hhj = hh_;
    for j = 2:length(DD)
        hhj = hh_*hhj;
        hhDD{j} = hhj*DD{j}; % (h_*D)^(2*j)
    end
    % hhDD{j} = (h_*D)^(2*j).
    %
    % Separately evaluate even and odd parts of P polynomial (PEven, POdd),
    % Eq's 55. Initially omit a factor of h_*D in POdd.
    PEven = polyval_(cEven,hhDD); % local function
    POdd = polyval_(cOdd,hhDD);
else % n==1
    % P = P1 in Eq 12.
    PEven = eye(size1,size1);
    if ~ismatrix(D)
        PEven = repmat(PEven,[1,1,size_]);
    end
    POdd = PEven; % without h_*D factor
end
% Put missing factor of h_*D in POdd.
POdd = pagemtimes(h_*D,POdd);

% Calculate Phi from Eq 17. Phi_ is Phi-diag(d).
Phi_ = pagemldivide(PEven-POdd,2*POdd);
% Eq 18:
if isempty(size_) % ismatrix(D)
    d = ones(size1,1); % Phi = Phi_+diag(d)
    for j = 1:jPwr
        d_ = diag(Phi_)+d;
        Phi_ = Phi_-diag(d_-d);
        d = d_;
        Phi_ = Phi_^2+d.*Phi_+Phi_.*d.'; % Phi_^2+diag(d)*Phi_+Phi_*diag(d)
        d = d.^2;
    end
    if minusI
        d = d-1;
    end
    E = Phi_+diag(d);
else
    d = ones([size1,1,size_]);
    idiag = (1:size1+1:size1^2).'+(0:prod(size_)-1)*size1^2;
    idiag = idiag(:);
    for j = 1:jPwr
        d_ = reshape(Phi_(idiag)+d(:),[size1,1,size_]);
        Phi_(idiag) = Phi_(idiag)-(d_(:)-d(:));
        d = d_;
        Phi_ = pagemtimes(Phi_,Phi_)+d.*Phi_+Phi_.*pagetranspose(d);
        d = d.^2;
    end
    if minusI
        d = d-1;
    end
    E = Phi_;
    E(idiag) = E(idiag)+d(:);
end

    function jPwr = calc_jPwr
        % Determine the scaling power (number of squaring operations, p in
        % Eq 3).
        %
        % The following shared variables are used by calc_jPwr: n, D, DD,
        % normD, normDpwr, cEven, cOdd, RelTol. D and DD are only used if
        % normDpwr is empty, in which case normDpwr will be initialized by
        % calc_jPwr.
        %
        jPwr = 0;
        log1p_RelTol = log1p(RelTol);
        % Increase jPwr, if necessary, until the two inequalities in Eq 46
        % hold.
        if in_tolerance(jPwr) % nested function
            return
        end
        % in_tolerance(jPwr)==false
        jPwr = [jPwr,jPwr+1];
        flag = false; % will be true when in_tolerance(jPwr(2))==true
        while true
            % Current state:
            % * jPwr(1)<jPwr(2).
            % * in_tolerance(jPwr(1))==false.
            % * If flag, then in_tolerance(jPwr(2))==true; otherwise
            % in_tolerance(jPwr(2)) is unknown.
            if flag
                % in_tolerance(jPwr(1))==false, in_tolerance(jPwr(2))==true
                if jPwr(2)==jPwr(1)+1
                    jPwr = jPwr(2);
                    return
                end
                jPwr_ = round((jPwr(1)+jPwr(2))/2);
                if in_tolerance(jPwr_)
                    jPwr(2) = jPwr_; % in_tolerance(jPwr(2))==true
                else
                    jPwr(1) = jPwr_; % in_tolerance(jPwr(1))==false
                end
            elseif in_tolerance(jPwr(2))
                flag = true;
            else % in_tolerance(jPwr(2))==false; extend search range.
                jPwr = [jPwr(2),jPwr(2)+2*(jPwr(2)-jPwr(1))];
                % in_tolerance(jPwr(1))==false.
            end
        end
        function test = in_tolerance(jPwr)
            h = 2^(-jPwr-1);
            norm_hhDD = h^2*normDD; % size-[1,1,...]
            norm_hD = abs(h*normD); % size-[1,1,...]
            % norm_hhDD = norm_hD.^2
            % Calculate the denominator factor in Eq 46. With X = norm_hD,
            %   den = 2-PEven(1i*X)^2+POdd(1i*X)^2
            % Eqs 55 with X replaced by 1i*X, leading i*X factor in POdd_
            % omitted:
            %   PEven_ = polyval(cEven(end:-1:1),-norm_hhDD);
            %   POdd_ = polyval(cOdd(end:-1:1),-norm_hhDD);
            % Don't use polyval; this is faster:
            PEven_ = [ones(size(norm_hhDD)), ...
                cumprod(repmat(-norm_hhDD,1,numel(cEven)-1),2)]; ...
                % [1,(1i*X)^2,(1i*X)^4,...]
            POdd_ = pagemtimes(PEven_,cOdd(:)); ...
                % cOdd(1)+cOdd(2)*(1i*X)^2+cOdd(3)*(1i*X)^4+...
            PEven_ = pagemtimes(PEven_,cEven(:)); ...
                % cEven(1)+cEven(2)*(1i*X)^2+cEven(3)*(1i*X)^4+...
            % The -norm_hhDD argument above is (1i*X)^2. POdd_ is missing a
            % leading factor of 1i*X, which is accounted for by the
            % -norm_hhDD factor in the following line.
            den = 2-PEven_.^2-norm_hhDD.*POdd_.^2;
            % The first inequality in Eq 46 requires den>0, but use a
            % tighter limit den>=0.1 to ensure that the reciprocal value is
            % not very large.
            if any(den<0.1)
                test = false;
                return
            end
            if isempty(normDpwr)
                % Calculate normDpwr = norm(D^(2*n+1)) (bounding estimate)
                % via Eq 60. Use pre-calculated powers DD{j} = (D^2)^j;
                % use the following identity:
                %   norm(D^(j1+j2+...))<=norm(D^j1)*norm(D^j2)*...,
                %   j1+j2+... = 2*n+1.
                % Initially calculate a bound for norm((D/normD)^(2*n+1))
                % to guard against numeric overflow, then multiply by
                % normD^(2*n+1).
                normDpwr = 1;
                pwr = 0;
                for j_= length(DD):-1:1
                    norm_ = [];
                    while pwr+j_<=2*n
                        if isempty(norm_)
                            norm_ = DD{j_}./normDD.^j_; ...
                                % possible divide-by-zero
                            norm_ = sqrt(sum( ...
                                real(norm_).^2+imag(norm_).^2,[1,2])); ...
                                % norm(DD{j_}/normDD^j_,'fro')
                        end
                        normDpwr = normDpwr.*norm_;
                        pwr = pwr+2*j_;
                    end
                end
                norm_ = D./normD;
                norm_ = sqrt(sum( ...
                    real(norm_).^2+imag(norm_).^2,[1,2])); ...
                    % norm(D/normD,'fro')
                normDpwr = normDpwr.*norm_;
                normDpwr = normDpwr.*normD.^(2*n+1);
                normDpwr(~isfinite(normDpwr(:))) = 0; % for case D==0
                % norm(D^(2*n+1),'fro')<=normDpwr.
            end
            % Calculate normDelta limit (Eq 50) and norm_err (Eq 46).
            % Eqs 55:
            %   PEven_ = polyval(cEven(end:-1:1),norm_hD^2);
            %   POdd_ = norm_hD*polyval(cOdd(end:-1:1),norm_hD^2);
            % Don't use polyval; this is faster:
            PEven_ = [ones(size(norm_hD)), ...
                cumprod(repmat(norm_hhDD,1,numel(cEven)-1),2)]; ...
                % [1,norm_hD^2,norm_hD^4,...]
            POdd_ = norm_hD.*pagemtimes(PEven_,cOdd(:)); ...
                % norm_hD*(cOdd(1)+cOdd(2)*norm_hD^2+cOdd(3)*norm_hD^4+...)
            PEven_ = pagemtimes(PEven_,cEven(:)); ...
                % cEven(1)+cEven(2)*norm_hD^2+cEven(3)*norm_hD^4+...
            % Eq 50:
            cosh_ = cosh(norm_hD);
            sinh_ = sinh(norm_hD);
            normDelta = 2*abs(h).^(2*n+1).*normDpwr.*cosh_/ ...
                ((2*n+1)*prod(2*n-1:-2:3)^2);
            % Eq 46:
            norm_err = 0.5*(1+den.\(1+ ...
                (cosh_-PEven_).^2+(sinh_-POdd_).^2+normDelta)).*normDelta;
            % Eqs 22, 24:
            test = all(norm_err(:)<=2^-jPwr*log1p_RelTol);
        end % in_tolerance
    end % calc_jPwr

end % pageexpm

function [cEven,cOdd,count] = PadeCoef(n)
% For Padé order n (must be odd), return the even-order and odd-order
% polynomial coefficients cEven and cOdd (Eq 53), and the numbber of matrix
% multiplies (count) required to evaluate the even and odd Padé polynomials
% (Eq 59).
%
% Use Eq 12 to get coefficients of the Padé polynomial P in Eq 10. The
% polynomial coefficients for Padé order n are stored in polynomial(n).c.
% For n odd, the corresponding even- and odd-order polynomial coefficients
% are in polynomial(n).cEven and polynomial(n).cOdd, and are zero-padded
% and reshaped to size-[N,M] (Eqs 57). The number of matrix multiplies
% required for evaluation of P is polynomial(n).count (Eq 59). (For n even,
% polynomial(n).cEven, polynomial(n).cOdd, and polynomial(n).count are all
% [].)
if mod(n,2)~=1
    % n must be odd (checked in input validation).
    error('pageexpm:logic','Internal logic error (1).');
end
persistent polynomial % Use "clear all" to clear polynomial.
if isempty(polynomial)
    % P coef's for Padé order n==1 and n==2 are [1;1] and [1;1;1/3],
    % respectively (Eq 12 recursion basis, don't use order 0).
    polynomial = struct( ...
        'c',{[1;1],[1;1;1/3]}, ...
        'cEven',{1,[]}, ...
        'cOdd',{1,[]}, ...
        'count',{1,[]});
end
n_ = length(polynomial);
if n_>=n
    cEven = polynomial(n).cEven;
    cOdd = polynomial(n).cOdd;
    count = polynomial(n).count;
    return
end
% n>n_>=2
polynomial(n).c = []; % cell allocation
while n_<n
    % Get coefficients for Padé order n_+1 from orders n_ and n_-1 (Eq 12
    % recursion).
    c = [polynomial(n_).c;0]+ ...
        [0;0;polynomial(n_-1).c/(4*n_^2-1)];
    n_ = n_+1;
    polynomial(n_).c = c;
    if mod(n_,2)==0
        continue
    end
    m = (n_-1)/2; % Eq 54
    % Separate c into even and odd polynomial coefficients.
    cEven = c(1:2:end);
    cOdd = c(2:2:end);
    % Determine N and M in Eq 57. Try all possible [N,M] combinations to
    % minimize the number of matrix multiplies.
    % Eqs 58:
    N = 1:m+1;
    M = ceil((m+1)./N);
    % Eq 59:
    count = 1 ...
        +max(0,N-2) ...
        +(M>1) ...
        +2*max(0,M-1) ...
        -2*(M>1 & N.*(M-1)==m) ...
        +(m>0);
    j = find(count==min(count),1,'last');
    N = N(j);
    M = M(j);
    count = count(j);
    NM = N*M;
    cEven(end+1:NM) = 0;
    cEven = reshape(cEven,N,M);
    cOdd(end+1:NM) = 0;
    cOdd = reshape(cOdd,N,M);
    polynomial(n_).cEven = cEven;
    polynomial(n_).cOdd = cOdd;
    polynomial(n_).count = count;
end
% n_==n, cEven==polynomial(n).cEven, cOdd==polynomial(n).cOdd,
% count==polynomial(n).count;
end % PadeCoef

function Y = polyval_(c,pwrX)
% c is a non-empty, non-scalar, size-[N,M] polynomial coefficient array
% (cEven or cOdd from PadeCoef). X is a multidimensional array with
% size(X,1)==size(X,2). pwrX is a non-empty cell vector storing X powers
% pwrX{j}(:,:,i) = X(:,:,i)^j, j = 1:N (or j = 1:N-1 if M==1). Evaluate
%   Y(:,:,i) = sum|k (sum|j c(j,k)*X(:,:,i)^(j-1)) * (X(:,:,i)^N)^(k-1);
%   j = 1:N, k = 1:M
% (Eqs 57). The inner j sum uses the precomputed X powers X(:,:,i)^(j-1);
% the outer sum is evaluated by Horner's method with precomputed
% X(:,:,i)^N. The sum is equivalent to
%   Y(:,:,i) = sum|j c(j)*X(:,:,i)^(j-1); j = 1:N*M
% but with fewer matrix multiplies when M>1.
[N,M] = size(c);
size_ = size(pwrX{1});
NX = size_(1); % size_ = [NX,NX,...]
if size_(2)~=NX
    error('pageexpm:logic','Internal logic error (2).');
end
size_(1:2) = 1; % size_ = [1,1,...]
Y = []; % partial sum over k ([] implicitly means all-zero)
for k = M:-1:1
    % Evaluate Y_(:,:,i) = sum|j c(j,k)*X(:,:,i)^(j-1).
    Y_ = []; % partial sum over j
    for j = N:-1:2
        c_ = c(j,k);
        if c_==0
            continue
        end
        if isempty(Y_)
            Y_ = c_*pwrX{j-1};
        else
            Y_ = Y_+c_*pwrX{j-1};
        end
    end
    c_ = c(1,k);
    % Increment Y_ by c_*eye(NX,NX).
    if c_~=0
        if isempty(Y_)
            % Y_ = c_*eye(NX,NX):
            Y_ = diag(repmat(c_,NX,1));
            if length(size_)>2
                Y_ = repmat(Y_,size_);
            end
        else
            % Y_ = Y_+c_*eye(NX,NX):
            Y_ = Y_+diag(repmat(c_,NX,1));
        end
    end
    % Increment Y by Y_.
    if isempty(Y)
        Y = Y_;
    elseif ~isempty(Y_)
        Y = Y+Y_;
    end
    if k>1
        % Multiply Y by X^N (Horner's method for sum|k ...).
        if ~isempty(Y)
            Y = pagemtimes(Y,pwrX{N});
        end
    end
end
if isempty(Y)
    Y = zeros([NX,NX,size_(3:end)]);
end
end % polyval_

