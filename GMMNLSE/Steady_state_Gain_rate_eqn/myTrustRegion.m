function xk = myTrustRegion(F,J,x0,N,e,M,x_total,gpu_yes)
%MYTRUSTREGION Trust-region method for solving the (num_x*num_y*M)'s 
%optimization problems in parallel. This parallelization is crucial in
%multimode rate-equation gain computations. MATLAB's fsolve() can solve
%only one problem, so I need to develop this code myself to solve multiple
%optimization problems in parallel.
%
%   F: coupled equation; it outupts a (n,1,num_x,num_y,1,1,1,M) array
%   J: Jacobian matrix; it outupts a (n,m,num_x,num_y,1,1,1,M) array
%   x0: initial guess for m's variables; (m,1,num_x,num_y,1,1,1,M) array
%   N: parameter for non-monotone function evaluations during optimization; a scalar
%   e: threshold value of stopping optimization processes; a scalar
%   M: the number of parallelizations in MPA
%   x_total: the total population; (1,1,num_x,num_y)
%   gpu_yes: use GPU or not
%
% Dimension explanation:
%   n: the number of coupled equations
%   m: the number of variables, which corresponds to the number of computed populations, N0, N1, N2, ...
%   num_x: the size of the "x" spatial dimension
%   num_y: the size of the "y" spatial dimension
%
% Please refer to
% 1. Hamid Esmaeili and Morteza Kimiaei, "A new adaptive trust-region
%    method for system of nonlinear equations," Appl. Math. Modeling 38,
%    3003-3015 (2014)
% 2. "Trust-Region Methods," Applied & Computational Mathematics Emphasis (ACME), Brigham Young University
%    https://acme.byu.edu/0000017a-1bb8-db63-a97e-7bfa0bd80000/vol2lab19trustregion-pdf
%
% Notations in this function are consistent with the first reference paper.
%
% Developed by Yi-Hao Chen from Frank Wise's group at Cornell University @6/5/2024

% Find the dimension of the computation
num_x = size(x0,3); % the size of the "x" spatial dimension
num_y = size(x0,4); % the size of the "y" spatial dimension
single_mode_yes = (num_x == 1 & num_y == 1); % whether the simulation is single-mode or not

if ~single_mode_yes
    MATLAB_version = version('-release'); MATLAB_version = str2double(MATLAB_version(1:4));
    if ~gpu_yes
        if MATLAB_version < 2021
            error('myTrustRegion:pageXXXXXXError',...
                  ['pagetranspose() and pagemtimes() exist only after R2020b.\n',...
                   'For multimode cases, you must use GPU if your MATLAB version is before 2021.']);
        elseif MATLAB_version == 2021
            error('myTrustRegion:pagemldivideError',...
                  ['pagemldivide() exists only after R2022a.\n',...
                   'For multimode cases, you must use GPU if your MATLAB version is before 2022.']);
        end
    end
end

% Trust-region method can fail with a bad initial guess.
% Because the solver includes parameters dependent on trust-region internal
% iterations, it's useful to re-run the solver with an updated initial guess
% obtained from the previous solver. This resets some parameters and can
% lead to the solution.
max_iterations = 5;
for i = 1:max_iterations
    [xk,success] = run_myTrustRegion(F,J,x0,N,e,M,x_total,gpu_yes,single_mode_yes);

    if success
        break;
    else
        x0 = xk; % prepare to run the next Trust-region method with the updated initial guess, from the previous solver
    end
end

end

%% Main function for the trust-region method
function [xk,success] = run_myTrustRegion(F,J,x0,N,e,M,x_total,gpu_yes,single_mode_yes)

x_total = repmat(x_total,1,1,1,1,1,1,1,M);
max_x_total = max(x_total(:));

normF2 = @(x) sum(abs( F(x) ).^2,1); % objective function

eta0 = 0.2; % parameter for finding the trust-region radius
mu = 0.2; % threshold of optimization
max_iterations = 10; % maximum iterations for the trust-region solver

% Find the dimension of the computation
num_x = size(x0,3); % the size of the "x" spatial dimension
num_y = size(x0,4); % the size of the "y" spatial dimension

%%
m = zeros(max_iterations,1); % parameter as a non-monotone modification of Armijo-type optimization
xk = x0;
xk1 = x0;
etak = [eta0,0,0]; % container for computing eta's recursion relation: eta(k+2)=( eta(k+1)+eta(k) )/2
Flk = sqrt(normF2(x0));
Rk = Flk;
dk = Rk;
normF_all = zeros(max_iterations,1,num_x,num_y,1,1,1,M);
normF_all(1,:,:,:,:,:,:,:) = sqrt(normF2(x0));
success = false;
for k = 1:max_iterations
    % Stopping criterions include
    % 1. All of the norms of the output vector (for each parallelization) is smaller than "e"
    % 2. The step size "dk" of the trust-region method becomes too small (< max_x_total/1e4)
    %
    %    In general, there might be no solution to the coupled equation, so
    %    having this step-size check is important to avoid this solver from
    %    running forever.
    if any(normF_all(k,:,:,:,:,:,:,:) >= e, 'all') && any(dk > max_x_total/1e4, 'all')
        ratio_k = zeros(1,1,num_x,num_y,1,1,1,M);
        while any(ratio_k(:) < mu) && any(sqrt(normF2(xk1)) >= e, 'all') && any(dk > max_x_total/1e4, 'all')
            [xk1,mk] = dogleg(F,J,xk,dk,num_x,num_y,M,x_total,gpu_yes,single_mode_yes);
            ratio_k = (normF_all(k,:,:,:,:,:,:,:).^2-normF2(xk1))./(normF_all(k,:,:,:,:,:,:,:).^2-mk);
            if any(ratio_k(:) < mu)
                dk = 0.5*dk;
            end
        end
        
        xk = xk1;
        normF_all(k+1,:,:,:,:,:,:,:) = sqrt(normF2(xk1));
        if k == 1
            etak(2) = etak(1)/2;
        elseif k == 2
            etak(3) = (etak(1)+etak(2))/2;
        else % k >= 3
            etak(1) = etak(2);
            etak(2) = etak(3);
            etak(3) = (etak(1)+etak(2))/2;
        end
        m(k+1) = min(m(k)+1,N);
        
        Flk = max(normF_all(k+1-m(k):k+1,:,:,:,:,:,:,:),[],1);
        Rk = etak(3)*Flk + (1-etak(3))*normF_all(k+1,:,:,:,:,:,:,:);
        dk = max(cat(1,Rk,dk),[],1);
    else % stopping criterion is satisfied; time to stop
        success = true;
        break;
    end
end

end

%% Helper function for the trust-region method
function [x,mk] = dogleg(F,J,x0,R,num_x,num_y,M,x_total,gpu_yes,single_mode_yes)
%DOGLEG It's a numerical method for solving the sub-problem of a
%trust-region optimization.
%
%   F: coupled equation; it outupts a (n,1,num_x,num_y,1,1,1,M) array
%   J: Jacobian matrix; it outupts a (n,m,num_x,num_y,1,1,1,M) array
%   x0: initial guess for m's variables; (m,1,num_x,num_y,1,1,1,M)
%   R: trust-region radius; (1,1,num_x,num_y,1,1,1,M)

Jx0 = J(x0); % Jacobian matrix at x=x0

% x0 can have less dimensions than what is needed, so extension of
% dimensions is required for correct computation below, especially when
% using idx1, idx2, and idx3.
if size(x0,8) ~= M % the number of parallelizations
    x0 = repmat(x0,1,1,1,1,1,1,1,M);
end

if gpu_yes
    g = pagefun(@mtimes,pagefun(@transpose,Jx0),F(x0)); % gradient array for the model function in the sub-problem of the trust-region method
    H = pagefun(@mtimes,pagefun(@transpose,Jx0),Jx0); % Hessian matrix for the model function in the sub-problem of the trust-region method
else
    if single_mode_yes
        g = (Jx0.')*F(x0);
        H = (Jx0.')*Jx0;
    else
        g = pagemtimes(pagetranspose(Jx0),F(x0)); % gradient array for the model function in the sub-problem of the trust-region method
        H = pagemtimes(pagetranspose(Jx0),Jx0); % Hessian matrix for the model function in the sub-problem of the trust-region method
    end
end

dx = zeros(size(x0,1),1,num_x,num_y,1,1,1,M);
mk = zeros(1,         1,num_x,num_y,1,1,1,M); % model function value at x=x+dx

% Core of the dogleg method
% pB: uncontrained minimizer of the model function
% pU: steepest descent for the model function
%
% Dogleg path follows
%   path(tau) = tau*pU,             0<=tau<=1
%               pU+(tau-1)*(pB-pU), 1<=tau<=2
%
% When computing pB, under situations of a weakly-pumped system, stimulated
% terms are all close to zero, leading to a Hessian matrix that is close to
% singular. This will trigger MATLAB's warning:
%   Warning: Matrix is close to singular or badly scaled. Results may be inaccurate.
% To resolve this (not really an) issue, it's better to scale the Hessian
% first when computing pB and then scale/recover it back.
% This is achievable because, in general, for a invertible matrix H,
%   (H*A)^(-1) = A^(-1)*H^(-1),
% where A = [c1  0  0  0 ...
%             0 c2  0  0 ...
%             0  0 c3  0 ...
%             0  0  0 c4 ...
%             .  .  .  . ...] is the identity matrix with each of its
% entries replaced with the corresponding scaling factor for each column.
scaling_H = max(H,[],1);
if single_mode_yes % no page-wise computation; I added this check because MATLAB adds pagemldivide only after R2022a
    pB = -mldivide(H./scaling_H,g)./scaling_H.';
else
    if gpu_yes
        transpose_scaling_H = pagefun(@transpose,scaling_H);
        pB = -pagefun(@mldivide,H./scaling_H,g)./transpose_scaling_H;
    else
        transpose_scaling_H = pagetranspose(scaling_H);
        pB = -pagemldivide(H./scaling_H,g)./transpose_scaling_H;
    end
    
    % Multimode computation code for MATLAB without GPU and MATLAB's version before 2022.
    % I disabled this for maximum performance. User can enable it if
    % necessary and under old MATLAB. Just remember to disable the previous
    % few computation lines.
    %pB = zeros(size(x0,1),1,num_x,num_y,1,1,1,M);
    %for i_x = 1:num_x
    %    for i_y = 1:num_y
    %        for i_M = 1:M
    %            pB(:,:,i_x,i_y,1,1,1,i_M) = -mldivide(H(:,:,i_x,i_y,1,1,1,i_M)./scaling_H(:,:,i_x,i_y,1,1,1,i_M),g(:,:,i_x,i_y,1,1,1,i_M))./scaling_H(:,:,i_x,i_y,1,1,1,i_M).';
    %        end
    %    end
    %end
end
norm_pB = sqrt(sum(abs(pB).^2,1));
if gpu_yes
    pU = -sum(g.^2,1)./sum(g.*pagefun(@mtimes,H,g),1).*g;
else
    if single_mode_yes
        pU = -sum(g.^2,1)./sum(g.*(H*g),1).*g;
    else
        pU = -sum(g.^2,1)./sum(g.*pagemtimes(H,g),1).*g;
    end
end
norm_pU = sqrt(sum(abs(pU).^2,1));
% if norm_pB <= R; in principle, this corresponds to tau>2 if it extends to R. It's stopped at tau=2, which is dx=pB.
idx1 = shiftdim(norm_pB <= R,2); % remain as a boolean array for more "if-else" computation below (idx2 and idx3)
if any(idx1,'all')
    dx(:,:,idx1) = pB(:,:,idx1);
    mk(:,:,idx1) = 0;
end
% elseif norm_pU >= R; 0<=tau<=1
idx2 = (~idx1 & shiftdim(norm_pU >= R,2)); % 0<=tau<=1; remain as a boolean array for more "if-else" computation below (idx3)
if any(idx2,'all')
    dx(:,:,idx2) = R(:,:,idx2).*pU(:,:,idx2)./norm_pU(:,:,idx2);
    if gpu_yes
        tmp = sum(abs(F(x0)+pagefun(@mtimes,Jx0,dx)).^2,1);
    else
        if single_mode_yes
            tmp = sum(abs(F(x0)+(Jx0*dx)).^2,1);
        else
            tmp = sum(abs(F(x0)+pagemtimes(Jx0,dx)).^2,1);
        end
    end
    mk(:,:,idx2) = tmp(:,:,idx2);
end
% else; 1<=tau<=2
idx3 = (~idx1 & ~idx2); % 1<=tau<=2
if any(idx3,'all')
    a = sum((pB(:,:,idx3)-pU(:,:,idx3)).^2,1);
    b = 2*sum(pU(:,:,idx3).*(pB(:,:,idx3)-pU(:,:,idx3)),1);
    c = sum(pU(:,:,idx3).^2,1)-R(:,:,idx3).^2;
    % (tau-1) is the "r"-solution of a*r^2+b*r+c=0
    tau = (-b+sqrt(b.^2-4*a.*c))/2./a + 1;
    
    dx(:,:,idx3) = pU(:,:,idx3) + (tau-1).*(pB(:,:,idx3)-pU(:,:,idx3));
    if gpu_yes
        tmp = sum(abs(F(x0)+pagefun(@mtimes,Jx0,dx)).^2,1);
    else
        if single_mode_yes
            tmp = sum(abs(F(x0)+(Jx0*dx)).^2,1);
        else
            tmp = sum(abs(F(x0)+pagemtimes(Jx0,dx)).^2,1);
        end
    end
    mk(:,:,idx3) = tmp(:,:,idx3);
end

x = x0 + dx;
x(x<0) = 0; % population can't be smaller than zero

% Total population of the m's population shouldn't be larger than x_total.
% If this happens, I normalize them back to x_total.
sum_x = sum(x,1);
idx = (sum_x > x_total);
x(idx) = x(idx).*x_total(idx)./sum_x(idx);

end