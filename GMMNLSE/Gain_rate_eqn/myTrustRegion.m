function xk = myTrustRegion(F,J,x0,N,e,M,x_total)
%MYTRUSTREGION Trust-region method for solving the (num_x*num_y*M)'s 
%optimization problems in parallel. This parallelization is crucial in
%multimode rate-equation gain computations.
%
%   F: coupled equation; it outupts a (n,1,num_x,num_y,1,1,1,M) array
%   J: Jacobian matrix; it outupts a (n,m,num_x,num_y,1,1,1,M) array
%   x0: initial guess for m's variables; (m,1,num_x,num_y,1,1,1,M) array
%   N: parameter for non-monotone function evaluations during optimization; a scalar
%   e: threshold value of stopping optimization processes; a scalar
%   M: the number of parallelizations in MPA
%   x_total: the total population; (1,1,num_x,num_y)
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

x_total = repmat(x_total,1,1,1,1,1,1,1,M);

normF2 = @(x) sum(abs( F(x) ).^2,1); % objective function

eta0 = 0.2; % parameter for finding the trust-region radius
mu = 0.2; % threshold of optimization
max_iterations = 1000; % maximum iterations for the trust-region solver

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
for k = 1:max_iterations
    % Stopping criterions include
    % 1. All of the norms of the output vector (for each parallelization) is smaller than "e"
    % 2. The step size "dk" of the trust-region method becomes too small
    %    In general, there might be no solution to the coupled equation, so
    %    having this step-size check is important to avoid this solver to
    %    run forever.
    if any(normF_all(k,:,:,:,:,:,:,:) >= e, 'all') && any(dk > max(x0(:))/1e4, 'all')
        ratio_k = zeros(1,1,num_x,num_y);
        while any(ratio_k(:) < mu) && any(sqrt(normF2(xk1)) >= e, 'all') && any(dk > max(x0(:))/1e4, 'all')
            [xk1,mk] = dogleg(F,J,xk,dk,num_x,num_y,M,x_total);
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
        break;
    end
end

end

%% Helper function
function [x,mk] = dogleg(F,J,x0,R,num_x,num_y,M,x_total)
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

g = pagemtimes(pagetranspose(Jx0),F(x0)); % gradient array for the model function in the sub-problem of the trust-region method
H = pagemtimes(pagetranspose(Jx0),Jx0); % Hessian matrix for the model function in the sub-problem of the trust-region method

dx = zeros(size(x0,1),1,num_x,num_y,1,1,1,M);
mk = zeros(1,         1,num_x,num_y,1,1,1,M); % model function value at x=x+dx

% Core of the dogleg method
% pB: uncontrained minimizer of the model function
% pU: steepest descent for the model function
%
% Dogleg path follows
%   path(tau) = tau*pU,             0<=tau<=1
%               pU+(tau-1)*(pB-pU), 1<=tau<=2
if num_x == 1 && num_y == 1 % no page-wise computation; I added this check because MATLAB adds pagemldivide only after R2022a
    pB = -mldivide(H,g);
else
    MATLAB_version = version('-release'); MATLAB_version = str2double(MATLAB_version(1:4));
    if MATLAB_version >= 2022
        try
            pB = -pagemldivide(H,g);
        catch
            error('myTrustRegion:pagemldivideError',...
                  'You must set "sim.gpu_yes=true" for the multi-level gain medium.');
        end
    else
        error('myTrustRegion:pagemldivideError',...
              ['pagemldivide() exists only after R2022a.\n',...
               'For multimode cases, you must use GPU if your MATLAB version is before 2022.']);
    end
    
    % Code for MATLAB without GPU and MATLAB's version before 2022.
    %pB = zeros(size(x0,1),1,num_x,num_y,1,1,1,M);
    %for i_x = 1:num_x
    %    for i_y = 1:num_y
    %        for i_M = 1:M
    %            pB(:,:,i_x,i_y,1,1,1,i_M) = -mldivide(H(:,:,i_x,i_y,1,1,1,i_M),g(:,:,i_x,i_y,1,1,1,i_M));
    %        end
    %    end
    %end
end
norm_pB = sqrt(sum(abs(pB).^2,1));
pU = -sum(g.^2,1)./sum(g.*pagemtimes(H,g),1).*g;
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
    tmp = sum(abs(F(x0)+pagemtimes(Jx0,dx)).^2,1);
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
    tmp = sum(abs(F(x0)+pagemtimes(Jx0,dx)).^2,1);
    mk(:,:,idx3) = tmp(:,:,idx3);
end

x = x0 + dx;
x(x<0) = 0; % population can't be smaller than zero

% Total population of the m's population should be larger than x_total.
% If this happens, I normalize them back to x_total.
sum_x = sum(x,1);
idx = (sum_x > x_total);
x(idx) = x(idx).*x_total(idx)./sum_x(idx);

end