% This code tests my trust-region method tuned for rate equations.
% Compared to MATLAB's "fsolve()" which also rely on trust region in some
% cases, my trust-region method can solve for many optimization problems
% simulaneously, which is useful for multimode computations.

% Target objective coupled equations
r = @(x) [10*(x(2,:,:,:)-x(1,:,:,:).^3).^2+(1-x(1,:,:,:)).^2;...
          10*(x(1,:,:,:)-x(2,:,:,:).^3).^2+(1-x(2,:,:,:)).^2];
% and the Jacobian matrix
J = @(x) [-60*(x(2,:,:,:).*x(1,:,:,:).^2-x(1,:,:,:).^5)-2*(1-x(1,:,:,:)), 20*(x(2,:,:,:)-x(1,:,:,:).^3);...
           20*(x(1,:,:,:)-x(2,:,:,:).^3),   -60*(x(1,:,:,:).*x(2,:,:,:).^2-x(2,:,:,:).^5)-2*(1-x(2,:,:,:))];
% Initial conditions; size: (2,1,3) for simultaneously three optimizations for each (x1;x2) pair
% 1. multimode test
%x0 = cat(3,[50;10],[10;-100],[-1;-200]);
% 2. signle-mode test
x0 = cat(3,[50;10]);
% Parameter for non-monotone function evaluations in trust-region method
N = 5;
% Threshold value of stopping optimization processes
e = 1e-5;
% Unimportant parameters for this test
M = 1; % the number of parallelizations in MPA pulse-stepping algorithm; set it to 1 for it to run
x_total = 1e5; % the total population
               % Because this trust-region method is tuned to work for population evolutions,
               % which has a fixed number of total population,
               % we need it to be a huge number to make the current test work.
% Run the trust-region optimization
xk = myTrustRegion(r,J,x0,N,e,1,1e5);

% Show results
disp(xk); % the point that reaches optimum
disp(r(xk)); % the objective function at the optimum