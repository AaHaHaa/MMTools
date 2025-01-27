function [r,kr,...
          l0,exp_prefactor,...
          Q] = Hankel_info(Nr,r_max,varargin)
%HANKEL_INFO This is the precomputation for the FHATHA, fast Hankel
%transform of high accuracy
%   This numerical scheme for Hankel transform follows
%
%   Magni et al., "High-accuracy fast Hankel transform for optical beam
%   propagation," J. Opt. Soc. Am. A 9(11), 2031-2033 (1992)
%
% In 2D-UPPE code, the radial dimension is in the second dimension, so the
% parameters generated in this code are all row vectors [size: (1,Nr)].
%
% Input arguments:
%   Nr: the number of (radial) sampling points
%   r_max: half the size of the spatial window = maximum radius
%
% Optional input argument:
%   k_max: the maximum k-radius;
%          If it's not supplied, it's automatically determined by the largest radial sampling spacing in this code.

n = 0:(Nr-1);

% Finding alpha requires MATLAB's symbolic toolbox
syms a;
alpha = double(vpasolve(exp(-a*(Nr-1))==1-exp(-a),a,0.01)); % I pick 0.01 as initial guess

%% Function-evaluation points
zeta0 = (1+exp(alpha))/2*exp(-alpha*Nr); % the parameter that makes normalized function-evaluation point at the center of each sampling interval

zeta = zeta0*exp(alpha*n); % function-evaluation points

% I set the kr_max to be determined by the largest sampling spacing such
% that any k-vector smaller than kr_max can all be captured in the entire
% spatial window.
if isempty(varargin)
    kr_max = 2*pi/(1-exp(-alpha))/r_max/2; % this is the maximum "radius" of the k-space
else
    kr_max = varargin{1};
end

r = r_max*zeta; % sampling radius
kr = kr_max*zeta; % sampling k-radius

%% Precompute Q = J1(...)
n = 0:(2*Nr-1);
Q = besselj(1,kr_max*r_max*zeta0*exp(alpha*(n+1-Nr)));

%% R's prefactor
l0 = exp(alpha)*(2+exp(alpha))/(1+exp(alpha))^2/(1-exp(-2*alpha));

n = 0:(Nr-1);
exp_prefactor = exp(alpha*(n+1-Nr));

end