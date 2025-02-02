function [r,kr,...
          dr,dkr,...
          l0,exp_prefactor,r2_prefactor,...
          ifftQ] = Hankel_info(Nr,r_max,varargin)
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
%
% Output arguments:
%   r: sampling radius (m); (1,Nr)
%   kr: sampling k-vector radius (2*pi/m); (1,Nr)
%   dr: samppling radius spacing (m); (1,Nr);
%       dr includes the effect from r=0, so it is computed from
%           dr = diff([0,r]);
%   dkr: samppling k-vector radius (2*pi/m); (1,Nr);
%       dr includes the effect from kr=0, so it is computed from
%           dkr = diff([0,kr]);
%   exp_prefactor: exponential prefactor used in computing R in FHATHA prepared for cross correlation
%   r2_prefactor: (=r^2/2) prefactor used in computing the Hankel-transformed signal at kr=0
%   ifftQ: the inverse-Fourier-transformed Q in FHATHA prepared for cross correlation

n = 0:(Nr-1);

% Finding alpha requires MATLAB's symbolic toolbox
syms a;
alpha = double(vpasolve(exp(-a*(Nr-1))==(1-exp(-a)),a,0.01)); % I pick 0.01 as initial guess

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

 r = [0, r_max*zeta]; % sampling radius
kr = [0,kr_max*zeta]; % sampling k-radius

 dr = diff(r);
dkr = diff(kr);

%% Precompute Q = J1(...)
n = 0:(2*Nr-1);
Q = besselj(1,kr_max*r_max*zeta0*exp(alpha*(n+1-Nr)));
ifftQ = ifft(Q,[],2);

%% R's prefactor
n = 0:(Nr-1);
exp_prefactor = exp(alpha*(n+1-Nr));

l0 = exp(alpha)*(2+exp(alpha))/(1+exp(alpha))^2/(1-exp(-2*alpha));

%% prefactor in A_H(k=0)'s computation
r2_prefactor = exp_prefactor.^2/2;

end