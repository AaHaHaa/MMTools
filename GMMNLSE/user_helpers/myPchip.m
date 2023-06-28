function y = myPchip( x0,y0,x,smooth_order,cuda_dir_path )
%MYPCHIP It computes pchip interpolation up to a maximum of the 6th
%order with cuda.
%   Since MATLAB supports only pchip, which considers only up to 2nd
%   order, smooth higher-order derivatives are sometimes desirable.
%
% Input:
%   x0: the input x; an monotonically-increasing array
%   y0: the input y; an array
%   x: the desired interpolated x points; an array
%   smooth_order: the derivative order that users want the interpolated y
%                 to be smooth up to.
%                 If it's 1, then this is simply linear interpolation.
%                 If it's 2, then this is simply cubic spline interpolation.
%
% Output:
%   y: the desired interpolated y at user's x

num_x0 = length(x0);
num_x = length(x);
num_lines = size(y0,2);

%% Interpolation
% This code calculates each derivative modified with harmonic mean.
% See F. N. Fritsch and J. Butland, 
% "A method for constructing local monotone piecewise cubic interpolants," 
% SIAM J. Sci. Comput., 5(2), 300-304 (1984). DOI:10.1137/0905021.
% Unlike spline interpolation, pchip has a "shape-preserving" feature
% obtained from the following derivative computation.
% The harmonic mean G:
% G(S1,S2) = 2*S1*S2/(S1+S2), if S1*S2>0
%            0,               otherwise
y0_full = cell(smooth_order+1, 1);
y0_full{1} = y0;
for i = 1:smooth_order
    S = (y0_full{i}(2:end,:) - y0_full{i}(1:end-1,:))./(x0(2:end) - x0(1:end-1));
    S1 = S(1:end-1,:);
    S2 = S(2:end,:);
    reG = 2*real(S1).*real(S2)./real(S1+S2); reG(real(S1).*real(S2)<=0) = 0;
    imG = 2*imag(S1).*imag(S2)./imag(S1+S2); imG(imag(S1).*imag(S2)<=0) = 0;
    y0_full{i+1} = [S1(1,:); reG+1i*imG; S1(end,:)];
end

addpath([cuda_dir_path,'/../GMMNLSE algorithm']); % to use "setup_kernel_simple()"
cuda_mySpline = setup_kernel_simple('mySpline',cuda_dir_path,num_x);
y = complex(zeros(num_x,num_lines,'gpuArray'));
for lidx = 1:num_lines
    y_lidx = complex(zeros(num_x0,smooth_order+1,'gpuArray'));
    for i = 1:smooth_order+1
        y_lidx(:,i) = y0_full{i}(:,lidx);
    end
    y(:,lidx) = feval(cuda_mySpline,...
                                    x0,complex(y_lidx),uint32(num_x0),...
                                    x,y(:,lidx),uint32(num_x),...
                                    uint32(smooth_order));
end
y = gather(y);

%% Extrapolation
extrap_idx = x<x0(1) | x>x0(end);
y(extrap_idx,:) = interp1(x0,y0,x(extrap_idx),'pchip','extrap');

end