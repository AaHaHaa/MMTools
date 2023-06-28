function y = mySpline( x0,y0,x,smooth_order,cuda_dir_path )
%MYSPLINE It computes spline interpolation up to a maximum of the 6th
%order with cuda.
%   Since MATLAB supports only cubic spline, which considers only up to 2nd
%   order, smooth higher-order derivatives are sometimes desirable.
%
% Input:
%   x0: the input x; an equally-spaced monotonically-increasing array
%   y0: the input y; an array
%   x: the desired interpolated x points; an array
%   smooth_order: the derivative order that users want the interpolated y
%                 to be smooth up to.
%                 If it's 1, then this is simply linear interpolation.
%                 If it's 2, then this is simply cubic spline interpolation.
%
% Output:
%   y: the desired interpolated y at user's x

diff_x0 = diff(x0);
if abs(max(diff_x0)/min(diff_x0)) > 1.01
    error('mySpline:x0Error',...
          'x0 needs to be an equally-spaced monotonically-varying array.');
end

dx0 = mean(diff(x0));
num_x0 = length(x0);
num_x = length(x);
num_lines = size(y0,2);

% This code calculates each derivative with central difference.
% For the 1st and 2nd order derivatives, two edge points, one on each edge, are discarded.
% For the 3rd and 4th order ones, four edge points, two on each edge, are discarded.
% For the 5rd and 6th order ones, six edge points, three on each edge, are discarded.
switch smooth_order
    case {1,2}
        x0_interp = x0(2:end-1);
        num_x0_interp = num_x0 - 2;
    case {3,4}
        x0_interp = x0(3:end-2);
        num_x0_interp = num_x0 - 4;
    case {5,6}
        x0_interp = x0(4:end-3);
        num_x0_interp = num_x0 - 6;
end

% Calculate the derivatives based on central difference
central_diff_coeff = {[   0,    0, -1/2,   0, 1/2,   0,   0],...
                      [   0,    0,    1,  -2,   1,   0,   0],...
                      [   0, -1/2,    1,   0,  -1, 1/2,   0],...
                      [   0,    1,   -4,   6,  -4,   1,   0],...
                      [-1/2,    2, -5/2,   0, 5/2,  -2, 1/2],...
                      [   1,   -6,   15, -20,  15,  -6,   1]};

y0_full = cell(smooth_order+1, 1);
y0_full{1} = y0;
for i = 1:smooth_order
    switch i
        case {1,2}
            y0_full{i+1} = central_diff_coeff{i}(3)*y0(1:end-2,:) + ...
                           central_diff_coeff{i}(4)*y0(2:end-1,:) + ...
                           central_diff_coeff{i}(5)*y0(3:end,:);
        case {3,4}
            y0_full{i+1} = central_diff_coeff{i}(2)*y0(1:end-4,:) + ...
                           central_diff_coeff{i}(3)*y0(2:end-3,:) + ...
                           central_diff_coeff{i}(4)*y0(3:end-2,:) + ...
                           central_diff_coeff{i}(5)*y0(4:end-1,:) + ...
                           central_diff_coeff{i}(6)*y0(5:end,:);
        case {5,6}
            y0_full{i+1} = central_diff_coeff{i}(1)*y0(1:end-6,:) + ...
                           central_diff_coeff{i}(2)*y0(2:end-5,:) + ...
                           central_diff_coeff{i}(3)*y0(3:end-4,:) + ...
                           central_diff_coeff{i}(4)*y0(4:end-3,:) + ...
                           central_diff_coeff{i}(5)*y0(5:end-2,:) + ...
                           central_diff_coeff{i}(6)*y0(6:end-1,:) + ...
                           central_diff_coeff{i}(7)*y0(7:end,:);
    end
    y0_full{i+1} = y0_full{i+1}/dx0^double(i);
end

% Make sure all derivatives have the same length of arrays.
switch smooth_order
    case {1,2}
        y0_full{1} = y0_full{1}(2:end-1,:);
    case {3,4}
        y0_full{1} = y0_full{1}(3:end-2,:);
        y0_full{2} = y0_full{2}(2:end-1,:);
        y0_full{3} = y0_full{3}(2:end-1,:);
    case {5,6}
        y0_full{1} = y0_full{1}(4:end-3,:);
        y0_full{2} = y0_full{2}(3:end-2,:);
        y0_full{3} = y0_full{3}(3:end-2,:);
        y0_full{4} = y0_full{4}(2:end-1,:);
        y0_full{5} = y0_full{5}(2:end-1,:);
end

addpath([cuda_dir_path,'/../GMMNLSE algorithm']);
cuda_mySpline = setup_kernel_simple('mySpline',cuda_dir_path,num_x);
y = complex(zeros(num_x,num_lines,'gpuArray'));
for lidx = 1:num_lines
    y_lidx = complex(zeros(num_x0_interp,smooth_order+1,'gpuArray'));
    for i = 1:smooth_order+1
        y_lidx(:,i) = y0_full{i}(:,lidx);
    end
    y(:,lidx) = feval(cuda_mySpline,...
                                    x0_interp,complex(y_lidx),uint32(num_x0_interp),...
                                    x,y(:,lidx),uint32(num_x),...
                                    uint32(smooth_order));
end
y = gather(y);

end