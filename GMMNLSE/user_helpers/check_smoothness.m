function roughness = check_smoothness( x,smooth_n )
%CHECK_SMOOTHNESS
%   It calculates the difference between the smooth data and the original
%   one.
%   If there are oscillatory structures, the root-sum-of-square will be
%   large; otherwise, it'll be small.
%   From some comparisons, I found out that if it's smooth (converged
%   pulses), roughness ~ 5e-4. If it starts to become structured, it's
%   about 0.001. Really structured pulse can go up to 0.005.
%
%   x: (N,...); the pulse shape, that is, abs(E) or abs(E)^2

sx = size(x);

if ~exist('smooth_n','var') || smooth_n < 1
    width = calc_RMS(1:size(x,1),x)*(2*sqrt(log(2)))*sqrt(2); % t0 -> tfwhm
    smooth_n = ceil(width/20);
else
    smooth_n = ceil(smooth_n);
end

roughness = zeros([1,sx(2:end)]);
for i = 1:prod(sx(2:end))
    roughness(i) = rssq( smooth(x(:,i),smooth_n(i))-x(:,i) )./trapz(x(:,i));  
end

end

