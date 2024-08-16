function coeff = MPA_AM_coeff( M )
%MPA_AM_COEFF It calculates the coefficients required for computing with
%MPA based on the Adams-Moulton method
%
% Input:
%   M: the number of parallelization in MPA which corresponds to the order
%      of the Adams-Moulton method
%
% For example:
%   M=1:
%
%   y(n+1) = y(n) + h*(    1/ 2*f(t(n+1),y(n+1)) 
%                       +  1/ 2*f(t(n  ),y(n  )) )
%
%   M=2:
%
%   y(n+2) = y(n) + h*(    5/12*f(t(n+2),y(n+2))
%                       + 14/12*f(t(n+1),y(n+1))
%                       +  5/12*f(t(n  ),y(n  )) )
%
%   M=3:
%
%   y(n+3) = y(n) + h*(    9/24*f(t(n+3),y(n+3))
%                       + 29/24*f(t(n+2),y(n+2))
%                       + 23/24*f(t(n+1),y(n+1))
%                       + 11/24*f(t(n  ),y(n  )) )

coeff = zeros(M+1);
for s = 1:M
    b = Adams_Moulton_coeff(s);
    
    coeff(s+1,1:s+1) = coeff(s,1:s+1) + b;
end
coeff = coeff(2:end,:);

end

%% helper function
function b = Adams_Moulton_coeff( s )
%ADAMS_MOULTON_COEFF It calculates the coefficients required for computing
%with the Adams-Moulton method
%
% Reference:
%   Wikipedia: Linear multistep method
%   https://en.wikipedia.org/wiki/Linear_multistep_method
% 
% y(n+s) -y(n+s-1) = h*(   b(s  )*f(t(n+s  ),y(n+s)  )
%                        + b(s-1)*f(t(n+s-1),y(n+s-1))
%                        + b(s-2)*f(t(n+s-2),y(n+s-2))
%                        ...
%                        + b(0  )*f(t(n    ),y(n    )) )
%
% For example:
%   y(n  ) = y(n-1) + h*          f(t(n),y(n))
%
%   y(n+1) = y(n  ) + h*(    1/ 2*f(t(n+1),y(n+1)) 
%                         +  1/ 2*f(t(n  ),y(n  )) )
%
%   y(n+2) = y(n+1) + h*(    5/12*f(t(n+2),y(n+2))
%                         +  8/12*f(t(n+1),y(n+1))
%                         -  1/12*f(t(n  ),y(n  )) )
%
%   y(n+3) = y(n+2) + h*(    9/24*f(t(n+3),y(n+3))
%                         + 19/24*f(t(n+2),y(n+2))
%                         -  5/24*f(t(n+1),y(n+1))
%                         +  1/24*f(t(n  ),y(n  )) )

x = linspace(0,1,1000)';

b = zeros(1,s+1);
for j = 0:s
    i = 0:s; i(i==j) = []; % i ~= j
    integrand = prod(x+i-1,2);
    b(s-j+1) = (-1)^j/factorial(j)/factorial(s-j)*trapz(x,integrand);
end

end

