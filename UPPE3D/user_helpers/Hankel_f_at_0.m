function A0 = Hankel_f_at_0(A,l0)
%HANKEL_F_F0 It approximates the field at r=0
%
% In Hankel transform, the field is not uniformly sampled. In particular,
% its smallest sampling radius is not zero. This code extrapolates the the
% field at r=0 with existing field A by assuming a parabola with zero 
% derivative at r=0.

A0 = l0*(A(:,1,:,:) - A(:,2,:,:)) + A(:,2,:,:);

end

