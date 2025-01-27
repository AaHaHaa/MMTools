function A0 = Hankel_f_at_0(A,l0)
%HANKEL_F_F0 It approximates the field at r=0
%
% In Hankel transform, the field is not uniformly sampled. In particular,
% its smallest sampling radius is not zero (at the spatial center).
% This code extrapolates the field at 
%     r'_1 = r_max*xi_1
% with existing field A by assuming a parabola with zero derivative at r=0.
% xi_n is the (non-uniform) sampling points in FHATHA.
% In principle, we can compute the field at r=0. However, FHATHA I
% implement doesn't output the "alpha" value (in principle, it can but I
% don't want to complicate the operation). Since my Hankel_info() outputs
% l0, which can be used to compute A(r'_1), this code gives the field at
% r'_1 instead.
%
% Just for the record, below is the field at r=0, A(r=0):
%
%    A0 = (A(:,1,:,:)*exp(2*alpha) - A(:,2,:,:)) / (exp(2*alpha)-1);
%

A0 = l0*(A(:,1,:,:) - A(:,2,:,:)) + A(:,2,:,:);

end

