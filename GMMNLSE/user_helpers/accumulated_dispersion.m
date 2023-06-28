function dispersion = accumulated_dispersion( betas,L )
%ACCUMULATED_DISPERSION Calculates the total accumulated dispersion
%
%   betas: a (m,num_modes,num_fiber) matrix, m:         the order of betas
%                                            num_modes: the number of modes
%                                            num_fiber: the number of fibers
%          Only the beta2 of each segment is extracted with "betas(3,:)".
%   L: a (1,num_fiber) array; the length of each fiber segment
%
% Use this code as:
%   dispersion = accumulated_dispersion(cat(3,fiber_all.betas),[fiber_all.L0]); % if fiber_all(i) represents the ith fiber, where fiber_all is a structure array

dispersion = sum(betas(3,:,:).*permute(L,[1 3 2]),3); % sum( beta2_i * L_i )

end

