function stdQ = rmc_std_Q( fiber,sim,num_modes,mode_profiles )
%RMC_STD_Q It calculates the coupling strength between or inside mode 
%groups, which will be multiplied with the corresponding random matrix with
%a standard deviation 1.
%
%   num_modes: the number of modes including polarization modes
%
%   stdQ: a matrix of the size (num_modes, num_modes); coupling strength

if sim.scalar
    num_spatial_modes = num_modes;
else
    num_spatial_modes = num_modes/2;
end

k0 = 2*pi/sim.rmc.lambda0;
neff = mean(fiber.betas(1,:))/k0;
if isfield(fiber,'MM_folder')
    % Find the variance of integral2( randn*F1*F2 ):
    % Its mean = 0
    % Its variance = integral2( (F1*F2)^2 )
    var_randn_F2 = permute(...
                           sum((             mode_profiles.mode_profiles.*...
                                permute(conj(mode_profiles.mode_profiles),[1,2,4,3])...
                               ).^2,...
                              [1,2])*mode_profiles.dx^4,...
                          [3,4,1,2]); % make it (num_modes,num_modes)
    % Generate the matrix of coupling strength stdQ
    stdQ = k0/2/neff*abs(sim.rmc.varn)*sqrt(var_randn_F2);
end

if ~sim.scalar
    % Extend the matrix to include the polarization modes
    tmp = mat2cell(ones(num_modes,num_modes),2*ones(1,num_spatial_modes),2*ones(1,num_spatial_modes));
    for i = 1:num_spatial_modes^2
        tmp{i} = tmp{i}*i;
    end
    stdQ_idx = cell2mat(tmp);
    stdQ = stdQ(stdQ_idx);

    % Include the polarization-mode coupling to the right elements of the coupling-strength matrix
    tmp = mat2cell(sim.rmc.stdQ_polarizedmode*randn(2,2,num_spatial_modes),2,2,ones(1,num_spatial_modes));
    stdQ_polarizedmode = blkdiag(tmp{:});
    stdQ(stdQ_polarizedmode~=0) = stdQ_polarizedmode(stdQ_polarizedmode~=0);
end

% The diagonal part should be zero because this is the cross-coupling matrix
% The diagonal part is actually the dispersion term and should be zero.
diag_idx = (1:num_modes) + (0:num_modes-1)*num_modes;
stdQ(diag_idx) = 0;

% The coupling between polarization modes should be determined only by stdQ_polarizedmode,
% so I kill the coupling between different spatial modes with different polarizations here.
if ~sim.scalar
    p1_idx = zeros(num_modes,1); p1_idx(1:2:end) = 1;
    p2_idx = zeros(num_modes,1); p2_idx(2:2:end) = 1;
    coupling_idx = p1_idx*p1_idx' + p2_idx*p2_idx'; % spatial-mode coupling
    coupling_idx(diag_idx) = 0;
    for midx = 1:2:num_modes % add polarization-mode coupling
        coupling_idx(midx,midx+1) = 1;
        coupling_idx(midx+1,midx) = 1;
    end
    stdQ = stdQ.*coupling_idx;
end

stdQ = (stdQ + stdQ')/2;

end