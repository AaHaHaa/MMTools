function L = calc_effective_rmc_matrix(fiber,sim,Nt,dt,...
                                       rmc_matrices,...
                                       varargin)
%CALC_EFFECTIVE_RMC_MATRIX It calculates the effective matrix of the entire
%linear operations, including dispersion and random mode coupling

%% Default optional input arguments
% Accept only 4 optional inputs at most
numvarargs = length(varargin);
if numvarargs > 1
    error('calc_effective_rmc_matrix:TooManyInputs', ...
        'It takes only at most 1 optional input');
end

% Set defaults for optional inputs
fields = [];
optargs = {fields};

% Now put these defaults into the valuesToUse cell array, 
% and overwrite the ones specified in varargin.
optargs(1:numvarargs) = varargin;

% Place optional args in memorable variable names
fields = optargs{:};

%% Compute the dispersion w.r.t. frequency
f = (-floor(Nt/2):floor((Nt-1)/2))'/Nt/dt + sim.f0; % THz
lambda = 299792458*1e-12./f; % m
idx0 = find(lambda>=sim.rmc.lambda0,1,'last');
if idx0 > floor(Nt/2)
    idx0 = idx0 - floor(Nt/2);
else
    idx0 = idx0 + ceil(Nt/2);
end

num_modes = size(fiber.betas,2)*(sim.scalar+(~sim.scalar)*2);
fiber = betas_expansion_including_polarization_modes(sim,fiber,num_modes);

D_op = calc_D_op(fiber,sim,Nt,dt,2*pi*ifftshift(f-sim.f0,1),fields);
D_op = 1i*imag(D_op(idx0,:)); % take the dispersion only at one wavelength

%% Compute the effective linear operation of the entire fiber
D_op_expand = diag(D_op);

% Entire linear operator
L_op = D_op_expand + rmc_matrices;

L = eye(size(D_op_expand,1));
for i = 1:size(rmc_matrices,3)
    L = expm(L_op(:,:,i).*sim.dz)*L;
end

end