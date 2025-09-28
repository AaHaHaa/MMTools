function varargout = solve_gain_rate_eqn(direction,...
                                         sim,gain_rate_eqn,...
                                         N,...
                                         signal_fields,signal_fields_backward,...
                                         Power_pump_forward,Power_pump_backward,...
                                         Power_ASE_forward,Power_ASE_backward,...
                                         Omega,dt,...
                                         first_backward_before_iterations)
%SOLVE_GAIN_RATE_EQN Solves the power of pump, ASE, and signal at
%z+-dz, where +- depends on the propagation direction, forward or
%backward.
% It first solves Ni's, the ion density of each energy level, and uses this
% to calculate the power.
%
% computational dimension: (Nx,Nx,num_spatial_modes,num_spatial_modes,Nt,M,num_polarization,num_cross_sections/num_levels), M: parallelization in MPA

% Load required helper functions
func = solve_gain_rate_eqn_helpers();

Nt = length(Omega); % the number of time/frequency points
time_window = Nt*(dt*1e-12); % unit: s
df = 1/(Nt*dt); % unit: THz
hbar = 6.62607015e-34/(2*pi);
E_photon = (hbar*1e12)*permute(Omega + 2*pi*sim.f0, [2 3 4 5 1]); % size: (1,1,1,1,Nt), unit: J

M = size(signal_fields,3); % the number of parallelization in the MPA-stepping method

if sim.scalar
    num_spatial_modes = size(signal_fields,2);
else % polarized fields
    num_spatial_modes = size(signal_fields,2)/2;
end

E_forward  = signal_fields         *sqrt(time_window/gain_rate_eqn.t_rep);
E_backward = signal_fields_backward*sqrt(time_window/gain_rate_eqn.t_rep); % change it to the correct physical unit of W
                                                                   % AmAn*time_window = pulse energy
                                                                   % pulse energy/t_rep = W
                                                                   % t_rep: the time for a pulse to finish one round trip
Pdf_ASE_forward  = Power_ASE_forward *df; % W
Pdf_ASE_backward = Power_ASE_backward*df; % W

% AmAn field tensor
diag_idx = shiftdim(permute(0:num_spatial_modes^2*Nt:(num_spatial_modes^2*Nt*(M-1)),[1 3 4 2]) + ...
                    permute(0:num_spatial_modes^2:(num_spatial_modes^2*(Nt-1)),[1 3 2]) + ...
                    (1:(num_spatial_modes+1):num_spatial_modes^2)',-2); % find out the indices of the diagonal elements
if ~first_backward_before_iterations
    if sim.scalar
        % Calculate AmAn
        AmAn_forward  = permute(E_forward ,[4 5 2 6 1 3]).*permute(conj(E_forward ),[4 5 6 2 1 3]); % size: (1,1,num_spatial_modes,num_spatial_modes,Nt,M)
        AmAn_backward = permute(E_backward,[4 5 2 6 1 3]).*permute(conj(E_backward),[4 5 6 2 1 3]); % size: (1,1,num_spatial_modes,num_spatial_modes,Nt,M)
    else
        % Calculate AmAn
        polarized_E_forward  = cat(4,E_forward(:,1:2:end-1,:),E_forward(:,2:2:end,:)); % separate the polarization modes; size: (Nt,num_spatial_modes,M,num_polarizations)
        AmAn_forward_p1  = permute(polarized_E_forward(:,:,:,1),[4 5 2 6 1 3]).*permute(conj(polarized_E_forward(:,:,:,1)),[4 5 6 2 1 3]);
        AmAn_forward_p2  = permute(polarized_E_forward(:,:,:,2),[4 5 2 6 1 3]).*permute(conj(polarized_E_forward(:,:,:,2)),[4 5 6 2 1 3]);
        AmAn_forward  = AmAn_forward_p1 + AmAn_forward_p2; % size: (1,1,num_modes,num_modes,Nt,M)
        polarized_E_backward = cat(4,E_backward(:,1:2:end-1,:),E_backward(:,2:2:end,:)); % separate the polarization modes; size: (Nt,num_spatial_modes,M,num_polarizations)
        AmAn_backward_p1 = permute(polarized_E_backward(:,:,:,1),[4 5 2 6 1 3]).*permute(conj(polarized_E_backward(:,:,:,1)),[4 5 6 2 1 3]);
        AmAn_backward_p2 = permute(polarized_E_backward(:,:,:,2),[4 5 2 6 1 3]).*permute(conj(polarized_E_backward(:,:,:,2)),[4 5 6 2 1 3]);
        AmAn_backward = AmAn_backward_p1 + AmAn_backward_p2; % size: (1,1,num_modes,num_modes,Nt,M)
        
        clear AmAn_forward_p1 AmAn_forward_p2 AmAn_backward_p1 AmAn_backward_p2;
    end
    if isscalar(signal_fields_backward)
        AmAn = AmAn_forward; % No backward fields
    else
        AmAn = AmAn_forward + AmAn_backward;
    end
    clear AmAn_forward AmAn_backward
else
    if sim.gpu_yes
        AmAn = zeros(1,1,num_spatial_modes,num_spatial_modes,Nt,M,'gpuArray');
    else
        AmAn = zeros(1,1,num_spatial_modes,num_spatial_modes,Nt,M);
    end
end

% -------------------------------------------------------------------------
% --------------------- Rate equation to get N ----------------------------
% -------------------------------------------------------------------------

% Ion density in various levels
N = func.solve_N(sim,gain_rate_eqn,...
                 N,gain_rate_eqn.N_total,...
                 E_photon,diag_idx,...
                 gain_rate_eqn.overlap_factor,gain_rate_eqn.cross_sections_pump,gain_rate_eqn.cross_sections,...
                 Power_pump_forward,Power_pump_backward,...
                 Pdf_ASE_forward,Pdf_ASE_backward,...
                 AmAn,...
                 first_backward_before_iterations); % unit: 1/um^3

% -------------------------------------------------------------------------
% -------------- Power equation to get Power_next_step --------------------
% -------------------------------------------------------------------------
if isequal(sim.step_method,'RK4IP') % single mode
    dx = [];
    A_core = pi*(gain_rate_eqn.core_diameter/2)^2;
    dz = sim.dz; % the step size in RK4IP
else % use MPA for multimode
    dx = gain_rate_eqn.mode_profile_dx; % unit: um
    A_core = [];
    switch direction
        case 'forward'
            dz = sim.small_dz; % the small step size in MPA
        case 'backward'
            % During 'backward' propagation, only pump and ASE powers are updated,
            % which don't require MPA.
            dz = sim.dz; % the large step size in MPA
    end
end

% Power
Power_SE = []; % spontaneous emission
switch direction
    case 'forward'
        % at z+dz
        Power_pump_forward = func.solve_Power('pump',...
                                              sim.scalar,...
                                              dz*1e6,dx,A_core,...
                                              num_spatial_modes,gain_rate_eqn.sponASE_spatial_modes,...
                                              gain_rate_eqn.overlap_factor.pump,gain_rate_eqn.cross_sections_pump,...
                                              N,gain_rate_eqn.N_total,...
                                              Power_pump_forward,[],[],[],...
                                              gain_rate_eqn.GammaN,[],...
                                              length(gain_rate_eqn.energy_levels),...
                                              gain_rate_eqn.plusminus,gain_rate_eqn.N_idx); % no spontaneous term for pump
        Power_pump_forward(Power_pump_forward<0) = 0; % pump power cannot be negative
        if gain_rate_eqn.include_ASE
            [Power_ASE_forward,Power_SE] = func.solve_Power('ASE',...
                                                            sim.scalar,...
                                                            dz*1e6,dx,A_core,...
                                                            num_spatial_modes,gain_rate_eqn.sponASE_spatial_modes,...
                                                            gain_rate_eqn.overlap_factor.signal,gain_rate_eqn.cross_sections,...
                                                            N,gain_rate_eqn.N_total,...
                                                            Power_ASE_forward,E_photon,sim.cs.cs,[],...
                                                            [],gain_rate_eqn.FmFnN,...
                                                            length(gain_rate_eqn.energy_levels),...
                                                            gain_rate_eqn.plusminus,gain_rate_eqn.N_idx);
            Power_ASE_forward(Power_ASE_forward<0) = 0; % ASE power cannot be negative
        end
        
        if isequal(sim.step_method,'RK4IP') % RK4IP treats the gain as dispersion, so there is no need to compute g*Am
            field_input = [];
        else % multimode
            if sim.scalar
                field_input = permute(signal_fields,[4 5 2 6 1 3]);
            else % polarized fields
                polarized_fields = cat(4,signal_fields(:,1:2:end-1,:),signal_fields(:,2:2:end,:)); % separate the polarization modes
                field_input = permute(polarized_fields,[5 6 2 7 1 3 4]);
            end
        end
        [~,~,G] = func.solve_Power('signal',...
                                   sim.scalar,...
                                   dz*1e6,dx,A_core,...
                                   num_spatial_modes,gain_rate_eqn.sponASE_spatial_modes,...
                                   gain_rate_eqn.overlap_factor.signal,gain_rate_eqn.cross_sections,...
                                   N,gain_rate_eqn.N_total,...
                                   [],[],[],field_input,...
                                   [],gain_rate_eqn.FmFnN,...
                                   length(gain_rate_eqn.energy_levels),...
                                   gain_rate_eqn.plusminus,gain_rate_eqn.N_idx); % no spontaneous term for signal
    case 'backward'
        % at z-dz
        Power_pump_backward = func.solve_Power('pump',...
                                               sim.scalar,...
                                               dz*1e6,dx,A_core,...
                                               num_spatial_modes,gain_rate_eqn.sponASE_spatial_modes,...
                                               gain_rate_eqn.overlap_factor.pump,gain_rate_eqn.cross_sections_pump,...
                                               N,gain_rate_eqn.N_total,...
                                               Power_pump_backward,[],[],[],...
                                               gain_rate_eqn.GammaN,[],...
                                               length(gain_rate_eqn.energy_levels),...
                                               gain_rate_eqn.plusminus,gain_rate_eqn.N_idx); % no spontaneous term for pump
        Power_pump_backward(Power_pump_backward<0) = 0; % pump power cannot be negative
        if gain_rate_eqn.include_ASE
            Power_ASE_backward = func.solve_Power('ASE',...
                                                  sim.scalar,...
                                                  dz*1e6,dx,A_core,...
                                                  num_spatial_modes,gain_rate_eqn.sponASE_spatial_modes,...
                                                  gain_rate_eqn.overlap_factor.signal,gain_rate_eqn.cross_sections,...
                                                  N,gain_rate_eqn.N_total,...
                                                  Power_ASE_backward,E_photon,sim.cs.cs,[],...
                                                  [],gain_rate_eqn.FmFnN,...
                                                  length(gain_rate_eqn.energy_levels),...
                                                  gain_rate_eqn.plusminus,gain_rate_eqn.N_idx);
            Power_ASE_backward(Power_ASE_backward<0) = 0; % ASE power cannot be negative
        end
end

% Change the size back to (N,num_modes)
if isequal(direction,'forward')
    if sim.scalar || isequal(sim.step_method,'RK4IP')
        G = permute(G,[5 3 6 1 2 4]);
    else % polarized fields
        % "recovery_idx" is used to put the separated polarization modes
        % back into the same dimension of array.
        recovery_idx = [(1:num_spatial_modes);(1:num_spatial_modes)+num_spatial_modes];
        recovery_idx = recovery_idx(:);
        G = cat(3,G(:,:,:,:,:,:,1),G(:,:,:,:,:,:,2));
        G = permute(G(:,:,recovery_idx,:,:,:),[5 3 6 1 2 4]);
    end
end

switch direction
    case 'forward'
        varargout = {Power_pump_forward, Power_ASE_forward,Power_SE,G,N};
    case 'backward'
        varargout = {Power_pump_backward,Power_ASE_backward, N};
end

end