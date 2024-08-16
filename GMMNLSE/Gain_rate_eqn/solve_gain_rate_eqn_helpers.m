function func = solve_gain_rate_eqn_helpers
%SOLVE_GAIN_RATE_EQN_HELPERS containers for the rate-equation solvers

func.solve_N = @solve_N; % solve the population
func.solve_Power = @solve_Power; % solve the pump and ASE powers, and the signal gain

end

%% Contained function: solve_N
function N = solve_N(sim,gain_rate_eqn,...
                     N,N_total,...
                     E_photon,diag_idx,...
                     overlap_factor,cross_sections_pump,cross_sections,...
                     Power_pump_forward,Power_pump_backward,...
                     Pdf_ASE_forward,Pdf_ASE_backward,...
                     AmAn,...
                     first_backward_before_iterations)
%SOLVE_N It solves the steady-state populations among various energy
%levels.
%
% This function is a container of the population solver.
% The solver is split into two helper functions, one for only-two-level
% systems and the other for the general multi-level systems. Because
% two-level systems don't have complicated nonlinear terms, such as
% cross-relaxation effects in Er, Tm, etc., their population can be easily
% solved with
%   N2 = absorption_term / (absorption_term + emission_term + 1/tau), tau: upper-state lifetime
%
% However, for multi-level systems, a numerical solver is required. I
% implemented the trust-region method with a dogleg sub-problem solver. In
% addition, to solve it fast, the population from the previous propagation
% step is sent as an initial guess.
% Due to these two processes, it's helpful to split the population solver 
% into two functions based on whether it's a two-level system or not.

if length(gain_rate_eqn.energy_levels) == 2 % two-level system
    N = solve_N_2levels(sim,gain_rate_eqn,...
                        N_total,...
                        E_photon,diag_idx,...
                        overlap_factor,cross_sections_pump,cross_sections,...
                        Power_pump_forward,Power_pump_backward,...
                        Pdf_ASE_forward,Pdf_ASE_backward,...
                        AmAn,...
                        first_backward_before_iterations);
else % multi-level system
    N = solve_N_levels(sim,gain_rate_eqn,...
                       N,N_total,...
                       E_photon,diag_idx,...
                       overlap_factor,cross_sections_pump,cross_sections,...
                       Power_pump_forward,Power_pump_backward,...
                       Pdf_ASE_forward,Pdf_ASE_backward,...
                       AmAn,...
                       first_backward_before_iterations);
end

end

function N1 = solve_N_2levels(sim,gain_rate_eqn,...
                              N_total,...
                              E_photon,diag_idx,...
                              overlap_factor,cross_sections_pump,cross_sections,...
                              Power_pump_forward,Power_pump_backward,...
                              Pdf_ASE_forward,Pdf_ASE_backward,...
                              AmAn,...
                              first_backward_before_iterations)
%SOLVE_N_2LEVELS solves the ion density in the upper state under the steady 
%state of the rate equation
%
% computational dimension here: (Nx,Nx,num_spatial_modes,num_spatial_modes,Nt,M,num_polarizations,num_cross_sections/num_levels), M: parallelization in MPA

h = 6.62607004e-34;
c = 299792458;

% Rate/photon_energy:
% size: (Nx,Nx,num_modes,Nt), unit: W
% pump
sN = size(N_total);
sP = size(Power_pump_forward,1:7); % the 1st-6th dimensions are the same as defined;
                                   % the 7th dimension is left empty;
                                   % the 8th dimension is to store various cross-section data
num_cross_sections = length(cross_sections_pump);
if sim.gpu_yes
    R_forward  = zeros([sN,sP(3:end),num_cross_sections],'gpuArray');
    R_backward = zeros([sN,sP(3:end),num_cross_sections],'gpuArray');
else
    R_forward  = zeros([sN,sP(3:end),num_cross_sections]);
    R_backward = zeros([sN,sP(3:end),num_cross_sections]);
end
if ismember(gain_rate_eqn.pump_direction,{'co','bi'})
    R_forward  = overlap_factor.pump.*cross_sections_pump.*Power_pump_forward /h*(gain_rate_eqn.pump_wavelength*1e-9)/c; % 1/s
end
if ismember(gain_rate_eqn.pump_direction,{'counter','bi'})
    R_backward = overlap_factor.pump.*cross_sections_pump.*Power_pump_backward/h*(gain_rate_eqn.pump_wavelength*1e-9)/c; % 1/s
end
R.pump = R_forward + R_backward; % 1/s

% ASE and signal
if gain_rate_eqn.include_ASE
    if sim.scalar
        Power_ASE = Pdf_ASE_forward + Pdf_ASE_backward;
    else % polarized fields
        Power_ASE = Pdf_ASE_forward(:,:,1:2:end-1,:,:,:)  + Pdf_ASE_forward(:,:,2:2:end,:,:,:) + ... % forward
                    Pdf_ASE_backward(:,:,1:2:end-1,:,:,:) + Pdf_ASE_backward(:,:,2:2:end,:,:,:);     % backward
    end
    if first_backward_before_iterations % For bi/counter-pumping cases, the first backward propagation doesn't consider the signal fields.
        AmAn(diag_idx) = Power_ASE; % AmAn(:,:,n,n,:,:) = Power_ASE(:,:,n,:,:,:)
    else
        AmAn(diag_idx) = AmAn(diag_idx) + Power_ASE; % AmAn(:,:,n,n,:,:) = AmAn(:,:,n,n,:,:) + Power_ASE(:,:,n,:,:,:)
    end
end

integral_mn = sum(AmAn.*(cross_sections./E_photon),5); % size: (1,1,num_modes,num_modes,1,M,1,[GSA01,emi10])

% For SMF (single mode only), the computations below all have the length 1 
% or M, and thus can be faster with CPU, instead of GPU.
if sim.gpu_yes && length(sim.midx) == 1 % single mode
    integral_mn = gather(integral_mn);
end

% It's real! Use "real" to save the memory.
R.ASE_signal = real(sum(integral_mn.*overlap_factor.signal,[3,4])); % 1/s; size: (Nx,Nx,1,1,1,M,1,[GSA01,emi10])

% ion density in the upper state
% N2 = N_total*sum(all the absorption terms)/
%      ( sum(all the absorption and emission terms) + upper-state-lifetime term )
total_absorption = R.pump(:,:,:,:,:,:,:,1) + R.ASE_signal(:,:,:,:,:,:,:,1);
total_emission   = R.pump(:,:,:,:,:,:,:,2) + R.ASE_signal(:,:,:,:,:,:,:,2);
N1 = N_total.*total_absorption./... % size: (Nx,Nx,1,1,1,M)
     (total_absorption + total_emission + ... % absorption and emission terms
      1/gain_rate_eqn.N.eqn.tau);                   % upper-state-lifetime term

end

function N = solve_N_levels(sim,gain_rate_eqn,...
                            N,N_total,...
                            E_photon,diag_idx,...
                            overlap_factor,cross_sections_pump,cross_sections,...
                            Power_pump_forward,Power_pump_backward,...
                            Pdf_ASE_forward,Pdf_ASE_backward,...
                            AmAn,...
                            first_backward_before_iterations)
%SOLVE_N_LEVELS It solves the steady-state populations among various energy
%levels.
%
% computational dimension here: (Nx,Nx,num_spatial_modes,num_spatial_modes,Nt,M,num_polarizations,num_cross_sections/num_levels), M: parallelization in MPA

h = 6.62607004e-34;
c = 299792458;

% Rate/photon_energy:
% size: (Nx,Nx,num_modes,Nt), unit: W
% pump
sN = size(N_total);
sP = size(Power_pump_forward,1:7); % the 1st-6th dimensions are the same as defined;
                                   % the 7th dimension is left empty;
                                   % the 8th dimension is to store various cross-section data
num_cross_sections = length(cross_sections_pump);
if sim.gpu_yes
    R_forward  = zeros([sN,sP(3:end),num_cross_sections],'gpuArray');
    R_backward = zeros([sN,sP(3:end),num_cross_sections],'gpuArray');
else
    R_forward  = zeros([sN,sP(3:end),num_cross_sections]);
    R_backward = zeros([sN,sP(3:end),num_cross_sections]);
end
if ismember(gain_rate_eqn.pump_direction,{'co','bi'})
    R_forward  = overlap_factor.pump.*cross_sections_pump.*Power_pump_forward /h*(gain_rate_eqn.pump_wavelength*1e-9)/c; % 1/s
end
if ismember(gain_rate_eqn.pump_direction,{'counter','bi'})
    R_backward = overlap_factor.pump.*cross_sections_pump.*Power_pump_backward/h*(gain_rate_eqn.pump_wavelength*1e-9)/c; % 1/s
end
R.pump = R_forward + R_backward;

% ASE and signal
if gain_rate_eqn.include_ASE
    if sim.scalar
        Power_ASE = Pdf_ASE_forward + Pdf_ASE_backward;
    else % polarized fields
        Power_ASE = Pdf_ASE_forward(:,:,1:2:end-1,:,:,:)  + Pdf_ASE_forward(:,:,2:2:end,:,:,:) + ... % forward
                    Pdf_ASE_backward(:,:,1:2:end-1,:,:,:) + Pdf_ASE_backward(:,:,2:2:end,:,:,:);     % backward
    end
    if first_backward_before_iterations % For bi/counter-pumping cases, the first backward propagation doesn't consider the signal fields.
        AmAn(diag_idx) = Power_ASE; % AmAn(:,:,n,n,:,:) = Power_ASE(:,:,n,:,:,:)
    else
        AmAn(diag_idx) = AmAn(diag_idx) + Power_ASE; % AmAn(:,:,n,n,:,:) = AmAn(:,:,n,n,:,:) + Power_ASE(:,:,n,:,:,:)
    end
end

integral_mn = sum(AmAn.*(cross_sections./E_photon),5); % size: (1,1,num_modes,num_modes,1,M,1,num_cross_sections)

% For SMF (single mode only), the computations below all have the length 1 
% or M, and thus can be faster with CPU, instead of GPU.
if sim.gpu_yes && length(sim.midx) == 1 % single mode
    integral_mn = gather(integral_mn);
end

% It's real! Use "real" to save the memory.
R.ASE_signal = real(sum(integral_mn.*overlap_factor.signal,[3,4])); % 1/s; size: (Nx,Nx,1,1,1,M,1,num_cross_sections)

% Compute the steady-state populations
R = shiftdim(R.pump + R.ASE_signal,-2); % shift the dimension to prepare for the solver
N_total = shiftdim(N_total,-2);
N = permute(N,[8,9,1,2,3,4,5,6,7]); % It's important to have N, in principle from the previous propagation step, as an initial guess so that the trust-region solver for the coupled equation of the population can find a solution fast.
F = @(N) gain_rate_eqn.N.eqn.ss.F(N,N_total,gain_rate_eqn.N.eqn.Arad,gain_rate_eqn.N.eqn.Gammai,gain_rate_eqn.N.eqn.kijkl,R);
J = @(N) gain_rate_eqn.N.eqn.ss.J(N,N_total,gain_rate_eqn.N.eqn.Arad,gain_rate_eqn.N.eqn.Gammai,gain_rate_eqn.N.eqn.kijkl,R);
N = myTrustRegion(F,J,N,10,1e-4,size(R,8),N_total,sim.gpu_yes);
N = permute(N,[3,4,5,6,7,8,9,1,2]); % change it to the original dimension, with the 8th dimension storing various N, population of each level

end

%% Contained function: solve_Power
function [Pnext,P_spon,signal_out] = solve_Power(field_type,...
                                                 isscalar,...
                                                 dz,dx,A_core,...
                                                 num_spatial_modes,sponASE_spatial_modes,...
                                                 overlap_factor,cross_sections,...
                                                 N,N_total,...
                                                 P0,E_photon,Am,...
                                                 GammaN,FmFnN,...
                                                 num_levels,...
                                                 plusminus,N_idx)
%SOLVE_POWER This is a container for the functions of solving the pump and 
%ASE powers, and the signal gain in the current propagating step.

if num_levels == 2 % two-level system
    [Pnext,P_spon,signal_out] = solve_Power_2levels(field_type,...
                                                    isscalar,...
                                                    dz,dx,A_core,...
                                                    num_spatial_modes,sponASE_spatial_modes,...
                                                    overlap_factor,cross_sections,...
                                                    N,N_total,...
                                                    P0,E_photon,Am,...
                                                    GammaN,FmFnN);
else % multi-level system
    [Pnext,P_spon,signal_out] = solve_Power_levels(field_type,...
                                                   isscalar,...
                                                   dz,dx,A_core,...
                                                   num_spatial_modes,sponASE_spatial_modes,...
                                                   overlap_factor,cross_sections,...
                                                   N,N_total,...
                                                   P0,E_photon,Am,...
                                                   plusminus,N_idx);
end

end

function [Pnext,P_spon,signal_out] = solve_Power_2levels(field_type,...
                                                         isscalar,...
                                                         dz,dx,A_core,...
                                                         num_spatial_modes,sponASE_spatial_modes,...
                                                         overlap_factor,cross_sections,...
                                                         N1,N_total,...
                                                         P0,E_photon,Am,...
                                                         GammaN,FmFnN)
%SOLVE_POWER_2LEVELS solves Power(z+dz) for pump, ASE, and signal.
%
%   dz: um
%

signal_out = [];
P_spon = [];

% Photon energy in the spontaneous emission
if isequal(field_type,'ASE') % ASE
    if isscalar
        % Even in scalar computations, the gain is saturated by the spontaneous emission of two polarizations
        E_photon = 2*E_photon*1e12*sponASE_spatial_modes; % W/THz; 1e12 is to transform Hz into THz (J=W/Hz)
    else
        E_photon = E_photon*1e12*sponASE_spatial_modes; % W/THz; 1e12 is to transform Hz into THz (J=W/Hz)
    end
end

cross_section_all = sum(cross_sections,8);

if isempty(dx) % single mode; with RK4IP
    overlap_factor = overlap_factor*A_core; % For the single mode, the integral w.r.t. x and y can be done first, which becomes overlap_factor here.
    
    switch field_type
        case 'signal' % amplification factor
            tmp = 1 + overlap_factor*(cross_section_all*N1 - cross_sections(:,:,:,:,:,:,:,1)*N_total)*dz;
            tmp(tmp < 0) = 0; % Sometimes, if the factor is too close zero, it can be negative due to the numerically precision.
            signal_out = sqrt(tmp);
        case 'ASE'
            fz_spon = ( overlap_factor*cross_sections(:,:,:,:,:,:,:,2)*N1.*E_photon )*dz; % spontaneous emission; unit: W/THz
            fz = ( overlap_factor*(cross_section_all*N1.*P0 - cross_sections(:,:,:,:,:,:,:,1)*N_total.*P0) )*dz + fz_spon; % unit: W/THz
        case 'pump' % no spontaneous term
            fz = overlap_factor*(cross_section_all*N1 - cross_sections(:,:,:,:,:,:,:,1)*N_total).*P0*dz; % unit: W
    end
    
else % multimode; with MPA
    trapz2 = @(x) sum(x,[1,2])*dx^2; % take the integral w.r.t. the x-y plane
    
    switch field_type
        case 'signal' % ignore spontaneous term
            FmFnN1 = trapz2(overlap_factor.*N1);
            
            % Calculate gA*dz after passing through the gain.
            signal_out = permute( sum( dz/2.*(cross_section_all.*FmFnN1 - cross_sections(:,:,:,:,:,:,:,1).*FmFnN).*Am ,3) ,[1 2 4 3 5 6 7]); % Am_after_gain - Am
        case 'ASE'
            if ~isscalar % polarized fields
                P0 = cat(7,P0(:,:,1:2:end-1,:,:,:),P0(:,:,2:2:end,:,:,:));
            end
            diag_idx = sub2ind([num_spatial_modes num_spatial_modes],1:num_spatial_modes,1:num_spatial_modes);
            GammaN1 = trapz2(overlap_factor(:,:,diag_idx).*N1); % overlap_factor*N2
            GammaN = FmFnN(:,:,diag_idx); % overlap_factor*N_total
            
            fz_spon = real( GammaN1.*cross_sections(:,:,:,:,:,:,:,2).*E_photon ).*dz; % spontaneous emission; W/THz
            fz = real( GammaN1.*cross_section_all.*P0 - GammaN.*cross_sections(:,:,:,:,:,:,:,1).*P0 ).*dz + fz_spon; % W/THz
        case 'pump' % no spontaneous term
            GammaN1 = trapz2(overlap_factor.*N1);
            fz = ( GammaN1*cross_section_all - GammaN*cross_sections(:,:,:,:,:,:,:,1) ).*P0.*dz; % W
    end
end

switch field_type
    case 'signal'
        Pnext = [];
    case 'ASE'
        P_spon = calc_Pnext(zeros(size(P0)),fz_spon,num_spatial_modes);
        Pnext = calc_Pnext(P0,fz,num_spatial_modes);
    case 'pump'
        Pnext = calc_Pnext(P0,fz,num_spatial_modes);
end

end

function [Pnext,P_spon,signal_out] = solve_Power_levels(field_type,...
                                                        isscalar,...
                                                        dz,dx,A_core,...
                                                        num_spatial_modes,sponASE_spatial_modes,...
                                                        overlap_factor,cross_sections,...
                                                        N,N_total,...
                                                        P0,E_photon,Am,...
                                                        plusminus,N_idx)
%SOLVE_POWER_LEVELS solves Power(z+dz) for pump, ASE, and signal.
%
%   dz: um
%

signal_out = [];
P_spon = [];

N = cat(8,N_total-sum(N,8),N); % include the ground-state population

% Photon energy in the spontaneous emission
if isequal(field_type,'ASE') % ASE
    if isscalar
        % Even in scalar computations, the gain is saturated by the spontaneous emission of two polarizations
        E_photon = 2*E_photon*1e12*sponASE_spatial_modes; % W/THz; 1e12 is to transform Hz into THz (J=W/Hz)
    else
        E_photon = E_photon*1e12*sponASE_spatial_modes; % W/THz; 1e12 is to transform Hz into THz (J=W/Hz)
    end
end

if isempty(dx) % single mode; with RK4IP
    overlap_factor = overlap_factor*A_core; % For the single mode, the integral w.r.t. x and y can be done first, which becomes overlap_factor here.
    
    switch field_type
        case 'signal' % amplification factor
            tmp = 1 + overlap_factor*sum(plusminus.*cross_sections.*N(:,:,:,:,:,:,:,N_idx),8)*dz;
            tmp(tmp < 0) = 0; % Sometimes, if the factor is too close zero, it can be negative due to the numerically precision.
            signal_out = sqrt(tmp);
        case 'ASE'
            spon_idx = (squeeze(plusminus)>0);
            fz_spon = overlap_factor*sum(plusminus(spon_idx).*cross_sections(:,:,:,:,:,:,:,spon_idx).*N(:,:,:,:,:,:,:,N_idx(spon_idx)),8).*E_photon*dz; % spontaneous emission; unit: W/THz
            fz = overlap_factor*sum(plusminus.*cross_sections.*N(:,:,:,:,:,:,:,N_idx),8).*P0*dz + fz_spon; % unit: W/THz
        case 'pump' % no spontaneous term
            fz = overlap_factor*sum(plusminus.*cross_sections.*N(:,:,:,:,:,:,:,N_idx),8).*P0*dz; % unit: W
    end
    
else % multimode; with MPA
    trapz2 = @(x) sum(x,[1,2])*dx^2; % take the integral w.r.t. the x-y plane
    
    switch field_type
        case 'signal' % ignore spontaneous term
            FmFnN = trapz2(overlap_factor.*N(:,:,:,:,:,:,:,N_idx));
            
            % Calculate gA*dz after passing through the gain.
            signal_out = permute( sum( dz/2.*sum(plusminus.*cross_sections.*FmFnN,8).*Am ,3) ,[1 2 4 3 5 6 7]); % Am_after_gain - Am
        case 'ASE'
            if ~isscalar % polarized fields
                P0 = cat(7,P0(:,:,1:2:end-1,:,:,:),P0(:,:,2:2:end,:,:,:));
            end
            diag_idx = sub2ind([num_spatial_modes num_spatial_modes],1:num_spatial_modes,1:num_spatial_modes);
            GammaN = trapz2(overlap_factor(:,:,diag_idx).*N(:,:,:,:,:,:,:,N_idx)); % overlap_factor*N2
            
            spon_idx = (squeeze(plusminus)>0);
            fz_spon = sum(plusminus(spon_idx).*cross_sections(:,:,:,:,:,:,:,spon_idx).*GammaN(:,:,:,:,:,:,:,spon_idx),8).*E_photon*dz; % spontaneous emission; unit: W/THz
            fz = sum(plusminus.*cross_sections.*GammaN,8).*P0*dz + fz_spon; % unit: W/THz
        case 'pump' % no spontaneous term
            GammaN = trapz2(overlap_factor.*N(:,:,:,:,:,:,:,N_idx));
            fz = sum(plusminus.*cross_sections.*GammaN,8).*P0*dz; % unit: W
    end
end

switch field_type
    case 'signal'
        Pnext = [];
    case 'ASE'
        P_spon = calc_Pnext(zeros(size(P0)),fz_spon,num_spatial_modes);
        Pnext = calc_Pnext(P0,fz,num_spatial_modes);
    case 'pump'
        Pnext = calc_Pnext(P0,fz,num_spatial_modes);
end

end

%% Helper function for solve_Power()
function Pnext = calc_Pnext(P0,fz,num_spatial_modes)
%CALC_PNEXT

if size(fz,6) == 1 % single mode with RK4IP (no parallelization; M=1) or
                   % 'backward' propagation which updates only pump and ASE powers
    Pnext = P0 + fz;
else % 'forward' propagation in multimode (with MPA)
    % Pump and ASE powers are implemented with MPA as well which relies
    % on parallelization.
    %
    % Pnext = trapz(z=(0:M)*dz,fz) + P0
    fz(:,:,:,:,:,1,:) = fz(:,:,:,:,:,1,:)/2;
    if size(P0,6) == 1
        Pnext = repmat(P0,1,1,1,1,1,size(fz,6),1);
    else
        Pnext = P0;
    end
    Pnext(:,:,:,:,:,2:end,:) = P0(:,:,:,:,:,1,:) + cumsum(fz(:,:,:,:,:,1:end-1,:),6) + fz(:,:,:,:,:,2:end,:)/2;
end

if size(P0,7) == 2 % polarized multimode ASE
    % "recovery_idx" is used to put the separated polarization modes
    % back into the same dimension of array.
    recovery_idx = [(1:num_spatial_modes);(1:num_spatial_modes)+num_spatial_modes];
    recovery_idx = recovery_idx(:);

    Pnext = cat(3,Pnext(:,:,:,:,:,:,1),Pnext(:,:,:,:,:,:,2));
    Pnext = Pnext(:,:,recovery_idx,:,:,:);
end

end