function [E1,a5,...
          opt_deltaZ,success] = stepping_RK4IP_nogain_adaptive(E0,...
                                                               sim, prefactor,...
                                                               D_op, W_op,...
                                                               fr, haw, hbw, sponRS_prefactor, dt,...
                                                               a5_1)
%STEPPING_RK4IP_NOGAIN_ADAPTIVE Take one step with RK4IP without gain
%
% Input:
%    A0 - initial field in the frequency domain (N, num_modes); sqrt(W)
%    dummy_dt - time grid point spacing; ps; but this isn't used in this function
%
%    sim.scalar - scalar or polarized fields
%    sim.gpu_yes - true = GPU, false = CPU
%
%    sim.cuda_SRSK - the cuda for computing SR and SK values
%
%    sim.Raman_model - which Raman model is used
%    sim.Raman_sponRS - consider spontaneous Raman or not
%
%    sim.adaptive_deltaZ.threshold - the accuracy used to determined whether to increase or decrease the step size
%
%    prefactor - 1i*n2*omega/c; m/W
%
%    SRa_info.SRa - SRa tensor; m^-2
%    SRa_info.nonzero_midx1234s - required SRa indices in total
%    SRa_info.nonzero_midx34s - required (SRa) indices for partial Raman term (only for CPU computation)
%    SRb_info.SRb - SRb tensor; m^-2
%    SRb_info.nonzero_midx1234s - required SRb indices in total
%    SRb_info.nonzero_midx34s - required (SRb) indices for partial Raman term (only for CPU computation)
%    SK_info.SK - SK tensor; m^2 (unempty if considering polarizaton modes)
%    SK_info.nonzero_midx1234s - required SK indices in total (unempty if considering polarizaton modes)
%
%    D_op - dispersion term Dz/2 (N, num_modes)
%    haw - isotropic Raman response in the frequency domain
%    hbw - anisotropic Raman response in the frequency domain
%
%    sponRS_prefactor - prefactor for the spontaneous Raman scattering
%
%    a5_1 - the previous RK4 term that can be reused
%
%    dummy_G - gain term (N, 1); but this isn't used in this function
%    dummy_saturation_intensity - saturation intensity; J/m^2; but this isn't used in this function
%    dummy_fr - Raman fraction; but this isn't used in this function
%
% Output:
%    A1 - the field (in the frequency domain) after one step size (N, num_modes)
%    a5 - the RK4 term that can be reused in the next step
%    opt_deltaZ - recommended step size
%    success - whether the current step size is sufficiently small for the required tolerance

anisotropic_Raman_included = ~sim.scalar & sim.Raman_model==2;

% Spontaneous Raman scattering
if sim.Raman_model ~= 0 && sim.Raman_sponRS
    sponRS = ifft(abs(fft(sponRS_prefactor{1}.*randn(size(sponRS_prefactor{1})).*exp(1i*2*pi*rand(size(sponRS_prefactor{1}))))).^2).*sponRS_prefactor{2};
    sponRS_Gamma = fft(haw.*sponRS);
else
    sponRS_Gamma = 0;
end

% Represented under the interaction picture
E_IP = add_DW_RK4IP(E0,D_op,W_op,sim.deltaZ/2,sim.adaptive_deltaZ.DW_threshold);

% Propagate through the nonlinearity
if isempty(a5_1)
    a5_1 = N_op(       E0,...
                sim, prefactor, fr, haw, hbw, sponRS_Gamma, anisotropic_Raman_included, dt);
end
a1 = add_DW_RK4IP(a5_1,D_op,W_op,sim.deltaZ/2,sim.adaptive_deltaZ.DW_threshold);
a2 =       N_op(       E_IP+a1*(sim.deltaZ/2),...
                sim, prefactor, fr, haw, hbw, sponRS_Gamma, anisotropic_Raman_included, dt);
a3 =       N_op(       E_IP+a2*(sim.deltaZ/2),...
                sim, prefactor, fr, haw, hbw, sponRS_Gamma, anisotropic_Raman_included, dt);
a4 =       N_op(add_DW_RK4IP(E_IP+a3*(sim.deltaZ),D_op,W_op,sim.deltaZ/2,sim.adaptive_deltaZ.DW_threshold),...
                sim, prefactor, fr, haw, hbw, sponRS_Gamma, anisotropic_Raman_included, dt);

E1 = add_DW_RK4IP(E_IP + (a1+2*a2+2*a3)*(sim.deltaZ/6),D_op,W_op,sim.deltaZ/2,sim.adaptive_deltaZ.DW_threshold) + a4*(sim.deltaZ/6);

% Local error estimate
a5 =       N_op(       E1,...
                sim, prefactor, fr, haw, hbw, sponRS_Gamma, anisotropic_Raman_included, dt);
err = sum(abs((a4(:)-a5(:))*(sim.deltaZ/10)).^2);

% Stepsize control
normE2 = sum(abs(E1(:)).^2);
err = sqrt(err/normE2);
if isnan(normE2) % the computation is just so wrong, so we reduce the step size and do it again
    opt_deltaZ = 0.5*sim.deltaZ;
    success = false;
elseif normE2 == 0 % all-zero field
    opt_deltaZ = 2*sim.deltaZ;
    success = true;
else
    opt_deltaZ = max(0.5,min(2,0.8*(sim.adaptive_deltaZ.threshold/err)^(1/4)))*sim.deltaZ;

    success = err < sim.adaptive_deltaZ.threshold;
end

end

function dEdz = N_op(E_wk,...
                     sim, prefactor,...
                     fr, haw, hbw, sponRS_Gamma, anisotropic_Raman_included, dt)
%N_op Calculate dEdz

if any(prefactor(:)) % Compute the nonlinearity only when n2 isn't zero
    if sim.scalar
        E_tr = fftn(E_wk);
    else % don't do fft for the polarization dimension
        E_tr = fft(fft(fft(E_wk,[],1),[],2),[],3);
    end

    % Kerr term
    if sim.scalar
        Kerr = (1-fr)*3*E_tr.*abs(E_tr).^2;
    else
        Kerr = (1-fr)*(conj(E_tr).*sum(E_tr.^2,5) + 2*E_tr.*sum(abs(E_tr).^2,5));
    end

    % Raman term
    if sim.Raman_model ~= 0
        if sim.scalar
            if ~anisotropic_Raman_included
                R = dt*fft(haw.*ifft(abs(E_tr).^2));
            else % polarized field with an anisotropic Raman
                R = dt*fft((haw+hbw).*ifft(abs(E_tr).^2));
            end
            nonlinear_tr = Kerr + 2*R.*E_tr; clearvars Kerr R;
        else
            Ra = dt*fft(haw.*ifft(sum(abs(E_tr).^2,5)));

            if ~anisotropic_Raman_included
                nonlinear_tr = Kerr + 2*Ra.*E_tr; clearvars Kerr Ra;
            else % polarized field with an anisotropic Raman
                Rb = dt*fft(hbw.*ifft(E_tr.*permute(conj(E_tr),[1,2,3,4,6,5]) + conj(E_tr).*permute(E_tr,[1,2,3,4,6,5])));

                nonlinear_tr = Kerr + ( 2*Ra.*E_tr + sum(Rb.*permute(E_tr,[1,2,3,4,6,5]),6) ); clearvars Kerr Ra Rb;
            end
        end
        nonlinear_tr = nonlinear_tr + sponRS_Gamma.*E_tr;
    else
        nonlinear_tr = Kerr; clearvars Kerr;
    end

    % Finish adding the prefactor
    dEdz = ifft(ifft(prefactor.*ifft(nonlinear_tr,[],1),[],2),[],3); % nonlinear polarization
else
    dEdz = 0;
end

end