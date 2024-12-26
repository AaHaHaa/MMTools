function [E1,a5,...
          opt_dz,success] = stepping_RK4IP_nogain_adaptive(E0,...
                                                           sim, n2_prefactor,...
                                                           D_op, W_op,...
                                                           fr, haw, hbw,...
                                                           E_tr_noise,...
                                                           a5_1)
%STEPPING_RK4IP_NOGAIN_ADAPTIVE Take one step with RK4IP without gain

% Represented under the interaction picture
E_IP = add_DW_RK4IP(E0,D_op,W_op,sim.dz/2,sim.adaptive_dz.DW_threshold);

% Propagate through the nonlinearity
if isempty(a5_1)
    a5_1 = N_op(       E0,...
                sim, n2_prefactor,...
                fr, haw, hbw,...
                E_tr_noise);
end
a1 = add_DW_RK4IP(a5_1,D_op,W_op,sim.dz/2,sim.adaptive_dz.DW_threshold);
a2 =       N_op(       E_IP+a1*(sim.dz/2),...
                sim, n2_prefactor,...
                fr, haw, hbw,...
                E_tr_noise);
a3 =       N_op(       E_IP+a2*(sim.dz/2),...
                sim, n2_prefactor,...
                fr, haw, hbw,...
                E_tr_noise);
a4 =       N_op(add_DW_RK4IP(E_IP+a3*(sim.dz),D_op,W_op,sim.dz/2,sim.adaptive_dz.DW_threshold),...
                sim, n2_prefactor,...
                fr, haw, hbw,...
                E_tr_noise);

E1 = add_DW_RK4IP(E_IP + (a1+2*a2+2*a3)*(sim.dz/6),D_op,W_op,sim.dz/2,sim.adaptive_dz.DW_threshold) + a4*(sim.dz/6);

% Local error estimate
a5 =       N_op(       E1,...
                sim, n2_prefactor,...
                fr, haw, hbw,...
                E_tr_noise);
err = sum(abs((a4(:)-a5(:))*(sim.dz/10)).^2);

% Stepsize control
normE2 = sum(abs(E1(:)).^2);
err = sqrt(err/normE2);
if isnan(normE2) % the computation is just so wrong, so we reduce the step size and do it again
    opt_dz = 0.5*sim.dz;
    success = false;
elseif normE2 == 0 % all-zero field
    opt_dz = 2*sim.dz;
    success = true;
else
    opt_dz = max(0.5,min(2,0.8*(sim.adaptive_dz.threshold/err)^(1/4)))*sim.dz;

    success = err < sim.adaptive_dz.threshold;
end

end

function dEdz = N_op(E_wk,...
                     sim, n2_prefactor,...
                     fr, haw, hbw,...
                     E_tr_noise)
%N_op Calculate dEdz

if any(n2_prefactor(:)) % Compute the nonlinearity only when n2 isn't zero
    E_tr = ifft(ifft(fft(E_wk,[],1),[],2),[],3);
    E_tr_wNoise = E_tr + E_tr_noise;

    % Kerr term
    if sim.scalar
        Kerr = (1-fr)*3*E_tr_wNoise.*abs(E_tr_wNoise).^2;
    else
        Kerr = (1-fr)*(conj(E_tr_wNoise).*sum(E_tr_wNoise.^2,5) + 2*E_tr_wNoise.*sum(abs(E_tr_wNoise).^2,5));
    end

    % Raman term
    if sim.include_Raman
        if sim.scalar
            Ra = fft(haw.*ifft(abs(E_tr_wNoise).^2));
            
            nonlinear_tr = Kerr + 3*Ra.*E_tr_wNoise;
        else
            Ra = fft(haw.*ifft(sum(abs(E_tr_wNoise).^2,5)));

            if isempty(hbw)
                nonlinear_tr = Kerr + 3*Ra.*E_tr_wNoise;
            else % polarized field with an anisotropic Raman
                Rb = fft(hbw.*ifft(E_tr_wNoise.*permute(conj(E_tr_wNoise),[1,2,3,4,6,5]) + conj(E_tr_wNoise).*permute(E_tr_wNoise,[1,2,3,4,6,5])));

                nonlinear_tr = Kerr + ( 3*Ra.*E_tr_wNoise + 3/2*sum(Rb.*permute(E_tr_wNoise,[1,2,3,4,6,5]),6) );
            end
        end
    else
        nonlinear_tr = Kerr;
    end

    % Finish adding the prefactor
    dEdz = fft(fft(n2_prefactor.*ifft(nonlinear_tr,[],1),[],2),[],3); % nonlinear polarization
else
    dEdz = 0;
end

end