function [E1,a5,...
          opt_dz,success] = stepping_RK4IP_nogain_adaptive(E0,...
                                                           sim, n2_prefactor,...
                                                           D_op, W_op,...
                                                           fr, haw, hbw, sponRS_prefactor, dt,...
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
                E_tr_noise, sponRS_prefactor,...
                dt);
end
a1 = add_DW_RK4IP(a5_1,D_op,W_op,sim.dz/2,sim.adaptive_dz.DW_threshold);
a2 =       N_op(       E_IP+a1*(sim.dz/2),...
                sim, n2_prefactor,...
                fr, haw, hbw,...
                E_tr_noise, sponRS_prefactor,...
                dt);
a3 =       N_op(       E_IP+a2*(sim.dz/2),...
                sim, n2_prefactor,...
                fr, haw, hbw,...
                E_tr_noise, sponRS_prefactor,...
                dt);
a4 =       N_op(add_DW_RK4IP(E_IP+a3*(sim.dz),D_op,W_op,sim.dz/2,sim.adaptive_dz.DW_threshold),...
                sim, n2_prefactor,...
                fr, haw, hbw,...
                E_tr_noise, sponRS_prefactor,...
                dt);

E1 = add_DW_RK4IP(E_IP + (a1+2*a2+2*a3)*(sim.dz/6),D_op,W_op,sim.dz/2,sim.adaptive_dz.DW_threshold) + a4*(sim.dz/6);

% Local error estimate
a5 =       N_op(       E1,...
                sim, n2_prefactor,...
                fr, haw, hbw,...
                E_tr_noise, sponRS_prefactor,...
                dt);
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
                     E_tr_noise, sponRS_prefactor,...
                     dt)
%N_op Calculate dEdz

if any(n2_prefactor(:)) % Compute the nonlinearity only when n2 isn't zero
    if sim.scalar
        E_tr = fftn(E_wk);
    else % don't do fft for the polarization dimension
        E_tr = fft(fft(fft(E_wk,[],1),[],2),[],3);
    end
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
            Ra = dt*fft(haw.*ifft(abs(E_tr).^2));
            Ra_sponRS = dt*fft(haw.*ifft(      E_tr.*conj(E_tr_noise)+...
                                         E_tr_noise.*conj(E_tr)+...
                                         E_tr_noise.*conj(E_tr_noise)).*sponRS_prefactor{2});

            nonlinear_tr = Kerr + 2*(Ra+Ra_sponRS).*E_tr;
        else
            Ra = dt*fft(haw.*ifft(sum(abs(E_tr).^2,5)));
            Ra_sponRS = dt*fft(haw.*ifft(sum(      E_tr.*conj(E_tr_noise)+...
                                             E_tr_noise.*conj(E_tr)+...
                                             E_tr_noise.*conj(E_tr_noise),5)).*sponRS_prefactor{2});

            if isempty(hbw)
                nonlinear_tr = Kerr + 2*(Ra+Ra_sponRS).*E_tr;
            else % polarized field with an anisotropic Raman
                Rb = dt*fft(hbw.*ifft(E_tr.*permute(conj(E_tr),[1,2,3,4,6,5]) + conj(E_tr).*permute(E_tr,[1,2,3,4,6,5])));
                Rb_sponRS = dt*fft(hbw.*ifft(      E_tr.*permute(conj(E_tr_noise),[1,2,3,4,6,5]) + conj(E_tr      ).*permute(E_tr_noise,[1,2,3,4,6,5])+...
                                             E_tr_noise.*permute(conj(E_tr      ),[1,2,3,4,6,5]) + conj(E_tr_noise).*permute(E_tr      ,[1,2,3,4,6,5])+...
                                             E_tr_noise.*permute(conj(E_tr_noise),[1,2,3,4,6,5]) + conj(E_tr_noise).*permute(E_tr_noise,[1,2,3,4,6,5])));

                nonlinear_tr = Kerr + ( 2*(Ra+Ra_sponRS).*E_tr + sum((Rb+Rb_sponRS).*permute(E_tr,[1,2,3,4,6,5]),6) );
            end
        end
    else
        nonlinear_tr = Kerr;
    end

    % Finish adding the prefactor
    dEdz = ifft(ifft(n2_prefactor.*ifft(nonlinear_tr,[],1),[],2),[],3); % nonlinear polarization
else
    dEdz = 0;
end

end