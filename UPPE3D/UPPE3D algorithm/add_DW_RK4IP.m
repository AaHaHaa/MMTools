function [last_E,...
          opt_dz] = add_DW_RK4IP(E0,...
                                 F_op,...
                                 D_op,W_op,...
                                 dz,...
                                 L0,...
                                 Kz_correction,...
                                 DW_threshold,...
                                 r)
%ADD_DW_RK4IP This computes the exp(D_op+W_op)*E0 with the adaptive-step RK4IP

if any(W_op(:))
    apply_W = @(x) Kz_correction.*F_op.Fk(W_op.*F_op.iFk(x)); % apply waveguide effect in position space
else
    apply_W = @(x) 0; % no waveguide effect
end

last_E = E0;
z = 0;
a5 = [];
dz = min(dz,L0/2);
while z+eps(z) < L0 % eps(z) here is necessary due to the numerical error
    ever_fail = false;
    previous_E = last_E;
    previous_a5 = a5;

    success = false;
    while ~success
        if ever_fail
            last_E = previous_E;
            a5 = previous_a5;
        end

        [last_E,a5,...
         opt_dz, success] = DW_stepping(last_E,a5,...
                                        D_op, apply_W,...
                                        dz,...
                                        DW_threshold,...
                                        r);

        if ~success
            ever_fail = true;

            dz = opt_dz;
        end
    end

    % Update z
    z = z + dz;
    dz = min([opt_dz,L0-z]);
end

end

function [E1,a5,...
          opt_dz, success] = DW_stepping(E0,a5_1,...
                                         D_op, apply_W,...
                                         dz,...
                                         DW_threshold,...
                                         r)

D = exp(D_op*dz/2);

% Represented under the interaction picture
E_IP = D.*E0;

% Propagate through the nonlinearity
if isempty(a5_1)
    a5_1 = apply_W(E0);
end
a1 = D.*a5_1;
a2 = apply_W(    E_IP+a1*(dz/2));
a3 = apply_W(    E_IP+a2*(dz/2));
a4 = apply_W(D.*(E_IP+a3*(dz)) );

E1 = D.*(E_IP + (a1+2*a2+2*a3)*(dz/6)) + a4*(dz/6);


% Normalization
% Dispersion and waveguide operations shouldn't induce loss. The artificial
% loss results from approximation error of sequentially computing two
% operations rather than only one. Baker–Campbell–Hausdorff formula shows
% that such sequential operation is simply an approximation with an error
% corresponding to their commutator [D,W].
% Therefore, below I make E1 have the same energy as E0 artificially.
%
% Because radially-symmetric scheme relies on the exponential sampling
% strategy of FHATHA, an integral with trapz with non-uniform sampling
% points is used, rather than a simple summation in the xy scheme.
norm2_E = sum(abs(E0(:)).^2);
if norm2_E ~= 0 % all-zero field
    if isempty(r)
        E1 = E1*sqrt(norm2_E/sum(abs(E1(:)).^2));
    else % use the radially-symmetric scheme
        norm2_E = sum(trapz(r,abs(E0).^2.*r,2),[1,3,4]);
        E1 = E1*sqrt(norm2_E/sum(trapz(r,abs(E1).^2.*r,2),[1,3,4]));
    end
end    

% Local error estimate
a5 = apply_W(    E1);

if isempty(r)
    err = sum(abs((a4(:)-a5(:))*(dz/10)).^2);
else % use the radially-symmetric scheme
    err = sum(trapz(r,abs((a4-a5)*(dz/10)).^2.*r,2),[1,3,4]);
end
err = sqrt(err/norm2_E);

% Stepsize control
if isnan(norm2_E) % the computation is just so wrong, so we reduce the step size and do it again
    opt_dz = 0.5*dz;
    success = false;
elseif norm2_E == 0 % all-zero field
    opt_dz = 2*dz;
    success = true;
else
    opt_dz = max(0.5,min(2,0.8*(DW_threshold/err)^(1/4)))*dz;

    success = err < DW_threshold;
end

end