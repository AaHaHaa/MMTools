function A_H = FHATHA(A0,...
                      r_max,...
                      r,kr,...
                      dr,dkr,...
                      l0,exp_prefactor,...
                      Q)
%FHATHA fast Hankel transform of high accuracy
%   This numerical scheme for Hankel transform follows
%
%   Magni et al., "High-accuracy fast Hankel transform for optical beam
%   propagation," J. Opt. Soc. Am. A 9(11), 2031-2033 (1992)
%
% To use this function, please call Hankel_info() first, which generates
% the required parameters to apply FHATHA.
%
% For inverse Hankel transform, simply supply k_max and r in the input
% arguments r_max and k, respectively. l0, exp_prefactor, and Q remains the
% same whether it's Hankel or inverse Hankel transform.
%
% In 2D-UPPE code, the radial dimension is in the second dimension, so the
% Hankel transform is applied only to the second dimension. The field
% should be up to four dimensions (Nt,Nr,Nz,Np) in 2D-UPPE.

if any(A0(:))
    sA = size(A0);
    Nr = sA(2); % the number of radial sampling points
    if isequal(class(A0),'gpuArray')
        A = cat(2,A0,zeros(sA,'gpuArray')); % zero-padding
    else
        A = cat(2,A0,zeros(sA)); % zero-padding
    end
    % Compute R for FHATHA
    R = (A(:,1:Nr,:,:) - A(:,2:Nr+1,:,:)).*exp_prefactor;
    R(:,1,:,:) = R(:,1,:,:)*l0;
    if isequal(class(R),'gpuArray')
        R = cat(2,R,zeros(size(R),'gpuArray')); % zero-padding
    else
        R = cat(2,R,zeros(size(R))); % zero-padding
    end
    
    % Compute the Hankel transform with FHATHA
    A_H = fft(fft(R,[],2).*ifft(Q,[],2),[],2);
    A_H = r_max./kr.*A_H(:,1:Nr,:,:);
    
    % Calibrate A_H by fixing its energy. Energy should be conserved.
    %E0 = repmat(trapz(r,abs(A0).^2.*r,2),1,Nr,1,1);
    %E_H = repmat(trapz(kr,abs(A_H).^2.*kr,2),1,Nr,1,1);
    %idx = E0 > max(E0(:))/1e5;
    %A_H(idx) = A_H(idx).*sqrt(E0(idx)./E_H(idx));
    E0 = sum(mytrapz(dr,abs(A0).^2.*r,Nr),[1,3,4]);
    E_H = sum(mytrapz(dkr,abs(A_H).^2.*kr,Nr),[1,3,4]);
    A_H = A_H.*sqrt(E0./E_H);
else
    if isequal(class(A0),'gpuArray')
        A_H = zeros(size(A0),'gpuArray');
    else
        A_H = zeros(size(A0));
    end
end

end

function z = mytrapz(dx,y,Nr)

z = sum(dx.*(y(:,1:Nr-1,:,:) + y(:,2:Nr,:,:))/2,2);

end