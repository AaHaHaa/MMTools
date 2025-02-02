function A_H = FHATHA(A,...
                      r_max,...
                      r,kr,...
                      dr,dkr,...
                      l0,exp_prefactor,r2_prefactor,...
                      ifftQ,...
                      varargin)
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
% In radially-symmetic 3D-UPPE code, the radial dimension is in the second 
% dimension, so the Hankel transform is applied only to the second 
% dimension. The field should be up to four dimensions (Nt,Nr,Nz,Np) in 
% radially-symmetric 3D-UPPE.
%
% Different from the refernce Magni's paper, I modified such that it
% considers A(r=0) and A_H(kr=0). This improves the FHATHA's performance
% significantly. Otherwise, FHATHA's failing to satisfy discrete Parseval's
% theorem creates significant deviation after many iterations of (Hankel
% and then inverse Hankel transforms) sequence.
%
% Input arguments:
%   A: input signal; (Nt,Nr,Nz,Np)
%   r_max: half the size of the spatial window = maximum radius
%
% [Below should be the outputs from the pre-computed Hankel_info();
%  * represents parameters used only for conserving energy before and after the Hankel transform]
%   r*: sampling radius (m); (1,Nr)
%   kr: sampling k-vector radius (2*pi/m); (1,Nr)
%   dr*: samppling radius spacing (m); (1,Nr);
%        dr includes the effect from r=0, so it is computed from
%            dr = diff([0,r]);
%   dkr*: samppling k-vector radius (2*pi/m); (1,Nr);
%         dkr includes the effect from kr=0, so it is computed from
%            dkr = diff([0,kr]);
%   l0: prefactor used in computing R in FHATHA prepared for cross correlation
%   exp_prefactor: exponential prefactor used in computing R in FHATHA prepared for cross correlation
%   r2_prefactor: (=r^2/2) prefactor used in computing the Hankel-transformed signal at kr=0
%   ifftQ: the inverse-Fourier-transformed Q in FHATHA prepared for cross correlation
%   energy_conserved_yes: whether to ensure the energy conservation or not (default: false)
%                         It's default to "false" because FHATHA is more accurate with energy-conservation operation through simple scaling.
%
% Output arguments:
%   A_H: Hankel-transformed signal; (Nt,Nr,Nz,Np)

%% Default optional input arguments
% Accept only 1 optional inputs at most
numvarargs = length(varargin);
if numvarargs > 1
    error('FHATHA:TooManyInputs', ...
          'It takes only at most 1 optional input.');
end

% Set defaults for optional inputs
energy_conserved_yes = false;
optargs = {energy_conserved_yes};

% Now put these defaults into the valuesToUse cell array, 
% and overwrite the ones specified in varargin.
optargs(1:numvarargs) = varargin;

% Place optional args in memorable variable names
energy_conserved_yes = optargs{:};

%% main() of FHATHA
if any(A(:))
    sA = size(A,1:4);
    Nr = sA(2)-1; % the number of radial sampling points excluding r=0

    if isequal(class(A),'gpuArray')
        A_middle = cat(2,A(:,2:end,:,:),zeros([sA(1),1,sA(3),sA(4)],'gpuArray')); % zero-padding
    else
        A_middle = cat(2,A(:,2:end,:,:),zeros([sA(1),1,sA(3),sA(4)])); % zero-padding
    end
    A_middle(:,1,:,:) = (A(:,1,:,:) + A(:,2,:,:))/2                *0.5 + ...
                        (l0*(A(:,2,:,:) - A(:,3,:,:)) + A(:,3,:,:))*0.5;
    A_factor = A_middle(:,1:Nr,:,:) - A_middle(:,2:Nr+1,:,:);

    % Compute R for FHATHA
    R = A_factor.*exp_prefactor;
    if isequal(class(R),'gpuArray')
        R = cat(2,R,zeros(size(R),'gpuArray')); % zero-padding
    else
        R = cat(2,R,zeros(size(R))); % zero-padding
    end
    
    % Compute the Hankel transform with FHATHA
    A_H = fft(fft(R,[],2).*ifftQ,[],2);
    A_H = r_max./kr(2:end).*A_H(:,1:Nr,:,:);

    A_H = cat(2,sum(A_factor.*(r_max^2*r2_prefactor),2),A_H);
    
    % Calibrate A_H by fixing its energy. Energy should be conserved.
    if energy_conserved_yes
        E0 = sum(mytrapz(dr,abs(A).^2.*r,Nr+1),[1,3,4]);
        E_H = sum(mytrapz(dkr,abs(A_H).^2.*kr,Nr+1),[1,3,4]);
        %idx = abs(A_H).^2 > max(abs(A_H).^2,[],2)*1e-2; % scale only for strong signal
        %A_H(idx) = A_H(idx).*sqrt(E0./E_H);
        A_H = A_H.*sqrt(E0./E_H);
    end
else
    if isequal(class(A),'gpuArray')
        A_H = zeros(size(A),'gpuArray');
    else
        A_H = zeros(size(A));
    end
end

end

function z = mytrapz(dx,y,Nr)

z = sum(dx.*(y(:,1:Nr-1,:,:) + y(:,2:Nr,:,:))/2,2);

end