function A_H = FHATHA(A,...
                      r_max,kr,...
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

sA = size(A);
Nr = sA(2);
A = cat(2,A,zeros(sA));

R = (A(:,1:Nr,:,:) - A(:,2:Nr+1,:,:)).*exp_prefactor;
R(:,1,:,:) = R(:,1,:,:)*l0;
R = cat(2,R,zeros(size(R)));

A_H = fft(fft(R,[],2).*ifft(Q,[],2),[],2);
A_H = r_max./kr.*A_H(:,1:Nr,:,:);

end