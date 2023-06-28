function [phi,theta,psi,chi] = calc_ellipticity( fields,ellipticity_of_birefringence )
%CALC_ELLIPTICITY It calculates the information about the polarization of
%the fields
%   
%   Input:
%       fields: a (N,2,num_fields) array for (field_e1,field_e2), where
%               (e1,e2) are orthogonal basis for the elliptical polarizations
%       ellipticity_of_birefringence: as words suggest (default to 0: linear polarization) and
%           e1 = (x+iry)/sqrt(1+r^2)
%           e2 = (rx-iy)/sqrt(1+r^2), r = ellipticity (Nonlinear Fiber Optics, Agrawal)
%
%   Output:
%       - Ex = |E|*cos(theta)*exp(i*phi_x)
%       - Ey = |E|*sin(theta)*exp(i*phi_y)
%
%       phi = phi_y - phi_x: the phase difference (in deg)
%       theta: the rotation of the polarization ellipse (in deg)
%
%       psi:
%       chi: 

if size(fields,2) ~= 2
    error('calc_ellipticity:inputFieldsError',...
        '"fields" should be of the size (N,2).');
end

r = ellipticity_of_birefringence;

e1 = fields(:,1,:);
e2 = fields(:,2,:);

if r == 0 % linear polarization
    x = e1;
    y = e2;
else
    x = (e1 + r*e2)/sqrt(1+r^2);
    y = -1i*(r*e1 - e2)/sqrt(1+r^2);
end

% Stokes parameters: S = [I;Q;U;V]
%I = abs(x).^2 + abs(y).^2;
Q = abs(x).^2 - abs(y).^2;
U = 2*real(x.*conj(y));
V = -2*imag(x.*conj(y));
psi = atan2(U,Q)/2;
chi = atan2(V,sqrt(Q.^2+U.^2))/2;

% phase difference ( restrict it in [0,2*pi) )
phi = mod(angle(x) - angle(y),2*pi);
% orientation of the polarization ellipse
abs_x = abs(x);
abs_y = abs(y);
theta = atan2(abs_y,abs_x);

% transform into degree
phi = phi*180/pi;
theta = theta*180/pi;
psi = psi*180/pi;
chi = chi*180/pi;

end