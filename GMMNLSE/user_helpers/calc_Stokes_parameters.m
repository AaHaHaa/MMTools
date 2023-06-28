function Stokes_parameters = calc_Stokes_parameters( t,fields,ellipticity_of_birefringence )
%CALC_STOKES_PARAMETERS It calculates the Stokes vector of the fields
%   
%   Input:
%       t: (N,1); time (ps)
%       fields: a (N,2,num_fields) array for (field_e1,field_e2), where
%               (e1,e2) are orthogonal basis for the elliptical polarizations (sqrt(W))
%       ellipticity_of_birefringence: as words suggest (default to 0: linear polarization) and
%           e1 = (x+iry)/sqrt(1+r^2)
%           e2 = (rx-iy)/sqrt(1+r^2), r = ellipticity (Nonlinear Fiber Optics, Agrawal)
%
%   Output:
%       Stokes_parameters: a (1,4,num_fields) array

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
if size(fields,1) == 1
    I = abs(x).^2 + abs(y).^2;
    Q = abs(x).^2 - abs(y).^2;
    U = 2*real(x.*conj(y));
    V = -2*imag(x.*conj(y));
else
    I = trapz(t,abs(x).^2 + abs(y).^2);
    Q = trapz(t,abs(x).^2 - abs(y).^2);
    U = trapz(t,2*real(x.*conj(y)));
    V = trapz(t,-2*imag(x.*conj(y)));
end

Stokes_parameters = cat(2,I,Q,U,V)/1e3; % pJ to nJ

end