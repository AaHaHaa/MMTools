function V = calc_circular_convergence( t,fields,ellipticity_of_birefringence )
%CALC_CIRCULAR_CONVERGENCE It calculates the degree of convergence from the
%3rd component of the Stokes vector
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
%       V: a (num_fields,1) array; the integrated abs(S3),
%            where S3 is the 3rd component of the Stokes vector

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
I = trapz(t,abs(x).^2 + abs(y).^2);
V = squeeze(trapz(t,abs(-2*imag(x.*conj(y))))./I);

end