function field_with_SPM = add_nonlinear_phase(field,n2,f0,Aeff,L)
%ADD_NONLINEAR_PHASE adds the nonlinear phase to a narrowband field
%
%   field: (Nt,1); the electric field with the unit of sqrt(W)
%   n2: a scalar; nonlinear refractive index (m^2/W)
%   f0: a scalar; pulse center frequency (Hz)
%   Aeff: a scalar; effective mode-field area (m^2)
%   L: a scalar; propagation length (m)

c = 299792458; % m/s
omega0 = 2*pi*f0; % 2*pi*Hz
gamma = n2*omega0/c/Aeff;

nonlinear_phase = gamma*abs(field).^2*L;

field_with_SPM = field.*exp(1i*nonlinear_phase);

end