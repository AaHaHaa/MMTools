function output = redefine_field_time_window( input,initial_center_wavelength,...
                                                     final_center_wavelength,final_dt,final_Nt)
%REDEFINE_FIELD_TIME_WINDOW It re-computes the fields based on the time window.
%
% Note that this works only for single-peak/narrowband field in the 
% frequency domain. If there are multiple wavelengths, the interpolation of
% the phase at different time points can be totally wrong.
%
% Input arguments:
%   input - a structure with the fieldnames, "dt" and "fields"
%   initial_center_wavelength - center wavelength for the input field (m)
%   final_center_wavelength - center wavelength for the final/targeted field (m)
%   final_dt - the final dt (ps)
%   final_Nt - the final Nt, the number of points
%
% Output arguments:
%   output - a structure with the fieldnames, final "dt" and "fields"
%            besides other fieldnames inherited from the input

initial_Nt = size(input.fields,1);
initial_dt = input.dt;

abs_fields = abs(input.fields);
phase_fields = unwrap(angle(input.fields));

initial_t = (-floor(initial_Nt/2):floor((initial_Nt-1)/2))'*initial_dt;
final_t = (-floor(final_Nt/2):floor((final_Nt-1)/2))'*final_dt;

abs_fields = interp1(initial_t,abs_fields,final_t,'pchip',0);
phase_fields = interp1(initial_t,phase_fields,final_t,'pchip',0);

c = 299792458*1e-12; % m/ps
output_field = abs_fields.*exp(1i*(phase_fields+2*pi*(c./final_center_wavelength-c./initial_center_wavelength).*final_t));

output = input;
output.fields = output_field;
output.dt = final_dt;

end