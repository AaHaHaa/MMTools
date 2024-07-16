function output = saturable_absorber_action_NPE(input, M1, M2, P1, P2)
% SATURABLE_ABSORBER_ACTION_SIMPLE Apply an nonlinear loss function as a saturable absorber effect to each mode individually
% input - either a struct or a matrix:
%   input.fields - a (N, num_modes, m) matrix with the fields for each mode, in the time domain
%
% References:
% 1. Li et al., "Geometrical description of the onset of multi-pulsing in 
%    mode-locked laser cavities," J. Opt. Soc. Am. B, 27(10), 2068-2077 
%    (2010)
% 2. Wei et al., "General description and understanding of the nonlinear
%    dynamics of mode-locked fiber lasers," Sci. Rep. 7(1), 1292 (2017)

if isstruct(input)
    input_fields = input.fields(:, :, end);
else
    error('saturable_sbsorber_action_NPE:inputError',...
          '"input" must be a structure with "dt" and "fields".');
end

powers = abs(input_fields).^2;
output_fields = input_fields.*sqrt(M1+M2*(1-cos(2*pi*(powers-P1)/P2)));

output = struct('fields', output_fields);
if isfield(input,'dt')
    output.dt = input.dt;
end

end