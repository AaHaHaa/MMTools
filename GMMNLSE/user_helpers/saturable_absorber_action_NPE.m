function output = saturable_absorber_action_NPE(input, M1, M2, P1, P2)
% SATURABLE_ABSORBER_ACTION_SIMPLE Apply an nonlinear loss function as a saturable absorber effect to each mode individually
% input - either a struct or a matrix:
%   input.fields - a (N, num_modes, m) matrix with the fields for each mode, in the time domain
%  OR
%   input - a (N, num_mdoes, m) matrix with the fields for each mode, in the time domain
%
% References:
% 1. Li et al., "Geometrical description of the onset of multi-pulsing in 
%    mode-locked laser cavities," J. Opt. Soc. Am. B, 27(10), 2068-2077 
%    (2010)
% 2. Wei et al., "General description and understanding of the nonlinear
%    dynamics of mode-locked fiber lasers," Sci. Rep. 7(1), 1292 (2017)

if isstruct(input)
    input_field = input.fields(:, :, end);
else
    input_field = input(:, :, end);
end

powers = abs(input_field).^2;
output_field = input_field.*sqrt(M1+M2*(1-cos(2*pi*(powers-P1)/P2)));

if isstruct(input)
    output = input;
    output.fields = output_field;
else
    output = output_field;
end