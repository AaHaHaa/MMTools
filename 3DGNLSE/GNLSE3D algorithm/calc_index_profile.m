function [index_profile,dx] = calc_index_profile(fiber,Nx,spatial_window,f)
%CALC_INDEX_PROFILE It calculates the index profile of the fiber.
%
%   fiber: a structure with
%
%       Below you can use the fiber in the repository (in "fiber_collections.m") 
%       or set all parameters yourself
%       (1) Set it yourself
%           fiber_type - a string with either "step', 'Nufern', or 'GRIN'
%           material - fiber material; usually 'silica'
%           diameter - fiber core diameter (um)
%           NA - numerical aperture
%           extra_params - extra parameters for different fiber types
%
%       (2) Use fiber repository
%           fiber - the name of the fiber
%                   Currently I have '1060XP',
%                                    'YB1200-4_125',
%                                    'YB1200-6_125DC', 'YB1200-6_125DC-PM',
%                                    'YB1200-10_125DC', 'YB1200-10_125DC-PM',
%                                    'YB1200-20_400DC',
%                                    'YB1200-25_250DC', 'YB1200-25_250DC-PM',
%                                    'ER30-4_125', 'ER110-4_125',
%                                    'ER16-8_125', 'ER80-8_125',
%                                    'M5-980-125', 'M12-980-125',
%                                    'OM1','OM2','OM3','OM4',
%                                    'FUD-7005',
%                                    'PLMA-YDF-30_400-VIII'.
%
% -------------------------------------------------------------------------
% (1) step fiber:
%     extra_params = [];
% (2) GRIN fiber:
%     extra_params.alpha = 2.08; % Shape parameter
% -------------------------------------------------------------------------
%          
%   Nx - number of spatial grid points
%   spatial_window - full spatial window size (um) (usually set to 100 um)
%   f - frequency points (from small to large) (THz)

if ~isfield(fiber,'alpha')
    fiber.alpha = 2;
end
if isfield(fiber,'fiber')
    [fiber.fiber_type,fiber.core_diameter,fiber.clad_diameter,fiber.core_NA,fiber.clad_NA,fiber.alpha] = fiber_collections(fiber.fiber);
    fiber.material = 'fused silica';
end

% Sellmeier coefficients
[a,b] = Sellmeier_coefficients(fiber.material);

%% Calculate the modes
core_radius = fiber.core_diameter/2; % core radius of fiber, in um
clad_radius = fiber.clad_diameter/2; % cladding radius, in um

% Set the range in frequency space, which is more objective
c = 299.792458; % um/ps
wavelength = c./f; % um

% Calculate the index difference using the Sellmeier equation to generate n(lambda)
Sellmeier_terms = @(lambda,a,b) a.*lambda.^2./(lambda.^2 - b.^2);
n_from_Sellmeier = @(lambda) sqrt(1+sum(Sellmeier_terms(lambda,a,b),2));

n_core = n_from_Sellmeier(wavelength);
n_clad = sqrt(n_core.^2-fiber.core_NA^2);
n_coat = sqrt(n_clad.^2-fiber.clad_NA^2);

%% Calculate at the max frequency first to determine the number of modes
n = {n_core,n_clad,n_coat};
[index_profile,...
 ~, dx] = run_profile_function(fiber.fiber_type,...
                               Nx, spatial_window,...
                               core_radius, clad_radius,...
                               n, fiber.alpha);

end

%% helper function
function [epsilon, x, dx] = run_profile_function(fiber_type, Nx, spatial_window, core_radius, clad_radius, n, alpha)

profile_function = str2func(sprintf('build_%s',fiber_type)); % function that builds the fiber
switch fiber_type
    case 'step'
        [epsilon, x, dx] = profile_function(Nx, spatial_window, core_radius, clad_radius, n);
    case 'GRIN'
        [epsilon, x, dx] = profile_function(Nx, spatial_window, core_radius, clad_radius, n, alpha);
end

end