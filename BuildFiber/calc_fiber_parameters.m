clearvars; close all;

addpath('helpers/');

%% Set parameters (users modify only this block)
modes_used = 1; % an array; in general the modes do not have to be consecutive
Nx = 400; % number of spatial grid points for each mode
gpu_yes = false; % true = run on GPU, false = run on CPU
folder_name = '1060XP_wavelength1060nm'; % folder containing the calculated modes

% File name parameters:
Nf = 10;
wavelength0 = 1060e-9; % center wavelength in nm; an "integer"
freq_range = 50; % THz; frequency range. Should be the same as the range used to calculate the modes. Usually 100 THz.
bandwidth = 'narrowband'; % 'narrowband' or 'broadband'

polynomial_fit_order = 7;
num_disp_orders = 5; % i.e. if this is 3, four coefficients will be calculated, including the 0th order

%%
calc_dispersion(modes_used,folder_name,Nf,wavelength0,freq_range,bandwidth,polynomial_fit_order,num_disp_orders);
calc_SRSK_tensors(modes_used,Nx,gpu_yes,folder_name,wavelength0);
delete_modeprofiles(modes_used,folder_name,Nf,wavelength0,freq_range)