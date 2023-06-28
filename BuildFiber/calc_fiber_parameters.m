clearvars; close all;

%% Set parameters (users modify only this block)
modes_used = 1:32; % In general the modes do not have to be consecutive
Nx = 800; % number of spatial grid points for each mode
gpu_yes = true; % true = run on GPU, false = run on CPU
folder_name = 'GRIN-62.5_245_wavelength1030nm'; % folder containing the calculated modes

% File name parameters:
Nf = 100;
lambda0 = 1030e-9; % center wavelength in nm; an "integer"
freq_range = 100; % THz; frequency range. Should be the same as the range used to calculate the modes. Usually 100 THz.
bandwidth = 'broadband'; % 'narrowband' or 'broadband'

polynomial_fit_order = 7;
num_disp_orders = 5; % i.e. if this is 3, four coefficients will be calculated, including the 0th order

%%
calc_dispersion(modes_used,folder_name,Nf,lambda0,freq_range,bandwidth,polynomial_fit_order,num_disp_orders);
calc_SRSK_tensors(modes_used,Nx,gpu_yes,folder_name,lambda0);
delete_modeprofiles(modes_used,folder_name,Nf,lambda0,freq_range)