%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script uses the propagation constants obtained from the mode
% calculations to plot the dispersion parameters of each mode
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('helpers/');

%% Set parameters
clearvars; %close all;

Nf = 10;
lambda0 = 1030e-9; % m, the center wavelength used to solve for the modes
freq_range = 100; % THz; frequency range. Should be the same as the range used to calculate the modes
modes_list = 1; % usually just 1:num_modes, but it can be different

polynomial_fit_order = 7;
num_disp_orders = 3; % i.e. if this is 3, 4 coefficients will be calculated, including the 0th order

diameter = 5.8; % Only needed for the file name
radius = diameter/2;
saved_folder = '1060XP_wavelength1030nm';

folder_name = ['Fibers/' saved_folder];

%% Load in the calculated effective n values
if freq_range == 0
    Nf = 1;
end

% Set the range in frequency space, which is more objective
c = 2.99792458e-4; % speed of ligth m/ps
if freq_range == 0
    error('Cannot calculate dispersion with only one frequency point');
else
    freq1 = c/lambda0 - freq_range/2;
    freq2 = c/lambda0 + freq_range/2;

    f = linspace(freq1,freq2,Nf)'; % THz
    wavelength = c./f*1e6; % um
end

num_modes = length(modes_list);

n_calc = zeros(Nf, num_modes);
for kk = 1:Nf
    lambda_kk = wavelength(kk);
    for ii = 1:num_modes
        fname = [folder_name '/mode' num2str(modes_list(ii)) 'wavelength' num2str(round(lambda_kk*10000))];
        load([fname '.mat'])
        n_calc(kk, ii) = neff;
    end
    fprintf('Loading lambda = %d um\n', round(lambda_kk*1000));
end

%% Calculate the propagation constants
w = 2*pi*f; % angular frequencies in 1/ps
beta = n_calc.*w/c; % beta in 1/m

%% Calculate the propagation constants
f_calc = c./wavelength*1e6; % THz; loaded from the file

f = linspace(f_calc(end),f_calc(1),Nf)'; % THz; resample for applying Taylor series expansion
abs_beta = interp1(f_calc,abs(beta),f,'pchip');
ang_beta = interp1(f_calc,unwrap(angle(beta),[],1),f,'pchip');
beta = abs_beta.*exp(1i*ang_beta);

omega = 2*pi*f; % angular frequencies in 1/ps
df = f(2)-f(1);
domega = 2*pi*df;
beta = real(beta); % beta in 1/m

%% Calculate the derivatives

% We need to use cell arrays because the higher orders are calculated from
% finite differences. This means that each order has one less data point
% than the previous one.
omega_full = cell(num_disp_orders+1, 1); % 2*pi*THz
wavelength_full = cell(num_disp_orders+1, 1); % um

central_diff_coeff = {[   0,    0, -1/2,   0, 1/2,   0,   0],...
                      [   0,    0,    1,  -2,   1,   0,   0],...
                      [   0, -1/2,    1,   0,  -1, 1/2,   0],...
                      [   0,    1,   -4,   6,  -4,   1,   0],...
                      [-1/2,    2, -5/2,   0, 5/2,  -2, 1/2],...
                      [   1,   -6,   15, -20,  15,  -6,   1]};

omega_full{1} = omega;
wavelength_full{1} = 2*pi*c./omega_full{1}*1e6;
for disp_order = 1:num_disp_orders
    switch disp_order
        case {1,2}
            omega_full{disp_order+1} = omega(2:end-1); % 2*pi*THz
        case {3,4}
            omega_full{disp_order+1} = omega(3:end-2);
        case {5,6}
            omega_full{disp_order+1} = omega(4:end-3);
    end
        
    wavelength_full{disp_order+1} = 2*pi*c./omega_full{disp_order+1}*1e6; % in um
end

% beta_full will have all of the orders, for each mode, as a function of wavlength
beta_full = cell(num_disp_orders+1, 1);
beta_full{1} = beta;
for disp_order = 1:num_disp_orders
    switch disp_order
        case {1,2}
            beta_full{disp_order+1} = central_diff_coeff{disp_order}(3)*beta(1:end-2,:) + ...
                                      central_diff_coeff{disp_order}(4)*beta(2:end-1,:) + ...
                                      central_diff_coeff{disp_order}(5)*beta(3:end,:);
        case {3,4}
            beta_full{disp_order+1} = central_diff_coeff{disp_order}(2)*beta(1:end-4,:) + ...
                                      central_diff_coeff{disp_order}(3)*beta(2:end-3,:) + ...
                                      central_diff_coeff{disp_order}(4)*beta(3:end-2,:) + ...
                                      central_diff_coeff{disp_order}(5)*beta(4:end-1,:) + ...
                                      central_diff_coeff{disp_order}(6)*beta(5:end,:);
        case {5,6}
            beta_full{disp_order+1} = central_diff_coeff{disp_order}(1)*beta(1:end-6,:) + ...
                                      central_diff_coeff{disp_order}(2)*beta(2:end-5,:) + ...
                                      central_diff_coeff{disp_order}(3)*beta(3:end-4,:) + ...
                                      central_diff_coeff{disp_order}(4)*beta(4:end-3,:) + ...
                                      central_diff_coeff{disp_order}(5)*beta(5:end-2,:) + ...
                                      central_diff_coeff{disp_order}(6)*beta(6:end-1,:) + ...
                                      central_diff_coeff{disp_order}(7)*beta(7:end,:);
    end
    beta_full{disp_order+1} = beta_full{disp_order+1}/domega^disp_order;
end

% ps^n/m to fs^n/mm
for disp_order = 1:num_disp_orders+1
    beta_full{disp_order} = beta_full{disp_order}*10^(3*disp_order-6);
end

%% Display the results
wavelength_min = c/freq2*1e6; % um
wavelength_max = c/freq1*1e6; % um

coo = distinguishable_colors(num_modes);

ylabels = cell(num_disp_orders+1, 1);
ylabels{1} = '1/mm';
ylabels{2} = 'fs/mm';
for disp_order = 2:num_disp_orders
    ylabels{disp_order+1} = ['fs^' num2str(disp_order) '/mm'];
end

% Plot the results
figure;
lowest_order = 1;
for disp_order = lowest_order+1:num_disp_orders+1
    subplot(1,num_disp_orders+1-lowest_order,disp_order-lowest_order)
    hold on
    for midx = 1:num_modes
        plot(wavelength_full{disp_order}, beta_full{disp_order}(:, midx), 'Color', coo(midx,:), 'linewidth', 2);
        xlim([wavelength_min,wavelength_max]);
    end
    hold off
    set(gca,'fontsize',20);
    ylabel(ylabels{disp_order});
    xlabel('\mum')
    title(['\beta_' num2str(disp_order-1)]);
end