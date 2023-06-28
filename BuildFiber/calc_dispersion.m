%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script uses the propagation constants obtained from the mode
% calculations to approximate the dispersion parameters of each mode
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function calc_dispersion(modes_used,folder_name,Nf,lambda0,freq_range,bandwidth,polynomial_fit_order,num_disp_orders)

    % File name parameters:
    dir_prefix = ['Fibers/', folder_name]; % folder containing the calculated modes

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
        lambda = c./f*1e6; % um
    end

    num_modes = length(modes_used);

    n_calc = zeros(Nf, num_modes);
    for kk = 1:Nf
        lambda_kk = lambda(kk);
        for ii = 1:num_modes
            fname = [dir_prefix '/mode' num2str(modes_used(ii)) 'wavelength' num2str(round(lambda_kk*10000))];
            load([fname '.mat'])
            n_calc(kk, ii) = neff;
        end
        fprintf('Loading lambda = %d um\n', round(lambda_kk*1000));
    end

    %% Calculate the propagation constants
    w = 2*pi*f; % angular frequencies in 1/ps
    beta_f = n_calc.*w/c; % beta in 1/m

    %% Fit the propagation constants to a polynomial and save the appropriate derivatives
    switch bandwidth
        case 'narrowband'
            w_dispersion = 2*pi*c/lambda0; % angular frequency at which dispersion is calculated, in 1/ps

            betas = zeros(num_disp_orders+1,num_modes); % the dispersion coefficients
            for midx = 1:num_modes
                beta_f_i = beta_f(:, midx);

                [beta_fit_last,~,mu] = polyfit(w, beta_f_i, polynomial_fit_order); % the fit coefficients
                betas(1,midx) = polyval(beta_fit_last, w_dispersion,[],mu)/1000; % Save beta_0 in 1/mm
                for disp_order = 1:num_disp_orders
                    % The derivatives can be calculated exactly from the coefficients
                    beta_fit_last = ((polynomial_fit_order-(disp_order-1)):-1:1)/mu(2).*beta_fit_last(1:(polynomial_fit_order-(disp_order-1)));
                    betas(disp_order+1,midx) = polyval(beta_fit_last, w_dispersion,[],mu)*(1000)^disp_order/1000; % beta_n in fs^n/mm
                end
            end
            save([dir_prefix '/betas'], 'betas');
        case 'broadband'
            betas = beta_f;
            save([dir_prefix '/betas'], 'betas','f');
    end
end