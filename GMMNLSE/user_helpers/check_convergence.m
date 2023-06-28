function [converged_yes,fig] = check_convergence( energy,field,freq,t,varargin)
%CHECK_CONVERGENCE Check the convergence of energies
%
%   energy - a (1,rt_max)-array of energies
%   field - the field for the last entry in "energy" array
%   freq - an array of frequency of the field (THz)
%   t - an array of time of the field (ps)
%   fig - figure handles for "energy" and "pulse, spectrum" plots
%   tol_convergence - tolerance for the convergence; (default: 1e-5)
%   plot_yes - true(1) or false(0); whether to plot the field or not (default: true)
%   reset_idx

converged_yes = false;

least_rt = 5; % the least number of roundtrips necessary to taken. It should be at least larger than 3.

%% Default optional input arguments
% Accept only 4 optional inputs at most
numvarargs = length(varargin);
if numvarargs > 4
    error('check_convergence:TooManyInputs', ...
        'It takes only at most 4 optional inputs');
end

% Set defaults for optional inputs
tol_convergence = 1e-5;
plot_yes = true;
reset_idx = 1;
optargs = {tol_convergence,plot_yes,reset_idx};

% Now put these defaults into the valuesToUse cell array, 
% and overwrite the ones specified in varargin.
optargs(1:numvarargs) = varargin;

% Place optional args in memorable variable names
[tol_convergence,plot_yes,reset_idx] = optargs{:};

diff_num = 3; % the size of the considered difference elements
rt_num0 = find(energy(1:end)~=0,1,'last');
rt_num = rt_num0-reset_idx +1;
energy0 = energy(1:rt_num0); % save the original energy vector
energy = energy(reset_idx:rt_num0);

%% Decide whether to stop or rerun the iteration
if rt_num > least_rt
    difference = relative_diff(energy);
    
    if rt_num <= diff_num+2
        considered_difference = difference(1:rt_num-2);
    else
        considered_difference = difference( (rt_num-diff_num+1)-2:rt_num-2 );
    end
    
    % Determine if it's going to stop.
    if ~any(abs(considered_difference) > tol_convergence)
        % Converged! Iteration finished!
        converged_yes = true;
    end
end

%% Parameters for plotting
if plot_yes
    field = field(:,:,end);

    N = size(field,1);
    dt = t(2)-t(1); % ps
    factor_correct_unit = (N*dt)^2/1e3; % to make the spectrum of the correct unit "nJ/THz"
                                        % "/1e3" is to make pJ into nJ
    I_freq = abs(fftshift(ifft(field),1)).^2*factor_correct_unit; % in frequency domain

    fig(1) = figure('Name','Output energy');
    h = plot(energy0);
    ylabel('Energy (nJ)')
    xlabel('Round trip')
    set(h,'linewidth',2); set(gca,'fontsize',14);

    % Spectrum
    fig(2) = figure('Name','Output field');
    subplot(1,2,1);
    h1 = plot(freq,I_freq);
    xlim([min(freq) max(freq)]);
    xlabel('Frequency (THz)');
    ylabel('Intensity (nJ/THz)');
    grid on;
    set(h1,'linewidth',2); set(gca,'fontsize',14);

    % Pulse
    subplot(1,2,2);
    h2 = plot(t,abs(field).^2);
    xlim([min(t),max(t)]);
    xlabel('Time (ps)');
    ylabel('Intensity (W)');
    grid on;
    set(h2,'linewidth',2); set(gca,'fontsize',14);

    drawnow;
end

end

%% relative_diff
function avg_diff = relative_diff(array)
%DIFF_ARRAY calculate the relative slope from central difference and
% divide it by the center value.

% difference
diff_array = diff(array);
diff_array = diff_array(2:end)./array(2:end-1);
% central difference
central_diff_array = (array(3:end)-array(1:end-2))/2./array(2:end-1);

% average difference
avg_diff = (diff_array + central_diff_array)/2;

end