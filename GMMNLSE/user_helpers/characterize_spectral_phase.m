function [quardratic_phase,cubic_phase,fitted_param,quintic_phase,pulse_inst_freq] = characterize_spectral_phase( f,At,fitted_order,verbose )
%CHARATERIZE_SPECTRAL_PHASE It fits a polynomial to the spectral phase and
%finds the GVD and TOD.
%
% Input:
%   f  - (Nt,1); frequency (THz); from small to large, usually generated with
%
%        f = sim.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
%
%   At - (Nt,1); the time-domain electric field
%   fitted_order - a scalar; the order of the fitted polynomial (default: 7)
%   verbose - 1(true) or 0(false); plot the fitted curve and print the
%              results (quardratic_phase and cubic_phase) (default: false)
%
% Output:
%   quardratic_phase
%   cubic_phase
%   fitted_param - {p,s,mu}, the fitted parameters for the phase
%   quintic_phase
%   pulse_inst_freq - (Nt,1); pulse's instantaneous frequency in time
%
%   Please check "p.248, Ch.6.2.1, Applications of Nonlinear Fiber Optics"
%   for details of the convention of phase.
%       phi = phi0 + phi1*(w-w0) + 1/2*phi2*(w-w0)^2 + 1/6*phi3*(w-w0)^3 + ......

switch nargin
    case 2
        fitted_order = 7;
        verbose = false;
    case 3
        verbose = false;
end

Nt = size(At,1);

if numel(At) ~= Nt || size(At,2) ~= 1
    error('characterize_spectral_phase:EfSizeError',...
          'At can only be (Nt,1) column vector.');
end

f0 = feval(@(x)x(1),ifftshift(f,1)); % center frequency; THz
frequency_window = max(f) - min(f); % THz
dt = 1/frequency_window; % ps
t = (-floor(Nt/2):floor((Nt-1)/2))'*dt; % ps

%% First remove the spectral phase due to temporal offset for correct phase unwrapping
[~,peak_position] = max(sum(abs(At).^2,2),[],1);
avg_center = sum((1:Nt)'.*abs(At).^2,[1,2])/sum(abs(At(:)).^2);
% If two positions are too far away, then it might indicate that the pulse
% is center at the left edge of the time window such that its other half is
% at the right edge of the window, due to periodic boundary condition under
% numerical DFT.
% If two positions are close, then we pick the avg_center as the pulse
% center.
if abs(peak_position - avg_center) > Nt/4
    fftshift_At = fftshift(At,1);
    avg_center = sum((1:Nt)'.*abs(fftshift_At).^2,[1,2])/sum(abs(fftshift_At(:)).^2);
    if avg_center >= floor(Nt/2)+1
        avg_center = avg_center - floor(Nt/2);
    else
        avg_center = avg_center + floor(Nt/2);
    end
end
pulse_center = avg_center;

phase_shift = ifftshift(2*pi*(1:Nt)'/Nt*pulse_center,1); % phase shift due to temporal offset

Aw_noTOffset = fftshift(ifft(At,[],1).*exp(-1i*phase_shift),1);

% At centered at the time window
At = fftshift(fft(ifftshift(Aw_noTOffset,1),[],1),1); % for temporal computation later to show the instantaneous frequency
At_phase = unwrap(angle(At));

% instantaneous frequency
pulse_inst_freq = -(At_phase(3:end)-At_phase(1:end-2))/(2*dt)/(2*pi)+f0;
pulse_inst_freq = interp1(t(2:end-1),pulse_inst_freq,t,'linear');

%% Consider only the central (non-zero) part of the spectrum
% Remove the weak spectral signal for correct bandwidth computation with the RMS scheme
spectrum = abs(Aw_noTOffset).^2;
threshold_factor = 20;
spectrum(spectrum<max(spectrum)/threshold_factor) = 0;

% RMS computation
[bandwidth,center_f] = calc_RMS(f,spectrum); bandwidth = bandwidth*2;
[~,center_idx] = min(abs(f-center_f));
% Find the left and right edges of the spectrum
bandwidth_idx = ceil(bandwidth/((f(end)-f(1))/(length(f)-1)));
span_idx = ceil(bandwidth_idx*1.5);
left_idx = max(1,center_idx - span_idx);
right_idx = min(length(f),center_idx + span_idx);

f = f(left_idx:right_idx);
omega = 2*pi*f;
Aw_noTOffset = Aw_noTOffset(left_idx:right_idx);

% Phase calculation region should be small, only around the center frequency
left_bandwidth_idx = max(1,ceil(length(f)/2) - floor(bandwidth_idx*0.7));
right_bandwidth_idx = min(length(f),ceil(length(f)/2) + floor(bandwidth_idx*0.7));
phase_calculation_region = left_bandwidth_idx:right_bandwidth_idx;

%% Phase fitting with a polynomial
% Fit the spectral phase with a polynomial of order "fitted_order"
% First subtract the zeroth and first order variation so that the fitting
% around the center frequency can be more accurate. This preprocessing is
% important when the phase is quite messy, such as those affected by Raman.
Aw_phase = unwrap(angle(Aw_noTOffset));
[p,s,mu] = polyfit(omega(phase_calculation_region),Aw_phase(phase_calculation_region),2);
fitted_phase = polyval(p(end-1:end),omega,s,mu);
Aw_phase = Aw_phase - fitted_phase;
[p,s,mu] = polyfit(omega(phase_calculation_region),Aw_phase(phase_calculation_region),fitted_order);
[fitted_phase,delta] = polyval(p,omega,s,mu);

fitted_param = {p,s,mu};

quardratic_phase = p(fitted_order-1)/mu(2)^2*2*1e6; % ps^2 to fs^2
cubic_phase = p(fitted_order-2)/mu(2)^3*6*1e9; % ps^3 to fs^3
if fitted_order >= 4
    quintic_phase = p(fitted_order-3)/mu(2)^4*24*1e12; % ps^4 to fs^4
else
    quintic_phase = 0;
end

%% Plot
if verbose
    spectrum = abs(Aw_noTOffset).^2;
    spectrum = spectrum/max(spectrum);
    % Don't show phases when the spectrum is too weak
    weak_idx = spectrum < 1/threshold_factor/5; % avoid strong clipping in visualization (by adding an extra factor of 5) when spectrum contains weak parts
    Aw_phase(weak_idx) = NaN;
    fitted_phase(weak_idx) = NaN;

    % Plot spectrum and phases
    figure;
    yyaxis left
    hI = plot(f,spectrum,'b');
    ylim([0,max(spectrum)*1.5]);
    ylabel('PSD (norm.)');
    set(gca,'YColor','b');
    yyaxis right
    hp = plot(f,Aw_phase,'k');
    hold on;
    hpf = plot(f,fitted_phase,'r');
    hpfr = plot(f,fitted_phase-2*delta,'m--',f,fitted_phase+2*delta,'m--');
    hold off;
    ylabel('Phase (rad)');
    if f(1) > f(end)
        min_f = f(end);
        max_f = f(1);
    else
        min_f = f(1);
        max_f = f(end);
    end
    xlim([min_f,max_f]);
    xlabel('Frequency (THz)');
    legend('PSD','Phase','Fitted Phase','95% Prediction Interval');
    set(hI,'linewidth',2);set(hp,'linewidth',2);set(hpf,'linewidth',2);set(hpfr,'linewidth',2);
    set(gca,'fontsize',14);

    % Find temporal plot range for later use
    intensity = abs(At).^2;
    intensity_plot = intensity;
    threshold_factor = 100;
    intensity_plot(intensity<max(intensity)/threshold_factor) = 0;
    left = find(intensity_plot~=0,1);
    right = find(intensity_plot~=0,1,'last');
    center = floor((left+right)/2);
    span_factor = 2;
    span = floor((right-left)/2)*span_factor;
    left = floor(center-span);
    right = ceil(center+span);
    if left < 1, left = 1; end
    if right > Nt, right = Nt; end

    % Plot instantaneous frequency
    figure;
    yyaxis left
    hI = plot(t,abs(At).^2/max(abs(At).^2),'b');
    ylim([0,1.5]);
    ylabel('Power (norm.)');
    set(gca,'YColor','b');
    yyaxis right
    hc = plot(t,pulse_inst_freq,'Color',[0.8510,0.3255,0.0980]);
    ylabel('\Delta\nu (THz)');
    xlim([min(t(left:right)),max(t(left:right))]);
    xlabel('Time (ps)');
    set(hI,'linewidth',2);set(hc,'linewidth',2);
    set(gca,'fontsize',14);
    
    c = 299792.458; % nm/ps
    
    % Print the result
    fprintf('Fitted f0:      %6.4f(THz)\n',mu(1)/2/pi);
    fprintf('     = lambda0: %6.4f(nm)\n',c/(mu(1)/2/pi))
    fprintf('quardratic_phase: %8.6g(fs^2)\n',quardratic_phase);
    fprintf('cubic_phase: %8.6g(fs^3)\n',cubic_phase);
    if fitted_order >= 4
        fprintf('quintic_phase: %8.6g(fs^4)\n',quintic_phase);
    end
end

end

