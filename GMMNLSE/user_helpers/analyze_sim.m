function func = analyze_sim
%ANALYZE_SIM It contains several functions used for analyzing the fields
%along the fibers, mostly used for those inside the cavity
%
% This is the main function. Call this function first and then use it to
% call the other functions.
% Please refer to the functions below for the information about how to use
% them.
%
% Use:
%   func = analyze_sim;
%
%   [saved_fields,saved_z] = func.extract_saved_fields(fields,current_z,save_period,num_max_save );
%   fig = func.analyze_fields(t,f,fields,saved_z,splice_z,midx);
%   fig = func.analyze_ellipticity(fields,saved_z,ellipticity);
%   fig = analyze_gain(saved_z,splice_z,pump,N2);
%   fig = analyze_ASE(f,ASE,saved_z,splice_z,midx);
%   fig = func.animation(t,f,fields,saved_z,varargin);

func.extract_saved_field = @extract_saved_fields;
func.analyze_fields      = @analyze_fields;
func.analyze_ellipticity = @analyze_ellipticity;
func.analyze_gain        = @analyze_gain;
func.analyze_ASE         = @analyze_ASE;
func.animation           = @animation;

end

%%
function [saved_fields,saved_z] = extract_saved_fields( fields,num_max_save,current_z,z_info )
%EXTRACT_SAVED_FIELDS It extracts "num_max_save" fields out of "fields".
%
% Input:
%   fields: (N,num_modes,num_saved_fields); the saved fields along the fibers
%   num_max_save: if save_period is too small, that is, "num_saved_fields" 
%                 is too high, but you actually don't want too many saved 
%                 fields, this variable will choose fields to be saved with
%                 an almost equally distributed distance.
%   current_z: the z position of the input end of the fiber. This is used to
%              calculate saved_z.
%   save_period: the z interval between each saved fields
%
% Output:
%   saved_fields: the output saved fields
%   saved_z: the z position of each saved fields

if length(z_info) == 1
    save_period = z_info;
else
    z = z_info;        
end

num_fields = size(fields,3);
if num_fields > num_max_save
    save_idx = [1 ceil((1:num_max_save-1)*num_fields/(num_max_save-1))];
else
    save_idx = 1:num_fields;
end
if length(z_info) == 1
    saved_z = current_z + save_period*(save_idx-1);
else
    saved_z = current_z + z(save_idx);
end
saved_fields = fields(:,:,save_idx);

end

%%
function fig = analyze_fields(t,f,fields,saved_z,splice_z,midx)
%ANALYZE_FIELDS_WITHIN_CAVITY This function generates two plots of the
%information of the fields along the fibers.
% If the number of modes considered here is larger than one, midx>1, it'll
% plot only the total energy along the fibers; otherwise, it'll plot two
% figures, one with energy, bandwidth, and duration, the other with two
% colormaps of fields and spectra along the fibers.
%
% Input:
%   t: (N,1); time (ps)
%   f: (N,1); frequency (THz)
%   fields: (N,num_modes,num_fields); the fields along the fibers
%   saved_z: (1,num_saved_z); the z position of each field
%   splice_z: (1,num_splice); the z position of the splice (default: [])
%   midx: the mode index (default: all the modes = 1:num_modes)
%   **splice_z and midx are optional arguments.
%
% Output:
%   fig: the figure handle of the plots
%
% -------------------------------------------------------------------------
% Unit explanation:
%   intensity = abs(field).^2;
%   energy = trapz(t,intensity) = trapz(intensity)*dt;       % pJ
%   
%   spectrum_unknown_unit = abs(fftshift(ifft(field),1)).^2;
%
%   Parseval's theorem: sum(intensity) = sum(spectrum_unknown_unit)*N;
%                       * Note that spectrum_unknown_unit is from "ifft".
%   therefore sum(intensity)*dt = sum(spectrum_unknown_unit)*N*dt
%                               = sum(spectrum_unknown_unit)*(N*dt)^2/(N*dt)
%                               = sum(spectrum_unknown_unit)*(N*dt)^2*df
%                               = sum(spectrum_f)*df
%
%   spectrum_f = spectrum_unknown_unit*(N*dt)^2;
%   energy = trapz(f,spectrum_f) = trapz(spectrum_f)*df      % pJ
%                                = trapz(spectrum_f)/(N*dt);
%
%   c = 299792.458;     % nm/ps
%   wavelength = c./f;  % nm
%   spectrum_wavelength = spectrum_f.*(c./wavelength.^2);
%   energy = -trapz(wavelength,spectrum_wavelength);         % pJ
% -------------------------------------------------------------------------

sf = size(fields);

Nt = sf(1);
num_modes = sf(2);

if ~exist('midx','var')
    midx = 1:num_modes;
end
if ~exist('splice_z','var')
    splice_z = [];
end

dt = t(2) - t(1); % ps
c = 299792.458; % nm/ps
wavelength = c./f; % nm
positive_wavelength = wavelength>0; % exclude those negative wavelengths due to a large frequency window
intensity = abs(fields).^2; % W
%peak_power = permute(max(intensity),[3 2 1]);
energy = permute(trapz(intensity),[3 2 1])*dt/1e3; % nJ
factor_correct_unit = (Nt*dt)^2/1e3; % to make the spectrum of the correct unit "nJ/THz"
                                     % "/1e3" is to make pJ into nJ
spectrum = abs(fftshift(ifft(fields,[],1),1)).^2*factor_correct_unit; % in frequency domain

%1.  "real" is just to suppress the warning if calc_RMS gives a warning about
%    its imaginary output.
%    If the input is really noisy, the calculated RMS can be imaginary.
% 2. Cleaning the background noise is important for RMS calculations.
intensity_clean_background = intensity; spectrum_clean_background = spectrum;
for mi = 1:num_modes
    for zidx = 1:size(fields,3)
        intensity_clean_background(intensity(:,mi,zidx)<max(intensity(:,mi,zidx))*1e-3,mi,zidx) = 0;
        spectrum_clean_background(spectrum(:,mi,zidx)<max(spectrum(:,mi,zidx))*1e-3,mi,zidx) = 0;
    end
end
duration = real(permute(calc_RMS(t,intensity_clean_background),[3 2 1])); % ps
bandwidth = real(permute(calc_RMS(wavelength(positive_wavelength),spectrum_clean_background(positive_wavelength,:,:)),[3 2 1])); % nm

factor = c./wavelength(positive_wavelength).^2; % change the spectrum from frequency domain into wavelength domain

%{
% Calculate the Strehl ratio
grating_incident_angle = pi/6;
grating_spacing = 1e-6; % assume a 1000 line/mm grating here
wavelength0 = feval(@(x)x(1), ifftshift(wavelength,1));
num_interp = 5;
dechirped_field = zeros(sf);
for i = 1:sf(3)
    [~,~,dechirped_field(:,1,i)] = pulse_compressor('t',grating_incident_angle,wavelength0,t,fields(:,midx,i),grating_spacing,false);
    [quardratic_phase(i),cubic_phase(i)] = characterize_spectral_phase( f,fields(:,midx,i));
end
transform_limited_field = calc_transform_limited(fields(:,midx,:),num_interp,t);
Strehl_ratio = squeeze(max(abs(dechirped_field).^2)./max(abs(transform_limited_field).^2));
%}

% Plot
if length(midx) == 1
    fig(1) = figure('Name','Cavity: multiple info');
    fp = get(fig(1),'position');
    screen_size = get(0,'ScreenSize');
    original_top = screen_size(4)-fp(2)-fp(4);
    set(fig(1),'position',[fp(1) screen_size(4)-original_top-fp(4)*7/4 fp(3)*5/4 fp(4)*7/4]);
    
    subplot(3,1,1); h = plot(saved_z,energy(:,midx)); xlim([saved_z(1) saved_z(end)]);
    xlabel('Propagation length (m)'); ylabel('Pulse enegy (nJ)');
    title('Energy within the cavity');
    set(h,'linewidth',2); set(gca,'fontsize',14);
    hold on;
    for i = 1:length(splice_z)
        h = plot([splice_z(i);splice_z(i)],[min(energy(:,midx));max(energy(:,midx))]);
        set(h,'linewidth',2,'linestyle','--','color','black','marker','none');
    end
    hold off;
    %{
    subplot(3,1,2); h = plot(saved_z,peak_power(:,midx)); xlim([saved_z(1) saved_z(end)]);
    xlabel('Propagation length (m)'); ylabel('Peak power (W)');
    title('Peak power within the cavity');
    set(h,'linewidth',2); set(gca,'fontsize',14);
    hold on;
    for i = 1:length(splice_z)
        h = plot([splice_z(i);splice_z(i)],[min(peak_power(:,midx));max(peak_power(:,midx))]);
        set(h,'linewidth',2,'linestyle','--','color','black','marker','none');
    end
    hold off;
    %}
    subplot(3,1,[2 3]);
    yyaxis left;
    h = plot(saved_z,bandwidth(:,midx)); xlim([saved_z(1) saved_z(end)]);
    xlabel('Propagation length (m)'); ylabel('Bandwidth (nm)');
    title('RMS info within the cavity');
    set(h,'linewidth',2); set(gca,'fontsize',14);
    hold on;
    for i = 1:length(splice_z)
        h = plot([splice_z(i);splice_z(i)],[min(bandwidth(:,midx));max(bandwidth(:,midx))]);
        set(h,'linewidth',2,'linestyle','--','color','black','marker','none');
    end
    hold off;
    yyaxis right;
    h = plot(saved_z,duration(:,midx)); xlim([saved_z(1) saved_z(end)]);
    xlabel('Propagation length (m)'); ylabel('Duration (ps)');
    set(h,'linewidth',2); set(gca,'fontsize',14);
    %{
    hold on;
    for i = 1:length(splice_z)
        h = plot([splice_z(i);splice_z(i)],[min(duration(:,midx));max(duration(:,midx))]);
        set(h,'linewidth',2,'linestyle','--','color','black','marker','none');
    end
    hold off;
    %}
    
    fig(2) = figure('Name','Cavity: field evolution');
    fp = get(fig(2),'position');
    original_top = screen_size(4)-fp(2)-fp(4);
    set(fig(2),'position',[fp(1) screen_size(4)-original_top-fp(4)*7/4 fp(3)*5/4 fp(4)*7/4]);
    subplot(2,1,1); pcolor(saved_z,t,permute(intensity(:,midx,:),[1 3 2])); shading interp; cb = colorbar; ylim([t(1) t(end)]);
    xlabel('Propagation length (m)'); ylabel('Time (ps)');
    title('Pulse within the cavity');
    set(cb,'Location','southoutside');
    set(gca,'fontsize',14);
    hold on;
    for i = 1:length(splice_z)
        h = plot(splice_z(i)*ones(Nt,1),t);
        set(h,'linewidth',2,'linestyle','--','color','white','marker','none');
    end
    hold off;
    subplot(2,1,2); pcolor(saved_z,wavelength(positive_wavelength),permute(spectrum(positive_wavelength,midx,:),[1 3 2]).*factor); shading interp; cb = colorbar; ylim([min(wavelength(positive_wavelength)) max(wavelength(positive_wavelength))]);
    xlabel('Propagation length (m)'); ylabel('Wavelength (nm)');
    title('Spectrum within the cavity');
    set(cb,'Location','southoutside');
    set(gca,'fontsize',14);
    hold on;
    for i = 1:length(splice_z)
        h = plot(splice_z(i)*ones(sum(positive_wavelength),1),wavelength(positive_wavelength));
        set(h,'linewidth',2,'linestyle','--','color','white','marker','none');
    end
    hold off;
    
    %{
    fig(3) = figure('Name','Cavity: Strehl ratio');
    h = plot(saved_z,Strehl_ratio);
    xlabel('Propagation length (m)'); ylabel('Strehl ratio');
    title('Strehl ratio within the cavity');
    set(h,'linewidth',2); set(gca,'fontsize',20);
    hold on;
    for i = 1:length(splice_z)
        h = plot([splice_z(i);splice_z(i)],[min(Strehl_ratio);max(Strehl_ratio)]);
        set(h,'linewidth',2,'linestyle','--','color','black','marker','none');
    end
    hold off;
    %}
else
    fig = figure('Name','Cavity: energy');
    total_energy = sum(energy,2);
    h = plot(saved_z,total_energy); xlim([saved_z(1) saved_z(end)]);
    xlabel('Propagation length (m)'); ylabel('Pulse enegy (nJ)');
    title('Energy within the cavity');
    set(h,'linewidth',2); set(gca,'fontsize',20);
    hold on;
    for i = 1:length(splice_z)
        h = plot([splice_z(i);splice_z(i)],[min(total_energy);max(total_energy)]);
        set(h,'linewidth',2,'linestyle','--','color','black','marker','none');
    end
    hold off;
end
drawnow;

end

%%
function fig = analyze_ellipticity(fields,saved_z,ellipticity)
%ANALYZE_ELLIPTICITY_WITHIN_CAVITY It calculates the ellipticity of the
%fields within the cavity.
%
% Input:
%   fields: (Nt,2,num_fields); the fields along the fibers
%           Because it calculates the ellipticity, so of course it needs 2
%           orthogonal polarization modes.
%   saved_z: (1,num_saved_z); the z position of each field
%   ellipticity: a scalar; the ellipticity of the orthogonal basis
%                  e.g. 0 = linear polarization
%                       1 = circular polarization
% Output:
%   fig: the figure handle of the plot

Nt = size(fields,1);

% Ellipticity within the cavity
[ellipticity_phi,ellipticity_theta] = calc_ellipticity( fields,ellipticity );
[~,peak_idx] = max(abs(fields(:,1,:)).^2);
peak_idx = peak_idx + permute(0:Nt:(Nt*length(peak_idx)-1),[1 3 2]);
peak_phi = squeeze(ellipticity_phi(peak_idx));
peak_theta = squeeze(ellipticity_theta(peak_idx));

fig = figure('Name','Cavity: ellipticity');
yyaxis left
h1 = plot(saved_z,peak_phi);
ylim([-5,365]); % 0~360 (deg)
xlabel('Propagation length (m)');  ylabel('\phi (deg)');
yyaxis right
h2 = plot(saved_z,peak_theta);
ylim([-5,95]); % 0~90 (deg)
ylabel('\theta (deg)');
xlim([saved_z(1),saved_z(end)]);
title('Ellipticity within the cavity');

set(h1,'linewidth',2); set(h2,'linewidth',2);
set(gca,'fontsize',20);

end

%%
function fig = analyze_gain(saved_z,splice_z,pump,population)
%ANALYZE_GAIN_WITHIN_CAVITY This function generates two plots of the
%pump power and the ion density of the upper state for the gain rate-eqn
%model.
%
%   saved_z - (1,num_z); the z coordinate of each point of the fiber
%   splice_z - (1,num_splice); the z coordinate of each splice
%   pump - (1,1,num_z); the pump power along the gain fiber
%   population - (1,1,num_z,num_population);
%        the ion density of energy states
%        for the computation of the fundamental-mode gain rate equation
%        For multimode, population is of the size (Nx,Nx,num_z), where Nx 
%        is the number of points of the fiber cross section, so it's not 
%        possible to analyze it with a 1D plot.

pump_forward = permute(pump.forward,[3,1,2]);
pump_backward = permute(pump.backward,[3,1,2]);

fig = figure('Name','Gain info');
fp = get(fig,'position');
screen_size = get(0,'ScreenSize');
original_top = screen_size(4)-fp(2)-fp(4);
set(gcf,'position',[fp(1) screen_size(4)-original_top-fp(4)*7/4 fp(3)*5.5/4 fp(4)*7/4]);
if nargin > 3
    subplot(2,1,1);
end
h = plot(saved_z,[pump_forward,pump_backward]);
xlim([min(saved_z) max(saved_z)]);
xlabel('Propagation distance (m)'); ylabel('Pump power (W)');
title('Pump power within the cavity');
set(h,'linewidth',2); set(gca,'fontsize',20);
hold on;
for i = 1:length(splice_z)
    h = plot(splice_z(i)*[1 1],[min([pump_forward;pump_backward]) max([pump_forward;pump_backward])]);
    set(h,'linewidth',2,'linestyle','--','color','black','marker','none');
end
hold off;
legend('Forward pump','Backward pump');

if nargin > 3
    population = cat(4,1-sum(population,4),population);
    population = population*100; % in "%"
    subplot(2,1,2);
    h = plot(saved_z,squeeze(population));
    xlim([min(saved_z) max(saved_z)]);
    xlabel('Propagation distance (m)'); ylabel('Doped ion excitation (%)');
    title('Doped ion excitation within the cavity');
    set(h,'linewidth',2); set(gca,'fontsize',20);
    hold on;
    for i = 1:length(splice_z)
        h = plot(splice_z(i)*[1 1],[min(population(:)) max(population(:))]);
        set(h,'linewidth',2,'linestyle','--','color','black','marker','none');
    end
    hold off;
    
    Nname = cell(1,size(population,4));
    for Ni = 1:size(population,4)
        Nname{Ni} = sprintf('N%u',Ni-1);
    end
    legend(Nname);
end

end

%%
function fig = analyze_ASE(f,ASE,saved_z,splice_z,midx)
%ANALYZE_ASE_WITHIN_CAVITY This function plots the information of the 
% forward and backward ASE along the fibers.
% If the number of modes considered here is larger than one, midx>1, it'll
% plot only the total energy along the fibers; otherwise, it'll plot two
% figures, one with energy, the other with a colormap of the spectra along 
% the fibers.
%
% Input:
%   f: (Nt,1); frequency (THz)
%   ASE.forward:  (Nt,num_modes,num_fields); the forward ASE along the fibers
%   ASE.backward: (Nt,num_modes,num_fields); the backward ASE along the fibers
%   saved_z: (1,num_saved_z); the z position of each field
%   splice_z: (1,num_splice); the z position of the splice (default: [])
%   midx: the mode index (default: all the modes = 1:num_modes)
%   **splice_z and midx are optional arguments.

sASE = size(ASE.forward);

Nt = sASE(1);
num_modes = sASE(2);

if ~exist('midx','var')
    midx = 1:num_modes;
end
if ~exist('splice_z','var')
    splice_z = [];
end

power.forward = permute(trapz(f,ASE.forward,1),[3 2 1])*1e3; % mW
power.backward = permute(trapz(f,ASE.backward,1),[3 2 1])*1e3;

if length(midx) == 1
    fig(1) = figure('Name','Cavity: ASE power');
    h  = plot(saved_z,power.forward(:,midx),'linewidth',2,'Color','b'); hold on;
    h2 = plot(saved_z,power.backward(:,midx),'linewidth',2,'Color','r');
    xlim([saved_z(1) saved_z(end)]);
    xlabel('Propagation length (m)'); ylabel('ASE power (mW)');
    set(gca,'fontsize',20);
    for i = 1:length(splice_z)
        h = plot([splice_z(i);splice_z(i)],[min([power.forward(:,midx);power.backward(:,midx)]);max([power.forward(:,midx);power.backward(:,midx)])]);
        set(h,'linewidth',2,'linestyle','--','color','black','marker','none');
    end
    hold off;
    legend('Forward ASE','Backward ASE');
    
    fig(2) = figure('Name','Cavity: ASE power evolution');
    fp = get(fig(2),'position');
    screen_size = get(0,'ScreenSize');
    original_top = screen_size(4)-fp(2)-fp(4);
    set(fig(2),'position',[fp(1) screen_size(4)-original_top-fp(4)*7/4 fp(3)*5/4 fp(4)*7/4]);
    subplot(2,1,1);
    pcolor(saved_z,f,permute(ASE.forward(:,midx,:),[1 3 2])); shading interp; cb = colorbar; ylim([f(1) f(end)]);
    xlabel('Propagation length (m)'); ylabel('Frequency (THz)');
    title('Forward ASE within the cavity');
    set(cb,'Location','southoutside');
    set(gca,'fontsize',14);
    hold on;
    for i = 1:length(splice_z)
        h = plot(splice_z(i)*ones(Nt,1),f);
        set(h,'linewidth',2,'linestyle','--','color','white','marker','none');
    end
    hold off;
    subplot(2,1,2);
    pcolor(saved_z,f,permute(ASE.backward(:,midx,:),[1 3 2])); shading interp; cb = colorbar; ylim([f(1) f(end)]);
    xlabel('Propagation length (m)'); ylabel('Frequency (THz)');
    title('Backward ASE within the cavity');
    set(cb,'Location','southoutside');
    set(gca,'fontsize',14);
    hold on;
    for i = 1:length(splice_z)
        h = plot(splice_z(i)*ones(Nt,1),f);
        set(h,'linewidth',2,'linestyle','--','color','white','marker','none');
    end
    hold off;
else
    fig = figure('Name','Cavity: ASE power');
    total_power.forward = sum(power.forward,2); total_power.backward = sum(power.backward,2);
    h  = plot(saved_z,total_power.forward,'linewidth',2,'Color','b'); hold on;
    h2 = plot(saved_z,total_power.backward,'linewidth',2,'Color','r');
    xlim([saved_z(1) saved_z(end)]);
    xlabel('Propagation length (m)'); ylabel('ASE power (mW)');
    set(gca,'fontsize',20);
    hold on;
    for i = 1:length(splice_z)
        h = plot([splice_z(i);splice_z(i)],[min([total_power.forward;total_power.backward]);max([total_power.forward;total_power.backward])]);
        set(h,'linewidth',2,'linestyle','--','color','black','marker','none');
    end
    hold off;
    legend('Forward ASE','Backward ASE');
end
drawnow;

end

%%
function fig = animation(t,f,fields,saved_z,varargin)
%ANIMATION_WITHIN_CAVITY It shows the evolution of the pulse by plotting
%its spectrogram, field, and spectrum along the fibers.
%
%   t - (Nt,1); time (ps)
%   f - (Nt,1); frequency (THz)
%   fields - (Nt,num_modes,num_fields_in_fibers); the pulse field
%   saved_z - (1,num_fields_in_fibers); the propagation length of each saved field
%
%   Optional arguments:
%
%       splice_z - (1,num_splice); the location of the splice (default: [])
%       midx - a scalar; the mode to plot (default: 1)
%       fps - a scalar; the frame rate (default: 20)
%       unit - 'norm' or 'real'; the unit for the spectrum. (default: 'norm')
%              'norm' means normalized data and 'real' means demonstrating
%              the spectrum under the physical unit
%       tlim - (1,2) matrix; the range of the time to plot (ps) (default: [])
%       wavelengthlim - (1,2) matrix; the range of the wavelength to plot (nm) (default: [])
%       filename - the filename for the exported movie file (default: 'no_specified_name')
%       export_format - 'gif' or 'avi'; the exported movie file format (default: '')
%       **filename and export_format should be specified simultaneously.

numvarargin = nargin - 4;
% Set defaults for optional inputs
optargs = {[],1,20,'norm',[],[],'no_specified_name',''};
% Now put these defaults into the valuesToUse cell array, 
% and overwrite the ones specified in varargin.
optargs(1:numvarargin) = varargin;
% Place optional args in memorable variable names
[splice_z,midx,fps,unit,tlim,wavelengthlim,filename,export_format] = optargs{:};
if nargin == 11
    error('analyze_sim:InputArgumentsError',...
          '"filename" and "export_format" should both be specified.');
elseif nargin < 4 || nargin > 12
    error('analyze_sim:InputArgumentsError',...
          'The number of input arguments doesn''t match the requirements. The minimum of 4 inputs are required.');
end
if ~any(strcmp(unit,{'norm','real'}))
    error('analyze_sim:UnitError',...
          'The unit of the spectrum is either "real" or "norm".');
end

save_point = size(fields,3);

c = 299792.458; % nm/ps
wavelength = c./f; % nm
Nt = size(fields,1);
dt = t(2)-t(1); % ps
factor_correct_unit = (Nt*dt)^2/1e3; % to make the spectrum of the correct unit "nJ/THz"
                                     % "/1e3" is to make pJ into nJ
factor = c./wavelength.^2; % change the spectrum from frequency domain into wavelength domain
spectrum = abs(fftshift(ifft(fields(:,midx,:)),1)).^2*factor_correct_unit.*factor; % in wavelength domain
intensity = abs(fields(:,midx,:)).^2;

f0 = feval(@(x)x(1), ifftshift(f,1));
pulse_phase = unwrap(angle(fields(:,midx,:)));
smooth_range = 50;
for i = 1:save_point
    pulse_phase(:,:,i) = smooth(pulse_phase(:,:,i),smooth_range);
end
pulse_inst_freq = -(pulse_phase(3:end,1,:)-pulse_phase(1:end-2,1,:))/(2*dt)/(2*pi)+f0; % THz; I use "central difference" to calculate the slope here
if isempty(wavelengthlim)
    max_pulse_inst_freq = max(pulse_inst_freq(:));
    min_pulse_inst_freq = min(pulse_inst_freq(:));
else
    max_pulse_inst_freq = c/wavelengthlim(1);
    min_pulse_inst_freq = c/wavelengthlim(2);
end

% -------------------------------------------------------------------------
% Unit explanation:
%   intensity = abs(field).^2;
%   energy = trapz(t,intensity) = trapz(intensity)*dt;       % pJ
%   
%   spectrum_unknown_unit = abs(fftshift(ifft(field),1)).^2;
%
%   Parseval's theorem: sum(intensity) = sum(spectrum_unknown_unit)*Nt;
%                       * Note that spectrum_unknown_unit is from "ifft".
%   therefore sum(intensity)*dt = sum(spectrum_unknown_unit)*Nt*dt
%                               = sum(spectrum_unknown_unit)*(Nt*dt)^2/(Nt*dt)
%                               = sum(spectrum_unknown_unit)*(Nt*dt)^2*df
%                               = sum(spectrum_f)*df
%
%   spectrum_f = spectrum_unknown_unit*(Nt*dt)^2;
%   energy = trapz(f,spectrum_f) = trapz(spectrum_f)*df      % pJ
%                                = trapz(spectrum_f)/(Nt*dt);
%
%   c = 299792.458;     % nm/ps
%   wavelength = c./f;  % nm
%   spectrum_wavelength = spectrum_f.*(c./wavelength.^2);
%   energy = -trapz(wavelength,spectrum_wavelength);         % pJ
% -------------------------------------------------------------------------

% the upper and lower limit of the plots
max_intensity = max(intensity(:));
min_intensity = min(intensity(:));

if isequal(unit,'norm')
    max_sp = max(spectrum);
    min_sp = min(spectrum);
    spectrum = (spectrum-min_sp)./(max_sp-min_sp);
    
    spectrum_figure_unit = 'norm.';
else
    spectrum_figure_unit = 'nJ/nm';
end
max_spectrum = max(spectrum(:));
min_spectrum = min(spectrum(:));

% Remove useless inst. frequency when the intensity is too weak
for i = 1:save_point
    pulse_inst_freq(intensity(2:end-1,1,i)<max(intensity(2:end-1,1,i))/100,1,i) = NaN;
end

% start the figure and set its size
fig = figure('Name','Spectrogram');
fp = get(fig,'position');
screen_size = get(0,'ScreenSize');
original_top = screen_size(4)-fp(2)-fp(4);
set(gcf,'position',[fp(1) screen_size(4)-original_top-fp(4)*7/4 fp(3)*5.5/4 fp(4)*7/4]);

% Start plotting
ax = gca;
ax.NextPlot = 'replaceChildren';
Frame(save_point) = struct('cdata',[],'colormap',[]);
for i = 1:save_point
    [psd,t_spectrogram,f_spectrogram] = calc_spectrogram(t,f,fields(:, midx, i),true,tlim,wavelengthlim,100,100,false);
    
    wavelength_spectrogram = c./f_spectrogram; % nm
    
    figure(fig);
    
    % spectrogram
    ax = subplot(4,3,[4,5,7,8]);
    positive_wavelength = find(wavelength_spectrogram>0,1);
    pcolor(t_spectrogram,wavelength_spectrogram(positive_wavelength:end),psd(positive_wavelength:end,:))
    shading interp; colormap(jet);
    colorbar('location','south','Color','[1 1 1]');
    %cb = colorbar('location','south','Color','[1 1 1]');
    %cb.Label.String = 'PSD (nJ/nm)'; % the unit of the spectrogram colorbar
    xlabel('Time (ps)');
    ylabel('Wavelength (nm)');
    title('Spectrogram');
    set(gca,'fontsize',14);
    % field
    ax2 = subplot(4,3,[1,2]);
    yyaxis left
    h1_1 = plot(t,intensity(:,1,i)); % intensity (W)
    ylim([min_intensity max_intensity]);
    ylabel('Power (W)');
    yyaxis right
    h1_2 = plot(t(2:end-1),pulse_inst_freq(:,1,i)); % instantaneous frequency (THz)
    ylabel('Inst. freq. (THz)');
    if ~isempty(tlim)
        xlim(tlim);
    end
    ylim([min_pulse_inst_freq,max_pulse_inst_freq]);
    title('Field');
    set(h1_1,'linewidth',2); set(h1_2,'linewidth',2); set(gca,'fontsize',14);
    % spectrum
    ax3 = subplot(4,3,[6,9]);
    h2 = plot(spectrum(:,1,i),wavelength);
    xlim([min_spectrum max_spectrum]);
    if ~isempty(wavelengthlim)
        ylim(wavelengthlim);
    else
        ylim([min(wavelength) max(wavelength)]);
    end
    xlabel(['PSD (' spectrum_figure_unit ')']);
    set(h2,'linewidth',2); set(gca,'fontsize',14);
    title('Spectrum');
    
    % link their axes so that they have the same range of x,y-axes
    linkaxes([ax,ax2],'x'); linkaxes([ax,ax3],'y');
    set(ax2,'XTickLabel',''); set(ax3,'YTickLabel','');
    
    % propagation
    ax4 = subplot(5,3,[13,14,15]);
    ax4p = get(ax4,'position');
    set(ax4,'position',[ax4p(1:2) ax4p(3) ax4p(4)/1.5]);
    energy = permute(trapz(intensity),[3 2 1])*dt/1e3; % nJ
    h3 = plot(saved_z,energy); xlim([saved_z(1) saved_z(end)]);
    xlabel('Propagation length (m)'); ylabel('Pulse enegy (nJ)');
    title('Energy within the cavity');
    set(h3,'linewidth',2); set(ax4,'fontsize',14);
    hold on;
    for j = 1:length(splice_z)
        h = plot([splice_z(j);splice_z(j)],[min(energy);max(energy)]);
        set(h,'linewidth',2,'linestyle','--','color','black','marker','none');
    end
    h4 = plot([saved_z(i);saved_z(i)],[min(energy);max(energy)]);
    set(h4,'linewidth',3);
    hc = get(h4,'Color'); set(h4,'Color',[hc,0.8]);
    hold off;
    
    Frame(i) = getframe(fig);
    
    % Export as "gif"
    if isequal(export_format,'gif')
        im = frame2im(Frame(i));
        [imind,cm] = rgb2ind(im,256);
        if i == 1
            imwrite(imind,cm,[filename '.gif'],'gif','Loopcount',inf,'DelayTime',1/fps);
        else
            imwrite(imind,cm,[filename '.gif'],'gif','WriteMode','append','DelayTime',1/fps);
        end
    end
end
close(fig);

% Movie
fig_movie = implay(Frame,fps);
set(fig_movie.Parent,'position',[fp(1) screen_size(4)-original_top-fp(4)*7/4 fp(3)*5.5/4 fp(4)*7/4]); % enlarge the figure size to fit in so many subplots

% Export as "avi"
if isequal(export_format,'avi')
    exportVideo = VideoWriter(filename);
    exportVideo.FrameRate = fps;
    open(exportVideo);
    writeVideo(exportVideo, Frame);
    close(exportVideo);
end

end
