function [Strehl_ratio,dechirped_FWHM,transform_limited_FWHM,peak_power,fig] = analyze_field( t,f,field,compressor_type,varargin )
%ANALYZE_FIELD It plots the field and the spectrum, as well as their
%dechirped and transform-limited counterparts.
%
% Input:
%   t: (N,1); time (ps)
%   f: (N,1); frequency (THz)
%   field: (Nt,?); the field to be analyzed
%                  If ? isn't 1, it'll choose the most energetic field
%   compressor_type: 'Treacy-r': reflection grating pair,
%                    'Treacy-t': transmission grating pair,
%                    'Treacy-beta2': consider only "beta2" term of a reflection grating pair
%                             (Please refer to Ch.6.2 in Applications of Nonlinear Fiber Optics (2ed), Agrawal)
%                    'Offner1': single-grating Offner compressor
%                    'Offner2': double-grating Offner compressor
%                               (Assume the off-center distance of the first grating is known and fixed.
%                                The grating separation is varied for pulse compression here.)
%
%   Required arguments for grating-based compressor:
%
%      grating_incident_angle: a scalar; the incident angle of light toward the grating
%      grating_spacing: a scalar; the line spacing of the grating
%
%   Required arguments for prism compressor:
%
%      alpha: the apex angle of prisms in a prism compressor (rad)
%      prism_material: material in Sellmeier_coefficients.m in GMMNLSE_algorithm/
%
%   Extra required arguments for the Offner compressor:
%
%       R: radius of curvature of the convex mirror (R and 2R for two mirrors; 2R for the concave mirror) (m)
%       offcenter: the off-center distance of the first grating (m)
%
% Optional input argument:
%   verbose: 1(true) or 0(false); whether to plot and display the results or not (default: true)
%   global_opt: 1(true) or 0(false); whether to use global optimization for the compressor (default: false)
%   ASE: ASE information which includes
%           ASE.t_rep: the repetition rate of the pulse in gain-rate-eqn model (s).
%                      This is used to compute the correct unit for the ASE spectrum.
%           ASE.spectrum: the ASE spectrum; a column vector
%   assumed_dechirped_duration: assumed pulse duration of the dechirped pulse.
%                               This is used to help interpolate the field such that the temporal sampling is high enough to resolute the duration after dechirped.

%% Move the required input arguments out of the optional input arguments, varargin
switch compressor_type
    case {'Treacy-r','Treacy-t','Treacy-beta2'}
        grating_incident_angle = varargin{1};
        grating_spacing = varargin{2};
        
        n = 2;
    case 'prism'
        alpha = varargin{1};
        prism_material = varargin{2};
        
        n = 2;
    case {'grism1','grism2'}
        grating_incident_angle = varargin{1};
        grating_spacing = varargin{2};
        alpha = varargin{3};
        prism_material = varargin{4};
        
        n = 4;
    case 'Offner1'
        grating_incident_angle = varargin{1};
        grating_spacing = varargin{2};
        R = varargin{3};
        
        n = 3;
    case 'Offner2'
        grating_incident_angle = varargin{1};
        grating_spacing = varargin{2};
        R = varargin{3};
        offcenter = varargin{4};

        n = 4;
    otherwise
        error('The value of compressor_type is wrong.');
end
if length(varargin) > n % optional input arguments
    varargin = varargin(n+1:end);
else
    varargin = {};
end

%% Default optional input arguments
% Accept only 4 optional inputs at most
numvarargs = length(varargin);
if numvarargs > 3
    error('analyze_field:TooManyInputs', ...
          'It takes only at most 3 optional inputs');
end

% Set defaults for optional inputs
verbose = true;
global_opt = false;
ASE = [];
optargs = {verbose,global_opt,ASE};

% Now put these defaults into the valuesToUse cell array, 
% and overwrite the ones specified in varargin.
optargs(1:numvarargs) = varargin;

% Place optional args in memorable variable names
[verbose,global_opt,ASE] = optargs{:};

if ~verbose
    fig = [];
end

%%
% Pick the strongest field only
[~,mi] = max(sum(abs(field).^2,1));
field = field(:,mi(1));

c = 299792.458; % nm/ps
wavelength = c./f; % nm
Nt = size(field,1);
dt = t(2)-t(1); % ps
factor_correct_unit = (Nt*dt)^2/1e3; % to make the spectrum of the correct unit "nJ/THz"
                                     % "/1e3" is to make pJ into nJ
spectrum = abs(fftshift(ifft(field),1)).^2*factor_correct_unit; % in frequency domain

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

%% Plot the field and the spectrum
if verbose
    [t_rms,t_c] = calc_RMS(t,abs(field).^2);
    factor = c./wavelength.^2; % change the spectrum from frequency domain into wavelength domain
    [wavelength_rms,wavelength_c] = calc_RMS(wavelength,spectrum.*factor);
    
    fig(1) = figure('Name','Field and Spectrum');
    fp = get(gcf,'position');
    screen_size = get(0,'ScreenSize');
    original_top = screen_size(4)-fp(2)-fp(4);
    set(gcf,'position',[fp(1) screen_size(4)-original_top-fp(4)*7/4 fp(3)*7/4 fp(4)]);
    subplot(1,2,1);
    h1 = plot(t,abs(field).^2);
    xlabel('Time (ps)'); ylabel('Power (W)');
    xlim(t_c+[-1,1]*t_rms*6);
    ax1 = gca;
    subplot(1,2,2);
    h2 = plot(wavelength,spectrum.*factor);
    xlabel('Wavelength (nm)'); ylabel('Spectrum (nJ/nm)');
    xlim(wavelength_c+[-1,1]*wavelength_rms*6);
    ax2 = gca;
    set(h1,'linewidth',2); set(h2,'linewidth',2);
    set(ax1,'fontsize',16); set(ax2,'fontsize',16);
end

%% Dechirped and Transform-limited
increase_sampling = true;
num_interp = 5;
FWHM_sampling_num = 50;
% Guess the required number of sampling points from the transform-limited pulse duration
while increase_sampling
    insert_idx = [linspace(1,Nt,(num_interp+1)*(Nt-1)+1)'; Nt+(1:num_interp)'/(num_interp+1)];
    t_interp = interp1(t,insert_idx,'linear','extrap');

    field_f = ifft(field);
    field_f = cat(1,field_f(1:ceil(Nt/2),:),zeros(Nt*num_interp,1),field_f(ceil(Nt/2)+1:end,:));
    f_interp = (-floor(Nt/2):ceil(Nt/2)-1)'/(Nt*dt) + f(floor(Nt/2)+1);
    field_interp = fft(field_f);
    
    [transform_limited_field,~,transform_limited_FWHM,pulse_FWHM] = calc_transform_limited( field_interp,0,t_interp );
    if transform_limited_FWHM/(dt*1e3/(1+num_interp)) < FWHM_sampling_num
        num_interp = max(5,ceil(dt*1e3/(transform_limited_FWHM/FWHM_sampling_num) - 1)); % "FWHM_sampling_num" sampling points for FWHM; lower-bound this value to 5
    else
        increase_sampling = false;
    end
end

% Dechirp the pulse
switch compressor_type
    case 'prism'
        [~,dechirped_FWHM,dechirped_field] = pulse_compressor(compressor_type,                    [],feval(@(x)x(1),ifftshift(c./f_interp,1)),t_interp,field_interp,alpha,prism_material,false,global_opt);
    case {'grism1','grism2'}
        [~,dechirped_FWHM,dechirped_field] = pulse_compressor(compressor_type,grating_incident_angle,feval(@(x)x(1),ifftshift(c./f_interp,1)),t_interp,field_interp,grating_spacing,alpha,prism_material,false,global_opt,-1);
    case {'Treacy-t','Treacy-r','Treacy-beta2'}
        [~,dechirped_FWHM,dechirped_field] = pulse_compressor(compressor_type,grating_incident_angle,feval(@(x)x(1),ifftshift(c./f_interp,1)),t_interp,field_interp,grating_spacing,false,global_opt,-1);
    case 'Offner1'
        [~,dechirped_FWHM,dechirped_field] = pulse_compressor(compressor_type,grating_incident_angle,feval(@(x)x(1),ifftshift(c./f_interp,1)),t_interp,field_interp,grating_spacing,R,false,global_opt,-1);
    case 'Offner2'
        [~,dechirped_FWHM,dechirped_field] = pulse_compressor(compressor_type,grating_incident_angle,feval(@(x)x(1),ifftshift(c./f_interp,1)),t_interp,field_interp,grating_spacing,R,offcenter,false,global_opt,-1);
end

% -------------------------------------------------------------------------

% Strehl ratio
peak_power = max(abs(dechirped_field).^2);
Strehl_ratio = peak_power/max(abs(transform_limited_field).^2);

% Plot only the central part of the tme window of the dechirped and 
% transform-limited pulse because their duration are too small compared to
% the time window
if verbose
    intensity = abs(dechirped_field).^2;
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
    if right > Nt*(num_interp+1), right = Nt*(num_interp+1); end

    fig(2) = figure('Name','Transform-limited vs. Dechirped');
    fp = get(gcf,'position');
    screen_size = get(0,'ScreenSize');
    original_top = screen_size(4)-fp(2)-fp(4);
    set(gcf,'position',[fp(1) screen_size(4)-original_top-fp(4)*5/4 fp(3:4)*5/4]);
    h3 = plot(t_interp*1e3,abs(transform_limited_field).^2/1e3);
    hold on;
    h4 = plot(t_interp*1e3,abs(dechirped_field).^2/1e3);
    hold off;
    xlim([min(t_interp(left:right)) max(t_interp(left:right))]*1e3);
    xlabel('Time (fs)'); ylabel('Power (kW)');
    title('Transform-limited vs. Dechirped');
    legend('transform-limited','dechirped');
    set(h3,'linewidth',2); set(h4,'linewidth',2);
    set(gca,'fontsize',16);

    % Print the results
    if exist('cprintf','file')
        cprintf('blue','Pulse duration: %6.4f(fs)\n',pulse_FWHM);
        cprintf('blue','Dechirped duration: %6.4f(fs)\n',dechirped_FWHM);
        cprintf('blue','Transform-limited duration: %6.4f(fs)\n',transform_limited_FWHM);
        cprintf('red','--> Strehl ratio = %6.4f\n',Strehl_ratio);
    else
        fprintf('Pulse duration: %6.4f(fs)\n',pulse_FWHM);
        fprintf('Dechirped duration: %6.4f(fs)\n',dechirped_FWHM);
        fprintf('Transform-limited duration: %6.4f(fs)\n',transform_limited_FWHM);
        fprintf('--> Strehl ratio = %6.4f\n',Strehl_ratio);
    end
    fprintf('Peak power = %6.4f(kW)\n',peak_power/1e3);
    
    % Spectrogram
    [~,~,~,fig(3)] = calc_spectrogram(t,f,field);
end

% Plot the total spectrum including pulse and ASE.
% Unit of the ASE spectrum:
% It is W/THz. However, to plot it with the pulse spectrum, it needs to be 
% transformed into nJ/nm. Unlike the noisy ASE "field" added to the pulse
% during pulse evolution (see pulse-propagation codes), the spectrometer
% measures the "total" ASE power, irrelevant to the numerical time window.
% Therefore, the energy of each pulse should include the ASE energy 
% throughout one repetition time:
%    ASE energy (J/THz) = ASE power (W/THz) * t_rep (s)
if verbose && ~isempty(ASE)
    fig(4) = figure('Name','Spectrum');
    spectrum_total = spectrum.*factor + (ASE.spectrum.*ASE.t_rep*1e9).*factor; % 1e9 in the ASE term is to make it nJ
    plot(299792.458./f(f>0),  spectrum_total(f>0),'r','linewidth',2); hold on;
    plot(299792.458./f(f>0),spectrum(f>0).*factor(f>0),'b','linewidth',2); hold off;
    l = legend('Pulse+ASE','Pulse'); set(l,'fontsize',16);
    xlabel('Wavelength (nm)');
    ylabel('Spectrum (nJ/nm)');
    xlim([min(299792.458./f(f>0)),max(299792.458./f(f>0))]);
    set(gca,'fontsize',16);
end

end