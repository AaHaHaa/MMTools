function pulse_compressor_animator( compressor_type,separation,theta_in,wavelength0,time,field,varargin )
%PULSE_COMPRESSOR_ANIMATOR Make an animation of spectrograms of the fields 
% after differernt grating/prism separations of the compressor
%
%   compressor_type: 'Treacy-r': reflective grating pair,
%                    'Treacy-t': transmissive grating pair,
%                    'Treacy-beta2': consider only "beta2" term of a Treacy reflective grating compressor
%                                    (Please refer to Ch.6.2 in Applications of Nonlinear Fiber Optics (2ed), Agrawal)
%                    'prism': prism pair
%                             It has the opposite sign of TOD from the silica fiber.
%                             Dechirper based on a grating pair, despite beta2 compensation, only adds beta3 to the field, deteriorating the pulse quality.
%                             This can be solved by using a prism pair.
%                             However, because prism compressor has a smaller GDD, it is not suitable for huge dispersion compensation.
%                    'grism1': grism pair
%                              It has the opposite sign of TOD from the silica fiber as 'prism'.
%                              Due to the addition of gratings, it can compensate the phase with a smaller more-practical grism separation.
%                              This 'grism1' relies on a configuration whose grating is attached to the prism.
%                    'grism2': grism pair
%                              It has the opposite sign of TOD from the silica fiber as 'prism'.
%                              Due to the addition of gratings, it can compensate the phase with a smaller more-practical grism separation.
%                              This 'grism2' relies on a configuration whose grating is not attached to the prism.
%                              It looks like a combination of 'Treacy-t' and 'prism', where prism pair is added between two gratings.
%                    'Offner1': single-grating Offner compressor
%                               Offner1 is the Offner compressor with one off-centered grating
%                    'Offner2': double-grating Offner compressor
%                               (Assume the off-center distance of the first grating is known and fixed.
%                                The grating separation is varied for pulse compression here.)
%                               Offner2 is the Offner compressor with two gratings, one off-center and the other some distance away.
%                    'Offner3': true aberration-free double-grating Offner compressor
%                               Offner3 is the Offner compressor with two gratings, one on-center and the other some distance away.
%                    'Martinez': Martinez stretcher used as a compressor
%                        grating-to-lens distance < focal length: add normal dispersion
%                        grating-to-lens distance > focal length: add anomalous dispersion
%   separation: grating spacings which the pulse is dechirped with.
%                       This is an array that forms the movie.
%   theta_in: a scalar; For a grating pair, the incident angle on the grating (rad)
%             'prism' doesn't need this input since it's determined by the operation of minimum deviation. Just use [].
%   wavelength0: a scalar; central wavelength (nm)
%   time: (N,1); the time grid points (ps)
%   field: (N,1); the electric field in time domain
%   grating_spacing: a scalar; the spacing between each grating line (m)
%   alpha: the apex angle of prisms in a prism compressor (rad)
%   prism_material: material in Sellmeier_coefficients.m in GMMNLSE_algorithm/
%
%   Extra required arguments for the Offner compressor:
%
%       R: radius of curvature of the convex mirror (R and 2R for two mirrors; 2R for the concave mirror) (m)
%       offcenter: the off-center distance of the first grating (m);
%                  required only for 'Offner2'
%
%   Optional arguments:
%
%       tlim - (1,2) matrix; the range of the time to plot (ps) (default: [])
%       wavelengthlim - (1,2) matrix; the range of the wavelength to plot (nm) (default: [])
%       t_feature - a scalar; the ratio of the tiny pulse structure vs. pulse duration you want to resolve;
%                   the larger the number, the higher the time resolution (default: 50)
%       f_feature - a scalar; the ratio of the tiny spectral structure vs. pulse bandwidth you want to revolve
%                   the larger the number, the higher the frequency resolution (default: 50)
%       lambda_or_f - true (use wavelength) or false (use frequency);
%                     plot with frequency (THz) or wavelength (nm);
%                     If plot_yes=false, this acts nothing
%                    (default: true)
%       m: a scalar; the diffraction order (default: -1)
% =========================================================================
% Use:
%   Treacy:
%      pulse_compressor_animator( compressor_type,separation,theta_in,wavelength0,time,field,grating_spacing,tlim,wavelengthlim,t_feature,f_feature,lambda_or_f,m )
%
%       * compressor_type is either 'Treacy-t', 'Treacy-r', or 'Treacy-beta2'
%       * [tlim,wavelengthlim,t_feature,f_feature,lambda_or_f,m] are optional inputs
%
%   Prism:
%      pulse_compressor_animator( 'prism',separation,[],wavelength0,time,field,alpha,prism_material,tlim,wavelengthlim,t_feature,f_feature,lambda_or_f )
%
%       * [tlim,wavelengthlim,t_feature,f_feature,lambda_or_f,m] are optional inputs
%
%   Grism1:
%      pulse_compressor_animator( 'grism1',separation,theta_in,wavelength0,time,field,grating_spacing,alpha,prism_material,tlim,wavelengthlim,t_feature,f_feature,lambda_or_f,m )
%
%       * [tlim,wavelengthlim,t_feature,f_feature,lambda_or_f,m] are optional inputs
%
%   Grism2:
%      pulse_compressor_animator( 'grism2',separation,theta_in,wavelength0,time,field,grating_spacing,alpha,prism_material,tlim,wavelengthlim,t_feature,f_feature,lambda_or_f,m )
%
%      * [tlim,wavelengthlim,t_feature,f_feature,lambda_or_f,m] are optional inputs
%
%   Offner1:
%       pulse_compressor_animator( 'Offner1',separation,theta_in,wavelength0,time,field,grating_spacing,R,tlim,wavelengthlim,t_feature,f_feature,lambda_or_f,m )
%
%       * [tlim,wavelengthlim,t_feature,f_feature,lambda_or_f,m] are optional inputs
%
%   Offner2:
%       pulse_compressor_animator( 'Offner2',separation,theta_in,wavelength0,time,field,grating_spacing,R,offcenter,tlim,wavelengthlim,t_feature,f_feature,lambda_or_f,m )
%
%       * Compared to "Offner1", "offcenter" is required for an input here
%       * [tlim,wavelengthlim,t_feature,f_feature,lambda_or_f,m] are optional inputs
%
%   Offner3:
%       pulse_compressor_animator( 'Offner3',separation,theta_in,wavelength0,time,field,grating_spacing,R,tlim,wavelengthlim,t_feature,f_feature,lambda_or_f,m )
%
%       * [tlim,wavelengthlim,t_feature,f_feature,lambda_or_f,m] are optional inputs
%
%   Martinez:
%      pulse_compressor_animator( 'Martinez',separation,theta_in,wavelength0,time,field,grating_spacing,focal_length,tlim,wavelengthlim,t_feature,f_feature,lambda_or_f,m )
%
%       * [tlim,wavelengthlim,t_feature,f_feature,lambda_or_f,m] are optional inputs

switch compressor_type
    case {'Treacy-r','Treacy-t','Treacy-beta2'}
        grating_spacing = varargin{1};
        
        n = 1;
    case 'prism'
        alpha = varargin{1};
        prism_material = varargin{2};
        
        n = 2;
    case {'grism1','grism2'}
        grating_spacing = varargin{1};
        alpha = varargin{2};
        prism_material = varargin{3};
        
        n = 3;
    case {'Offner1','Offner3'}
        grating_spacing = varargin{1};
        R = varargin{2};
        
        n = 2;
    case 'Offner2'
        grating_spacing = varargin{1};
        R = varargin{2};
        offcenter = varargin{3};
        
        n = 3;
    case 'Martinez'
        grating_spacing = varargin{1};
        focal_length = varargin{2};
        
        n = 2;
    otherwise
        error('The value of compressor_type is wrong.');
end
if length(varargin) > n % optional input arguments
    varargin = varargin(n+1:end);
else
    varargin = {};
end

% Default parameters
optargs = {[],[],50,50,true,-1};
% Load paramters
optargs(1:length(varargin)) = varargin;
[tlim,wavelengthlim,t_feature,f_feature,lambda_or_f,m] = optargs{:};

switch compressor_type
    case {'Treacy-beta2','Treacy-r','Treacy-t'}
        Frame = Treacy(compressor_type,separation,theta_in,wavelength0,time,field,grating_spacing,                     m,tlim,wavelengthlim,t_feature,f_feature,lambda_or_f);
    case 'prism'
        Frame = prism(                 separation,         wavelength0,time,field,               alpha,prism_material,  tlim,wavelengthlim,t_feature,f_feature,lambda_or_f);
    case {'grism1','grism2'}
        Frame = grism(compressor_type,separation,theta_in,wavelength0,time,field,grating_spacing,alpha,prism_material,m,tlim,wavelengthlim,t_feature,f_feature,lambda_or_f);
    case 'Offner1'
        Frame = Offner(compressor_type,separation,theta_in,wavelength0,time,field,grating_spacing,R,0,                 m,tlim,wavelengthlim,t_feature,f_feature,lambda_or_f);
    case 'Offner2'
        Frame = Offner(compressor_type,separation,theta_in,wavelength0,time,field,grating_spacing,R,offcenter,         m,tlim,wavelengthlim,t_feature,f_feature,lambda_or_f);
    case 'Offner3'
        Frame = Offner(compressor_type,separation,theta_in,wavelength0,time,field,grating_spacing,R,0,                 m,tlim,wavelengthlim,t_feature,f_feature,lambda_or_f);
    case 'Martinez'
        Frame = Martinez(              separation,theta_in,wavelength0,time,field,grating_spacing,focal_length,          tlim,wavelengthlim,t_feature,f_feature,lambda_or_f);
    otherwise
        error('The value of compressor_type is wrong.');
end
fig_movie = implay(Frame,2);
fig = figure;
fp = get(fig,'position');
screen_size = get(0,'ScreenSize');
original_top = screen_size(4)-fp(2)-fp(4);
close(fig);
set(fig_movie.Parent,'position',[fp(1) screen_size(4)-original_top-fp(4)*7/4 fp(3)*5.5/4 fp(4)*7/4]); % enlarge the figure size to fit in so many subplots

end

%%
function Frame = Treacy(compressor_type,separation,theta_in,wavelength0,time,field,grating_spacing,m,tlim,wavelengthlim,t_feature,f_feature,lambda_or_f)
%GRATING_PAIR

N = length(time); % the number of time points
dt = abs(time(2)-time(1)); % ps
c = 299792.458; % nm/ps
wavelength = c./((-N/2:N/2-1)'/N/dt + c/wavelength0); % nm
if size(time,1) == 1
    time = time';
end

field_w = fftshift(ifft(field),1); % frequency domain

% Find the center wavelength of the pulse spectrum with its second moment
f = c./wavelength; % THz
f_c = calc_f_c(f,abs(field_w).^2);
wavelength_c = c/f_c;
wavelength = c./(f+f_c-c/wavelength0);
field = field.*exp(1i*2*pi*(f_c-c/wavelength0)*time);
field_w = fftshift(ifft(field),1);

if isequal(compressor_type,'Treacy-beta2')
    theta_out = asin( m*wavelength_c/(grating_spacing*1e9) + sin(theta_in) ); % the transmitted/reflected angle of the -1st order diffraction
else % 'r' and 't'
    theta_out = asin( m*wavelength/(grating_spacing*1e9) + sin(theta_in) ); % the transmitted/reflected angle of the m-th order diffraction
end

% The dechirped pulse
save_point = length(separation);
Frame(save_point) = struct('cdata',[],'colormap',[]);
for si = 1:save_point
    dechirped_field = Treacy_dechirping(compressor_type,separation(si),theta_in,theta_out,wavelength,wavelength_c,time,field_w,grating_spacing,m);
    
    [~,~,~,fig,ax,cb] = calc_spectrogram(time,c./wavelength,dechirped_field,true,tlim,wavelengthlim,t_feature,f_feature,true,lambda_or_f);
    colormap(whitejet_lower(512)); set(cb,'Color','[0 0 0]');
    ax.NextPlot = 'replaceChildren';
    Frame(si) = getframe(fig);
    close(fig);
end

end

%%
function [field,y,total_phase] = Treacy_dechirping(compressor_type,separation,theta_in,theta_out,wavelength,wavelength0,time,field_w,grating_spacing,m)
%TREACY_DECHIRPING It finds the field after propagating through the 
%compressor

if separation > 0
    % The light can go out of the grating pairs under certain angles or
    % wavelengths, so this part of light is rejected from the compressor.
    rejected_part_f = imag(theta_out)~=0 | wavelength<0;
    rejected_part = rejected_part_f | abs(field_w)<max(abs(field_w(:)))/1e3;
    if all(rejected_part)
        error('All the spectral components are rejected with the selected diffractive order.');
    end
    field_w(rejected_part) = 0;
    theta_out(rejected_part_f) = 0;

    if ~isequal(compressor_type,'Treacy-beta2')
        switch compressor_type
            case 'Treacy-r' % reflective grating pair
                propagation_distance = (separation*1e9)*sec(theta_out).*(1+cos(theta_in+theta_out)); % nm
                grating_phase = pi + m*2*pi*separation*tan(-theta_out)/grating_spacing;
            case 'Treacy-t' % transmissive grating pair
                propagation_distance = (separation*1e9)*sec(theta_out).*(1-cos(theta_in-theta_out)); % nm
                grating_phase = m*2*pi*separation*tan(-theta_out)/grating_spacing;
        end

        propagation_phase = 2*pi./wavelength.*propagation_distance;

        total_phase = zeros(size(propagation_phase));
        total_phase(~rejected_part_f) = mod(2*(propagation_phase(~rejected_part_f)+grating_phase(~rejected_part_f)),2*pi); % "2" accounts for back-and-forth propagation
                                                                                                                     % (~rejected_part) accounts for the light going into the compressor; if this is ignored, "mod" will give error when the input argument isn't real
    else % 'beta2'
        c = 299792.458; % nm/ps
        omega = 2*pi*c./wavelength; % 2*pi*THz
        omega0 = feval(@(x) x(1), ifftshift(omega,1));
        phi2 = -m^2*wavelength0^3*(separation*1e9)/(2*pi*c^2*(grating_spacing*1e9)^2)*(1-(-m*wavelength0/(grating_spacing*1e9)-sin(theta_in))^2)^(-1.5);
        total_phase = (phi2/2*(omega-omega0).^2)*2; % "2", at the end, accounts for back-and-forth propagation
    end

    % Propagate the light through the grating and transform it back to time domain
    field = fft( ifftshift(field_w.*exp(1i*total_phase),1) );

    y = separation*tan(-theta_out);
else
    field = fft( ifftshift(field_w,1) );
    y = zeros(size(theta_out));
    total_phase = zeros(size(theta_out));
end

% Shift the pulse to the center of the time window
[~,pulse_center] = max(sum(abs(field).^2,2));
index_shift = pulse_center-floor(length(time)/2)-1;
field = double(circshift(field,-index_shift,1));

end

%%
function Frame = prism(separation,wavelength0,time,field,alpha,prism_material,tlim,wavelengthlim,t_feature,f_feature,lambda_or_f)
%PRISM

N = length(time); % the number of time points
dt = abs(time(2)-time(1)); % ps
c = 299792.458; % nm/ps
wavelength = c./((-N/2:N/2-1)'/N/dt + c/wavelength0); % nm
if size(time,1) == 1
    time = time';
end

field_w = fftshift(ifft(field),1); % frequency domain

% Find the center wavelength of the pulse spectrum with its second moment
f = c./wavelength; % THz
f_c = calc_f_c(f,abs(field_w).^2);
wavelength_c = c/f_c;
wavelength = c./(f+f_c-c/wavelength0);
field = field.*exp(1i*2*pi*(f_c-c/wavelength0)*time);
field_w = fftshift(ifft(field),1);

[a,b] = Sellmeier_coefficients(prism_material);
Sellmeier_terms = @(lambda,a,b) a.*lambda.^2./(lambda.^2 - b.^2);
n_from_Sellmeier = @(lambda) sqrt(1+sum(Sellmeier_terms(lambda,a,b),2));
n = n_from_Sellmeier(wavelength/1e3); % wavelength: nm

% Find the angle of incidence satisfying minimum deviation after a prism
n_c = n_from_Sellmeier(wavelength_c/1e3);
theta_in = asin(n_c*sin(alpha/2));

% Find the angles related to the first prism
theta2 = asin(sin(theta_in)./n);
theta3 = alpha - theta2;
theta4 = asin(n.*sin(theta3));

% Determine the beam width after the first prism to find the minimum prism size
I = sum(abs(field_w).^2,2);
considered_regime = I > max(I)/1e4 & ~arrayfun(@(x)any(imag(x)),theta4) & wavelength>0;
max_theta4 = max(theta4(considered_regime));

% The dechirped pulse
save_point = length(separation);
Frame(save_point) = struct('cdata',[],'colormap',[]);
for si = 1:save_point
    dechirped_field = prism_dechirping(separation(si),max_theta4,alpha,theta_in,theta2,theta3,theta4,n,wavelength,time,field_w);
    
    [~,~,~,fig,ax,cb] = calc_spectrogram(time,c./wavelength,dechirped_field,true,tlim,wavelengthlim,t_feature,f_feature,true,lambda_or_f);
    colormap(whitejet_lower(512)); set(cb,'Color','[0 0 0]');
    ax.NextPlot = 'replaceChildren';
    Frame(si) = getframe(fig);
    close(fig);
end

end

%%
function [field,y,y2,total_phase] = prism_dechirping(separation,max_theta4,alpha,theta_in,theta2,theta3,theta4,n,wavelength,time,field_w)
%PRISM_DECHIRPING It finds the field after propagating through the 
%compressor

if separation > 0
    % The light can go out of the prism pairs under certain angles or
    % wavelengths, so this part of light is rejected from the compressor.
    rejected_part_f = imag(theta4)~=0 | wavelength<0;
    rejected_part = rejected_part_f | abs(field_w)<max(abs(field_w(:)))/1e3;
    if all(rejected_part)
        error('All the spectral components are rejected with the selected diffractive order.');
    end
    field_w(rejected_part) = 0;
    theta2(rejected_part_f) = 0;
    theta3(rejected_part_f) = 0;
    theta4(rejected_part_f) = 0;
    
    dh = separation*(tan(theta4)*cos(alpha/2) - sin(alpha/2));
    max_dh = separation*(tan(max_theta4)*cos(alpha/2) - sin(alpha/2));
    y = max_dh - dh;
    y2 = y*sec(alpha/2).*cos(theta3)./cos(theta2);
    
    ls = separation*sec(theta4);
    lp = 2*y.*sin(alpha/2)./cos(theta2);
    lM = -y.*(cos(alpha)+sin(alpha)*tan(theta2))./cos(alpha/2)*sin(theta_in);
    
    propagation_distance = (ls + lp.*n + lM)*1e9; % nm
    propagation_phase = 2*pi./wavelength.*propagation_distance;
    total_phase = zeros(size(propagation_phase));
    total_phase(~rejected_part_f) = mod(2*(propagation_phase(~rejected_part_f)),2*pi); % "2" accounts for back-and-forth propagation
                                                                                       % (~rejected_part) accounts for the light going into the compressor; if this is ignored, "mod" will give error when the input argument isn't real
    total_phase(y<0) = 0; % this is rejected as well. Those y<0 means that the beam doesn't hit the second prism, passing over its apex
    
    % Propagate the light through the grating and transform it back to time domain
    field = fft( ifftshift(field_w.*exp(1i*total_phase),1) );
else
    field = fft( ifftshift(field_w,1) );
    y = zeros(size(theta4));
    y2 = zeros(size(theta4));
    total_phase = zeros(size(theta4));
end

% Shift the pulse to the center of the time window
[~,pulse_center] = max(sum(abs(field).^2,2));
index_shift = pulse_center-floor(length(time)/2)-1;
field = double(circshift(field,-index_shift,1));

end

%%
function Frame = grism(compressor_type,separation,theta_in,wavelength0,time,field,grating_spacing,alpha,prism_material,m,tlim,wavelengthlim,t_feature,f_feature,lambda_or_f)
%GRISM

N = length(time); % the number of time points
dt = abs(time(2)-time(1)); % ps
c = 299792.458; % nm/ps
wavelength = c./((-N/2:N/2-1)'/N/dt + c/wavelength0); % nm
if size(time,1) == 1
    time = time';
end

field_w = fftshift(ifft(field),1); % frequency domain

% Find the center wavelength of the pulse spectrum with its second moment
f = c./wavelength; % THz
f_c = calc_f_c(f,abs(field_w).^2);
wavelength_c = c/f_c;
wavelength = c./(f+f_c-c/wavelength0);
field = field.*exp(1i*2*pi*(f_c-c/wavelength0)*time);
field_w = fftshift(ifft(field),1);

[a,b] = Sellmeier_coefficients(prism_material);
Sellmeier_terms = @(lambda,a,b) a.*lambda.^2./(lambda.^2 - b.^2);
n_from_Sellmeier = @(lambda) sqrt(1+sum(Sellmeier_terms(lambda,a,b),2));
n = n_from_Sellmeier(wavelength/1e3); % wavelength: nm

% Find the angles related to the first prism
switch compressor_type
    case 'grism1'
        theta2 = asin((m*wavelength*1e-9/grating_spacing + sin(theta_in))./n);
        theta3 = alpha + theta2;
        theta4 = asin(n.*sin(theta3));
        
        beta = []; % dummy variable for the code to run
        theta5 = [];
        theta6 = [];
    case 'grism2'
        % Find the grating tilt angle satisfying minimum deviation after a prism
        n_c = n_from_Sellmeier(wavelength_c/1e3);
        theta2_c = asin(m*wavelength_c*1e-9/grating_spacing + sin(theta_in));
        theta3_c = asin(n_c*sin(alpha/2));
        beta = theta2_c + theta3_c - alpha/2; % grating tile angle

        theta2 = asin(m*wavelength*1e-9/grating_spacing + sin(theta_in));
        theta3 = alpha/2 + beta - theta2;
        theta4 = asin(sin(theta3)./n);
        theta5 = alpha - theta4;
        theta6 = asin(n.*sin(theta5));
end

% Determine the beam width after the first prism to find the minimum prism size
I = sum(abs(field_w).^2,2);
switch compressor_type
    case 'grism1'
        considered_max_theta = theta4;
    case 'grism2'
        considered_max_theta = theta6;
end
considered_regime = I > max(I)/1e4 & ~arrayfun(@(x)any(imag(x)),considered_max_theta) & wavelength>0;
max_theta = max(considered_max_theta(considered_regime));

% The dechirped pulse
save_point = length(separation);
Frame(save_point) = struct('cdata',[],'colormap',[]);
for si = 1:save_point
    dechirped_field = grism_dechirping(compressor_type,separation(si),max_theta,considered_regime,alpha,beta,theta_in,theta2,theta3,theta4,theta5,theta6,n,wavelength,time,field_w,grating_spacing,m);
    
    [~,~,~,fig,ax,cb] = calc_spectrogram(time,c./wavelength,dechirped_field,true,tlim,wavelengthlim,t_feature,f_feature,true,lambda_or_f);
    colormap(whitejet_lower(512)); set(cb,'Color','[0 0 0]');
    ax.NextPlot = 'replaceChildren';
    Frame(si) = getframe(fig);
    close(fig);
end

end

%%
function [field,xdh,RMs,total_phase] = grism_dechirping(compressor_type,separation,max_theta,considered_regime,alpha,beta,theta_in,theta2,theta3,theta4,theta5,theta6,n,wavelength,time,field_w,grating_spacing,m)
%GRISM_DECHIRPING It finds the field after propagating through the 
%compressor

if separation > 0
    % The light can go out of the prism pairs under certain angles or
    % wavelengths, so this part of light is rejected from the compressor.
    if all(~considered_regime)
        error('All the spectral components are rejected with the selected diffractive order.');
    end
    field_w(~considered_regime) = 0;
    theta2(~considered_regime) = 0;
    theta3(~considered_regime) = 0;
    theta4(~considered_regime) = 0;
    if isequal(compressor_type,'grism2')
        theta5(~considered_regime) = 0;
        theta6(~considered_regime) = 0;
    end
    
    switch compressor_type
        case 'grism1'
            dh = separation*(tan(theta4)*cos(alpha) - sin(alpha));
            x = separation*(tan(max_theta)*cos(alpha) - sin(alpha)); % = maximum dh
            xdh = x - dh;
            y1 = xdh*sec(alpha);
            y2 = y1.*cos(theta3)./cos(theta2);

            ls = separation*sec(theta4);
            lp = xdh./cos(theta2)*tan(alpha);
            lM = xdh.*(cos(alpha)-sin(alpha)*tan(theta2))*sec(alpha)*sin(theta_in);

            propagation_distance = (ls + lp.*n + lM)*1e9; % nm
            grating_phase = m*2*pi*y2/grating_spacing;
        case 'grism2'
            dh = separation*(tan(theta6)*cos(alpha/2) - sin(alpha/2));
            hx = separation*(tan(max_theta)*cos(alpha/2) - sin(alpha/2));
            xdh = hx - dh;
            y1 = xdh*sec(alpha/2);
            y2 = y1.*cos(theta5)./cos(theta4);
            y3 = max(y2(considered_regime));
            y4 = (y3-y2)./cos(theta2).*cos(theta3);
            
            ls = separation*sec(theta6);
            lp = y1*sin(alpha)./cos(theta4);
            lpg = (y3-y2)./cos(theta2)*sin(alpha/2+beta);
            lM = -y4*sin(theta_in);
            
            propagation_distance = (ls + lp.*n + lpg + lM)*1e9; % nm
            grating_phase = m*2*pi*y4/grating_spacing;
    end
    propagation_phase = 2*pi./wavelength.*propagation_distance;
    total_phase = zeros(size(propagation_phase));
    total_phase(considered_regime) = 2*(propagation_phase(considered_regime) + grating_phase(considered_regime)); % "2" accounts for back-and-forth propagation
                                                                                                                         % (~rejected_part) accounts for the light going into the compressor; if this is ignored, "mod" will give error when the input argument isn't real
    total_phase( xdh<0 ) = 0; % this is rejected as well. Those y<0 means that the beam doesn't hit the second prism, passing over its apex
    
    % Propagate the light through the grating and transform it back to time domain
    field = fft( ifftshift(field_w.*exp(1i*total_phase),1) );
    
    % for finding the roof mirror size
    switch compressor_type
        case 'grism1'
            RMs = y2*cos(theta_in);
        case 'grism2'
            RMs = y4*cos(theta_in);
    end
else
    field = fft( ifftshift(field_w,1) );
    xdh = zeros(size(theta2));
    RMs = zeros(size(theta2));
    total_phase = zeros(size(theta2));
end

% Shift the pulse to the center of the time window
[~,pulse_center] = max(sum(abs(field).^2,2));
index_shift = pulse_center-floor(length(time)/2)-1;
field = double(circshift(field,-index_shift,1));

end

%%
function Frame = Offner(compressor_type,separation,theta_in,wavelength0,time,field,grating_spacing,R,offcenter,m,tlim,wavelengthlim,t_feature,f_feature,lambda_or_f)
%OFFNER

N = length(time); % the number of time points
dt = abs(time(2)-time(1)); % ps
c = 299792.458; % nm/ps
wavelength = c./((-N/2:N/2-1)'/N/dt + c/wavelength0); % nm
if size(time,1) == 1
    time = time';
end

field_w = fftshift(ifft(field),1); % frequency domain

% Find the center wavelength of the pulse spectrum with its second moment
f = c./wavelength; % THz
f_c = calc_f_c(f,abs(field_w).^2);
wavelength = c./(f+f_c-c/wavelength0);
field = field.*exp(1i*2*pi*(f_c-c/wavelength0)*time);
field_w = fftshift(ifft(field),1);

theta_out = asin( m*wavelength/(grating_spacing*1e9) + sin(theta_in) ); % the transmitted/reflected angle of the m-th order diffraction

% The dechirped pulse
save_point = length(separation);
Frame(save_point) = struct('cdata',[],'colormap',[]);
for si = 1:save_point
    switch compressor_type(end) 
        case '1' % single-grating Offner compressor
            dechirped_field = Offner_dechirping_with_aberration(0,separation(si),theta_in,theta_out,wavelength,time,field_w,grating_spacing,R,m);
        case '2' % double-grating version for the single-grating case
            dechirped_field = Offner_dechirping_with_aberration(separation(si),offcenter,theta_in,theta_out,wavelength,time,field_w,grating_spacing,R,m);
        case '3' % true aberration-free double-grating Offner compressor
            dechirped_field = Offner_dechirping_no_aberration(separation(si),theta_in,theta_out,wavelength,field_w,grating_spacing,R,m);
    end
    
    [~,~,~,fig,ax,cb] = calc_spectrogram(time,c./wavelength,dechirped_field,true,tlim,wavelengthlim,t_feature,f_feature,true,lambda_or_f);
    colormap(whitejet_lower(512)); set(cb,'Color','[0 0 0]');
    ax.NextPlot = 'replaceChildren';
    Frame(si) = getframe(fig);
    close(fig);
end

end

%%
function [field,y,concave_leftmost,concave_rightmost,convex_size,total_phase] = Offner_dechirping_with_aberration(separation,offcenter,theta_in,theta_out,wavelength,time,field_w,grating_spacing,R,m)
%OFFNER_DECHIRPING_WITH_ABERRATION It finds the field after propagating 
%through the Offner compressor

if offcenter > -R && offcenter < R
    % The light can go out of the grating pairs under certain angles or
    % wavelengths, so this part of light is rejected from the compressor.
    rejected_part_f = imag(theta_out)~=0 | wavelength<0;
    rejected_part = rejected_part_f | abs(field_w)<max(abs(field_w(:)))/1e3;
    if all(rejected_part)
        error('All the spectral components are rejected with the selected diffractive order.');
    end
    field_w(rejected_part) = 0;
    theta_out(rejected_part_f) = 0;

    theta_diff = theta_in - theta_out;
    theta = asin(-offcenter/2/R*cos(theta_diff));
    psi = asin(-offcenter/R*cos(theta_diff));
    delta = psi - 2*theta;

    l1 = 2*R*(cos(theta)+tan(theta_diff).*sin(theta));
    l2 = R*(2*cos(theta)-cos(psi));
    tmp = 2*(l2.*sin(psi)-l1.*sin(delta));
    x = tmp./(cos(delta).*cot(delta-theta_out)-sin(delta));
    y = tmp./(cos(delta-theta_out)-tan(delta).*sin(delta-theta_out));

    tmp2 = 2*(delta+theta_diff-pi/2);
    theta1 = 2*theta_in-theta_out-pi/2-tmp2;
    theta2 = theta_diff-tmp2;
    delta_l = separation*csc(theta1).*cos(theta2);

    propagation_distance = 2*(l1+l2)-x-y.*sin(theta_in)+delta_l; % m
    % The difference between the transmissive and reflective grating pairs,
    % the only difference is an extra "pi" phase which doesn't affect the analysis here.
    grating_phase = -m*2*pi*y/grating_spacing;

    propagation_phase = 2*pi./wavelength.*propagation_distance*1e9;

    total_phase = zeros(size(propagation_phase));
    total_phase(~rejected_part_f) = mod(2*(propagation_phase(~rejected_part_f)+grating_phase(~rejected_part_f)),2*pi); % "2" accounts for back-and-forth propagation

    % Propagate the light through the grating and transform it back to time domain
    field = fft( ifftshift(field_w.*exp(1i*total_phase),1) );

    % The leftmost and rightmost positions on the concave mirror
    % Below are to compute the minimum size required for a concave mirror
    phi = theta_diff-pi/2-theta;
    phi2 = psi-theta+phi;
    phi3 = (psi-theta)*2+phi;
    concave_leftmost = 2*R*phi;
    concave_rightmost = 2*R*phi3;
    % Below is to compute the minimum size required for a convex mirror
    convex_size = R*phi2;
else
    field = fft(ifftshift(field_w,1));
    y = zeros(size(theta_out));
    concave_leftmost = zeros(size(theta_out));
    concave_rightmost = zeros(size(theta_out));
    convex_size = zeros(size(theta_out));
    total_phase = zeros(size(theta_out));
end

% Shift the pulse to the center of the time window
[~,pulse_center] = max(sum(abs(field).^2,2));
index_shift = pulse_center-floor(length(time)/2)-1;
field = double(circshift(field,-index_shift,1));

end

%%
function [field,total_phase] = Offner_dechirping_no_aberration(separation,theta_in,theta_out,wavelength,field_w,grating_spacing,R,m)
%OFFNER_DECHIRPING_NO_ABERRATION It finds the field after propagating 
%through the aberration-free double-grating Offner compressor

if separation > -2*R && separation < 2*R
    % The light can go out of the grating pairs under certain angles or
    % wavelengths, so this part of light is rejected from the compressor.
    rejected_part = imag(theta_out)~=0 | wavelength<0;
    if all(rejected_part)
        error('All the spectral components are rejected with the selected diffractive order.');
    end
    
    original_energy = sum(abs(field_w).^2,1);
    field_w(rejected_part) = 0;
    updated_energy = sum(abs(field_w).^2,1);
    theta_out(rejected_part) = 0;
    
    if abs(updated_energy/original_energy - 1) < 0.05
        x = separation*sec(-theta_out);
        phi = -theta_out - (pi/2 - theta_in);

        propagation_distance = -x.*(1+sin(phi)); % m
        % The difference between the transmissive and reflective grating pairs,
        % the only difference is an extra "pi" phase which doesn't affect the analysis here.
        grating_phase = -m*2*pi*(separation*tan(-theta_out))/grating_spacing;

        propagation_phase = 2*pi./wavelength.*propagation_distance*1e9;

        total_phase = zeros(size(propagation_phase));
        total_phase(~rejected_part) = mod(2*(propagation_phase(~rejected_part)+grating_phase(~rejected_part)),2*pi); % "2" accounts for back-and-forth propagation

        % Propagate the light through the grating and transform it back to time domain
        field = fft( ifftshift(field_w.*exp(1i*total_phase),1) );
    else
        field = fft( ifftshift(field_w,1) );
        total_phase = zeros(size(theta_out));
    end
else
    field = fft( ifftshift(field_w,1) );
    total_phase = zeros(size(theta_out));
end

% Shift the pulse to the center of the time window
[~,pulse_center] = max(sum(abs(field).^2,2));
index_shift = pulse_center-floor(length(wavelength)/2)-1;
field = double(circshift(field,-index_shift,1));

end

%%
function Frame = Martinez(separation,theta_in,wavelength0,time,field,grating_spacing,focal_length,tlim,wavelengthlim,t_feature,f_feature,lambda_or_f)
%Martinez

N = length(time); % the number of time points
dt = abs(time(2)-time(1)); % ps
c = 299792.458; % nm/ps
wavelength = c./((-N/2:N/2-1)'/N/dt + c/wavelength0); % nm
if size(time,1) == 1
    time = time';
end

field_w = fftshift(ifft(field),1); % frequency domain

% Find the center wavelength of the pulse spectrum with its second moment
f = c./wavelength; % THz
f_c = calc_f_c(f,abs(field_w).^2);
%wavelength_c = c/f_c;
wavelength = c./(f+f_c-c/wavelength0);
field = field.*exp(1i*2*pi*(f_c-c/wavelength0)*time);
field_w = fftshift(ifft(field),1);

theta_out = asin( m*wavelength/(grating_spacing*1e9) + sin(theta_in) ); % the transmitted/reflected angle of the m-th order diffraction

% The dechirped pulse
save_point = length(separation);
Frame(save_point) = struct('cdata',[],'colormap',[]);
for si = 1:save_point
    dechirped_field = Martinez_dechirping(separation(si),theta_in,theta_out,wavelength,field_w,grating_spacing,focal_length,m);
    
    [~,~,~,fig,ax,cb] = calc_spectrogram(time,c./wavelength,dechirped_field,true,tlim,wavelengthlim,t_feature,f_feature,true,lambda_or_f);
    colormap(whitejet_lower(512)); set(cb,'Color','[0 0 0]');
    ax.NextPlot = 'replaceChildren';
    Frame(si) = getframe(fig);
    close(fig);
end

end

%%
function [field,h_lens1,h_lens2,y,total_phase] = Martinez_dechirping(grating_lens_distance,theta_in,theta_out,wavelength,field_w,grating_spacing,focal_length,m)
%MARTINEZ_DECHIRPING It finds the field after propagating through the 
%Martinez stretcher

if grating_lens_distance > 0
    % The light can go out of the grating pairs under certain angles or
    % wavelengths, so this part of light is rejected from the compressor.
    rejected_part = imag(theta_out)~=0 | abs(real(theta_out))>pi/2*0.9 | wavelength<0;
    if all(rejected_part)
        error('All the spectral components are rejected with the selected diffractive order.');
    end
    
    original_energy = sum(abs(field_w).^2,1);
    field_w(rejected_part) = 0;
    updated_energy = sum(abs(field_w).^2,1);
    theta_out(rejected_part) = 0;
    theta0 = sum(theta_out.*abs(field_w).^2)/sum(abs(field_w).^2);
    
    if abs(updated_energy/original_energy - 1) < 0.05
        phi = theta0 - theta_out;
        
        h_lens1 = grating_lens_distance*tan(phi);
        [theta_out_lens1,beyond_lens,lens1_path_length] = lens_op(h_lens1,phi,focal_length);
        
        rejected_part = rejected_part | beyond_lens;
        
        h_lens2 = h_lens1 + 2*focal_length*tan(theta_out_lens1);
        [theta_out_lens2,~,lens2_path_length] = lens_op(h_lens2,theta_out_lens1,focal_length);
        
        phi2 = -theta_out_lens2; % this should be "phi" under paraxial approximation
        
        ll1 = h_lens2.*cot(phi2); ll1(isnan(ll1)) = 2*focal_length-grating_lens_distance;
        x = (ll1-grating_lens_distance)./cos(theta_out).*cos(theta_out+phi2);
        y = (ll1-grating_lens_distance)./cos(theta_out).*sin(phi2);

        ll2 = h_lens2.*csc(phi2); ll2(isnan(ll2)) = 2*focal_length-grating_lens_distance;
        propagation_distance = grating_lens_distance*sec(phi) + 2*focal_length*sec(theta_out_lens1) + ll2-x - y*sin(theta_in) + lens1_path_length + lens2_path_length; % m
        % The difference between the transmissive and reflective grating pairs,
        % the only difference is an extra "pi" phase which doesn't affect the analysis here.
        grating_phase = -m*2*pi*y/grating_spacing;

        propagation_phase = 2*pi./wavelength.*propagation_distance*1e9;

        total_phase = zeros(size(propagation_phase));
        total_phase(~rejected_part) = mod(2*(propagation_phase(~rejected_part)+grating_phase(~rejected_part)),2*pi); % "2" accounts for back-and-forth propagation

        % Propagate the light through the grating and transform it back to time domain
        field = fft( ifftshift(field_w.*exp(1i*total_phase),1) );
        
        % The leftmost and rightmost positions on the telescope
        % Below are to compute the minimum size required for the first lens
        % For the first lens, it's h_lens1.
        % For the second lens, it's h_lens2.
        % For the grating size, it's y.
    else
        error('All the spectral components are rejected with the selected diffractive order.');
    end
else % l can't be negative
    field = fft( ifftshift(field_w,1) );
    h_lens1 = zeros(size(theta_out));
    h_lens2 = zeros(size(theta_out));
    y = zeros(size(theta_out));
    total_phase = zeros(size(theta_out));
end

% Shift the pulse to the center of the time window
[~,pulse_center] = max(sum(abs(field).^2,2));
index_shift = pulse_center-floor(length(wavelength)/2)-1;
field = double(circshift(field,-index_shift,1));

end

%%
function f_c = calc_f_c(f,A2_f)

idx = f>0;
f = f(idx);
A2_f = A2_f(idx);

A2_f(A2_f<max(A2_f)/10) = 0; % it's necessary to kill the background

area = trapz(f,A2_f);
f_c = trapz(f,f.*A2_f)./area;

end

function [theta_out,beyond_lens,added_path_length] = lens_op(h_in,theta_in,f)
%LENS_OP It calculates the beam output angle after propagating through the
%lens

% lens is circular, so h can't be too large unless it's a fancy flat lens
% such as metalens.
n_lens = 1.45;
n_air = 1;
R = (n_lens - n_air)*f; % from lensmaker's formula
beyond_lens = abs(h_in) > R*0.8;

theta_out = theta_in - atan(h_in/f); % with paraxial approximation

% For Martinez stretcher, the lens thickness needs to be considered.
% When l=f, it should provide no dispersion at all. However, differernt
% colors still propagate with different path lengths. This is compensated
% by the lens thickness.
added_path_length = f - sqrt(h_in.^2 + f^2);

end