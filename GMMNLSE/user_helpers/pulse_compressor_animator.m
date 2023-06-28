function pulse_compressor_animator( compressor_type,optimal_separation,theta_in,wavelength0,time,field,grating_spacing,varargin )
%PULSE_COMPRESSOR_ANIMATOR Make an animation of spectrograms of the fields 
% after differernt grating separations of the compressor from zero to input
% optimal separation where the pulse should be fully dechirped.
%
%   compressor_type: 'Treacy-r': reflective grating pair,
%                    'Treacy-t': transmissive grating pair,
%                    'Treacy-beta2': consider only "beta2" term of a Treacy reflective grating compressor
%                                    (Please refer to Ch.6.2 in Applications of Nonlinear Fiber Optics (2ed), Agrawal)
%                    'Offner1': single-grating Offner compressor
%                               Offner1 is the Offner stretcher with one off-centered grating
%                    'Offner2': double-grating Offner compressor
%                               (Assume the off-center distance of the first grating is known and fixed.
%                                The grating separation is varied for pulse compression here.)
%                               Offner2 is the Offner stretcher with two gratings, one off-center and the other some distance away.
%   optimal_separation: grating spacing where the pulse should be fully dechirped
%   theta_in: a scalar; For a grating pair, the incident angle on the grating (rad)
%   wavelength0: a scalar; central wavelength (nm)
%   time: (N,1); the time grid points (ps)
%   field: (N,1); the electric field in time domain
%   grating_spacing: a scalar; the spacing between each grating line (m)
%
%   Extra required arguments for the Offner compressor:
%
%       R: radius of curvature of the convex mirror (R and 2R for two mirrors; 2R for the concave mirror) (m)
%       offcenter: the off-center distance of the first grating (m)
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
%      pulse_compressor_animator( compressor_type,optimal_separation,theta_in,wavelength0,time,field,grating_spacing,tlim,wavelengthlim,t_feature,f_feature,lambda_or_f,m )
%
%       * compressor_type is either 'Treacy-t', 'Treacy-r', or 'Treacy-beta2'
%       * [tlim,wavelengthlim,t_feature,f_feature,lambda_or_f,m] are optional inputs
%
%   Offner1:
%       pulse_compressor_animator( 'Offner1',optimal_separation,theta_in,wavelength0,time,field,grating_spacing,tlim,wavelengthlim,t_feature,f_feature,lambda_or_f,R,m )
%
%       * [tlim,wavelengthlim,t_feature,f_feature,lambda_or_f,m] are optional inputs
%
%   Offner2:
%       pulse_compressor_animator( 'Offner2',optimal_separation,theta_in,wavelength0,time,field,grating_spacing,tlim,wavelengthlim,t_feature,f_feature,lambda_or_f,R,offcenter,m )
%
%       * Compared to "Offner1", "offcenter" is required for an input here
%       * [tlim,wavelengthlim,t_feature,f_feature,lambda_or_f,m] are optional inputs

switch compressor_type
    case 'Offner1'
        R = varargin{1};
        
        if length(varargin) > 1
            varargin = varargin(2:end);
        end
    case 'Offner2'
        R = varargin{1};
        offcenter = varargin{2};
        
        if length(varargin) > 2
            varargin = varargin(3:end);
        end
end

% Default parameters
optargs = {[],[],50,50,true,-1};
% Load paramters
optargs(1:length(varargin)) = varargin;
[tlim,wavelengthlim,t_feature,f_feature,lambda_or_f,m] = optargs{:};

separation = linspace(0,optimal_separation,20);
switch compressor_type
    case {'Treacy-beta2','Treacy-r','Treacy-t'}
        Frame = Treacy(compressor_type,separation,theta_in,wavelength0,time,field,grating_spacing,            m,tlim,wavelengthlim,t_feature,f_feature,lambda_or_f);
    case 'Offner1'
        Frame = Offner(compressor_type,separation,theta_in,wavelength0,time,field,grating_spacing,R,0,        m,tlim,wavelengthlim,t_feature,f_feature,lambda_or_f);
    case 'Offner2'
        Frame = Offner(compressor_type,separation,theta_in,wavelength0,time,field,grating_spacing,R,offcenter,m,tlim,wavelengthlim,t_feature,f_feature,lambda_or_f);
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
    
    [~,~,~,fig,ax,cb] = calc_spectrogram(time,c./wavelength,dechirped_field,tlim,wavelengthlim,t_feature,f_feature,true,lambda_or_f);
    colormap(whitejet_lower(512)); set(cb,'Color','[0 0 0]');
    ax.NextPlot = 'replaceChildren';
    Frame(si) = getframe(fig);
    close(fig);
end

end

%%
function [field,y] = Treacy_dechirping(compressor_type,separation,theta_in,theta_out,wavelength,wavelength0,time,field_w,grating_spacing,m)
%DECHIRPING It finds the field after propagating through the compressor

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
            grating_phase = pi - m*2*pi*separation*tan(theta_out)/grating_spacing;
        case 'Treacy-t' % transmissive grating pair
            propagation_distance = (separation*1e9)*sec(theta_out).*(1-cos(theta_in-theta_out)); % nm
            grating_phase = -m*2*pi*separation*tan(theta_out)/grating_spacing;
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

% Shift the pulse to the center of the time window
[~,pulse_center] = max(sum(abs(field).^2,2));
index_shift = pulse_center-floor(length(time)/2);
field = double(circshift(field,-index_shift,1));

y = -separation*tan(theta_out);

end

%%
function Frame = Offner(compressor_type,separation,theta_in,wavelength0,time,field,grating_spacing,R,offcenter,m,tlim,wavelengthlim,t_feature,f_feature,lambda_or_f)
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
wavelength = c./(f+f_c-c/wavelength0);
field = field.*exp(1i*2*pi*(f_c-c/wavelength0)*time);
field_w = fftshift(ifft(field),1);

theta_out = asin( m*wavelength/(grating_spacing*1e9) + sin(theta_in) ); % the transmitted/reflected angle of the m-th order diffraction

% The dechirped pulse
save_point = length(separation);
Frame(save_point) = struct('cdata',[],'colormap',[]);
for si = 1:save_point
    if compressor_type(end) == '1' % single-grating
        dechirped_field = Offner_dechirping(0,separation(si),theta_in,theta_out,wavelength,time,field_w,grating_spacing,R,m);
    else % '2', double-grating
        dechirped_field = Offner_dechirping(separation(si),offcenter,theta_in,theta_out,wavelength,time,field_w,grating_spacing,R,m);
    end

    [~,~,~,fig,ax,cb] = calc_spectrogram(time,c./wavelength,dechirped_field,tlim,wavelengthlim,t_feature,f_feature,true,lambda_or_f);
    colormap(whitejet_lower(512)); set(cb,'Color','[0 0 0]');
    ax.NextPlot = 'replaceChildren';
    Frame(si) = getframe(fig);
    close(fig);
end

end

%%
function [field,y,concave_leftmost,concave_rightmost,convex_size] = Offner_dechirping(separation,offcenter,theta_in,theta_out,wavelength,time,field_w,grating_spacing,R,m)
%STRETCHING It finds the field after propagating through the compressor

% The light can go out of the grating pairs under certain angles or
% wavelengths, so this part of light is rejected from the compressor.
rejected_part_f = arrayfun(@(x)any(imag(x)),theta_out) | wavelength<0;
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

% Shift the pulse to the center of the time window
[~,pulse_center] = max(sum(abs(field).^2,2));
index_shift = pulse_center-floor(length(time)/2);
field = double(circshift(field,-index_shift,1));

% The leftmost and rightmost positions on the concave mirror
% Below are to compute the minimum size required for a concave mirror
phi = theta_diff-pi/2-theta;
phi2 = psi-theta+phi;
phi3 = (psi-theta)*2+phi;
concave_leftmost = 2*R*phi;
concave_rightmost = 2*R*phi3;
% Below is to compute the minimum size required for a convex mirror
convex_size = R*phi2;

end

function f_c = calc_f_c(f,A2_f)

idx = f>0;
f = f(idx);
A2_f = A2_f(idx);

A2_f(A2_f<max(A2_f)/10) = 0; % it's necessary to kill the background

area = trapz(f,A2_f);
f_c = trapz(f,f.*A2_f)./area;

end