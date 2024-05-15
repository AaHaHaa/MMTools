function [optimal_offcenter_l,stretched_field,grating_size,roof_mirror_size,size1,size2,recover_info] = pulse_stretcher_addNormalDispersion( stretcher_type,desired_duration,theta_in,wavelength0,time,field,grating_spacing,R_or_f,varargin )
%PULSE_STRETCHER_ADDNORMALDISPERSION Find the optimum distances and the 
%corresponding stretched field after the stretcher adds normal
%dispersion
%   
%   stretcher_type: 'single-Offner', 'double-Offner', or 'Martinez';
%                    single-Offner is the Offner stretcher with one off-centered grating
%                    double-Offner is the Offner stretcher with two gratings, one on-center and the other some distance away.
%                        It features aberration-free.
%                    Martinez is the common Martinez stretcher with two lenses
%   desired_duration: desired pulse duration (ps)
%   theta_in: a scalar; For a grating pair, the incident angle on the grating (rad);
%   wavelength0: a scalar; central wavelength (nm)
%   time: (N,1); the time grid points (ps)
%   field: (N,1); the electric field in time domain
%   grating_spacing: a scalar; the spacing between each grating line (m)
%   R_or_f: For Offner type -> radius of curvature of the convex mirror (R and 2R for two mirrors) (m)
%           For Martinez -> focal length of the lenses
%
%   Optional arguments:
%
%       verbose: display the result in the command window;
%                true(1) or false(0) (default: false)
%       global_opt: turn on global optimization to find the optimal compression;
%                   This is necessary sometimes if the pulse is messy.
%                   true(1) or false(0) (default: false)
%       m: a scalar; the diffraction order (default: -1)
%
%   Output arguments:
%
%       optimal_offcenter_l: For Offner type -> the optimal offcenter distance of the grating
%                            For Martinez -> the distance between the grating and the lens
%s      stretched_field: the stretched field
%       grating_size: the size of the grating (m)
%       roof_mirror_size: the size of the roof_mirror (m)
%       size1: For Offner type -> the size of the concave mirror
%              For Martinez -> the size of the lens 1 (closer to the input)
%       size2: For Offner type -> the size of the convex mirror
%              For Martinez -> the size of the lens 2 (closer to the roof mirror)
% =========================================================================
% Use:
%    Offner type:
%      [optimal_offcenter,stretched_field,grating_size,roof_mirror_size,concave_size,convex_size,recover_info] = pulse_stretcher( stretcher_type,desired_duration,theta_in,wavelength0,time,field,grating_spacing,R )
%      [optimal_offcenter,stretched_field,grating_size,roof_mirror_size,concave_size,convex_size,recover_info] = pulse_stretcher( stretcher_type,desired_duration,theta_in,wavelength0,time,field,grating_spacing,R,verbose )
%      [optimal_offcenter,stretched_field,grating_size,roof_mirror_size,concave_size,convex_size,recover_info] = pulse_stretcher( stretcher_type,desired_duration,theta_in,wavelength0,time,field,grating_spacing,R,verbose,global_opt )
%      [optimal_offcenter,stretched_field,grating_size,roof_mirror_size,concave_size,convex_size,recover_info] = pulse_stretcher( stretcher_type,desired_duration,theta_in,wavelength0,time,field,grating_spacing,R,verbose,global_opt,m )
%    Martinez type:
%      [optimal_l,stretched_field,grating_size,roof_mirror_size,lens1_size,lens2_size,recover_info] = pulse_stretcher( 'Martinez',desired_duration,theta_in,wavelength0,time,field,grating_spacing,focal_length )
%      [optimal_l,stretched_field,grating_size,roof_mirror_size,lens1_size,lens2_size,recover_info] = pulse_stretcher( 'Martinez',desired_duration,theta_in,wavelength0,time,field,grating_spacing,focal_length,verbose )
%      [optimal_l,stretched_field,grating_size,roof_mirror_size,lens1_size,lens2_size,recover_info] = pulse_stretcher( 'Martinez',desired_duration,theta_in,wavelength0,time,field,grating_spacing,focal_length,verbose,global_opt )
%      [optimal_l,stretched_field,grating_size,roof_mirror_size,lens1_size,lens2_size,recover_info] = pulse_stretcher( 'Martinez',desired_duration,theta_in,wavelength0,time,field,grating_spacing,focal_length,verbose,global_opt,m )
%

if ~ismember(stretcher_type,{'single-Offner','double-Offner','Martinez'})
    error('The value of stretcher_type is wrong.');
end

% Default parameters
optargs = {false false -1};
% Load paramters
optargs(1:length(varargin)) = varargin;
[verbose,global_opt,m] = optargs{:};
[optimal_offcenter_l,stretched_field,grating_size,roof_mirror_size,size1,size2,recover_info] = grating_pair(stretcher_type,desired_duration,theta_in,wavelength0,time,field,grating_spacing,R_or_f,verbose,m,global_opt);

end

%%
function [optimal_value,stretched_field,grating_size,roof_mirror_size,size1,size2,recover_info] = grating_pair(stretcher_type,desired_duration,theta_in,wavelength0,time,field,grating_spacing,R_or_f,verbose,m,global_opt)
%GRATING_PAIR

N = length(time); % the number of time points
dt = abs(time(2)-time(1)); % ps
c = 299792.458; % nm/ps
wavelength = c./((-N/2:N/2-1)'/N/dt + c/wavelength0); % nm
if size(time,1) == 1
    time = time'; % make sure it's a column vector
end

field_w = fftshift(ifft(field),1); % frequency domain

% Find the center wavelength of the pulse spectrum with its second moment
freq = c./wavelength; % THz
freq_c = calc_f_c(freq,abs(field_w).^2);
wavelength_c = c/freq_c;
wavelength = c./(freq+freq_c-c/wavelength0);
field = field.*exp(1i*2*pi*(freq_c-c/wavelength0)*time);
field_w = fftshift(ifft(field),1); % recompute the spectrum

theta_out = asin( m*wavelength/(grating_spacing*1e9) + sin(theta_in) ); % the transmitted angle of the m-th order diffraction

find_FWHM = @(s) calc_FWHM(stretcher_type,s,theta_in,theta_out,wavelength,time,field_w,grating_spacing,R_or_f,m);
find_optimum_stretcher_distance = @(s) abs(find_FWHM(s) - desired_duration);

t_fwhm = calc_RMS(time,abs(field).^2)*(2*sqrt(log(2))); % assume Gaussian shape
omega_fwhm = calc_RMS(freq-freq_c,abs(field_w).^2)*2*pi*(2*sqrt(log(2))); % assume Gaussian shape
guess_GVD = min(-(desired_duration-t_fwhm)/omega_fwhm*1e6/2,0);

% Run the optimization process to find the optimal distance, offcenter for
% Offner stretchers or l for Martinez stretcher
% -------------------------------------------------------------------------
switch stretcher_type
    case 'single-Offner'
        min_distance = 0;
        max_distance = R_or_f*0.9; % if it's R, the transmissive grating hits the convex mirror
        initial_guess = min_distance; % no added dispersion at offcenter = 0
    case 'double-Offner'
        min_distance = 0;
        max_distance = R_or_f*2*0.9; % if it's 2*R, the transmissive grating hits the concave mirror
        initial_guess = min_distance; % no added dispersion at offcenter = 0
    case 'Martinez'
        min_distance = 0;
        max_distance = R_or_f;
        %initial_guess = max_distance; % no added dispersion at L=f
        initial_guess = guess_GVD/2/(m^2*wavelength_c^3/(2*pi*c^2*grating_spacing^2)*(1-(-m*(wavelength_c*1e-9)/grating_spacing-sin(theta_in))^2)^(-1.5)*1e-3) + R_or_f;
end

% This is the local optimization I used before which might be stuck at
% local minimum if the pulse isn't smooth enough.
option = optimset('TolX',1e-10);
%option = optimset('PlotFcns',@optimplotfval,'TolX',1e-10); % plot the process of optimization
[optimal_value1,feval1] = fminsearch(find_optimum_stretcher_distance,initial_guess,option);
if feval1/desired_duration > 0.001 % sometimes it fails to find the optimum, so try to run with random initial guess again
    [optimal_value1,feval1] = fminsearch(find_optimum_stretcher_distance,rand(1)*(max_distance-min_distance)+min_distance,option);
end
% Because fminsearch is an unconstrained method, I set the constrain below.
if optimal_value1 < min_distance
    optimal_value1 = min_distance;
elseif optimal_value1 > max_distance
    optimal_value1 = max_distance;
end
    
% Global optimization
if global_opt && abs(feval1/desired_duration) > 0.001
    problem = createOptimProblem('fmincon',...
        'objective',find_optimum_stretcher_distance,...
        'x0',initial_guess,... % try with a different initial value
        'lb',min_distance,'ub',max_distance,...
        'options',optimoptions(@fmincon,'Display','off','TolX',1e-10));
    gs = GlobalSearch('MaxTime',60,'NumTrialPoints',130,'NumStageOnePoints',20,'MaxWaitCycle',5,'Display','off');
    %gs = MultiStart('MaxTime',120,'Display','off','UseParallel',true,'StartPointsToRun','bounds-ineqs');
    %[optimal_offcenter2,feval2] = run(gs,problem,20);
    [optimal_value2,feval2] = run(gs,problem);
else
    optimal_value2 = 0;
    feval2 = inf;
end

if feval1 < feval2
    optimal_value = optimal_value1;
else
    optimal_value = optimal_value2;
end

% The final stretched pulse
switch stretcher_type
    case 'single-Offner'
        [stretched_field,y,concave_leftmost,concave_rightmost,convex_size,added_phase] = single_Offner(0,optimal_value,theta_in,theta_out,wavelength,field_w,grating_spacing,R_or_f,m);
    case 'double-Offner'
        [stretched_field,added_phase] = double_Offner(optimal_value,theta_in,theta_out,wavelength,field_w,grating_spacing,R_or_f,m);
    case 'Martinez'
        [stretched_field,lens1,lens2,grating,added_phase] = Martinez(optimal_value,theta_in,theta_out,wavelength,field_w,grating_spacing,R_or_f,m);
end

tol_range = 1e-3;
if optimal_value < (min_distance+tol_range) || optimal_value > (max_distance-tol_range)
    warning('The distance of the grating stretcher is too close to the edge of the distance range used in the optimization.');
end
if verbose
    stretched_FWHM = calc_FWHM(stretcher_type,optimal_value,theta_in,theta_out,wavelength,time,field_w,grating_spacing,R_or_f,m);
    fprintf('Pulse duration after the stretcher = %6.4f(ps)\n',stretched_FWHM);
end

stretched_field_w = fftshift(ifft(stretched_field),1);

considered_regime = abs(stretched_field_w) > max(abs(stretched_field_w))/1e2 & imag(theta_out)==0 & wavelength>0;
switch stretcher_type
    case 'single-Offner'
        % Minimum size required for the concave mirror
        concave1_leftmost = min(concave_leftmost(considered_regime));
        concave1_rightmost = max(concave_leftmost(considered_regime));
        concave2_leftmost = min(concave_rightmost(considered_regime));
        concave2_rightmost = max(concave_rightmost(considered_regime));
        concave_size = {[concave1_leftmost,concave1_rightmost],...
                        [concave2_leftmost,concave2_rightmost]};
        % Minimum size required for the concave mirror
        convex_size = max(convex_size(considered_regime)) - min(convex_size(considered_regime));
        % Minimum size required for the grating
        grating_leftmost = min(y(considered_regime));
        grating_rightmost = max(y(considered_regime));
        grating_size = [grating_leftmost,grating_rightmost];
        roof_mirror_size = abs(diff(grating_size))*cos(theta_in);
    case 'double-Offner'
        grating_size = max(optimal_value*tan(theta_out(considered_regime))) - min(optimal_value*tan(theta_out(considered_regime)));
        concave_size = max(2*R_or_f*theta_out(considered_regime)) - min(2*R_or_f*theta_out(considered_regime));
        convex_size = concave_size/2;
        roof_mirror_size = grating_size*cos(theta_in);
    case 'Martinez'
        lens1_size = max(lens1(considered_regime)) - min(lens1(considered_regime));
        lens2_size = max(lens2(considered_regime)) - min(lens2(considered_regime));
        grating_size = max(grating(considered_regime)) - min(grating(considered_regime));
        roof_mirror_size = grating_size*cos(theta_in);
end

% Plot the result, if "verbose" is true.
if verbose
    switch stretcher_type
        case {'single-Offner','double-Offner'}
            show_result(stretcher_type,time,stretched_field,optimal_value,grating_size,roof_mirror_size,concave_size,convex_size);
        case 'Martinez'
            show_result(stretcher_type,time,stretched_field,optimal_value,grating_size,roof_mirror_size,lens1_size,lens2_size);
    end
end

% Recover back to the input frequency range
stretched_field = stretched_field.*exp(-1i*2*pi*(freq_c-c/wavelength0)*time);

switch stretcher_type
    case {'single-Offner','double-Offner'}
        size1 = concave_size;
        size2 = convex_size;
    case 'Martinez'
        size1 = lens1_size;
        size2 = lens2_size;
end

% Information to recover the field
recover_info = {exp(1i*2*pi*(freq_c-c/wavelength0)*time),exp(-1i*added_phase)};

end

%%
function [field,y,concave_leftmost,concave_rightmost,convex_size,total_phase] = single_Offner(separation,offcenter,theta_in,theta_out,wavelength,field_w,grating_spacing,R,m)
%SINGLE_OFFNER It finds the field after propagating through the 
%single-grating Offner stretcher

if offcenter >0 && offcenter < R
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
        total_phase(~rejected_part) = mod(2*(propagation_phase(~rejected_part)+grating_phase(~rejected_part)),2*pi); % "2" accounts for back-and-forth propagation

        % Propagate the light through the grating and transform it back to time domain
        field = fft( ifftshift(field_w.*exp(1i*total_phase),1) );

        % Shift the pulse temporally
        field0 = fft( ifftshift(field_w,1) ); field0(abs(field0)<max(abs(field0))/3) = 0; % remove noise for the original field
        % 1. Shift the pulse to where the previous peak is first:
        % This initial step is important in case the stretched pulse goes
        % beyond the time window and fails the intensity-weighted
        % computation.
        [~,pulse_center0] = max(abs(field0).^2);
        [~,pulse_center] = max(abs(field).^2);
        index_shift = round(pulse_center - pulse_center0);
        field = double(circshift(field,-index_shift,1));
        % 2. Then shift the pulse according to intensity-weighted pulse
        % center.
        % This assumes that the input doesn't have a field that wrap around
        % the time window.
        field1 = field; field1(abs(field1)<max(abs(field1))/3) = 0; % remove noise for the stretched field
        [~,pulse_center0] = calc_RMS((1:length(field0))',abs(field0).^2);
        [~,pulse_center] = calc_RMS((1:length(field1))',abs(field1).^2);
        index_shift = round(pulse_center - pulse_center0);
        field = circshift(field,-index_shift,1);

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
        field = fft( ifftshift(field_w,1) );
        y = zeros(size(wavelength));
        concave_leftmost = y;
        concave_rightmost = y;
        convex_size = y;
        total_phase = y;
    end
else
    field = fft( ifftshift(field_w,1) );
    y = zeros(size(wavelength));
    concave_leftmost = y;
    concave_rightmost = y;
    convex_size = y;
    total_phase = y;
end

end

%%
function [field,total_phase] = double_Offner(separation,theta_in,theta_out,wavelength,field_w,grating_spacing,R,m)
%DOUBLE_OFFNER It finds the field after propagating through the 
%aberration-free double-grating Offner stretcher

if separation > 0 && separation < 2*R
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

        % Shift the pulse temporally
        field0 = fft( ifftshift(field_w,1) ); field0(abs(field0)<max(abs(field0))/3) = 0; % remove noise for the original field
        % 1. Shift the pulse to where the previous peak is first:
        % This initial step is important in case the stretched pulse goes
        % beyond the time window and fails the intensity-weighted
        % computation.
        [~,pulse_center0] = max(abs(field0).^2);
        [~,pulse_center] = max(abs(field).^2);
        index_shift = round(pulse_center - pulse_center0);
        field = double(circshift(field,-index_shift,1));
        % 2. Then shift the pulse according to intensity-weighted pulse
        % center.
        % This assumes that the input doesn't have a field that wrap around
        % the time window.
        field1 = field; field1(abs(field1)<max(abs(field1))/3) = 0; % remove noise for the stretched field
        [~,pulse_center0] = calc_RMS((1:length(field0))',abs(field0).^2);
        [~,pulse_center] = calc_RMS((1:length(field1))',abs(field1).^2);
        index_shift = round(pulse_center - pulse_center0);
        field = circshift(field,-index_shift,1);
    else
        field = fft( ifftshift(field_w,1) );
        total_phase = zeros(size(wavelength));
    end
else
    field = fft( ifftshift(field_w,1) );
    total_phase = zeros(size(wavelength));
end

end

%%
function [field,h_lens1,h_lens2,y,total_phase] = Martinez(grating_lens_distance,theta_in,theta_out,wavelength,field_w,grating_spacing,focal_length,m)
%MARTINEZ It finds the field after propagating through the Martinez 
%stretcher

if grating_lens_distance > 0 && grating_lens_distance < focal_length
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

        % Shift the pulse temporally
        field0 = fft( ifftshift(field_w,1) ); field0(abs(field0)<max(abs(field0))/3) = 0; % remove noise for the original field
        % 1. Shift the pulse to where the previous peak is first:
        % This initial step is important in case the stretched pulse goes
        % beyond the time window and fails the intensity-weighted
        % computation.
        [~,pulse_center0] = max(abs(field0).^2);
        [~,pulse_center] = max(abs(field).^2);
        index_shift = round(pulse_center - pulse_center0);
        field = double(circshift(field,-index_shift,1));
        % 2. Then shift the pulse according to intensity-weighted pulse
        % center.
        % This assumes that the input doesn't have a field that wrap around
        % the time window.
        field1 = field; field1(abs(field1)<max(abs(field1))/3) = 0; % remove noise for the stretched field
        [~,pulse_center0] = calc_RMS((1:length(field0))',abs(field0).^2);
        [~,pulse_center] = calc_RMS((1:length(field1))',abs(field1).^2);
        index_shift = round(pulse_center - pulse_center0);
        field = circshift(field,-index_shift,1);
        
        % The leftmost and rightmost positions on the telescope
        % Below are to compute the minimum size required for the first lens
        % For the first lens, it's h_lens1.
        % For the second lens, it's h_lens2.
        % For the grating size, it's y.
    else
        error('All the spectral components are rejected with the selected diffractive order.');
    end
else % grating_lens_distance can't be negative
    field = fft( ifftshift(field_w,1) );
    h_lens1 = [];
    h_lens2 = [];
    y = [];
    total_phase = zeros(size(wavelength));
end

end

%%
function FWHM = calc_FWHM(stretcher_type,offcenter,theta_in,theta_out,wavelength,time,field_w,grating_spacing,R,m)
%CALC_FWHM It finds the FWHM of the pulse

switch stretcher_type
    case 'single-Offner'
        field = single_Offner(0,offcenter,theta_in,theta_out,wavelength,field_w,grating_spacing,R,m);
    case 'double-Offner'
        field = double_Offner(offcenter,theta_in,theta_out,wavelength,field_w,grating_spacing,R,m);
    case 'Martinez'
        % For Martinez,
        % l = offcenter value here
        % f = R value here
        field = Martinez(offcenter,theta_in,theta_out,wavelength,field_w,grating_spacing,R,m);
end

threshold = max(abs(field).^2)/1.0001;
[~,~,FWHM,~] = findpeaks(abs(field).^2,time,'MinPeakHeight',threshold,'WidthReference','halfheight','MinPeakProminence',threshold/2);
try
FWHM = FWHM(1); % just in case the pulse is stretched so badly that it has many peaks and is thus outputed many FWHMs
catch
    disp('');
end
end

%%
function show_result(stretcher_type,time,stretched_field,optimal_value,grating_size,roof_mirror_size,size1,size2)
%SHOW_RESULT

switch stretcher_type
    case {'single-Offner','double-Offner'}
        concave_size = size1;
        convex_size = size2;
    case 'Martinez'
        lens1_size = size1;
        lens2_size = size2;
end

intensity = abs(stretched_field).^2;

% Plot only the central part of the time window
intensity_plot = intensity;
threshold_factor = 100;
intensity_plot(intensity<max(intensity)/threshold_factor) = 0;
left = find(intensity_plot~=0,1);
right = find(intensity_plot~=0,1,'last');
center = floor((left+right)/2);
span_factor = 2;
span = floor((right-left)/2)*span_factor;
left = max(1,floor(center-span));
right = min(length(time),ceil(center+span));

figure('Name','stretched pulse');
h = plot(time(left:right),intensity(left:right));
xlim([min(time(left:right)) max(time(left:right))]);
xlabel('Time (ps)');
ylabel('Power (W)');
title('Stretched pulse');

set(h,'linewidth',2);
set(gca,'fontsize',14);

switch stretcher_type
    case 'single-Offner'
        fprintf('Grating offcenter = %6.4f(cm)\n',optimal_value*100);
        fprintf('(1) Light position on the grating = %6.4f ~ %6.4f(cm)\n',grating_size(1)*100,grating_size(2)*100);
        fprintf('--> Minimum size = %6.4f(cm)\n',(grating_size(2)-grating_size(1))*100);
        fprintf('(2) Light position on the concave mirror 1 = %6.4f ~ %6.4f(cm)\n',concave_size{1}(1)*100,concave_size{1}(2)*100);
        fprintf('--> Minimum size = %6.4f(cm)\n',(concave_size{1}(2)-concave_size{1}(1))*100);
        fprintf('(3) Light position on the concave mirror 2 = %6.4f ~ %6.4f(cm)\n',concave_size{2}(1)*100,concave_size{2}(2)*100);
        fprintf('--> Minimum size = %6.4f(cm)\n',(concave_size{2}(2)-concave_size{2}(1))*100);
        fprintf('(4) Size of the convex mirror = %6.4f(cm)\n',convex_size*100);
        fprintf('(5) Size of the roof mirror = %6.4f(cm)\n',roof_mirror_size*100);
    case 'double-Offner'
        fprintf('Grating separation = %6.4f(cm)\n',optimal_value*100);
        fprintf('(1) Grating min. size = %6.4f(cm)\n',grating_size*100);
        fprintf('(2) Concave mirror min. size = %6.4f(cm)\n',concave_size*100);
        fprintf('(3) Convex mirror min. size = %6.4f(cm)\n',convex_size*100);
        fprintf('(4) Roof mirror min. size = %6.4f(cm)\n',roof_mirror_size*100);
    case 'Martinez'
        fprintf('Grating to lens = %6.4f(cm)\n',optimal_value*100);
        fprintf('(1) Grating min. size = %6.4f(cm)\n',grating_size*100);
        fprintf('(2) Lens 1 min. size = %6.4f(cm)\n',lens1_size*100);
        fprintf('(3) Lens 2 min. size = %6.4f(cm)\n',lens2_size*100);
        fprintf('(4) Roof mirror min. size = %6.4f(cm)\n',roof_mirror_size*100);
end

end

function f_c = calc_f_c(f,A2_f)
%CALC_F_C It finds the center frequency

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