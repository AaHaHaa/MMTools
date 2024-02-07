function [optimal_value,dechirped_FWHM,dechirped_field,varargout] = pulse_compressor( compressor_type,theta_in,wavelength0,time,field,varargin )
%PULSE_COMPRESSOR Find the optimum grating separation and the corresponding 
%dechirped field after different types of compressors
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
%   theta_in: a scalar; For a grating pair, the incident angle on the grating (rad);
%                       For a prism pair, it's also the Brewster angle (rad)
%   wavelength0: a scalar; central wavelength (nm)
%   time: (N,1); the time grid points (ps)
%   field: (N,1); the electric field in time domain
%   grating_spacing: a scalar; the spacing between each grating line (m)
%   alpha: the apex angle of prisms in a prism compressor (rad)
%   prism_material: 'N-SF10'
%
%   Extra required arguments for the Offner compressor:
%
%       R: radius of curvature of the convex mirror (R and 2R for two mirrors; 2R for the concave mirror) (m)
%       offcenter: the off-center distance of the first grating (m);
%                  required only for 'Offner2'
%
%   Extra required arguments for the Martinez stretcher(used as a compressor):
%
%       focal_length: focal length of the lenses (m)
%
%   Optional arguments:
%
%       verbose: display the result in the command window;
%                true(1) or false(0) (default: false)
%       global_opt: turn on global optimization to find the optimal compression;
%                   This is necessary sometimes if the pulse is messy.
%                   true(1) or false(0) (default: false)
%       m: a scalar; the diffraction order (default: -1)
% =========================================================================
% Use:
%   Treacy:
%      [optimal_separation,dechirped_FWHM,dechirped_field,grating_size,roof_mirror_size,recover_info] = pulse_compressor( compressor_type,theta_in,wavelength0,time,field,grating_spacing,verbose,global_opt,m )
%
%       * compressor_type is either 'Treacy-t', 'Treacy-r', or 'Treacy-beta2'
%       * grating_size and roof_mirror_size are optional outputs
%       * verbose, global_opt, and m are optional inputs
%
%   Prism:
%      [optimal_separation,dechirped_FWHM,dechirped_field,prism_height,roof_mirror_size,recover_info] = pulse_compressor( 'prism',[],wavelength0,time,field,alpha,prism_material,verbose,global_opt )
%
%       * prism_height and roof_mirror_size are optional outputs
%       * verbose and global_opt are optional inputs
%
%   Offner1:
%       [optimal_offcenter,dechirped_FWHM,dechirped_field,grating_size,roof_mirror_size,concave_size,convex_size,recover_info] = pulse_compressor( 'Offner1',theta_in,wavelength0,time,field,grating_spacing,R,verbose,global_opt,m )
%
%       * grating_size, roof_mirror_size, concave_size, convex_size are optional outputs
%       * verbose, global_opt, and m are optional inputs
%
%   Offner2:
%       [optimal_separation,dechirped_FWHM,dechirped_field,grating_size,roof_mirror_size,concave_size,convex_size,recover_info] = pulse_compressor( 'Offner2',theta_in,wavelength0,time,field,grating_spacing,R,offcenter,verbose,global_opt,m )
%
%       * Compared to "Offner1", "offcenter" is required for an input here
%       * grating_size, roof_mirror_size, concave_size, convex_size are optional outputs
%       * verbose, global_opt, and m are optional inputs
%
%   Offner3:
%       [optimal_separation,dechirped_FWHM,dechirped_field,grating_size,roof_mirror_size,concave_size,convex_size,recover_info] = pulse_compressor( 'Offner3',theta_in,wavelength0,time,field,grating_spacing,R,verbose,global_opt,m )
%
%       * grating_size, roof_mirror_size, concave_size, convex_size are optional outputs
%       * verbose, global_opt, and m are optional inputs
%
%   Martinez:
%       [optimal_separation,dechirped_FWHM,dechirped_field,grating_size,roof_mirror_size,lens1_size,lens2_size,recover_info] = pulse_compressor( 'Martinez',theta_in,wavelength0,time,field,grating_spacing,focal_length,verbose,global_opt,m )
%
%       * grating_size, roof_mirror_size, lens1_size, lens2_size are optional outputs
%       * verbose, global_opt, and m are optional inputs

switch compressor_type
    case {'Treacy-r','Treacy-t','Treacy-beta2'}
        grating_spacing = varargin{1};
        
        n = 1;
    case 'prism'
        alpha = varargin{1};
        prism_material = varargin{2};
        
        n = 2;
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
end
if length(varargin) > n % optional input arguments
    varargin = varargin(n+1:end);
else
    varargin = {};
end

% Default parameters
optargs = {false false -1};
% Load paramters
optargs(1:length(varargin)) = varargin;
[verbose,global_opt,m] = optargs{:};
switch compressor_type
    case {'Treacy-beta2','Treacy-r','Treacy-t'}
        [optimal_separation,dechirped_FWHM,dechirped_field,grating_size,roof_mirror_size,recover_info] = Treacy(compressor_type,theta_in,wavelength0,time,field,grating_spacing,verbose,m,global_opt);
        optimal_value = optimal_separation;
        varargout = {grating_size,roof_mirror_size,recover_info};
    case 'prism'
        [optimal_separation,dechirped_FWHM,dechirped_field,prism_height,roof_mirror_size,recover_info] = prism(wavelength0,time,field,alpha,prism_material,verbose,global_opt);
        optimal_value = optimal_separation;
        varargout = {prism_height,roof_mirror_size,recover_info};
    case 'Offner1'
        [optimal_offcenter,dechirped_FWHM,dechirped_field,grating_size,roof_mirror_size,concave_size,convex_size,recover_info] = Offner(compressor_type,theta_in,wavelength0,time,field,grating_spacing,R,0,         verbose,m,global_opt);
        optimal_value = optimal_offcenter;
        varargout = {grating_size,roof_mirror_size,concave_size,convex_size,recover_info};
    case 'Offner2'
        [optimal_separation,dechirped_FWHM,dechirped_field,grating_size,roof_mirror_size,concave_size,convex_size,recover_info] = Offner(compressor_type,theta_in,wavelength0,time,field,grating_spacing,R,offcenter,verbose,m,global_opt);
        optimal_value = optimal_separation;
        varargout = {grating_size,roof_mirror_size,concave_size,convex_size,recover_info};
    case 'Offner3'
        [optimal_separation,dechirped_FWHM,dechirped_field,grating_size,roof_mirror_size,concave_size,convex_size,recover_info] = Offner(compressor_type,theta_in,wavelength0,time,field,grating_spacing,R,0,        verbose,m,global_opt);
        optimal_value = optimal_separation;
        varargout = {grating_size,roof_mirror_size,concave_size,convex_size,recover_info};
    case 'Martinez'
        [optimal_separation,dechirped_FWHM,dechirped_field,grating_size,roof_mirror_size,lens1_size,lens2_size,recover_info] = Martinez(theta_in,wavelength0,time,field,grating_spacing,focal_length,verbose,m,global_opt);
        optimal_value = optimal_separation;
        varargout = {grating_size,roof_mirror_size,lens1_size,lens2_size,recover_info};
    otherwise
        error('The value of compressor_type is wrong.');
end

end

%%
function [optimal_separation,dechirped_FWHM,dechirped_field,grating_size,roof_mirror_size,recover_info] = Treacy(compressor_type,theta_in,wavelength0,time,field,grating_spacing,verbose,m,global_opt)
%TREACY

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

find_optimum_compressor_separation = @(d) -Treacy_calc_max_intensity(compressor_type,d,theta_in,theta_out,wavelength,wavelength_c,time,field_w,grating_spacing,m);
find_FWHM = @(d) Treacy_calc_FWHM(compressor_type,d,theta_in,theta_out,wavelength,wavelength_c,time,field_w,grating_spacing,m);

% Run the global optimization process to find the optimal grating separation
% -------------------------------------------------------------------------
min_separation = 0;
GVD_before_dechirping = characterize_spectral_phase(c./wavelength,fftshift(ifft(ifftshift(field(:,1),1)),1),3,false); % consider only the 1st mode
% Obtain the initial guess from the beta2 of a grating compressor.
% This number should be close to the answer already.
initial_guess = GVD_before_dechirping/2/(m^2*wavelength_c^3/(2*pi*c^2*grating_spacing^2)*(1-(-m*(wavelength_c*1e-9)/grating_spacing-sin(theta_in))^2)^(-1.5)*1e-3);
if initial_guess < min_separation
    initial_guess = 0;
elseif ~isreal(initial_guess) % in case that the phase is a mess
    initial_guess = 0.1;
end

% This is the local optimization I used before which might be stuck at
% local minimum if the pulse isn't smooth enough.
option = optimset('TolX',1e-20);
%option = optimset('PlotFcns',@optimplotfval,'TolX',1e-20); % plot the process of optimization
[optimal_separation1,feval1] = fminsearch(find_optimum_compressor_separation,initial_guess,option);
% Because fminsearch is unconstrained method, I set the constrain below.
if optimal_separation1 < min_separation
    optimal_separation1 = min_separation;
end
% Global optimization
optimal_separation2 = 0;    feval2 = 0;
if global_opt
    transform_limited_field = calc_transform_limited( field );
    transform_limited_peak_power = max(abs(transform_limited_field).^2,[],1);
    dechirped_field = Treacy_dechirping(compressor_type,optimal_separation1,theta_in,theta_out,wavelength,wavelength_c,time,field_w,grating_spacing,m);
    dechirped_peak_power = max(abs(dechirped_field).^2,[],1);
    
    % Apply global optimization only when the dechirped field isn't good enough
    if dechirped_peak_power < transform_limited_peak_power*0.5
        problem = createOptimProblem('fmincon',...
            'objective',find_optimum_compressor_separation,...
            'x0',initial_guess,...
            'lb',min_separation,...
            'options',optimoptions(@fmincon,'Display','off','TolX',1e-20));
        gs = GlobalSearch('MaxTime',60,'NumTrialPoints',130,'NumStageOnePoints',20,'MaxWaitCycle',5,'Display','off');
        %gs = MultiStart('MaxTime',120,'Display','off','UseParallel',true,'StartPointsToRun','bounds-ineqs');
        %[optimal_separation2,feval2] = run(gs,problem,20);
        [optimal_separation2,feval2] = run(gs,problem);
    end
end

if feval1 < feval2
    optimal_separation = optimal_separation1;
else
    optimal_separation = optimal_separation2;
end

% The final dechirped pulse
[dechirped_field,y,added_phase] = Treacy_dechirping(compressor_type,optimal_separation,theta_in,theta_out,wavelength,wavelength_c,time,field_w,grating_spacing,m);
dechirped_FWHM = find_FWHM(optimal_separation); % ps
dechirped_FWHM = dechirped_FWHM*1e3; % fs

tol_range = 1e-3; % 1 mm
if optimal_separation < (min_separation+tol_range)
    warning('The distance of the grating compressor is too small.');
    
    if ~global_opt
        fprintf('Try running with global optimization.\n');
    end
end

dechirped_field_w = fftshift(ifft(dechirped_field),1);

% Minimum size required for the grating and the mirror
I = sum(abs(dechirped_field_w).^2,2);
considered_regime = I > max(I)/1e4 & ~arrayfun(@(x)any(imag(x)),theta_out) & wavelength>0;
if ~isequal(compressor_type,'Treacy-beta2')
    grating_size = max(y(considered_regime)) - min(y(considered_regime));
    roof_mirror_size = grating_size*cos(theta_in);
else
    grating_size = 0;
    roof_mirror_size = 0;
end

% Plot the result, if "verbose" is true.
if verbose
    show_result(compressor_type,time,dechirped_field,optimal_separation,0,dechirped_FWHM,grating_size,roof_mirror_size);
end

% Recover back to the input frequency range
dechirped_field = dechirped_field.*exp(-1i*2*pi*(f_c-c/wavelength0)*time);

% Information to recover the field
recover_info = {exp(1i*2*pi*(f_c-c/wavelength0)*time),exp(-1i*added_phase)};

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
index_shift = pulse_center-floor(length(time)/2);
field = double(circshift(field,-index_shift,1));

end

%%
function [optimal_separation,dechirped_FWHM,dechirped_field,prism_height,roof_mirror_size,recover_info] = prism(wavelength0,time,field,alpha,prism_material,verbose,global_opt)
%TREACY

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

% Calculate the derivatives based on central difference
central_diff_coeff = {[-1/2;  0; 1/2],...
                      [   1; -2;   1]};
wavelength_halfrange = 1; % nm
n_find_derivatives = n_from_Sellmeier((wavelength_c+[-wavelength_halfrange,0,wavelength_halfrange]')/1e3); % find the derivatives Dn and D2n from +/-1 nm around pulse ceter wavelength
Dn  = sum(n_find_derivatives.*central_diff_coeff{1})/(wavelength_halfrange*1e-9);
Dn2 = sum(n_find_derivatives.*central_diff_coeff{2})/(wavelength_halfrange*1e-9)^2;
% Find theta4 for the center wavelength
theta2_c = asin(sin(theta_in)./n_c);
theta3_c = alpha - theta2_c;
theta4_c = asin(n_c.*sin(theta3_c));

find_optimum_compressor_separation1 = @(d) -prism_calc_max_intensity(d,[],max_theta4,alpha,theta_in,theta2,theta3,theta4,n,wavelength,time,field_w);
find_FWHM = @(d,dh) prism_calc_FWHM(d,dh,max_theta4,alpha,theta_in,theta2,theta3,theta4,n,wavelength,time,field_w);

% Run the global optimization process to find the optimal prism separation
% -------------------------------------------------------------------------
min_separation = 0;
GVD_before_dechirping = characterize_spectral_phase(c./wavelength,fftshift(ifft(ifftshift(field(:,1),1)),1),3,false); % consider only the 1st mode
initial_guess = -GVD_before_dechirping*1e-30/(wavelength_c*1e-9)^3*(2*pi*(c*1e3)^2)/(4*sec(max_theta4))/((Dn2+(2*n_c-1/n_c^3)*Dn^2)*sin(max_theta4-theta4_c)-2*Dn^2*cos(max_theta4-theta4_c));
if initial_guess < min_separation
    initial_guess = 0;
elseif ~isreal(initial_guess) % in case that the phase is a mess
    initial_guess = 0.1;
end

% This is the local optimization I used before which might be stuck at
% local minimum if the pulse isn't smooth enough.
option = optimset('TolX',1e-20);
%option = optimset('PlotFcns',@optimplotfval,'TolX',1e-20); % plot the process of optimization
[optimal_separation1,feval1] = fminsearch(find_optimum_compressor_separation1,initial_guess,option);
% Because fminsearch is unconstrained method, I set the constrain below.
if optimal_separation1 < min_separation
    optimal_separation1 = min_separation;
end
% Global optimization
optimal_separation2 = 0;    feval2 = 0;
if global_opt
    transform_limited_field = calc_transform_limited( field );
    transform_limited_peak_power = max(abs(transform_limited_field).^2,[],1);
    max_dh1 = optimal_separation1*(tan(max_theta4)*cos(alpha/2) - sin(alpha/2));
    dechirped_field = prism_dechirping(optimal_separation1,max_dh1,max_theta4,alpha,theta_in,theta2,theta3,theta4,n,wavelength,time,field_w);
    dechirped_peak_power = max(abs(dechirped_field).^2,[],1);
    
    % Apply global optimization only when the dechirped field isn't good enough
    if dechirped_peak_power < transform_limited_peak_power*0.5
        problem = createOptimProblem('fmincon',...
            'objective',find_optimum_compressor_separation1,...
            'x0',initial_guess,...
            'lb',min_separation,...
            'options',optimoptions(@fmincon,'Display','off','TolX',1e-20));
        gs = GlobalSearch('MaxTime',60,'NumTrialPoints',130,'NumStageOnePoints',20,'MaxWaitCycle',5,'Display','off');
        %gs = MultiStart('MaxTime',120,'Display','off','UseParallel',true,'StartPointsToRun','bounds-ineqs');
        %[optimal_separation2,feval2] = run(gs,problem,20);
        [optimal_separation2,feval2] = run(gs,problem);
    end
end

if feval1 < feval2
    optimal_separation = optimal_separation1;
else
    optimal_separation = optimal_separation2;
end

% The final dechirped pulse
max_dh = optimal_separation*(tan(max_theta4)*cos(alpha/2) - sin(alpha/2));
[dechirped_field,y,y2,added_phase] = prism_dechirping(optimal_separation,max_dh,max_theta4,alpha,theta_in,theta2,theta3,theta4,n,wavelength,time,field_w);
dechirped_FWHM = find_FWHM(optimal_separation,max_dh); % ps
dechirped_FWHM = dechirped_FWHM*1e3; % fs

tol_range = 1e-3; % 1 mm
if optimal_separation < (min_separation+tol_range)
    warning('The distance of the grating compressor is too small.');
    
    if ~global_opt
        fprintf('Try running with global optimization.\n');
    end
end

% Minimum size required for the prism and the mirror
prism_height = max(y(considered_regime)) - min(y(considered_regime));
roof_mirror_size = (max(y2(considered_regime)) - min(y2(considered_regime)))*cos(theta_in);

% Plot the result, if "verbose" is true.
if verbose
    show_result('prism',time,dechirped_field,optimal_separation,0,dechirped_FWHM,0,roof_mirror_size,prism_height);
end

% Recover back to the input frequency range
dechirped_field = dechirped_field.*exp(-1i*2*pi*(f_c-c/wavelength0)*time);

% Information to recover the field
recover_info = {exp(1i*2*pi*(f_c-c/wavelength0)*time),exp(-1i*added_phase)};

end

%%
function [field,y,y2,total_phase] = prism_dechirping(separation,max_dh,max_theta4,alpha,theta_in,theta2,theta3,theta4,n,wavelength,time,field_w)
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
    if isempty(max_dh)
        max_dh = separation*(tan(max_theta4)*cos(alpha/2) - sin(alpha/2));
    end
    y = max_dh - dh;
    y2 = y*sec(alpha/2).*cos(theta3)./cos(theta2);
    
    ls = separation*sec(theta4);
    lp = 2*y.*sin(alpha/2)./cos(theta2).*n;
    lM = -y.*(cos(alpha)+sin(alpha)*tan(theta2))./cos(alpha/2)*sin(theta_in);
    
    propagation_distance = (ls + lp + lM)*1e9; % nm
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
index_shift = pulse_center-floor(length(time)/2);
field = double(circshift(field,-index_shift,1));

end

%%
function [optimal_value,dechirped_FWHM,dechirped_field,grating_size,roof_mirror_size,concave_size,convex_size,recover_info] = Offner(compressor_type,theta_in,wavelength0,time,field,grating_spacing,R,offcenter,verbose,m,global_opt)
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

switch compressor_type(end)
    case '1' % single-grating Offner compressor
        find_optimum_compressor_value = @(o) -Offner_calc_max_intensity(compressor_type,0,o,theta_in,theta_out,wavelength,time,field_w,grating_spacing,R,m);
        find_FWHM = @(o) Offner_calc_FWHM(compressor_type,0,o,theta_in,theta_out,wavelength,time,field_w,grating_spacing,R,m);
    case '2' % double-grating version for the single-grating case
        find_optimum_compressor_value = @(s) -Offner_calc_max_intensity(compressor_type,s,offcenter,theta_in,theta_out,wavelength,time,field_w,grating_spacing,R,m);
        find_FWHM = @(s) Offner_calc_FWHM(compressor_type,s,offcenter,theta_in,theta_out,wavelength,time,field_w,grating_spacing,R,m);
    case '3' % true aberration-free double-grating Offner compressor
        find_optimum_compressor_value = @(s) -Offner_calc_max_intensity(compressor_type,s,0,theta_in,theta_out,wavelength,time,field_w,grating_spacing,R,m);
        find_FWHM = @(s) Offner_calc_FWHM(compressor_type,s,0,theta_in,theta_out,wavelength,time,field_w,grating_spacing,R,m);
end

% Run the global optimization process to find the optimal grating separation
% -------------------------------------------------------------------------
switch compressor_type(end)
    case {'1','2'} % Offner compressor with aberration
        min_value = -R*0.9;
        max_value =  R*0.9;
    case '3' % Offner compressor without aberration
        min_value = -2*R*0.9;
        max_value =  2*R*0.9;
end
initial_guess = 0;

% This is the local optimization I used before which might be stuck at
% local minimum if the pulse isn't smooth enough.
option = optimset('TolX',1e-20);
%option = optimset('PlotFcns',@optimplotfval,'TolX',1e-20); % plot the process of optimization
[optimal_value1,feval1] = fminsearch(find_optimum_compressor_value,0,option);
% Because fminsearch is unconstrained method, I set the constrain below.
if optimal_value1 < min_value
    optimal_value1 = min_value;
elseif optimal_value1 > max_value
    optimal_value1 = max_value;
end
% Global optimization
optimal_value2 = 0;    feval2 = 0;
if global_opt
    transform_limited_field = calc_transform_limited( field );
    transform_limited_peak_power = max(abs(transform_limited_field).^2,[],1);
    switch compressor_type(end) 
        case '1' % single-grating Offner compressor
            dechirped_field = Offner_dechirping_with_aberration(0,optimal_value1,theta_in,theta_out,wavelength,time,field_w,grating_spacing,R,m);
        case '2' % double-grating version for the single-grating case
            dechirped_field = Offner_dechirping_with_aberration(optimal_value1,offcenter,theta_in,theta_out,wavelength,time,field_w,grating_spacing,R,m);
        case '3' % true aberration-free double-grating Offner compressor
            dechirped_field = Offner_dechirping_no_aberration(optimal_value1,theta_in,theta_out,wavelength,field_w,grating_spacing,R,m);
    end
    dechirped_peak_power = max(abs(dechirped_field).^2,[],1);
    
    % Apply global optimization only when the dechirped field isn't good enough
    if dechirped_peak_power < transform_limited_peak_power*0.5
        problem = createOptimProblem('fmincon',...
            'objective',find_optimum_compressor_value,...
            'x0',initial_guess,...
            'lb',min_value,'ub',max_value,...
            'options',optimoptions(@fmincon,'Display','off','TolX',1e-20));
        gs = GlobalSearch('MaxTime',60,'NumTrialPoints',130,'NumStageOnePoints',20,'MaxWaitCycle',5,'Display','off');
        %gs = MultiStart('MaxTime',120,'Display','off','UseParallel',true,'StartPointsToRun','bounds-ineqs');
        %[optimal_separation2,feval2] = run(gs,problem,20);
        [optimal_value2,feval2] = run(gs,problem);
    end
end

if feval1 < feval2
    optimal_value = optimal_value1;
else
    optimal_value = optimal_value2;
end

% The final dechirped pulse
switch compressor_type(end) 
    case '1' % single-grating Offner compressor
        [dechirped_field,y,concave_leftmost,concave_rightmost,convex_size,added_phase] = Offner_dechirping_with_aberration(0,optimal_value,theta_in,theta_out,wavelength,time,field_w,grating_spacing,R,m);
    case '2' % double-grating version for the single-grating case
        [dechirped_field,y,concave_leftmost,concave_rightmost,convex_size,added_phase] = Offner_dechirping_with_aberration(optimal_value,offcenter,theta_in,theta_out,wavelength,time,field_w,grating_spacing,R,m);
    case '3' % true aberration-free double-grating Offner compressor
        [dechirped_field,added_phase] = Offner_dechirping_no_aberration(optimal_value,theta_in,theta_out,wavelength,field_w,grating_spacing,R,m);
end
dechirped_FWHM = find_FWHM(optimal_value); % ps
dechirped_FWHM = dechirped_FWHM*1e3; % fs

tol_range = 1e-3; % 1 mm
if optimal_value < (min_value+tol_range) || optimal_value > (max_value-tol_range)
    warning('The distance of the grating compressor is too close to the edge of the distance range used in the optimization.');
    
    if ~global_opt
        fprintf('Try running with global optimization.\n');
    end
end

dechirped_field_w = fftshift(ifft(dechirped_field),1);

% Minimum size required for the grating and the mirror
I = sum(abs(dechirped_field_w).^2,2);
considered_regime = I > max(I)/1e4 & ~arrayfun(@(x)any(imag(x)),theta_out) & wavelength>0;
switch compressor_type(end)
    case {'1','2'} % Offner compressor with aberration
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
        roof_mirror_size = grating_size*cos(theta_in);
    case '3' % Offner compressor without aberration
        grating_size = max(optimal_value*tan(theta_out(considered_regime))) - min(optimal_value*tan(theta_out(considered_regime)));
        concave_size = max(2*R*theta_out(considered_regime)) - min(2*R*theta_out(considered_regime));
        convex_size = concave_size/2;
        roof_mirror_size = grating_size*cos(theta_in);
end

% Plot the result, if "verbose" is true.
if verbose
    switch compressor_type(end)
        case '1' % single-grating Offner compressor
            show_result(compressor_type,time,dechirped_field,0,optimal_value,        dechirped_FWHM,grating_size,roof_mirror_size,concave_size,convex_size);
        case '2' % double-grating version for the single-grating case
            show_result(compressor_type,time,dechirped_field,optimal_value,offcenter,dechirped_FWHM,grating_size,roof_mirror_size,concave_size,convex_size);
        case '3' % true aberration-free double-grating Offner compressor
            show_result(compressor_type,time,dechirped_field,optimal_value,0,        dechirped_FWHM,grating_size,roof_mirror_size,concave_size,convex_size);
    end
end

% Recover back to the input frequency range
dechirped_field = dechirped_field.*exp(-1i*2*pi*(f_c-c/wavelength0)*time);

% Information to recover the field
recover_info = {exp(1i*2*pi*(f_c-c/wavelength0)*time),exp(-1i*added_phase)};

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
index_shift = pulse_center-floor(length(time)/2);
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
index_shift = pulse_center-floor(length(wavelength)/2);
field = double(circshift(field,-index_shift,1));

end

%%
function [optimal_separation,dechirped_FWHM,dechirped_field,grating_size,roof_mirror_size,lens1_size,lens2_size,recover_info] = Martinez(theta_in,wavelength0,time,field,grating_spacing,focal_length,verbose,m,global_opt)
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
wavelength_c = c/f_c;
wavelength = c./(f+f_c-c/wavelength0);
field = field.*exp(1i*2*pi*(f_c-c/wavelength0)*time);
field_w = fftshift(ifft(field),1);

theta_out = asin( m*wavelength/(grating_spacing*1e9) + sin(theta_in) ); % the transmitted/reflected angle of the m-th order diffraction

find_optimum_compressor_separation = @(d) -Martinez_calc_max_intensity(d,theta_in,theta_out,wavelength,field_w,grating_spacing,focal_length,m);
find_FWHM = @(d) Martinez_calc_FWHM(d,theta_in,theta_out,wavelength,time,field_w,grating_spacing,focal_length,m);

% Run the global optimization process to find the optimal grating separation
% -------------------------------------------------------------------------
min_separation = 0;
%initial_guess = focal_length;
GVD_before_dechirping = characterize_spectral_phase(c./wavelength,fftshift(ifft(ifftshift(field(:,1),1)),1),3,false); % consider only the 1st mode
% Obtain the initial guess from the beta2 of a grating compressor.
% This number should be close to the answer already.
initial_guess = GVD_before_dechirping/2/(m^2*wavelength_c^3/(2*pi*c^2*grating_spacing^2)*(1-(-m*(wavelength_c*1e-9)/grating_spacing-sin(theta_in))^2)^(-1.5)*1e-3) + focal_length;

% This is the local optimization I used before which might be stuck at
% local minimum if the pulse isn't smooth enough.
option = optimset('TolX',1e-20);
%option = optimset('PlotFcns',@optimplotfval,'TolX',1e-20); % plot the process of optimization
[optimal_separation1,feval1] = fminsearch(find_optimum_compressor_separation,initial_guess,option);
% Because fminsearch is unconstrained method, I set the constrain below.
if optimal_separation1 < min_separation
    optimal_separation1 = min_separation;
end
% Global optimization
optimal_separation2 = 0;    feval2 = 0;
if global_opt
    transform_limited_field = calc_transform_limited( field );
    transform_limited_peak_power = max(abs(transform_limited_field).^2,[],1);
    dechirped_field = Martinez_dechirping(optimal_separation1,theta_in,theta_out,wavelength,field_w,grating_spacing,focal_length,m);
    dechirped_peak_power = max(abs(dechirped_field).^2,[],1);
    
    % Apply global optimization only when the dechirped field isn't good enough
    if dechirped_peak_power < transform_limited_peak_power*0.5
        problem = createOptimProblem('fmincon',...
            'objective',find_optimum_compressor_separation,...
            'x0',initial_guess,...
            'lb',min_separation,...
            'options',optimoptions(@fmincon,'Display','off','TolX',1e-20));
        gs = GlobalSearch('MaxTime',60,'NumTrialPoints',130,'NumStageOnePoints',20,'MaxWaitCycle',5,'Display','off');
        %gs = MultiStart('MaxTime',120,'Display','off','UseParallel',true,'StartPointsToRun','bounds-ineqs');
        %[optimal_separation2,feval2] = run(gs,problem,20);
        [optimal_separation2,feval2] = run(gs,problem);
    end
end

if feval1 < feval2
    optimal_separation = optimal_separation1;
else
    optimal_separation = optimal_separation2;
end

% The final dechirped pulse
[dechirped_field,lens1,lens2,grating,added_phase] = Martinez_dechirping(optimal_separation,theta_in,theta_out,wavelength,field_w,grating_spacing,focal_length,m);
dechirped_FWHM = find_FWHM(optimal_separation); % ps
dechirped_FWHM = dechirped_FWHM*1e3; % fs

tol_range = 1e-3; % 1 mm
if optimal_separation < (min_separation+tol_range)
    warning('The distance between the grating and the lens is too small.');
    
    if ~global_opt
        fprintf('Try running with global optimization.\n');
    end
end

dechirped_field_w = fftshift(ifft(dechirped_field),1);

% Minimum size required for the grating and the mirror
I = sum(abs(dechirped_field_w).^2,2);
considered_regime = I > max(I)/1e4 & ~arrayfun(@(x)any(imag(x)),theta_out) & wavelength>0;
lens1_size = max(lens1(considered_regime)) - min(lens1(considered_regime));
lens2_size = max(lens2(considered_regime)) - min(lens2(considered_regime));
grating_size = max(grating(considered_regime)) - min(grating(considered_regime));
roof_mirror_size = grating_size*cos(theta_in);

% Plot the result, if "verbose" is true.
if verbose
    show_result('Martinez',time,dechirped_field,optimal_separation,0,dechirped_FWHM,grating_size,roof_mirror_size,lens1_size,lens2_size);
end

% Recover back to the input frequency range
dechirped_field = dechirped_field.*exp(-1i*2*pi*(f_c-c/wavelength0)*time);

% Information to recover the field
recover_info = {exp(1i*2*pi*(f_c-c/wavelength0)*time),exp(-1i*added_phase)};

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
index_shift = pulse_center-floor(length(wavelength)/2);
field = double(circshift(field,-index_shift,1));

end

%%
function max_intensity = Treacy_calc_max_intensity(compressor_type,separation,theta_in,theta_out,wavelength,wavelength0,time,field_w,grating_spacing,m)
%TREACY_CALC_MAX_INTENSITY

field = Treacy_dechirping(compressor_type,separation,theta_in,theta_out,wavelength,wavelength0,time,field_w,grating_spacing,m);

max_intensity = max(sum(abs(field).^2,2));

end

function max_intensity = prism_calc_max_intensity(separation,max_dh,max_theta4,alpha,theta_in,theta2,theta3,theta4,n,wavelength,time,field_w)
%PRISM_CALC_MAX_INTENSITY

field = prism_dechirping(separation,max_dh,max_theta4,alpha,theta_in,theta2,theta3,theta4,n,wavelength,time,field_w);

max_intensity = max(sum(abs(field).^2,2));

end

function max_intensity = Offner_calc_max_intensity(compressor_type,separation,offcenter,theta_in,theta_out,wavelength,time,field_w,grating_spacing,R,m)
%OFFNER_CALC_MAX_INTENSITY

switch compressor_type(end)
    case {'1','2'} % Offner compressor with aberration
        field = Offner_dechirping_with_aberration(separation,offcenter,theta_in,theta_out,wavelength,time,field_w,grating_spacing,R,m);
    case '3'  % Offner compressor without aberration
        field = Offner_dechirping_no_aberration(separation,theta_in,theta_out,wavelength,field_w,grating_spacing,R,m);
end
max_intensity = max(sum(abs(field).^2,2));

end

function max_intensity = Martinez_calc_max_intensity(grating_lens_distance,theta_in,theta_out,wavelength,field_w,grating_spacing,focal_length,m)
%MARTINEZ_CALC_MAX_INTENSITY

field = Martinez_dechirping(grating_lens_distance,theta_in,theta_out,wavelength,field_w,grating_spacing,focal_length,m);

max_intensity = max(sum(abs(field).^2,2));

end

%%
function FWHM = Treacy_calc_FWHM(compressor_type,separation,theta_in,theta_out,wavelength,wavelength0,time,field_w,grating_spacing,m)
%TREACY_CALC_FWHM It finds the FWHM

field = Treacy_dechirping(compressor_type,separation,theta_in,theta_out,wavelength,wavelength0,time,field_w,grating_spacing,m);

I = sum(abs(field).^2,2);
threshold = max(I)/1.0001;
[~,~,FWHM,~] = findpeaks(I,time,'MinPeakHeight',threshold,'WidthReference','halfheight','MinPeakProminence',threshold/2);
FWHM = FWHM(1); % just in case the pulse is dechirped so badly that it has many peaks and is thus outputed many FWHMs

end

function FWHM = prism_calc_FWHM(separation,max_dh,max_theta4,alpha,theta_in,theta2,theta3,theta4,n,wavelength,time,field_w)
%PRISM_CALC_FWHM It finds the FWHM

field = prism_dechirping(separation,max_dh,max_theta4,alpha,theta_in,theta2,theta3,theta4,n,wavelength,time,field_w);

I = sum(abs(field).^2,2);
threshold = max(I)/1.0001;
[~,~,FWHM,~] = findpeaks(I,time,'MinPeakHeight',threshold,'WidthReference','halfheight','MinPeakProminence',threshold/2);
FWHM = FWHM(1); % just in case the pulse is dechirped so badly that it has many peaks and is thus outputed many FWHMs

end

function FWHM = Offner_calc_FWHM(compressor_type,separation,offcenter,theta_in,theta_out,wavelength,time,field_w,grating_spacing,R,m)
%OFFNER_CALC_FWHM It finds the FWHM

switch compressor_type(end)
    case {'1','2'} % Offner compressor with aberration
        field = Offner_dechirping_with_aberration(separation,offcenter,theta_in,theta_out,wavelength,time,field_w,grating_spacing,R,m);
    case '3' % Offner compressor without aberration
        field = Offner_dechirping_no_aberration(separation,theta_in,theta_out,wavelength,field_w,grating_spacing,R,m);
end

I = sum(abs(field).^2,2);
threshold = max(I)/1.0001;
[~,~,FWHM,~] = findpeaks(I,time,'MinPeakHeight',threshold,'WidthReference','halfheight','MinPeakProminence',threshold/2);
FWHM = FWHM(1); % just in case the pulse is dechirped so badly that it has many peaks and is thus outputed many FWHMs

end

function FWHM = Martinez_calc_FWHM(grating_lens_distance,theta_in,theta_out,wavelength,time,field_w,grating_spacing,focal_length,m)
%MARTINEZ_CALC_FWHM It finds the FWHM

field = Martinez_dechirping(grating_lens_distance,theta_in,theta_out,wavelength,field_w,grating_spacing,focal_length,m);

I = sum(abs(field).^2,2);
threshold = max(I)/1.0001;
[~,~,FWHM,~] = findpeaks(I,time,'MinPeakHeight',threshold,'WidthReference','halfheight','MinPeakProminence',threshold/2);
FWHM = FWHM(1); % just in case the pulse is dechirped so badly that it has many peaks and is thus outputed many FWHMs

end

%%
function show_result(compressor_type,time,dechirped_field,separation,offcenter,dechirped_FWHM,grating_size,roof_mirror_size,varargin)
%SHOW_RESULT

switch compressor_type
    case 'prism'
        prism_height = varargin{1};
    case {'Offner1','Offner2','Offner3'}
        concave_size = varargin{1};
        convex_size = varargin{2};
    case 'Martinez'
        lens1_size = varargin{1};
        lens2_size = varargin{2};
end

time = time*1e3; % fs
intensity = sum(abs(dechirped_field).^2,2);

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

figure('Name','Dechirped pulse');
h = plot(time(left:right),intensity(left:right));
xlim([min(time(left:right)) max(time(left:right))]);
xlabel('Time (fs)');
ylabel('Power (W)');
title('Dechirped pulse');

set(h,'linewidth',2);
set(gca,'fontsize',14);

fprintf('the shortest FWHM = %6.4f(fs)\n',dechirped_FWHM);
switch compressor_type
    case {'Treacy-t','Treacy-r'} % Treacy compressor
        fprintf('Grating separation = %6.4f(cm)\n',separation*100);
        fprintf('Minimum size of the grating = %6.4f(cm)\n',grating_size*100);
        fprintf('Minimum size of the mirror = %6.4f(cm)\n',roof_mirror_size*100);
    case 'prism'
        fprintf('Prism separation = %6.4f(cm)\n',separation*100);
        fprintf('Minimum size of the prism = %6.4f(cm)\n',prism_height*100);
        fprintf('Minimum size of the mirror = %6.4f(cm)\n',roof_mirror_size*100);
    case {'Offner1','Offner2'} % Offner compressor with aberration
        if compressor_type(end) == '1' % single-grating Offner compressor
            fprintf('Grating offcenter = %6.4f(cm)\n',offcenter*100);
        else % '2', double-grating version for the single-grating case
            fprintf('Grating separation = %6.4f(cm)\n',separation*100);
        end
        fprintf('(1) Light position on the grating = %6.4f ~ %6.4f(cm)\n',grating_size(1)*100,grating_size(2)*100);
        fprintf('--> Minimum size = %6.4f(cm)\n',(grating_size(2)-grating_size(1))*100);
        fprintf('(2) Light position on the concave mirror 1 = %6.4f ~ %6.4f(cm)\n',concave_size{1}(1)*100,concave_size{1}(2)*100);
        fprintf('--> Minimum size = %6.4f(cm)\n',(concave_size{1}(2)-concave_size{1}(1))*100);
        fprintf('(3) Light position on the concave mirror 2 = %6.4f ~ %6.4f(cm)\n',concave_size{2}(1)*100,concave_size{2}(2)*100);
        fprintf('--> Minimum size = %6.4f(cm)\n',(concave_size{2}(2)-concave_size{2}(1))*100);
        fprintf('(4) Size of the convex mirror = %6.4f(cm)\n',convex_size*100);
    case 'Offner3' % true aberration-free double-grating Offner compressor
        fprintf('Grating separation = %6.4f(cm)\n',separation*100);
        fprintf('Minimum size of the grating = %6.4f(cm)\n',grating_size*100);
        fprintf('Minimum size of the mirror = %6.4f(cm)\n',roof_mirror_size*100);
        fprintf('Minimum size of the concave mirror = %6.4f(cm)\n',concave_size*100);
        fprintf('Minimum size of the convex mirror = %6.4f(cm)\n',convex_size*100);
    case 'Martinez'
        fprintf('Grating to lens = %6.4f(cm)\n',separation*100);
        fprintf('Grating min. size = %6.4f(cm)\n',grating_size*100);
        fprintf('Lens 1 min. size = %6.4f(cm)\n',lens1_size*100);
        fprintf('Lens 2 min. size = %6.4f(cm)\n',lens2_size*100);
        fprintf('Roof mirror min. size = %6.4f(cm)\n',roof_mirror_size*100);
end

end

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