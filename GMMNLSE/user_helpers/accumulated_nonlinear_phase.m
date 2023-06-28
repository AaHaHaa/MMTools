function nonlinear_phase = accumulated_nonlinear_phase( L0,Aeff,f0,field,z,dt )
%ACCUMULATED_NONLINEAR_PHASE Calculates the total accumulated nonlinear
%phase from the propagation.
%   For oscillators, there are usually more than one type of fibers, mostly
%   including passive and active fibers, so this code also takes this into
%   account and this is why there's a dimension called "num_fiber" below.
%   save_num is the total number of saved fields.
%   
%   L0:   (1,num_fiber); the length of each fiber (m)
%   Aeff: (1,num_spatial_modes,num_fiber); the 1/( (m,m,m,m) overlap integral) of the mth eigenmodes (m^2)
%         num_spatial_modes = num_modes or num_modes/2
%   f0:   a scalar; central frequency (THz)
%   
%   field: (N,num_modes,save_num); fields during propagations (sqrt(W))
%               or
%           a structure with each field corresponding to its saved electric fields.
%               e.g. field.smf1 = (N,num_modes,save_num1) array
%                    field.gain = (N,num_modes,save_num2) array
%                    field.smf2 = (N,num_modes,save_num3) array
%                    ==> send this "field" as an input argument
%                        Also, save_num = save_num1 + save_num2 + save_num3
%
%   z: (1,save_num); the z position of each saved field (m)
%   dt: a scalar; the difference between each time point (ps)

if isstruct(field)
    fibername = fieldnames(field); % the fiber type
    N = size(field.(fibername{1}),1); % the number of time points
    num_modes = size(field.(fibername{1}),2); % the number of modes
    field_isstruct = true;
else
    N = size(field,1); % the number of time points
    num_modes = size(field,2); % the number of modes
    field_isstruct = false;
end
num_spatial_modes = size(Aeff,2);

num_segments = length(L0);

switch num_spatial_modes
    case num_modes/2
        Aeff = reshape(repmat(Aeff,[2,1,1]),[1,num_modes,size(Aeff,3)]);
    case num_modes
    otherwise
        error('accumulated_nonlinear_phase:AeffError',...
              'Check the number of modes of Aeff.');
end

%% Extract Aeff and field/"last_field(->field)" from the input arguments
save_num = length(z);
each_save_num = zeros(1,num_segments);
each_save_num(1) = sum(z<=L0(1));
extend_idx_Aeff = cell(1,save_num);
extend_idx_Aeff{1} = ones(1,each_save_num(1));
for i = 2:num_segments
    each_save_num(i) = sum(z>L0(i-1) & z<=L0(i));
    extend_idx_Aeff{i} = i*ones(1,each_save_num(i));
end
extend_idx_Aeff = cell2mat(extend_idx_Aeff);

Aeff = Aeff(:,:,extend_idx_Aeff);

if field_isstruct
    last_field = zeros(N,num_modes,save_num);
    last_field(:,:,1:each_save_num(1)) = field.(fibername{1});
end
ii = each_save_num(1);
for i = 2:num_segments
    if field_isstruct
        last_field(:,:,ii+1:ii+each_save_num(i)) = field.(fibername{i})(:,:,2:end);
    end
    
    ii = ii+each_save_num(i);
end
if field_isstruct
    field = last_field;
end

%% Some parameters
c = 2.99792458e-4; % speed of ligth m/ps
w0 = 2*pi*f0; % angular frequency (THz)
n2 = 2.3e-20; % m^2 W^-1
nonlin_const = n2*w0/c; % W^-1 m

%% Raman response function for silica
fr = 0.18; % for silica
t1 = 12.2e-3; % raman parameter t1 [ps]
t2 = 32e-3; % raman parameter t2 [ps]
t_shifted = dt*(0:N-1)'; % time starting at 0
 
hr = ((t1^2+t2^2)/(t1*t2^2)).*exp(-t_shifted/t2).*sin(t_shifted/t1);

hrw = ifft(hr)*N; % The factor of N is needed and used later in the convolution theorem because of how ifft is defined.
                  % The convolution theorem is PQ=F[(p*q)]/N, where (p*q) is the discrete-form convolution, for discrete Fourier transforms.
                  % Normal (integral-form) convolution is p*q=(p*q)*dt.

%% Calculate accumulated nonlinear phase
I = abs(field).^2./Aeff; % intensity: |E|^2/Aeff (|E|^2 is normalized to P, as in Agrawal)

if num_spatial_modes == num_modes/2
    I = I(:,1:2:num_modes-1,:) + I(:,2:2:num_modes,:);
end

nonlinear_phase_Kerr  = (1-fr)*nonlin_const*trapz(z,I,3);
nonlinear_phase_Raman =    fr *nonlin_const*trapz(z, real(fft(ifft(I).*hrw)) ,3)*dt; % This calculation gives a real result, so "real()" is just to simplify the notation and for "max()" below.

nonlinear_phase = max(nonlinear_phase_Kerr + nonlinear_phase_Raman);

end

