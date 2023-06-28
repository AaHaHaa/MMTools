function [SK_info, SRa_info, SRb_info] = calc_SRSK(fiber,sim,num_spatial_modes)
%CALC_SRSK It computes the SR and SK values required for the GMMNLSE
%
%   mode_info.nonzero_midx1234s - a (4,?) matrix
%   mode_info.SRa - a (num_nonzeros_midx1234s,number_frequency_points) array
%
%   ellipticity - the ellipticity of the polarization modes; Please refer to "Nonlinear Fiber Optics, eq (6.1.18) Agrawal" for the equations.
%                 0: linear polarization   -> (+,-)=(x,y)
%                 1: circular polarization -> (+,-)=(right,left)
%
%   include_anisotropic_Raman - true or false

% Because of the symmetry, there are possibly zeros in SR. To improve the
% performance, store only the indices of each nonzero elements and their
% corresponding values.
nonzero_midx1234s = find(permute(fiber.SR,[4 3 2 1])); % The indices of nonzero elements in SR will be stored in nonzero_midx1234s.
nonzero_midx1234s = reverse_linear_indexing(nonzero_midx1234s,num_spatial_modes);
[midx1,midx2,midx3,midx4] = ind2sub(num_spatial_modes*ones(1,4),nonzero_midx1234s'); % restore linear indexing back to subscripts
SRa_info.nonzero_midx1234s = uint8([midx1;midx2;midx3;midx4]);
SRa_info.SRa = fiber.SR(nonzero_midx1234s); % the corresponding values of each nonzero elements
clear nonzero_midx1234s midx1 midx2

if sim.gpu_yes
   SRa_info.SRa = gpuArray(SRa_info.SRa);
else % If not using the GPU, we also need to calculate the indices that don't have all zero coefficients for any given last two indices
    tmp = cellfun(@(SR12) any(SR12(:)), mat2cell(fiber.SR,num_spatial_modes,num_spatial_modes,ones(num_spatial_modes,1),ones(num_spatial_modes,1)));
    nonzero_midx34s = find(squeeze(tmp));
    [midx3,midx4] = ind2sub(num_spatial_modes*ones(1,2),nonzero_midx34s');
    SRa_info.nonzero_midx34s = uint8([midx3;midx4]);
end
clear midx3 midx4

if sim.scalar
    [SRa_info.nonzero_midx1234s,SRa_info.SRa,...
     SRa_info.beginning_nonzero,SRa_info.ending_nonzero] = refine_scalar_S(SRa_info.nonzero_midx1234s,SRa_info.SRa,...
                                                                    num_spatial_modes,...
                                                                    sim.gpu_yes);
    
    switch sim.ellipticity
        case 0 % linear polarization
            SK_factor = 1;
        case 1 % circular polarization
            SK_factor = 2/3;
        otherwise
            error('calc_SRSK:ellipticityError',...
                'The scalar mode supports only linear and circular polarizations.');
    end
    SK_info = struct('SK',SK_factor*SRa_info.SRa,...
                     'nonzero_midx1234s', SRa_info.nonzero_midx1234s,...
                     'beginning_nonzero',SRa_info.beginning_nonzero,...
                     'ending_nonzero',SRa_info.ending_nonzero);
    SRb_info = []; % dummy output argument
else
    [SK_info, SRa_info, SRb_info] = calc_polarized_SRSK(SRa_info,sim.ellipticity,num_spatial_modes*2,sim.Raman_model==2,sim.gpu_yes);
end

% Incoporate (1-fiber.fr) and fiber.fr into SRa, SRb, and SK to save the computational time
if sim.Raman_model ~= 0
    SK_info.SK = (1-fiber.fr)*SK_info.SK;
    
    if sim.Raman_model ~= 0
        SRa_info.SRa = fiber.fr*SRa_info.SRa;
        
        if ~isempty(SRb_info)
            SRb_info.SRb = fiber.fr*SRb_info.SRb;
        end
    end
end

% Set up dummy outputs for the cuda file to run normally, but it won't be
% used during the computations
if isempty(SRb_info)
    SRb_info = struct('SRb',0,...
                      'nonzero_midx1234s', uint8(zeros(4,1)),...
                      'beginning_nonzero', uint32(0),...
                      'ending_nonzero', uint32(0));
end

if sim.gpu_yes
    SK_info = struct('SK',gpuArray(SK_info.SK),...
                     'nonzero_midx1234s', gpuArray(SK_info.nonzero_midx1234s),...
                     'beginning_nonzero', gpuArray(SK_info.beginning_nonzero),...
                     'ending_nonzero', gpuArray(SK_info.ending_nonzero));
    SRa_info = struct('SRa',gpuArray(SRa_info.SRa),...
                      'nonzero_midx1234s', gpuArray(SRa_info.nonzero_midx1234s),...
                      'beginning_nonzero', gpuArray(SRa_info.beginning_nonzero),...
                      'ending_nonzero', gpuArray(SRa_info.ending_nonzero));
    SRb_info = struct('SRb',gpuArray(SRb_info.SRb),...
                      'nonzero_midx1234s', gpuArray(SRb_info.nonzero_midx1234s),...
                      'beginning_nonzero', gpuArray(SRb_info.beginning_nonzero),...
                      'ending_nonzero', gpuArray(SRb_info.ending_nonzero));
end

end

%% Helper function
function [SK_info, SRa_info, SRb_info] = calc_polarized_SRSK(mode_info,ellipticity,num_modes,include_anisotropic_Raman,use_gpu)
%CALC_POLARIZED_SRSK It computes the SR and SK, when considering polarization
%modes, from "scalar SR values under orthogonal polarizations".
%
%   mode_info.nonzero_midx1234s - a (4,?) matrix
%   mode_info.SRa - a (num_nonzeros_midx1234s,number_frequency_points) array
%
%   ellipticity - the ellipticity of the polarization modes; Please refer to "Nonlinear Fiber Optics, eq (6.1.18) Agrawal" for the equations.
%                 0: linear polarization   -> (+,-)=(x,y)
%                 1: circular polarization -> (+,-)=(right,left)
%
%   num_modes - the number of modes, including polarization modes
%
%   include_anisotropic_Raman - true or false

sSRa = size(mode_info.SRa);

if sSRa(2) ~= 1
    mode_info.SRa = mode_info.SRa.'; % SR needs to be a column vector
end
if sSRa(1) == 1 && sSRa(2) == 1 % only one spatial mode
    single_mode = true;
else
    single_mode = false;
end

oddidx = @(x) 2*x-1;
evenidx = @(x) 2*x;

%% SRa: isotropic Raman term
odd12 = oddidx(mode_info.nonzero_midx1234s([1 2],:));
even12 = evenidx(mode_info.nonzero_midx1234s([1 2],:));
odd34 = oddidx(mode_info.nonzero_midx1234s([3 4],:));
even34 = evenidx(mode_info.nonzero_midx1234s([3 4],:));

SRa_midx = cat(2, [odd12; odd34],...
                  [odd12; even34],...
                  [even12; odd34],...
                  [even12; even34]);
SRa = repmat(mode_info.SRa,4,1,1);

SRa_info = mode_info;

% Sort SR indices
[sort_SRa_midx,sort_idx] = sortrows(SRa_midx.');
SRa_info.nonzero_midx1234s = sort_SRa_midx.';
SRa_info.SRa = SRa(sort_idx);
[SRa_info.nonzero_midx1234s,SRa_info.SRa,...
 SRa_info.beginning_nonzero,SRa_info.ending_nonzero] = refine_polarized_S(SRa_info.nonzero_midx1234s,SRa_info.SRa,...
                                                                          num_modes,...
                                                                          use_gpu);

%% SK
odd1 = oddidx(mode_info.nonzero_midx1234s(1,:)); even1 = evenidx(mode_info.nonzero_midx1234s(1,:));
odd2 = oddidx(mode_info.nonzero_midx1234s(2,:)); even2 = evenidx(mode_info.nonzero_midx1234s(2,:));
odd3 = oddidx(mode_info.nonzero_midx1234s(3,:)); even3 = evenidx(mode_info.nonzero_midx1234s(3,:));
odd4 = oddidx(mode_info.nonzero_midx1234s(4,:)); even4 = evenidx(mode_info.nonzero_midx1234s(4,:));

switch ellipticity
    case 0 % linear polarizations
        odd1 = oddidx(mode_info.nonzero_midx1234s(1,:));         odd4 = oddidx(mode_info.nonzero_midx1234s(4,:));
        even1 = evenidx(mode_info.nonzero_midx1234s(1,:));       even4 = evenidx(mode_info.nonzero_midx1234s(4,:));
        odd23 = oddidx(mode_info.nonzero_midx1234s([2 3],:));
        even23 = evenidx(mode_info.nonzero_midx1234s([2 3],:));

        % Sk_midx = cat(2, [odd1;   odd23;  odd4],...
        %                  [odd1;  even23;  odd4],...
        %                  [even1;  odd23; even4],...
        %                  [even1; even23; even4]);
        % SK = 2/3*SRa + 1/3*Sk
        %
        SK_midx = cat(2, [odd12; odd34],...        % SR
                         [even12; even34],...      % SR
                         [odd12; even34],...       % 2/3*SR
                         [even12; odd34],...       % 2/3*SR
                         [odd1; even23; odd4],...  % 1/3*SR
                         [even1; odd23; even4]);   % 1/3*SR
        SK = mode_info.SRa.*[1 1 2/3 2/3 1/3 1/3];
    case 1 % circular polarization
        odd1 = oddidx(mode_info.nonzero_midx1234s(1,:)); even1 = evenidx(mode_info.nonzero_midx1234s(1,:));
        odd2 = oddidx(mode_info.nonzero_midx1234s(2,:)); even2 = evenidx(mode_info.nonzero_midx1234s(2,:));
        odd3 = oddidx(mode_info.nonzero_midx1234s(3,:)); even3 = evenidx(mode_info.nonzero_midx1234s(3,:));
        odd4 = oddidx(mode_info.nonzero_midx1234s(4,:)); even4 = evenidx(mode_info.nonzero_midx1234s(4,:));

        SK_midx = cat(2, [ odd1;  odd2;  odd3;  odd4],...  % 0
                         [ odd1;  odd2; even3; even4],...  % 1
                         [ odd1; even2;  odd3; even4],...  % 1
                         [even1;  odd2; even3;  odd4],...  % 1
                         [even1; even2;  odd3;  odd4],...  % 1
                         [even1; even2; even3; even4]);    % 0
        % SK = 2/3*SRa + 1/3*Sk
        Sk_term = [0 1 1 1 1 0];
        SRa_term = [1 1 0 0 1 1];
        SK = mode_info.SRa.*(2/3*SRa_term+1/3*Sk_term);
    otherwise % elliptical polarization
        % basis_o = (x+iry)/sqrt(1+r^2)
        % basis_e = (rx-iy)/sqrt(1+r^2)
        %
        % Notice that r=0 corresponds to basis_e = -iy. Since I separate the
        % linear polarization above, which has (basis_o = x, basis_e = y), it doesnt' matter here.
        %
        r = ellipticity; % match the notation with the "Nonlinear Fiber Optics, Agrawal"
        oo = (1-r^2)/(1+r^2);
        ee = -oo;
        oe = 2*r/(1+r^2);

        SK_midx = cat(2, [ odd1;  odd2;  odd3;  odd4],...  % (oo)^2
                         [ odd1;  odd2;  odd3; even4],...  % (oo)(oe)
                         [ odd1;  odd2; even3;  odd4],...  % ...
                         [ odd1; even2;  odd3;  odd4],...  % ...
                         [even1;  odd2;  odd3;  odd4],...  % ...
                         [ odd1;  odd2; even3; even4],...  % (oe)^2
                         [ odd1; even2;  odd3; even4],...  % ...
                         [even1;  odd2;  odd3; even4],...  % (oo)(ee)
                         [ odd1; even2; even3;  odd4],...  % ...
                         [even1;  odd2; even3;  odd4],...  % (oe)^2
                         [even1; even2;  odd3;  odd4],...  % ...
                         [ odd1; even2; even3; even4],...  % (ee)(oe)
                         [even1;  odd2; even3; even4],...  % ...
                         [even1; even2;  odd3; even4],...  % ...
                         [even1; even2; even3;  odd4],...  % ...
                         [even1; even2; even3; even4]);    % (ee)^2
        % SK = 2/3*SRa + 1/3*Sk
        Sk_term = [oo^2 ...
                   oo*oe oo*oe oo*oe oo*oe ...
                   oe^2 oe^2 ...
                   oo*ee oo*ee ...
                   oe^2 oe^2 ...
                   ee*oe ee*oe ee*oe ee*oe ...
                   ee^2];
        SRa_term = [1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1];
        SK = mode_info.SRa.*(2/3*SRa_term+1/3*Sk_term);
end

% Sort SK indices
[sort_SK_midx,sort_idx] = sortrows(SK_midx.');
SK_info.nonzero_midx1234s = sort_SK_midx.';
SK_info.SK = SK(sort_idx);
[SK_info.nonzero_midx1234s,SK_info.SK,...
 SK_info.beginning_nonzero,SK_info.ending_nonzero] = refine_polarized_S(SK_info.nonzero_midx1234s,SK_info.SK,...
                                                                        num_modes,...
                                                                        use_gpu);

if single_mode % it needs to be a column vector
                % If it's a scalar field, SRa is a scalar such that SK,
                % calculated from SRa, is a row vector. Hence, SK needs a
                % transpose.
    SK_info.SK = SK_info.SK.';
end

%% SRb: anisotropic Raman term
if include_anisotropic_Raman
    switch ellipticity
        case 0 % linear polarizations

            % Srb_midx = cat(2, [ odd1;  odd2;  odd3;  odd4],...
            %                   [ odd1; even2;  odd3; even4],...
            %                   [even1;  odd2; even3;  odd4],...
            %                   [even1; even2; even3; even4]);
            %
            % Sk_midx = cat(2, [odd1;   odd23;  odd4],...
            %                  [odd1;  even23;  odd4],...
            %                  [even1;  odd23; even4],...
            %                  [even1; even23; even4]);
            %
            % SRb = 1/2*(Srb + Sk)
            %
            SRb_midx = cat(2, [ odd12;  odd34],...              % SR
                              [even12; even34],...              % SR
                              [ odd1; even2;  odd3; even4],...  % 1/2*SR
                              [even1;  odd2;  odd3; even4],...  % 1/2*SR
                              [ odd1; even2; even3;  odd4],...  % 1/2*SR
                              [even1;  odd2; even3;  odd4]);    % 1/2*SR
            SRb = mode_info.SRa.*[1 1 1/2 1/2 1/2 1/2];
        case 1 % circular polarization
            SRb_midx = cat(2, [ odd1;  odd2;  odd3;  odd4],...  % 1/2
                              [ odd1;  odd2; even3; even4],...  % 1/2
                              [ odd1; even2;  odd3; even4],...  % 1
                              [even1;  odd2; even3;  odd4],...  % 1
                              [even1; even2;  odd3;  odd4],...  % 1/2
                              [even1; even2; even3; even4]);    % 1/2
            % SRb = 1/2*(Srb + Sk)
            SRb = mode_info.SRa.*[1/2 1/2 1 1 1/2 1/2];
        otherwise % elliptical polarization
            % basis_o = (x+iry)/sqrt(1+r^2)
            % basis_e = (rx-iy)/sqrt(1+r^2)
            %
            % Notice that r=0 corresponds to basis_e = -iy. Since I separate the
            % linear polarization above, which has (basis_o = x, basis_e = y), it doesnt' matter here.
            %

            SRb_midx = SK_midx;
            
            % SRb = 1/2*(Srb + Sk)
            Srb_term = [1 0 0 0 0 0 1 0 0 1 0 0 0 0 0 1];
            SRb = mode_info.SRa.*(1/2*(Srb_term+Sk_term));
    end

    % Sort SRb indices
    [sort_SRb_midx,sort_idx] = sortrows(SRb_midx.');
    SRb_info.nonzero_midx1234s = sort_SRb_midx.';
    SRb_info.SRb = SRb(sort_idx); 
    [SRb_info.nonzero_midx1234s,SRb_info.SRb,...
     SRb_info.beginning_nonzero,SRb_info.ending_nonzero] = refine_polarized_S(SRb_info.nonzero_midx1234s,SRb_info.SRb,...
                                                                              num_modes,...
                                                                              use_gpu);
    
    if single_mode % it needs to be a column vector
                    % If it's a scalar field, SRa is a scalar such that SRb,
                    % calculated from SRa, is a row vector. Hence, SRb needs a
                    % transpose.
        SRb_info.SRb = SRb_info.SRb.';
    end
else
    SRb_info = []; % dummy output argument
end


%% midx_34s: Part of Raman terms calculations under CPU
if isfield(mode_info,'nonzero_midx34s')
    SRa_midx34s = [odd34 even34];
    SRa_info.nonzero_midx34s = sortrows(SRa_midx34s.').';
    if include_anisotropic_Raman
        SRb_midx34s = cat(2, odd34,         ...
                            even34,         ...
                            [ odd3; even4], ...
                            [even3;  odd4]);
        SRb_info.nonzero_midx34s = sortrows(SRb_midx34s.').';
    end
end

end

function rev_idx = reverse_linear_indexing(idx,n)
%REVERSE_LINEAR_INDEXING
%
%   MATLAB linear indexing starts from the first dimension and so on, e.g.,
%   A = [7 5; is expanded into [7 8 1 9 5 2 4 1] while its indices correspond to [1 5;
%        8 2;                                                                     2 6;
%        1 4;                                                                     3 7;
%        9 1]                                                                     4 8]
%
%   Please refer to MATLAB linear indexing for details.
%
%   However, in SR and SK summation calculation, "Ra" and "nonlinear" loop  
%   over the first two indices. If "nonzero_midx" can be sorted in a way 
%   that doesn't move back and forth the memory all the time, it should be
%   faster...probably. So I add this function to correct the linear 
%   indexing direction given from "find" function.
%
%   Take the above as an example, the result indexing should now be [1 2;  while the first two rows correspond to midx1 and midx2.
%                                                                    3 4;
%                                                                    5 6;
%                                                                    7 8]
% -------------------------------------------------------------------------
% Original linear indexing direction:
%   nonzero_midx1234s = find(fiber.SR);
%
% Corrected linear indexing direction:
%   nonzero_midx1234s = find(permute(fiber.SR,[4 3 2 1]));
%   nonzero_midx1234s = reverse_linear_indexing(nonzero_midx1234s,num_spatial_modes);

rev_matrix = permute(reshape(1:n^4,n,n,n,n),[4 3 2 1]);
rev_idx = rev_matrix(idx);

end

function [beginning_nonzero,ending_nonzero] = beginning_idx(nonzero_midx1234s,num_modes)
%BEGINNING_IDX It finds the first position of each (midx1,midx2) in
%nonzero_midx1234s

beginning_nonzero = zeros(num_modes,num_modes);
ending_nonzero = zeros(num_modes,num_modes);
for midx1 = 1:num_modes
    beginning_nonzero_midx1 = find(nonzero_midx1234s(1,:) == midx1,1);
    if ~isempty(beginning_nonzero_midx1)
        for midx2 = 1:num_modes
            beginning_nonzero_midx2 = find(nonzero_midx1234s(2,beginning_nonzero_midx1:end) == midx2,1);
            if ~isempty(beginning_nonzero_midx2) && nonzero_midx1234s(1,beginning_nonzero_midx1 + beginning_nonzero_midx2 -1) == midx1
                beginning_nonzero(midx1,midx2) = beginning_nonzero_midx1 + beginning_nonzero_midx2 -1;

                if midx1 == num_modes && midx2 == num_modes
                    ending_nonzero(midx1,midx2) = size(nonzero_midx1234s,2)+1;
                else
                    ending_nonzero_midx2 = find(nonzero_midx1234s(2,beginning_nonzero(midx1,midx2):end) ~= midx2,1);
                    ending_nonzero(midx1,midx2) = beginning_nonzero(midx1,midx2) + ending_nonzero_midx2 -1;
                end
            end
        end
    end
end
beginning_nonzero = uint32(beginning_nonzero');
ending_nonzero = uint32(ending_nonzero');

end

function [nonzero_midx1234s,S,...
          beginning_nonzero,ending_nonzero] = refine_scalar_S(nonzero_midx1234s,S,...
                                                              num_modes,...
                                                              use_gpu)
%REFINE_SCALAR_S It refine the nonzero indices by removing (midx1,midx2,midx4,midx3)
%when midx3~=midx4 in scalar simulations. The cuda file in GMMNLSE is 
%updated accordingly.

[beginning_nonzero,ending_nonzero] = beginning_idx(nonzero_midx1234s,num_modes);

if use_gpu
    for midx1 = 1:num_modes
        for midx2 = 1:num_modes
            if beginning_nonzero(midx1,midx2) ~= 0
                for i = beginning_nonzero(midx1,midx2):(ending_nonzero(midx1,midx2)-2)
                    if nonzero_midx1234s(3,i) < nonzero_midx1234s(4,i)
                        duplicate_idx = find(ismember(nonzero_midx1234s([3,4],(i+1):(ending_nonzero(midx1,midx2)-1))', [nonzero_midx1234s(4,i),nonzero_midx1234s(3,i)], 'rows'),1);
                        
                        S(duplicate_idx + i) = 0;
                    end
                end
            end
        end
    end

    to_be_remove = (S==0);
    nonzero_midx1234s = nonzero_midx1234s(:,~to_be_remove);
    S = S(~to_be_remove);

    [beginning_nonzero,ending_nonzero] = beginning_idx(nonzero_midx1234s,num_modes);
end

end

function [nonzero_midx1234s,S,...
          beginning_nonzero,ending_nonzero] = refine_polarized_S(nonzero_midx1234s,S,...
                                                                 num_modes,...
                                                                 use_gpu)
%REFINE_POLARIZED_S It refine the nonzero indices by removing (midx1,midx2,midx4,midx3)
%when midx3~=midx4 in polarized simulations. The cuda file in GMMNLSE is 
%updated accordingly.

[beginning_nonzero,ending_nonzero] = beginning_idx(nonzero_midx1234s,num_modes);

if use_gpu
    for midx1 = 1:num_modes
        for midx2 = 1:num_modes
            if beginning_nonzero(midx1,midx2) ~= 0
                for i = beginning_nonzero(midx1,midx2):(ending_nonzero(midx1,midx2)-2)
                    if nonzero_midx1234s(3,i) < nonzero_midx1234s(4,i) && mod(nonzero_midx1234s(3,i),2) == mod(nonzero_midx1234s(4,i),2)
                        duplicate_idx = find(ismember(nonzero_midx1234s([3,4],(i+1):(ending_nonzero(midx1,midx2)-1))', [nonzero_midx1234s(4,i),nonzero_midx1234s(3,i)], 'rows'),1);
                        
                        S(duplicate_idx + i) = 0;
                    end
                end
            end
        end
    end

    to_be_remove = (S==0);
    nonzero_midx1234s = nonzero_midx1234s(:,~to_be_remove);
    S = S(~to_be_remove);

    [beginning_nonzero,ending_nonzero] = beginning_idx(nonzero_midx1234s,num_modes);
end

end