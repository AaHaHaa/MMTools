function [Arad,Gamma,energy_levels] = calc_Judd_Ofelt(gain_rate_eqn)
%CALC_JUDD_OFELT Computes radiative and non-radiative transition rates, 
% Arad and Gamma

%% 
skip_yes = false;
switch gain_rate_eqn.gain_medium
    case 'Er'
        if ismember(lower(gain_rate_eqn.base_medium),{'silica','silicate','fluorophosphate','chalcohalide','zblan','fluoride','germanate','tellurite'})
            skip_yes = true;
        end
        
        % integral(sigma_abs,lambda) is applied for
        % 4I15/2 --> 2H9/2:            ~411 nm
        % 4I15/2 --> (4F3/2 +) 4F5/2:  ~449 nm
        % 4I15/2 --> 4F7/2:            ~487 nm
        % 4I15/2 --> 4S3/2 (+ 2H11/2): ~525 nm
        % 4I15/2 --> 4F9/2:            ~653 nm
        % 4I15/2 --> 4I9/2:            ~800 nm
        % 4I15/2 --> 4I11/2:           ~975 nm
        % 4I15/2 --> 4I13/2:          ~1529 nm
        doping_data = load('JO_erbium.mat');
        band_idx = [11,13,14,16:20]; % row indices of the corresponding transitions (GSA) in JO_erbium.mat
        G_J = 16; % 2*J+1=2*6+1=16 for the gound-state 4I15/2
        
        %num_levels = 9; % N0, N1,..., N8
        level_idx = fliplr([5,7,8,10:15]); % indices of the transitions among considered energy levels
        transition_idx = [ 0,105,104,102,99,95,84,77,60;...
                          105, 0,103,101,98,94,83,76,59;...
                          104,103, 0,100,97,93,82,75,58;...
                          102,101,100, 0,96,92,81,74,57;...
                           99, 98, 97,96, 0,91,80,73,56;...
                           95, 94, 93,92,91, 0,79,72,55;...
                           84, 83, 82,81,80,79, 0,70,53;...
                           77, 76, 75,74,73,72,70, 0,52;...
                           60, 59, 58,57,56,55,53,52, 0];
        transition_idx(transition_idx==0) = 1; % these entries won't be used, but need to be set for the code to run correctly
        
        [xx,yy] = meshgrid(level_idx); level_matrix_idx = min(cat(3,xx,yy),[],3);
        Arad_J = doping_data.num(level_matrix_idx);
    case 'Nd'
        if ismember(lower(gain_rate_eqn.base_medium),{'silica','silicate'})
            skip_yes = true;
        end
        
        % integral(sigma_abs,lambda) is applied for
        % 4I9/2 --> 2D5/2 + 2P1/2:                       ~434 nm
        % 4I9/2 --> 2G11/2 + 2(D,P)3/2 + 2G9/2 + 2K15/2: ~470 nm
        % 4I9/2 --> 4G9/2 + 4G7/2 + 2K13/2:              ~520 nm
        % 4I9/2 --> 2G7/2 + 4G5/2:                       ~575 nm
        % 4I9/2 --> 4F9/2:                               ~680 nm
        % 4I9/2 --> 4S3/2 + 4F7/2:                       ~740 nm
        % 4I9/2 --> 2H9/2 + 4F5/2:                       ~795 nm
        % 4I9/2 --> 4F3/2:                               ~865 nm
        doping_data = load('JO_neodymium.mat');
        band_idx = [2,3,6,11,12,14,16,17,18]; % row indices of the corresponding transitions (GSA) in JO_neodymium.mat
        G_J = 10; % 2*J+1=10 for the gound-state 4I9/2
        
        %num_levels = 12; % N0, N1,..., N11
        level_idx = fliplr([2,3,6,11,12,14,16:21]); % indices of the transitions among considered energy levels
        transition_idx = [  0 210,209,207,204,200,189,174,165,105, 57, 39;...
                          210,  0,208,206,203,199,188,173,164,104, 56, 38;...
                          209,208,  0,205,202,198,187,172,163,103, 55, 37;...
                          207,206,205,  0,201,197,186,171,162,102, 54, 36;...
                          204,203,202,201,  0,196,185,170,161,101, 53, 35;...
                          200,199,198,197,196,  0,184,169,160,100, 52, 34;...
                          189,188,187,186,185,184,  0,167,158, 98, 50, 32;...
                          174,173,172,171,170,169,167,  0,156, 96, 48, 30;...
                          165,164,163,162,161,160,158,156,  0, 95, 47, 29;...
                          105,104,103,102,101,100, 98, 96, 95,  0, 42, 24;...
                           57, 56, 55, 54, 53, 52, 50, 48, 47, 42,  0, 21;...
                           39, 38, 37, 36, 35, 34, 32, 30, 29, 24, 21,  0];
        transition_idx(transition_idx==0) = 1; % these entries won't be used, but need to be set for the code to run correctly
        
        [xx,yy] = meshgrid(level_idx); level_matrix_idx = min(cat(3,xx,yy),[],3);
        Arad_J = doping_data.num(level_matrix_idx);
    case 'Tm'
        if ismember(lower(gain_rate_eqn.base_medium),{'silica','silicate','fluorophosphate','chalcohalide','zblan','fluoride','germanate','tellurite'})
            skip_yes = true;
        end
        
        % integral(sigma_abs,lambda) is applied for
        % 3H6 --> 1G4:        ~443 nm
        % 3H6 --> 3F2 + 3F3:  ~680 nm
        % 3H6 --> 3H4:        ~786 nm
        % 3H6 --> 3H5:       ~1209 nm
        % 3H6 --> 3F4:       ~1630 nm
        doping_data = load('JO_thulium.mat');
        band_idx = [7,9,10,11,12]; % row indices of the corresponding transitions (GSA) in JO_thulium.mat
        G_J = 13; % 2*J+1=2*6+1=13 for the gound-state 3H6
        
        % For N4, I pick 3F3, rather than 3F2, because 3F2 should decay to 3F3 fast.
        %num_levels = 8; % N0, N1,..., N7
        level_idx = fliplr([5,6,7,9:13]); % indices of the transitions among considered energy levels
        transition_idx = [ 0,78,77,75,72,63,57,50;...
                          78, 0,76,74,71,62,56,49;...
                          77,76, 0,73,70,61,55,48;...
                          75,74,73, 0,69,60,54,47;...
                          72,71,70,69, 0,59,53,46;...
                          63,62,61,60,59, 0,51,44;...
                          57,56,55,54,53,51, 0,43;...
                          50,49,48,47,46,44,43, 0];
        transition_idx(transition_idx==0) = 1; % these entries won't be used, but need to be set for the code to run correctly
        
        [xx,yy] = meshgrid(level_idx); level_matrix_idx = min(cat(3,xx,yy),[],3);
        Arad_J = doping_data.num(level_matrix_idx);
    case 'Ho'
        if ismember(lower(gain_rate_eqn.base_medium),{'silica','silicate','fluorophosphate','chalcohalide','zblan','fluoride','germanate','tellurite'})
            skip_yes = true;
        end
        
        % integral(sigma_abs,lambda) is applied for
        % 5I8 --> 5I7: ~1950 nm
        % 5I8 --> 5I6: ~1150 nm
        % 5I8 --> 5I5:  ~883 nm
        % 5I8 --> 5I4:  ~642 nm
        doping_data = load('JO_holmium.mat');
        band_idx = 31:34; % row indices of the corresponding transitions (GSA) in JO_holmium.mat
        G_J = 17; % 2*J+1=2*8+1=17 for the gound-state 5I8
        
        %num_levels = 5; % N0, N1,..., N4
        level_idx = fliplr(11:15); % indices of the transitions among considered energy levels
        transition_idx = [ 0,91,90,88,85;...
                          91, 0,89,87,84;...
                          90,89, 0,86,83;...
                          88,87,86, 0,82;...
                          85,84,83,82, 0];
        transition_idx(transition_idx==0) = 1; % these entries won't be used, but need to be set for the code to run correctly
        
        [xx,yy] = meshgrid(level_idx); level_matrix_idx = min(cat(3,xx,yy),[],3);
        Arad_J = doping_data.num(level_matrix_idx);
end

%% Use the measured ground-state absorption cross sections to obtain Judd-Ofeld's Omega parameters
if ~skip_yes
    switch gain_rate_eqn.gain_medium
        case 'Tm'
            lambda = linspace(426,2500,1000)'*1e-9;
            [GSA05,GSA04,GSA03,GSA02,GSA01] = read_cross_sections_Tm(gain_rate_eqn.cross_section_filename, lambda);
            GSA = [GSA05,GSA04,GSA03,GSA02,GSA01];
    end
    
    band_lambda = zeros(size(GSA,2),1);
    band_sum = zeros(size(GSA,2),1);
    for i = 1:size(GSA,2)
        band_sum(i) = trapz(lambda,GSA(:,i)); % m^3
        band_lambda(i) = trapz(lambda,lambda.*GSA(:,i))/band_sum(i); % m
    end
    JO_GSA_of_U = doping_data.layerM(band_idx); %#ok not used; names of the transitions
    JO_U = doping_data.raree(band_idx,:); % reduced matrix elements for the unit tensor operator U(2), U(4), and U(6)
end

% Load base material of the fiber
% Check "base_data.material" below which contains the names.
% It's typically 'silica' or 'ZBLAN' for fiber lasers.
base_medium = gain_rate_eqn.base_medium;
if strcmpi(gain_rate_eqn.base_medium,'silica') % Sellmeier mat file has only 'sio2', so we need to change the name to 'silica'
    base_medium = 'sio2';
end
base_data = load('JO_Sellmeier_base_material.mat');
base_idx = strcmpi(base_data.material,base_medium);
JO_Sellmeier_coeff = base_data.sellmeier(base_idx,:);
if ~skip_yes
    band_n = sqrt(JO_Sellmeier_coeff(1) + ...
                  JO_Sellmeier_coeff(2)*(band_lambda*1e6).^2./((band_lambda*1e6).^2 - JO_Sellmeier_coeff(3)) + ...
                  JO_Sellmeier_coeff(4)*(band_lambda*1e6).^2./((band_lambda*1e6).^2 - JO_Sellmeier_coeff(5)));
    chi = ( (band_n.^2 + 2)/3 ).^2; % a factor in calculations

    % Calculation of the Judd-Ofeld's Omega and line strength
    Sm = band_sum*10.413487e33.*band_n./chi*G_J./(band_lambda*1e9); % note that: 3*h*c/(8*pi^3*e^2) = 10.413487)
    Omega = (JO_U'*JO_U)\(JO_U'*Sm); % size: (3,1); Judd-Ofelt least-square fit for absorption line strengths
else
    switch gain_rate_eqn.gain_medium
        case 'Er'
            switch lower(gain_rate_eqn.base_medium)
                case 'erzsg' % ErZSG
                    % Yamasaki et al., "Optical properties of Er^{3+} 
                    % heavily doped silica glass fabricated by zeolite 
                    % method," Journal of Non-Crystalline Solids 543 120149
                    % (2020)
                    Omega = [11.9;2.47;2.62];
                case 'silica'
                    % Bobu et al., "Optical spectroscopy, 1.5μm emission, 
                    % and upconversion properties of Er3+-doped 
                    % metaphosphate laser glasses," J. Opt. Soc. Am. B 24,
                    % 2218-2228 (2007)
                    Omega = [8.36;1.76;0.82];
                    % D. Ning, "High Quantum Efficiency and High
                    % Concentration Erbium-Doped Silica Glasses Fabricated
                    % by Sintering Nanoporous Glasses," Journal of Rare
                    % Earths 24, 761-764 (2006)
                   %Omega = [7.93;1.70;1.33];
                case 'silicate'
                    % Richard Quimby, "Excited state absorption at 980 nm
                    % in erbuim doped glass," Fiber Laser Sources and
                    % Amplifiers III 1581, 72-79 (1992)
                    Omega = [4.26;0.81;0.46];
                case 'fluorophosphate' % high F
                    % Richard Quimby, "Excited state absorption at 980 nm
                    % in erbuim doped glass," Fiber Laser Sources and
                    % Amplifiers III 1581, 72-79 (1992)
                    Omega = [2.75;1.69;1.15];
               %case 'fluorophosphate' % low F
               %    Omega = [3.39;1.79;1.20];
                case 'phosphate'
                    % Richard Quimby, "Excited state absorption at 980 nm
                    % in erbuim doped glass," Fiber Laser Sources and
                    % Amplifiers III 1581, 72-79 (1992)
                    Omega = [5.42;1.55;0.98];
                    % W. Krupke, "Induced-emission cross sections in
                    % neodymium laser glasses," IEEE Journal of Quantum
                    % Electronics 10, 450-457 (1974)
                   %Omega = [3.53;0.76;0.56];
                case {'zblan','fluoride'}
                    % Richard Quimby, "Excited state absorption at 980 nm
                    % in erbuim doped glass," Fiber Laser Sources and
                    % Amplifiers III 1581, 72-79 (1992)
                    Omega = [2.51;1.41;0.956];
                    % Tanabe et al., "Compositional dependence of Judd-Ofelt
                    % parameters of Er3+ ions in alkali-metal borate
                    % glasses," J. Phys. Rev. B 46(6) 3305-3310 (1992)
                   %Omega = [2.91;1.27;1.11];
                    % G. S. Ofelt, "Intensities of Crystal Spectra of
                    % Rare-Earth Ions," The Journal of Chemical Physics 37,
                    % S11-S20 (1962)
                   %Omega = [2.91;1.78;1.00];
                case 'zbla'
                    % Richard Quimby, "Excited state absorption at 980 nm
                    % in erbuim doped glass," Fiber Laser Sources and
                    % Amplifiers III 1581, 72-79 (1992)
                    Omega = [2.54;1.39;0.965];
                case 'tellurite'
                    % Gomes et al., "Energy level decay and excited state
                    % absorption processes in erbium-doped tellurite
                    % glass," J. Appl. Phys. 110, 083111 (2011)
                    Omega = [5.54;1.70;1.22]; % 1e-20*cm^2
                    %Lachheb et al., "Characterization of Tm3+ doped TNZL
                    %glass laser material," Journal of Luminescence 161
                    %281-287 (2015)
                   %Omega = [4.64;1.61;1.26];
            end
        case 'Nd'
            switch lower(gain_rate_eqn.base_medium)
                case 'silicate'
                    % William F. Krupke, "Induced-Emission Cross Sections in
                    % Neodymium Laser Glasses," IEEE J. Quantum Electron. 10(4), 
                    % 450-457 (1974)
                    Omega = [3.30;4.68;5.18]; % 1e-20*cm^2
                case 'silica'
                    % Yanbo et al., "Spectroscopic Properties of
                    % Nd^{3+}-Doped High Silica Glass Prepared
                    % by Sintering Porous Glass," J. Rare Earth, 24(6),
                    % 765-770 (2006)
                    Omega = [5.02,3.81,2.56]; % 1e-20*cm^2
                case 'fluorophosphate'
                    % Naftaly and Jha, "Nd^{3+}-doped fluoroaluminate
                    % glasses for a 1.3 um amplifier," J. Appl. Phys. 87,
                    % 2098-2104 (2000)
                    Omega = [2.7;3.2;5.1]; % 1e-20*cm^2
                case 'zblan'
                    % Naftaly and Jha, "Nd^{3+}-doped fluoroaluminate
                    % glasses for a 1.3 um amplifier," J. Appl. Phys. 87,
                    % 2098-2104 (2000)
                    Omega = [2.2;2.8;3.9]; % 1e-20*cm^2
            end
        case 'Tm'
            % Wang et al., "Spectroscopic and structural characterization
            % of barium tellurite glass fibers for mid-infrared ultra-broad
            % tunable fiber lasers," Opt. Mater. Express 6(6), 2095-2107
            % (2016)
            switch lower(gain_rate_eqn.base_medium)
                case 'silica'
                    Omega = [6.23;1.91;1.36];
                case 'silicate'
                    Omega = [3.08;0.99;0.40];
                case 'fluorophosphate'
                    Omega = [3.01;2.56;1.54];
                case 'chalcohalide'
                    Omega = [5.80;1.60;1.30];
                case {'zblan','fluoride'}
                    Omega = [1.96;1.36;1.16];
                case 'germanate'
                    Omega = [6.11;1.41;1.31];
                case 'tellurite'
                    Omega = [4.73;0.84;1.15];
            end
        case 'Ho'
            % Wang et al., "Spectroscopic and structural characterization
            % of barium tellurite glass fibers for mid-infrared ultra-broad
            % tunable fiber lasers," Opt. Mater. Express 6(6), 2095-2107
            % (2016)
            switch lower(gain_rate_eqn.base_medium)
                case 'silica'
                    Omega = [1.97;1.47;1.23];
                case 'silicate'
                    Omega = [3.14;3.04;0.94];
                case 'fluorophosphate'
                    Omega = [2.05;3.77;1.28];
                case 'chalcohalide'
                    Omega = [0.10;4.97;0.98];
                case {'zblan','fluoride'}
                    Omega = [1.86;1.90;1.32];
                case 'germanate'
                    Omega = [7.83;6.37;2.05];
                case 'tellurite'
                    Omega = [4.33;2.49;0.97];
            end
    end
end
%% Compute spontaneous transition rates Arad required for the latter rate-equation gain computations
transition_lambda = 1e4./abs(doping_data.energy(level_idx)'-doping_data.energy(level_idx)); % um
transition_n = sqrt(JO_Sellmeier_coeff(1) + ...
                    JO_Sellmeier_coeff(2)*(transition_lambda).^2./((transition_lambda).^2 - JO_Sellmeier_coeff(3)) + ...
                    JO_Sellmeier_coeff(4)*(transition_lambda).^2./((transition_lambda).^2 - JO_Sellmeier_coeff(5)));
transition_n(abs(imag(transition_n)) > 0) = 1.45; % For super-small transition wavelengths where two energy levels are too close to each other, transition_n can be weird; hence, I just set it to n=1.45 in silica

chi = ( (transition_n.^2 + 2)/3 ).^2; % a factor in calculations

Sfit = doping_data.s1(transition_idx).*Omega(1) + ...
       doping_data.s2(transition_idx).*Omega(2) + ...
       doping_data.s3(transition_idx).*Omega(3); % line strengths by the Judd-Ofelt equation
Aed = (7.235432e10*transition_n.*chi).*(Sfit*1e-20)./(Arad_J.*(transition_lambda*1e-4).^3); % electric-dipole transition probability
Amd = 2.697348e10*transition_n.^3./(Arad_J.*(transition_lambda*1e3).^3).*doping_data.s4(transition_idx); % magnetic-dipole transition probability
Arad = Aed + Amd; % 1/s; (num_levels,num_levels); total spontaneous transition rates

%% Nonradiative multiphonon decay rates
[Bel,alpha,hw] = multiphonon_decay_parameters(gain_rate_eqn.base_medium);
dE = doping_data.energy(level_idx(2:end)) - doping_data.energy(level_idx(1:end-1)); % cm^(-1)
Gamma = Bel.*exp(-alpha*(dE-2*hw)); % 1/s; multiphonon decay rates of each level (decaying to its nearest lower level)

%% Name of the considered energy levels (just for verbose)
energy_levels = doping_data.layere(level_idx)';

end