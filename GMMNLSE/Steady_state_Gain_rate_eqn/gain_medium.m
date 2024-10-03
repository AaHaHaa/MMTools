function func = gain_medium()
%GAIN_MEDIUM This loads relevant parameters for a specific gain medium.

func.load_medium_parameters = @load_medium_parameters;
func.load_N_related_parameters = @load_N_related_parameters;
func.load_cross_sections = @load_cross_sections;

end

%% Helper functions
function gain_rate_eqn = load_medium_parameters(gain_rate_eqn)
%LOAD_MEDIUM_PARAMETERS Load some settings and parameters

switch gain_rate_eqn.gain_medium
    case 'Yb'
        % Because I have two cross-sectional data for Yb, from Nufern and
        % from Liekki, I open the option for the user to select which file
        % to use. It's default to using Liekki's, which is consistent with
        % Thorlabs fibers our group mostly use.
        if isfield(gain_rate_eqn,'cross_section_filename') && ...
            contains(gain_rate_eqn.cross_section_filename,'Nufern')
            gain_rate_eqn.cross_section_filename = 'Yb_Gen_VIII_Cross_Section (Nufern).txt';
        else
            gain_rate_eqn.cross_section_filename = 'Liekki Yb_AV_20160530.txt';
        end
        gain_rate_eqn.N.eqn.tau = 840e-6; % lifetime of Yb in F_(5/2) state (Paschotta et al., "Lifetime quenching in Yb-doped fibers"); in "s"
        gain_rate_eqn.base_medium = 'silica'; % Yb fibers are basically made of silica
        gain_rate_eqn.energy_levels = {'2F7/2','2F5/2'}; % from low to high energy
    case 'Er'
        gain_rate_eqn.cross_section_filename = 'Er.txt';
        %gain_rate_eqn.tau = 9e-3; % lifetime of Er in (^4)I_(13/2) state; (Z. Y. Zhang et al., "Fluorescence decay-time characteristics of erbium-doped optical fiber at elevated temperatures")
        
        [gain_rate_eqn.N.eqn.Arad,gain_rate_eqn.N.eqn.Gammai,gain_rate_eqn.energy_levels] = calc_Judd_Ofelt(gain_rate_eqn);
        
        gain_rate_eqn.N.eqn.ss = FJ_Er();
    case 'Nd'
        gain_rate_eqn.cross_section_filename = 'Nd.txt';
        if ~isfield(gain_rate_eqn,'base_medium')
            gain_rate_eqn.base_medium = 'silica'; % Nd fibers are basically made of silica
        end
        
        [gain_rate_eqn.N.eqn.Arad,gain_rate_eqn.N.eqn.Gammai,gain_rate_eqn.energy_levels] = calc_Judd_Ofelt(gain_rate_eqn);
        
        gain_rate_eqn.N.eqn.ss = FJ_Nd();
    case 'Tm'
        gain_rate_eqn.cross_section_filename = 'Tm.txt';
        
        [gain_rate_eqn.N.eqn.Arad,gain_rate_eqn.N.eqn.Gammai,gain_rate_eqn.energy_levels] = calc_Judd_Ofelt(gain_rate_eqn);
         
        gain_rate_eqn.N.eqn.ss = FJ_Tm();
    case 'Ho'
        gain_rate_eqn.cross_section_filename = 'Ho.txt';
        %gain_rate_eqn.tau = 1.9e-3; % lifetime of Ho in 5I_7 state (Gouet et al., "Realization and simulation of high-power holmium doped fiber lasers for long-range transmission"); in "s"
        
        [gain_rate_eqn.N.eqn.Arad,gain_rate_eqn.N.eqn.Gammai,gain_rate_eqn.energy_levels] = calc_Judd_Ofelt(gain_rate_eqn);
        
        gain_rate_eqn.N.eqn.ss = FJ_Ho();
end

end

function gain_rate_eqn = load_N_related_parameters(gain_rate_eqn,N_total)
%LOAD_N_RELATED_PARAMETERS Load some settings and parameters

silica_density = 2.2e-12; % g/um^3; silica density
silica_atomic_mass = 60.0835; % a.u.; SiO2 = Si 28.0855 a.u. + (O 15.999 a.u.)*2
Avogadro_number = 6.02214076e23;
N_silica = silica_density/silica_atomic_mass*Avogadro_number; % 1/um^3; silica number density
ratio_Silica2DopedIon = N_silica/max(N_total(:)); % ratio of number density of doped rare ions to silica (in a silica-based fiber)

switch gain_rate_eqn.gain_medium
    case 'Yb'
        gain_rate_eqn.N.eqn.kijkl = [];
    case 'Er'
        gain_rate_eqn.N.eqn.kijkl = population_nonlinear_terms(gain_rate_eqn.gain_medium,[]);
    case 'Nd'
        gain_rate_eqn.N.eqn.kijkl = population_nonlinear_terms(gain_rate_eqn.gain_medium,[]);
    case 'Tm'
        % For Coherent/Nufern's LMA-TDF-25P/400, LMA-TDF-25P-250, and
        % SM-TDF-10P-130, the Tm is doped to be ~5 wt.%.
        % From:
        % Richard Tumminelli, Vincent Petit, Adrian Carter, Alexander Hemming, Nikita Simakov, John Haub, 
        % "Highly doped and highly efficient Tm doped fiber laser (Conference Presentation)," 
        % Proc. SPIE 10512, Fiber Lasers XV: Technology and Systems, 105120M (14 March 2018); 
        % https://doi.org/10.1117/12.2295936
        % (It has an online video, which I highly recommended to watch.)
        %
        % This conference talk reports the highest Tm doping up to 8.5
        % wt.% to achieve 93.2% quantum efficiency with the 74.5% slope 
        % efficiency through the cross-relaxation effect.
        doped_atomic_mass = 168.93421; % a.u.; Tm atomic mass
        % From
        % Jollivet et al, "HIGH PERFORMANCE LARGE-MODE AREA DOUBLE-CLAD FIBERS FOR KW POWER SCALING OF FIBER LASERS FROM 1 TO 2 MICRONS," 7th International Workshop on Specialty Optical Fibers (2022),
        % Nufern's 25P/400 HC should have 8 wt.% Tm in the core  with 5.2
        % dB/m cladding absorption at 793 nm. This factor calibrates my
        % calculation value to fit this number. I still don't understand my
        % calculation deviates from it.
        calibration_ratio = 1.523;
        doped_ion_wt = doped_atomic_mass/(doped_atomic_mass + silica_atomic_mass*ratio_Silica2DopedIon)*100*calibration_ratio; % weight percent; wt.%
        
        gain_rate_eqn.N.eqn.kijkl = population_nonlinear_terms(gain_rate_eqn.gain_medium,doped_ion_wt);
    case 'Ho'
        gain_rate_eqn.N.eqn.kijkl = population_nonlinear_terms(gain_rate_eqn.gain_medium,[]);
end

end

%% Read cross sections from the file
function [gain_rate_eqn,...
          cross_sections,cross_sections_pump,...
          GSA_find_Ntotal] = load_cross_sections(gain_rate_eqn,lambda)
%LOAD_CROSS_SECTIONS Load cross-section data and some relevant parameters

necessary_lambda = [gain_rate_eqn.absorption_wavelength_to_get_N_total*1e-9; ...
                    gain_rate_eqn.pump_wavelength*1e-9; ...
                    lambda*1e-9];
switch gain_rate_eqn.gain_medium
    case 'Yb'
        [GSA01,emi10] = read_cross_sections_2levels(gain_rate_eqn.cross_section_filename,necessary_lambda); % read the file
        GSA01 = GSA01*1e12; % change the unit to um^2
        emi10 = emi10*1e12;
        
        cross_sections_pump = struct('GSA01',GSA01(2),'emi10',emi10(2)); % pump
        cross_sections = struct('GSA01',GSA01(3:end),'emi10',emi10(3:end)); % signal, ASE
        
        gain_rate_eqn.plusminus = [-1,1];
        gain_rate_eqn.N_idx = [1,2];
        
        GSA_find_Ntotal = GSA01(1);
    case 'Er'
        [GSA08,GSA07,GSA06,GSA05,GSA04,GSA03,GSA02,GSA01,...
         ESA15,ESA26,ESA12,ESA23,...
         emi10,emi21,emi32] = read_cross_sections_Er(gain_rate_eqn.cross_section_filename,necessary_lambda); % read the file
        GSA08 = GSA08*1e12; % change the unit to um^2
        GSA07 = GSA07*1e12;
        GSA06 = GSA06*1e12;
        GSA05 = GSA05*1e12;
        GSA04 = GSA04*1e12;
        GSA03 = GSA03*1e12;
        GSA02 = GSA02*1e12;
        GSA01 = GSA01*1e12;
        ESA15 = ESA15*1e12;
        ESA26 = ESA26*1e12;
        ESA12 = ESA12*1e12;
        ESA23 = ESA23*1e12;
        emi10 = emi10*1e12;
        emi21 = emi21*1e12;
        emi32 = emi32*1e12;
        
        cross_sections_pump = struct('GSA08',GSA08(2),'GSA07',GSA07(2),'GSA06',GSA06(2),'GSA05',GSA05(2),'GSA04',GSA04(2),'GSA03',GSA03(2),'GSA02',GSA02(2),'GSA01',GSA01(2),...
                                     'ESA15',ESA15(2),'ESA26',ESA26(2),'ESA12',ESA12(2),'ESA23',ESA23(2),...
                                     'emi10',emi10(2),'emi21',emi21(2),'emi32',emi32(2)); % pump
        cross_sections = struct('GSA08',GSA08(3:end),'GSA07',GSA07(3:end),'GSA06',GSA06(3:end),'GSA05',GSA05(3:end),'GSA04',GSA04(3:end),'GSA03',GSA03(3:end),'GSA02',GSA02(3:end),'GSA01',GSA01(3:end),...
                                'ESA15',ESA15(3:end),'ESA26',ESA26(3:end),'ESA12',ESA12(3:end),'ESA23',ESA23(3:end),...
                                'emi10',emi10(3:end),'emi21',emi21(3:end),'emi32',emi32(3:end)); % signal, ASE
        
        % Signs for the gain calculations.
        % 1 for emission and -1 for absorption
        gain_rate_eqn.plusminus = [-1,-1,-1,-1,-1,-1,-1,-1,...
                                   -1,-1,-1,-1,...
                                    1, 1, 1];
        % Indices to select the correct population to multiply with the
        % cross sections during gain calculations.
        gain_rate_eqn.N_idx = [1,1,1,1,1,1,1,1,...
                               2,3,2,3,...
                               2,3,4];
        
        GSA_find_Ntotal = GSA01(1) + GSA02(1) + GSA03(1) + GSA04(1) + GSA05(1) + GSA06(1) + GSA07(1) + GSA08(1);
    case 'Nd'
        [GSA011,GSA010,GSA09,GSA08,GSA07,GSA06,GSA05,GSA04,GSA03,...
         ESA49,...
         emi40,emi41,emi42] = read_cross_sections_Nd(gain_rate_eqn.cross_section_filename,necessary_lambda); % read the file
        GSA011 = GSA011*1e12; % change the unit to um^2
        GSA010 = GSA010*1e12;
        GSA09 = GSA09*1e12;
        GSA08 = GSA08*1e12;
        GSA07 = GSA07*1e12;
        GSA06 = GSA06*1e12;
        GSA05 = GSA05*1e12;
        GSA04 = GSA04*1e12;
        GSA03 = GSA03*1e12;
        ESA49 = ESA49*1e12;
        emi40 = emi40*1e12;
        emi41 = emi41*1e12;
        emi42 = emi42*1e12;
        
        cross_sections_pump = struct('GSA011',GSA011(2),'GSA010',GSA010(2),'GSA09',GSA09(2),'GSA08',GSA08(2),'GSA07',GSA07(2),'GSA06',GSA06(2),'GSA05',GSA05(2),'GSA04',GSA04(2),'GSA03',GSA03(2),...
                                     'ESA49',ESA49(2),...
                                     'emi40',emi40(2),'emi41',emi41(2),'emi42',emi42(2)); % pump
        cross_sections = struct('GSA011',GSA011(3:end),'GSA010',GSA010(3:end),'GSA09',GSA09(3:end),'GSA08',GSA08(3:end),'GSA07',GSA07(3:end),'GSA06',GSA06(3:end),'GSA05',GSA05(3:end),'GSA04',GSA04(3:end),'GSA03',GSA03(3:end),...
                                'ESA49',ESA49(3:end),...
                                'emi40',emi40(3:end),'emi41',emi41(3:end),'emi42',emi42(3:end)); % signal, ASE

        % Signs for the gain calculations.
        % 1 for emission and -1 for absorption
        gain_rate_eqn.plusminus = [-1,-1,-1,-1,-1,-1,-1,-1,-1,...
                                   -1,...
                                    1, 1, 1];
        % Indices to select the correct population to multiply with the
        % cross sections during gain calculations.
        gain_rate_eqn.N_idx = [1,1,1,1,1,1,1,1,1,...
                               5,...
                               5,5,5];
        
        GSA_find_Ntotal = GSA03(1) + GSA04(1) + GSA05(1) + GSA06(1) + GSA07(1) + GSA08(1) + GSA09(1) + GSA010(1) + GSA011(1);
    case 'Tm'
        [GSA05,GSA04,GSA03,GSA02,GSA01,...
         ESA14,ESA35,ESA13,ESA23,...
         emi50,emi30,emi31,emi10,emi32] = read_cross_sections_Tm(gain_rate_eqn.cross_section_filename,necessary_lambda); % read the file
        GSA05 = GSA05*1e12; % change the unit to um^2
        GSA04 = GSA04*1e12;
        GSA03 = GSA03*1e12;
        GSA02 = GSA02*1e12;
        GSA01 = GSA01*1e12;
        ESA14 = ESA14*1e12;
        ESA35 = ESA35*1e12;
        ESA13 = ESA13*1e12;
        ESA23 = ESA23*1e12;
        emi50 = emi50*1e12;
        emi30 = emi30*1e12;
        emi31 = emi31*1e12;
        emi10 = emi10*1e12;
        emi32 = emi32*1e12;
        
        cross_sections_pump = struct('GSA05',GSA05(2),'GSA04',GSA04(2),'GSA03',GSA03(2),'GSA02',GSA02(2),'GSA01',GSA01(2),...
                                     'ESA14',ESA14(2),'ESA35',ESA35(2),'ESA13',ESA13(2),'ESA23',ESA23(2),...
                                     'emi50',emi50(2),'emi30',emi30(2),'emi31',emi31(2),'emi10',emi10(2),'emi32',emi32(2)); % pump
        cross_sections = struct('GSA05',GSA05(3:end),'GSA04',GSA04(3:end),'GSA03',GSA03(3:end),'GSA02',GSA02(3:end),'GSA01',GSA01(3:end),...
                                'ESA14',ESA14(3:end),'ESA35',ESA35(3:end),'ESA13',ESA13(3:end),'ESA23',ESA23(3:end),...
                                'emi50',emi50(3:end),'emi30',emi30(3:end),'emi31',emi31(3:end),'emi10',emi10(3:end),'emi32',emi32(3:end)); % signal, ASE
        
        % Signs for the gain calculations.
        % 1 for emission and -1 for absorption
        gain_rate_eqn.plusminus = [-1,-1,-1,-1,-1,...
                                   -1,-1,-1,-1,...
                                    1, 1, 1, 1, 1];
        % Indices to select the correct population to multiply with the
        % cross sections during gain calculations.
        gain_rate_eqn.N_idx = [1,1,1,1,1,...
                               2,4,2,3,...
                               6,4,4,2,4];
        
        GSA_find_Ntotal = GSA01(1) + GSA02(1) + GSA03(1) + GSA04(1) + GSA05(1);
    case 'Ho'
        [GSA04,GSA03,GSA02,GSA01,...
         ESA14,...
         emi10] = read_cross_sections_Ho(gain_rate_eqn.cross_section_filename,necessary_lambda); % read the file
        GSA04 = GSA04*1e12; % change the unit to um^2
        GSA03 = GSA03*1e12;
        GSA02 = GSA02*1e12;
        GSA01 = GSA01*1e12;
        ESA14 = ESA14*1e12;
        emi10 = emi10*1e12;
        
        cross_sections_pump = struct('GSA04',GSA04(2),'GSA03',GSA03(2),'GSA02',GSA02(2),'GSA01',GSA01(2),...
                                     'ESA14',ESA14(2),...
                                     'emi10',emi10(2)); % pump
        cross_sections = struct('GSA04',GSA04(3:end),'GSA03',GSA03(3:end),'GSA02',GSA02(3:end),'GSA01',GSA01(3:end),...
                                'ESA14',ESA14(3:end),...
                                'emi10',emi10(3:end)); % signal, ASE
        
        % Signs for the gain calculations.
        % 1 for emission and -1 for absorption
        gain_rate_eqn.plusminus = [-1,-1,-1,-1,...
                                   -1,...
                                    1];
        % Indices to select the correct population to multiply with the
        % cross sections during gain calculations.
        gain_rate_eqn.N_idx = [1,1,1,1,...
                               2,...
                               2];
        
        GSA_find_Ntotal = GSA01(1) + GSA02(1) + GSA03(1) + GSA04(1);
end

% Computational dimension:
% (Nx,Nx,num_spatial_modes,num_spatial_modes,Nt,M,num_polarization,num_cross_sections),
% M: parallelization in MPA (if used)
cross_sections = structfun(@(x) permute(x,[2 3 4 5 1]),cross_sections,'UniformOutput',false); % change it to the size (1,1,1,1,Nt)
if ~isempty(gain_rate_eqn.plusminus)
    gain_rate_eqn.plusminus = permute(gain_rate_eqn.plusminus,[1,3,4,5,6,7,8,2]); % this is for cross-sections computations, so put it to the 8th dimension
end

end