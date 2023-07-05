function [fiber_type,...
          core_diameter,clad_diameter,...
          core_NA,clad_NA,...
          fname_user_defined,...
          alpha] = fiber_collections(fiber,lambda0)
%FIBER_COLLECTIONS Stores several pre-defined fiber parameters from
%Thorlabs

alpha = 0; % Shape parameter; only for graded-index fibers

switch fiber
    case '1060XP' % Thorlabs SM fiber
        fiber_type = 'step';
        core_diameter = 5.8;
        clad_diameter = 125; % coating diameter : 245um
        core_NA = 0.14;
        clad_NA = 0.22;
    case 'YB1200-4_125' % Thorlabs core-pumped SM active fiber
        fiber_type = 'step';
        core_diameter = 4.4; % given only the MFD at 1060nm
        clad_diameter = 125; % coating diameter : 245um
        core_NA = 0.2;
        clad_NA = 0.22;
    case {'YB1200-6_125DC','YB1200-6_125DC-PM'} % Thorlabs cladding-pumped non-PM and PM active fibers
        fiber_type = 'step';
        core_diameter = 7;
        clad_diameter = 125; % coating diameter : 245um
        core_NA = 0.12;
        clad_NA = 0.22;
    case {'YB1200-10_125DC','YB1200-10_125DC-PM'} % Thorlabs cladding-pumped non-PM and PM active fiber
        fiber_type = 'step';
        core_diameter = 10;
        clad_diameter = 125; % coating diameter : 245um
        core_NA = 0.08;
        clad_NA = 0.22;
    case 'YB1200-20_400DC' % Thorlabs cladding-pumped active fiber
        fiber_type = 'step';
        core_diameter = 20;
        clad_diameter = 400; % coating diameter : 520um
        core_NA = 0.07;
        clad_NA = 0.22;
    case 'YB1200-25_250DC' % Thorlabs cladding-pumped active fiber
        fiber_type = 'step';
        core_diameter = 25;
        clad_diameter = 250; % coating diameter : 350um
        core_NA = 0.07;
        clad_NA = 0.22;
    case 'YB1200-25_250DC-PM' % Thorlabs cladding-pumped PM active fiber
        fiber_type = 'step';
        core_diameter = 25;
        clad_diameter = 250; % coating diameter : 350um
        core_NA = 0.062;
        clad_NA = 0.22;
    case {'ER30-4_125','ER110-4_125'} % Thorlabs Liekki SM cladding-pumped active fiber
                                      % Feature: 
                                      %    ER30-4_125:  Extremely high, >50% conversion efficiency in the L band
                                      %    ER110-4_125: Extremely high doping concentration for short device length and reduced nonlinearity
        fiber_type = 'step';
        core_diameter = 3.5; % given only the MFD at 1550nm
        clad_diameter = 125; % coating diameter : 245um
        core_NA = 0.2;
        clad_NA = 0.22;
    case {'ER16-8_125','ER80-8_125'} % Thorlabs Liekki LMA cladding-pumped active fiber
                                      % Feature: 
                                      %    ER16-8_125: Good spliceability, power conversion efficiency, and spectral reproducibility
                                      %    ER80-8_125: For 980 nm pumps with emission at 1550 nm. Large core and good spliceability.
        fiber_type = 'step';
        core_diameter = 8; % given only the MFD at 1550nm
        clad_diameter = 125; % coating diameter : 245um
        core_NA = 0.13;
        clad_NA = 0.22;
    case 'M5-980-125' % Thorlabs MetroGain SM cladding-pumped active fiber
                      % Feature: effective for high-power C-Band use (1530 - 1565 nm) when pumped at 1480 nm
                      % Emission wavelength: C band
        fiber_type = 'step';
        core_diameter = 4; % given only the MFD at 1550nm: 5.5-6.3um
        clad_diameter = 125; % coating diameter : 245um
        core_NA = 0.23; % 0.21-0.24
        clad_NA = 0.22;
    case 'M12-980-125' % Thorlabs MetroGain SM cladding-pumped active fiber
                       % Feature: optimized for L-band emission with a 980 nm pump source. Its high absorption allows for shorter active fiber lengths compared to conventional Er-doped fibers emitting in the L band.
                       % Emission wavelength: L band
        fiber_type = 'step';
        core_diameter = 4; % given only the MFD at 1550nm: 5.7-6.6um
        clad_diameter = 125; % coating diameter : 245um
        core_NA = 0.23; % 0.21-0.24
        clad_NA = 0.22;
    case 'OM1' % Thorlabs graded-index multimode fiber
        fiber_type = 'GRIN';
        core_diameter = 62.5;
        clad_diameter = 125; % coating diameter : 245um
        core_NA = 0.272;
        clad_NA = 0.22;
        alpha = 2.08; % Shape parameter
    case {'OM2','OM3','OM4'} % Thorlabs graded-index multimode fiber
        fiber_type = 'GRIN';
        core_diameter = 50;
        clad_diameter = 125; % coating diameter : 245um
        core_NA = 0.2;
        clad_NA = 0.22;
        alpha = 2.08; % Shape parameter
    case 'PLMA-YDF-30_400-VIII' % Nufern Yb-doped LMA fiber
        fiber_type = 'step';
        core_diameter = 30;
        clad_diameter = 400; % coating diameter: 550um
        core_NA = 0.06; % inner core NA
        clad_NA = 0.22;
    case 'PLMA-YDF-25_400-VIII' % Nufern Yb-doped LMA fiber
        fiber_type = 'step';
        core_diameter = 25;
        clad_diameter = 400; % coating diameter: 550um
        core_NA = 0.060; % inner core NA
        clad_NA = 0.22;
end

fname_user_defined = sprintf('%s_wavelength%4unm',fiber,lambda0*1e9);

end