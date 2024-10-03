function kijkl = population_nonlinear_terms(gain_medium,doped_wt)
%POPULATION_NONLINEAR_TERMS Loads nonlinear terms considered in population 
%evolutions

% Upconversion factors should be proportional to C^2/(C^2+C0^2), where C is
% the doping concentration and C0 is the characteristic concentration.
%
% Gunnar Rustad and Knut Stenersen, "Modeling of Laser-Pumped Tm and Ho
% lasers Accounting for Upconversion and Ground-State Depletion," IEEE
% Journal of Quantum Electronics 32(9), 1645-1656 (1996).

switch gain_medium
    case 'Er'
        % Henderson-Sapir et al., "New energy-transfer upconversion process
        % in Er3+:ZBLAN mid-infrared fiber lasers," Opt. Express 24,
        % 6869-6883 (2016)
        k1103 = 0.6e-23; % m^3/s
        k2206 = 0.2e-23; % m^3/s
        k5031 = 0.38e-23; % m^3/s
        k4251 = 1.55e-23; % m^3/s
        
        kijkl = [k1103,k2206,k5031,k4251]*1e18; % change them to "um^3/s" for latter computations
    case 'Nd'
        % Skrzypczak et al., "Comprehensive Rate Equation Analysis of
        % Upconversion Luminescence Enhancement Due to BaCl2 Nanocrystals
        % in Neodymium-Doped Fluorozirconate-Based Glass Ceramics," J.
        % Phys. Chem. C 118(24), 13087-13098 (2014)
        k4438 = 39.2e-24; % m^3/s
        k4429 = (53.2 + 15.1)*1e-24; % m^3/s
        k44110 = (46.1 + 16.3)*1e-24; % m^3/s
        k44011 = (59.8 + 1.3)*1e-24; % m^3/s
        k4033 = 5.6e-24; % m^3/s
        
        kijkl = [k4438,k4429,k44110,k44011,k4033]*1e18; % change them to "um^3/s" for latter computations
    case 'Tm'
        % [1]
        % Kamradek et al., "Energy transfer coefficients in thulium-doped
        % silica fibers," Opt. Mat. Express 11(6), 1805-1814 (2021)
        % [2]
        % Jackson and King, "Theoretical Modeling of Tm-Doped Silica Fiber
        % Lasers," Journal of Lightwave Technology 17(5), 948-956 (1999)
        % [3]
        % Tao et al., "Cross Relaxation in Tm-doped Fiber Lasers," 2nd
        % International Symposium on Laser Interaction with Matter (LIMIS
        % 2012) 8789 (2013)
        C0 = 4.3; % wt.%
        wt_factor = doped_wt^2/(doped_wt^2+C0^2);
        k3101 = 2.8e-22*wt_factor; % 1.8e-22 m^3/s at 5.75 wt.% [1-3]
        k1310 = 0.084*k3101; % m^3/s [1-3]
        k2101 = 1.089e-22*wt_factor; % 7e-23 m^3/s at 5.75 w.% [1]
        k1012 = 0.5*k2101; % m^3/s [1]
        k5203 = 1.56e-23*wt_factor; % 1e-23 m^3/s at 5.75wt.% [1]
        k5631 = 9.33e-23*wt_factor; % 6.0e-23 m^3/s at 5.75wt.% [1]
        k5751 = 1.87e-23*wt_factor; % 1.2e-23 m^3/s at 5.75wt.% [1]
        k5313 = 2.49e-23*wt_factor; % 1.6e-23 m^3/s at 5.75wt.% [1]
        
        kijkl = [k3101,k1310,k2101,k1012,k5203,k5631,k5751,k5313]*1e18; % change them to "um^3/s" for latter computations
    case 'Ho'
        % Huang et al., "Theoretical Modeling of Ho-Doped Fiber Lasers
        % Pumped by Laser-Diodes Around 1.125 um," J. Light. Technol.
        % 30(20), 3235-3240 (2012)
        k2101 = 2e-24; % m^3/s
        k1012 = 4e-23; % m^3/s
        k3101 = 7.8e-22; % m^3/s
        k1013 = 2.3e-24; % m^3/s
        
        kijkl = [k2101,k1012,k3101,k1013]*1e18; % change them to "um^3/s" for latter computations
end

