function [Bel,alpha,hw] = multiphonon_decay_parameters(material)
%MULTIPHONON_DECAY_PARAMETERS It computes the required parameters for
%finding the multiphonon decay rates of each energy levels used in
%rate-equation gain computations.
%   
% tellurite, phosphate, borate, silicate, germanate, als, gls, zbla are from
%
%   Reisfeld and Eyal, "POSSIBLE WAYS OF RELAXATIONS FOR EXCITED STATES OF 
%   RARE EARTH IONS IN AMORPHOUS MEDIA," Journal de Physique Colloques 46 
%   (C7), 349-355 (1985).
%
% sio2 and yag(Y2O3) are from
%
%   Yu et al., "The Temperature-Dependence of Multiphonon Relaxation of
%   Rare-Earth Ions in Solid-State Hosts," J. Phys. Chem. C Nanomater
%   Interfaces 120(18), 9958-9964 (2016).
%
% "find_parameters.m" in the "Multiphonon decay rate" folder also has some 
% fitting results from 
%
%   Wetenkamp, et al., "Optical properties of rare earth-doped ZBLAN 
%   glasses," Journal of Non-Crystalline Solids 140, 35-40 (1992).

switch lower(material)
    case 'tellurite'
        B = 6.3e10; % 1/s
        alpha = 4.7e-3; % cm
        hw = 700; % cm^(-1)
    case 'phosphate'
        B = 5.4e12; % 1/s
        alpha = 4.7e-3; % cm
        hw = 1200; % cm^(-1)
    case 'borate'
        B = 2.9e12; % 1/s
        alpha = 3.8e-3; % cm
        hw = 1400; % cm^(-1)
    case {'silicate'}
        B = 5.4e12; % 1/s
        alpha = 4.7e-3; % cm
        hw = 1200; % cm^(-1)
    case 'germanate'
        B = 3.4e10; % 1/s
        alpha = 4.9e-3; % cm
        hw = 900; % cm^(-1)
    case {'als','gls'}
        B = 1e10; % 1/s
        alpha = 2.9e-3; % cm
        hw = 350; % cm^(-1)
    case 'zbla'
        B = 1.88e10; % 1/s
        alpha = 5.77e-3; % cm
        hw = 460; % cm^(-1)
    case 'silica' % SiO2
        B = 1.4e12; % 1/s
        alpha = 3.8e-3; % cm
        hw = 1100; % cm^(-1)
    case 'zblan'
        B = 1.89e5; % 1/s
        alpha = 2.08e-3; % cm
        hw = 550; % cm^(-1)
    case 'yag'
        % phonon energy is taken from
        % Jeffrey Eldridge, "Luminescence decay-based Y2O3:Er phosphor
        % thermometry: Temperature sensitivity governed by multiphonon
        % emission with an effective phonon energy transition," Journal of
        % Luminescence 214, 116535 (2019)
        B = 1.204e8; % 1/s
        alpha = 3.53e-3; % cm
        hw = 600; % cm^(-1)
end

Bel = B*exp(-2*alpha*hw); % 1/s

end

