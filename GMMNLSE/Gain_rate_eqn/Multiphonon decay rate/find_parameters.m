% This script finds the required parameters for calculating the multiphonon
% decay rates.
%
% Its data is from
% Wetenkamp, et al., "Optical properties of rare earth-doped ZBLAN
% glasses," Journal of Non-Crystalline Solids 140, 35-40 (1992).

close all; clearvars;

% ZBLA
data_ZBLA = readmatrix('ZBLA.csv');
linefitcoeff_ZBLA = line_fitting( data_ZBLA(:,1),log(data_ZBLA(:,2)) );
fprintf('ZBLA:\n');
fprintf('B=%6.4e(1/s), alpha=%1.2e(cm)\n',exp(linefitcoeff_ZBLA(2)), -linefitcoeff_ZBLA(1));

% ZBLAN
data_ZBLAN = readmatrix('ZBLAN.csv');
linefitcoeff_ZBLAN = line_fitting( data_ZBLAN(:,1),log(data_ZBLAN(:,2)) );
fprintf('ZBLAN:\n');
fprintf('B=%6.4e(1/s), alpha=%1.2e(cm)\n',exp(linefitcoeff_ZBLAN(2)), -linefitcoeff_ZBLAN(1));

% YAG
data_YAG = readmatrix('YAG.csv');
linefitcoeff_YAG = line_fitting( data_YAG(:,1),log(data_YAG(:,2)) );
fprintf('YAG:\n');
fprintf('B=%6.4e(1/s), alpha=%1.2e(cm)\n',exp(linefitcoeff_YAG(2)), -linefitcoeff_YAG(1));

% SiO2
data_SiO2 = readmatrix('SiO2.csv');
linefitcoeff_SiO2 = line_fitting( data_SiO2(:,1),log(data_SiO2(:,2)) );
fprintf('SiO2/silica:\n');
fprintf('B=%6.4e(1/s), alpha=%1.2e(cm)\n',exp(linefitcoeff_SiO2(2)), -linefitcoeff_SiO2(1));