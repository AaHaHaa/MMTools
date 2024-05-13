function [a,b] = Sellmeier_coefficients(material)
%SELLMEIER_COEFFICIENTS It loads the Sellmeier coefficients based on the
%material.
%
%   Use:
%       Sellmeier_terms = @(lambda,a,b) a.*lambda.^2./(lambda.^2 - b.^2);
%       n_from_Sellmeier = @(lambda) sqrt(1+sum(Sellmeier_terms(lambda,a,b),2));
%       n = n_from_Sellmeier(lambda); % lambda: um

switch material
    case 'fused silica' % From "https://www.osapublishing.org/josa/abstract.cfm?uri=josa-55-10-1205",
                        % I. H. Malitson, "Interspecimen Comparison of the Refractive Index of Fused Silica"
                        % n^2 - 1 = ...
        a = [0.6961663, 0.4079426, 0.8974794];
        b = [0.0684043, 0.1162414, 9.896161];
    case 'As2S3' % From "Supercontinuum generation in photonic crystal fibers made from highly nonlinear glasses" by A.V. Husakou
                 % n^2 - 1 = ...
        a = [1.898367, 1.922297, 0.87651, 0.11887, 0.95699];
        b = [0.15,     0.25,     0.35,    0.45,    27.3861];
    case 'H2' % From E. R. Peck and S. Hung. Refractivity and dispersion of hydrogen in the visible and near infrared, J. Opt. Soc. Am. 67, 1550-1554 (1977)
              % n - 1 = ...
              % At 0 degree Celsius temperature and 101325 Pa
        a =      [0.0148956/180.7, 0.0049037/92];
        b = sqrt([        1/180.7,         1/92]);
    case 'Xe' % From A. Börzsönyi,Z. Heiner, M. P. Kalashnikov, A. P. Kovács, and K. Osvay, Dispersion measurement of inert gases and gas mixtures at 800 nm, Appl. Opt. 47, 4856-4863 (2008)
              % n^2 - 1 = ...
              % At 0 degree Celsius temperature and 101325 Pa
        a =      [103701.61e-8, 31228.61e-8];
        b = sqrt([12.75e-3,     0.561e-3]);
    case 'air' % Check "refractiveindex.info"
               % From P. E. Ciddor. Refractive index of air: new equations for the visible and near infrared, Appl. Optics 35, 1566-1573 (1996)
               % n - 1 = ...
        a =      [0.05792105/238.0185, 0.00167917/57.362];
        b = sqrt([         1/238.0185,          1/57.362]);
    case 'N2' % Check "refractiveindex.info"
              % From E. R. Peck and B. N. Khanna. Dispersion of nitrogen, J. Opt. Soc. Am. 56, 1059-1063 (1966)
              % n - 1 = ...
              % At 0 degree Celsius temperature and 101325 Pa
        a =      [6.8552e-5, 3.243157e-2/144];
        b = sqrt([        0,           1/144]);
    %{
    case 'N2' % Check "refractiveindex.info"
              % From A. Börzsönyi, Z. Heiner, M. P. Kalashnikov, A. P. Kovács, and K. Osvay, Dispersion measurement of inert gases and gas mixtures at 800 nm, Appl. Opt. 47, 4856-4863 (2008)
              % n^2 - 1 = ...
              % At 0 degree Celsius temperature and 101325 Pa
        a =      [39209.95e-8, 18806.48e-8];
        b = sqrt([ 1146.24e-6,   13.476e-3]);
    %}
    case 'O2' % Check "refractiveindex.info"
              % From J. Zhang, Z. H. Lu, and L. J. Wang. Precision refractive index measurements of air, N2, O2, Ar, and CO2 with a frequency comb, Appl. Opt. 47, 3143-3151 (2008)
              %     and
              %      Petr Křen. Comment on "Precision refractive index measurements of air, N2, O2, Ar, and CO2 with a frequency comb", Appl. Opt. 50, 6484-6485 (2011)
              % n - 1 = ...
              % The paper measured with 20 degree Celsius temperature and 101325 Pa, I recovered it back to 0 degree.
        a =      [1.181494e-4, 9.708931e-3/75.4]*293.15/273.15;
        b = sqrt([        0,             1/75.4]);
    case 'Ar' % From EDSON R. PECK AND DONALD J. FIISIER, "dispersion of Argon" (1964)
              % n-1 = ...
              % At 0 degree temperature and 101325 Pa
        a =      [6.786711e-5, 3.0182943e-2/144];
        b = sqrt([          0,            1/144]);
    case 'Kr' % From A. Börzsönyi,Z. Heiner, M. P. Kalashnikov, A. P. Kovács, and K. Osvay, Dispersion measurement of inert gases and gas mixtures at 800 nm, Appl. Opt. 47, 4856-4863 (2008)
              % n^2 - 1 = ...
              % At 0 degree Celsius temperature and 101325 Pa
        a =      [26102.88e-8, 56946.82e-8];
        b = sqrt([    2.01e-6,   10.043e-3]);
    case 'Ne' % From A. Börzsönyi,Z. Heiner, M. P. Kalashnikov, A. P. Kovács, and K. Osvay, Dispersion measurement of inert gases and gas mixtures at 800 nm, Appl. Opt. 47, 4856-4863 (2008)
              % n^2 - 1 = ...
              % At 0 degree Celsius temperature and 101325 Pa
        a =      [9154.48e-8, 4018.63e-8];
        b = sqrt([656.97e-6,    5.728e-3]);
    case 'He' % From A. Börzsönyi,Z. Heiner, M. P. Kalashnikov, A. P. Kovács, and K. Osvay, Dispersion measurement of inert gases and gas mixtures at 800 nm, Appl. Opt. 47, 4856-4863 (2008)
              % n^2 - 1 = ...
              % At 0 degree Celsius temperature and 101325 Pa
        a =      [4977.77e-8, 1856.94e-8];
        b = sqrt([28.54e-6,     7.760e-3]);
    case 'CH4' % From R. RoLLEFsoN AND R. HAVENs, Index of Refraction of Methane in the Infra-Red and the Diyole Momentof the CH Bond (1940)
               % n^2 - 1 = ...
               % At 0 degree Celsius temperature and 101325 Pa
               %
               % In the paper, it says that at an infinite wavelength, n-1 = 435.9e-6, so n^2-1 = 871.9e-6.
               % The resonance contributions are transformed from K/(vi^2-v), where K has the unit of cm^(-2) and v has cm^(-1).
               % I realized that the paper is wrong such that it forgot a factor of 2 in the resonance contributions. Its given values are for n-1; for n^2-1, they should be approximately twice as large.
        a = [871.9e-6, 2*16.50e-8/(3019.6e-4)^2, 2*7.5e-8./(1304e-4)^2];
        b = [        0,            1/3019.6e-4,           1/1304e-4];
    case 'sapphire' % From https://refractiveindex.info/?shelf=3d&book=crystals&page=sapphire
                    % n^2 - 1 = ...
                    % At 20 degree Celsius
        a = [1.4313493, 0.65054713, 5.3414021];
        b = [0.0726631,  0.1193242, 18.028251];
    case 'BK7' % From https://refractiveindex.info/?shelf=glass&book=BK7&page=SCHOTT
                    % n^2 - 1 = ...
                    % At 20 degree Celsius
        a =      [   1.03961212,  0.231792344, 1.01046945];
        b = sqrt([0.00600069867, 0.0200179144, 103.560653]);
    case 'N-SF1' % From https://refractiveindex.info/?shelf=glass&book=SCHOTT-SF&page=N-SF1
                    % n^2 - 1 = ...
                    % At 20 degree Celsius
        a =      [  1.60865158,  0.237725916, 1.51530653];
        b = sqrt([0.0119654879, 0.0590589722, 135.521676]);
    case 'N-SF10' % From https://refractiveindex.info/?shelf=glass&book=SF10&page=SCHOTT
                    % n^2 - 1 = ...
                    % At 20 degree Celsius
        a =      [  1.62153902,  0.256287842, 1.64447552];
        b = sqrt([0.0122241457, 0.0595736775, 147.468793]);
    case 'N-SF11' % From https://refractiveindex.info/?shelf=glass&book=SF11&page=SCHOTT
                    % n^2 - 1 = ...
                    % At 20 degree Celsius
        a =      [ 1.73759695,   0.313747346, 1.89878101];
        b = sqrt([0.013188707,  0.0623068142,  155.23629]);
    case 'N-SF14' % From https://refractiveindex.info/?shelf=glass&book=SCHOTT-SF&page=N-SF14
                    % n^2 - 1 = ...
                    % At 20 degree Celsius
        a =      [  1.69022361,  0.288870052, 1.7045187];
        b = sqrt([0.0130512113, 0.061369188, 149.517689]);
    case 'CaF2' % From https://refractiveindex.info/?shelf=main&book=CaF2&page=Malitson
                    % n^2 - 1 = ...
                    % At 20 degree Celsius
        a = [  0.5675888, 0.4710914, 3.8484723];
        b = [0.050263605, 0.1003909, 34.649040];
    case 'N-LAK21' % From https://refractiveindex.info/?shelf=glass&book=SCHOTT-LaK&page=N-LAK21
                    % n^2 - 1 = ...
                    % At 20 degree Celsius
        a =      [   1.22718116,  0.420783743, 1.01284843];
        b = sqrt([0.00602075682, 0.0196862889, 88.4370099]);
    otherwise
        error('Material isn''t in this repository yet.');
end

end