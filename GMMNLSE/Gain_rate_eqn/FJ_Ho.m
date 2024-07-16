function func = FJ_Ho
%FJ_HO container for the coupled equations of Ho population evolutions and
%their Jacobian matrix
%
% 1. Huang et al., "Theoretical Modeling of Ho-Doped Fiber Lasers Pumped
%    by Laser-Diodes Around 1.125 um," J. Light. Technol. 30(20), 3235-3240
%    (2012)
% 2. Wang et al., "Numerical Modeling of in-Band Pumped Ho-Doped Silica
%    Fiber Lasers," J. Light. Technol. 36(24), 5863-5880 (2018)
% 3. Wang et al., "Effects of ion clustering and excited state absorption
%    on the performance of Ho-doped fiber lasers," Opt. Express 27(10),
%    14283-14297 (2019)
% 4. Khamis and Alabbasi, "Influence of upconversion transition process on
%    the performance of inband pumped holmium-doped silica fiber 
%    amplifiers," Eur. Phys. J. D 75(1), 1-12 (2021)
% 5. Alharbi et al., "Performance Optimization of Holmium Doped Fiber
%    Amplifiers for Optical Communication Applications in 2–2.15 um
%    Wavelength Range," Photonics 9(4) (2022)

func.F = @F;
func.J = @J;

end

%% Contained functions
function dNdt = F(N,N_total,A,Gamma,kijkl,R)
%F Population coupled equations among Ho's energy levels
%
% Input arguments:
%   N: population; (4,1,num_x,num_y,1,1,1,M)
%   N_total: total population; a (1,1,num_x,num_y) array
%   A: spontaneous transition rate/probabilities; a (5,5) matrix
%   Gamma: multiphonon decay rates; a (4,1) vector
%   kijkl: nonlinear coupling coefficients for upconversion or cross relaxation effects; a (1,4) vector
%   R: photon rate; (1,1,num_x,num_y,1,1,1,M,1,num_cross_sections)
%
% Output argument:
%   dNdt: variation of populations; a (4,1,num_x,num_y,1,1,1,M) array

% Populations of each energy level
N0 = N(1,:,:,:,:,:,:,:); % ground state, 5I8
N1 = N(2,:,:,:,:,:,:,:); % 5I7
N2 = N(3,:,:,:,:,:,:,:); % 5I6
N3 = N(4,:,:,:,:,:,:,:); % 5I5
N4 = N_total - sum(N,1); % 5I4

% Nonlinear coupling coefficients
k2101 = kijkl(1);
k1012 = kijkl(2);
k3101 = kijkl(3);
k1013 = kijkl(4);

% Stimulated photon rates (absorption and emission)
wa04 = R(:,:,:,:,:,:,:,:,:,1);
wa03 = R(:,:,:,:,:,:,:,:,:,2);
wa02 = R(:,:,:,:,:,:,:,:,:,3);
wa01 = R(:,:,:,:,:,:,:,:,:,4);
wa14 = R(:,:,:,:,:,:,:,:,:,5);
we10 = R(:,:,:,:,:,:,:,:,:,6);

% Variation of populations
% dN0 + dN1 + dN2 + dN3 + dN4 = dN_total = 0
dN0dt = A(2,1)*N1+A(3,1)*N2+A(4,1)*N3+A(5,1)*N4 ...
        + Gamma(1)*N1 ...
        + we10.*N1 ...
        - (wa01 + wa02 + wa03 + wa04).*N0 ...
        - k2101*N2.*N0 + k1012*N1.^2 - k3101*N3.*N0 + k1013*N1.^2;
dN1dt = A(3,2)*N2+A(4,2)*N3+A(5,2)*N4 ...
        + Gamma(2)*N2 ...
        - (A(2,1) + Gamma(1))*N1 ...
        - we10.*N1 ...
        + wa01.*N0 - wa14.*N1 ...
        + 2*k2101*N2.*N0 - 2*k1012*N1.^2 + 2*k3101*N3.*N0 - 2*k1013*N1.^2;
dN2dt = A(4,3)*N3+A(5,3)*N4 ...
        + Gamma(3)*N3 ...
        - (A(3,1) + A(3,2) + Gamma(2))*N2 ...
        + wa02.*N0 ...
        - k2101*N2.*N0 + k1012*N1.^2;
dN3dt = A(5,4)*N4 ...
        + Gamma(4)*N4 ...
        - (A(4,1) + A(4,2) + A(4,3) + Gamma(3))*N3 ...
        + wa03.*N0 ...
        - k3101*N3.*N0 + k1013*N1.^2;
%dN4dt = - (A(5,1) + A(5,2) + A(5,3) + A(5,4) + Gamma(4))*N4 ...
%        + wa04.*N0 ...
%        + wa14.*N1;

dNdt = cat(1,dN0dt,dN1dt,dN2dt,dN3dt);

end

function Jacobian_Ho = J(N,A,Gamma,kijkl,R)
%J Jacobian matrix of the population coupled equations among Ho's energy
%levels
%
% Input arguments:
%   N: population, but it's a ratio here over N_total; (4,1,num_x,num_y,1,1,1,M)
%   A: spontaneous transition rate/probabilities; a (5,5) matrix
%   Gamma: multiphonon decay rates; a (4,1) vector
%   kijkl: nonlinear coupling coefficients for upconversion or cross relaxation effects; a (1,4) vector
%   R: photon rate; (1,1,num_x,num_y,1,1,1,M,1,num_cross_sections)
%
% Output argument:
%   Jacobian_Ho: Jacobian matrix; a (4,4,num_x,num_y,1,1,1,M) array

% Populations of each energy level
N0 = N(1,:,:,:,:,:,:,:); % ground state, 5I8
N1 = N(2,:,:,:,:,:,:,:); % 5I7
N2 = N(3,:,:,:,:,:,:,:); % 5I6
N3 = N(4,:,:,:,:,:,:,:); % 5I5
%N4 = N_total - sum(N,1); % 5I4

% Nonlinear coupling coefficients
k2101 = kijkl(1);
k1012 = kijkl(2);
k3101 = kijkl(3);
k1013 = kijkl(4);

% Stimulated photon rates (absorption and emission)
wa04 = R(:,:,:,:,:,:,:,:,:,1);
wa03 = R(:,:,:,:,:,:,:,:,:,2);
wa02 = R(:,:,:,:,:,:,:,:,:,3);
wa01 = R(:,:,:,:,:,:,:,:,:,4);
wa14 = R(:,:,:,:,:,:,:,:,:,5);
we10 = R(:,:,:,:,:,:,:,:,:,6);

% Jacobian matrix
dN0dt_N0 = -A(5,1) ...
           - (wa01 + wa02 + wa03 + wa04) ...
           - k2101*N2 - k3101*N3;
dN0dt_N1 = A(2,1)-A(5,1) ...
           + Gamma(1) ...
           + we10 ...
           + 2*k1012*N1 + 2*k1013*N1;
dN0dt_N2 = A(3,1)-A(5,1) ...
           - k2101*N0;
dN0dt_N3 = A(4,1)-A(5,1) ...
           - k3101*N0;
dN1dt_N0 = -A(5,2) ...
           + wa01 ...
           + 2*k2101*N2 + 2*k3101*N3;
dN1dt_N1 = -A(5,2) ...
           - (A(2,1) + Gamma(1)) ...
           - we10 ...
           - wa14 ...
           - 4*k1012*N1 - 4*k1013*N1;
dN1dt_N2 = A(3,2)-A(5,2) ...
           + Gamma(2) ...
           + 2*k2101*N0;
dN1dt_N3 = A(4,2)-A(5,2) ...
           + 2*k3101*N0;
dN2dt_N0 = -A(5,3) ...
           + wa02 ...
           - k2101*N2;
dN2dt_N1 = -A(5,3) ...
           + 2*k1012*N1;
dN2dt_N2 = -A(5,3) ...
           - (A(3,1) + A(3,2) + Gamma(2)) ...
           - k2101*N0;
dN2dt_N3 = A(4,3)-A(5,3) ...
           + Gamma(3);
dN3dt_N0 = -A(5,4) ...
           - Gamma(4) ...
           + wa03 ...
           - k2101*N0;
dN3dt_N1 = -A(5,4) ...
           - Gamma(4) ...
           + 2*k1013*2*N1;
dN3dt_N2 = -A(5,4) ...
           - Gamma(4);
dN3dt_N3 = -A(5,4) ...
           - Gamma(4) ...
           - (A(4,1) + A(4,2) + A(4,3) + Gamma(3)) ...
           - k3101*N0;

% Some computations don't involve photon rates, so they have
% less-dimensional arrays, which makes array catenation into one Jacobian
% matrix fail. Extension of their dimensions is thus required.
num_x = size(R,3);
num_y = size(R,4);
M = size(R,8);
if size(N,8) ~= M % the number of parallelizations
    % Before extension, they're (1,1,num_x,num_y,1,1,1,1).
    dN0dt_N2 = repmat(dN0dt_N2,1,1,1,1,1,1,1,M);
    dN0dt_N3 = repmat(dN0dt_N3,1,1,1,1,1,1,1,M);
    dN1dt_N2 = repmat(dN1dt_N2,1,1,1,1,1,1,1,M);
    dN1dt_N3 = repmat(dN1dt_N3,1,1,1,1,1,1,1,M);
    dN2dt_N1 = repmat(dN2dt_N1,1,1,1,1,1,1,1,M);
    dN2dt_N2 = repmat(dN2dt_N2,1,1,1,1,1,1,1,M);
    dN3dt_N1 = repmat(dN3dt_N1,1,1,1,1,1,1,1,M);
    dN3dt_N3 = repmat(dN3dt_N3,1,1,1,1,1,1,1,M);
end
% Before extension, they're scalars with the dimension (1,1,1,1,1,1,1,1).
dN2dt_N3 = repmat(dN2dt_N3,1,1,num_x,num_y,1,1,1,M);
dN3dt_N2 = repmat(dN3dt_N2,1,1,num_x,num_y,1,1,1,M);

Jacobian_Ho = cat(1,cat(2,dN0dt_N0,dN0dt_N1,dN0dt_N2,dN0dt_N3),...
                    cat(2,dN1dt_N0,dN1dt_N1,dN1dt_N2,dN1dt_N3),...
                    cat(2,dN2dt_N0,dN2dt_N1,dN2dt_N2,dN2dt_N3),...
                    cat(2,dN3dt_N0,dN3dt_N1,dN3dt_N2,dN3dt_N3));

end