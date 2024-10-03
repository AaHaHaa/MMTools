function func = FJ_Tm()
%FJ_TM container for the coupled equations of Tm population evolutions and
%their Jacobian matrix
%
% 1. Rustad and Stenersen, "Modeling of laser-pumped Tm and Ho lasers
%    accounting for upconversion and ground-state depletion," IEEE J.
%    Quantum Electron. 32(9), 1645-1656 (1996)
% 2. Jackson and King, "Efficient gain-switched operation of a Tm-doped
%    silica fiber laser," IEEE J. Quantum Electron. 34(5), 779-789 (1998)
% 3. Jackson and King, "Theoretical modeling of Tm-doped silica fiber
%    lasers," J. Light. Technol. 17(5), 948-956 (1999)
% 4. Jackson and Mossman, "Efficiency dependence on the Tm^3+ and Al^3+
%    concentrations for Tm^3+-doped silica double-clad fiber lasers," Appl.
%    Opt. 42(15), 2702-2707 (2003)
% 5. Zhao et al., "The slope efficiency of 2 μm thulium doped fiber laser,"
%    Proc.SPIE 7843, 784306 (2010)
% 6. Peterka et al., "Theoretical modeling of fiber laser at 810 nm based
%    on thulium-doped silica fibers with enhanced ^3H_4 level lifetime,"
%    Opt. Express 19(3), 2773-2781 (2011)
% 7. Kamradek et al., "Energy transfer coefficients in thulium-doped silica
%    fibers," Opt. Mat. Express 11(6), 1805-1814 (2021)
% 8. Louot et al., "Emission Wavelength Limits of a Continuous-Wave
%    Thulium-Doped Fiber Laser Source Operating at 1.94 μm, 2.09 μm or 2.12
%    μm," Photonics 11(3) (2024)

func.F = @F;
func.J = @J;

end

%% Contained functions
function dNdt = F(N,N_total,A,Gamma,kijkl,R)
%F Population coupled equations among Tm's energy levels
%
% Input arguments:
%   N: population; (7,1,num_x,num_y,1,1,1,M)
%   N_total: total population; a (1,1,num_x,num_y) array
%   A: spontaneous transition rate/probabilities; a (8,8) matrix
%   Gamma: multiphonon decay rates; a (7,1) vector
%   kijkl: nonlinear coupling coefficients for upconversion or cross relaxation effects; a (1,8) vector
%   R: photon rate; (1,1,num_x,num_y,1,1,1,M,1,num_cross_sections)
%
% Output argument:
%   dNdt: variation of populations; a (7,1,num_x,num_y,1,1,1,M) array

% Populations of each energy level
N0 = N_total - sum(N,1); % ground state, 3H6
N1 = N(1,:,:,:,:,:,:,:); % 3F4
N2 = N(2,:,:,:,:,:,:,:); % 3H5
N3 = N(3,:,:,:,:,:,:,:); % 3H4
N4 = N(4,:,:,:,:,:,:,:); % 3F2 + 3F3
N5 = N(5,:,:,:,:,:,:,:); % 1G4
N6 = N(6,:,:,:,:,:,:,:); % 1D2
N7 = N(7,:,:,:,:,:,:,:); % 1I6

% Nonlinear coupling coefficients
k3101 = kijkl(1);
k1310 = kijkl(2);
k2101 = kijkl(3);
k1012 = kijkl(4);
k5203 = kijkl(5);
k5631 = kijkl(6);
k5751 = kijkl(7);
k5313 = kijkl(8);

% Stimulated photon rates (absorption and emission)
wa05 = R(:,:,:,:,:,:,:,:,:,1);
wa04 = R(:,:,:,:,:,:,:,:,:,2);
wa03 = R(:,:,:,:,:,:,:,:,:,3);
wa02 = R(:,:,:,:,:,:,:,:,:,4);
wa01 = R(:,:,:,:,:,:,:,:,:,5);
wa14 = R(:,:,:,:,:,:,:,:,:,6);
wa35 = R(:,:,:,:,:,:,:,:,:,7);
wa13 = R(:,:,:,:,:,:,:,:,:,8);
wa23 = R(:,:,:,:,:,:,:,:,:,9);
we50 = R(:,:,:,:,:,:,:,:,:,10);
we30 = R(:,:,:,:,:,:,:,:,:,11);
we31 = R(:,:,:,:,:,:,:,:,:,12);
we10 = R(:,:,:,:,:,:,:,:,:,13);
we32 = R(:,:,:,:,:,:,:,:,:,14);

% Variation of populations
% dN0 + dN1 + dN2 + dN3 + dN4 + dN5 + dN6 + dN7 = dN_total = 0
%dN0dt = A(2,1)*N1+A(3,1)*N2+A(4,1)*N3+A(5,1)*N4+A(6,1)*N5+A(7,1)*N6+A(8,1)*N7 ...
%        + Gamma(1)*N1 ...
%        + we10.*N1 + we30.*N3 + we50.*N5 ...
%        - (wa01 + wa02 + wa03 + wa04 + wa05).*N0 ...
%        - k3101*N3.*N0 + k1310*N1.^2 - k2101*N2.*N0 + k1012*N1.^2 - k5203*N5.*N0;
dN1dt = A(3,2)*N2+A(4,2)*N3+A(5,2)*N4+A(6,2)*N5+A(7,2)*N6+A(8,2)*N7 ...
        + Gamma(2)*N2 ...
        - (A(2,1) + Gamma(1))*N1 ...
        - we10.*N1 + we31.*N3 ...
        + wa01.*N0 - (wa13 + wa14).*N1 ...
        + 2*k3101*N3.*N0 - 2*k1310*N1.^2 + 2*k2101*N2.*N0 - 2*k1012*N1.^2 + k5631*N5.*N3 + k5751*N5.^2 - k5313*N5.*N1;
dN2dt = A(4,3)*N3+A(5,3)*N4+A(6,3)*N5+A(7,3)*N6+A(8,3)*N7 ...
      + Gamma(3)*N3 ...
      - (A(3,1) + A(3,2) + Gamma(2))*N2 ...
      + we32.*N3 ...
      + wa02.*N0 - wa23.*N2 ...
      - k2101*N2.*N0 + k1012*N1.^2 + k5203*N5.*N0;
dN3dt = A(5,4)*N4+A(6,4)*N5+A(7,4)*N6+A(8,4)*N7 ...
        + Gamma(4)*N4 ...
        - (A(4,1) + A(4,2) + A(4,3) + Gamma(3))*N3 ...
        - (we30 + we31 + we32).*N3 ...
        + wa03.*N0 + wa13.*N1 + wa23.*N2 - wa35.*N3 ...
        - k3101*N3.*N0 + k1310*N1.^2 - k5631*N5.*N3 + k5313*N5.*N1 + k5203*N5.*N0;
dN4dt = A(6,5)*N5+A(7,5)*N6+A(8,5)*N7 ...
      + Gamma(5)*N5 ...
      - (A(5,1) + A(5,2) + A(5,3) + A(5,4) + Gamma(4))*N4 ...
      + wa04.*N0 + wa14.*N1 ...
      + k5313*N5.*N1;
dN5dt = A(7,6)*N6+A(8,6)*N7 ...
        + Gamma(6)*N6 ...
        - (A(6,1) + A(6,2) + A(6,3) + A(6,4) + A(6,5) + Gamma(5))*N5 ...
        - we50.*N5 ...
        + wa05.*N0 + wa35.*N3 ...
        - k5631*N5.*N3 - 2*k5751*N5.^2 - k5313*N5.*N1 - k5203*N5.*N0;
dN6dt = A(8,7)*N7 ...
        + Gamma(7)*N7 ...
        - (A(7,1) + A(7,2) + A(7,3) + A(7,4) + A(7,5) + A(7,6) + Gamma(6))*N6 ...
        + k5631*N5.*N3;
dN7dt = - (A(8,1) + A(8,2) + A(8,3) + A(8,4) + A(8,5) + A(8,6) + A(8,7) + Gamma(7))*N7 ...
        + k5751*N5.^2;

% Some computations don't involve photon rates, so they have
% less-dimensional arrays, which makes array catenation into one dNdt
% vector fail. Extension of their dimensions is thus required.
if size(N,8) ~= size(R,8) % the number of parallelizations
    dN6dt = repmat(dN6dt,1,1,1,1,1,1,1,size(R,8));
    dN7dt = repmat(dN7dt,1,1,1,1,1,1,1,size(R,8));
end

dNdt = cat(1,dN1dt,dN2dt,dN3dt,dN4dt,dN5dt,dN6dt,dN7dt);

end

function Jacobian_Tm = J(N,N_total,A,Gamma,kijkl,R)
%J Jacobian matrix of the population coupled equations among Tm's energy
%levels
%
% Input arguments:
%   N: population, but it's a ratio here over N_total; (7,1,num_x,num_y,1,1,1,M)
%   A: spontaneous transition rate/probabilities; a (8,8) matrix
%   Gamma: multiphonon decay rates; a (7,1) vector
%   kijkl: nonlinear coupling coefficients for upconversion or cross relaxation effects; a (1,8) vector
%   R: photon rate; (1,1,num_x,num_y,1,1,1,M,1,num_cross_sections)
%
% Output argument:
%   Jacobian_Tm: Jacobian matrix; a (7,7,num_x,num_y,1,1,1,M) array

% Populations of each energy level
N0 = N_total - sum(N,1); % ground state, 3H6
N1 = N(1,:,:,:,:,:,:,:); % 3F4
N2 = N(2,:,:,:,:,:,:,:); % 3H5
N3 = N(3,:,:,:,:,:,:,:); % 3H4
%N4 = N(4,:,:,:,:,:,:,:); % 3F2 + 3F3
N5 = N(5,:,:,:,:,:,:,:); % 1G4
%N6 = N(6,:,:,:,:,:,:,:); % 1D2
%N7 = N(7,:,:,:,:,:,:,:); % 1I6

% Nonlinear coupling coefficients
k3101 = kijkl(1);
k1310 = kijkl(2);
k2101 = kijkl(3);
k1012 = kijkl(4);
k5203 = kijkl(5);
k5631 = kijkl(6);
k5751 = kijkl(7);
k5313 = kijkl(8);

% Stimulated photon rates (absorption and emission)
wa05 = R(:,:,:,:,:,:,:,:,:,1);
wa04 = R(:,:,:,:,:,:,:,:,:,2);
wa03 = R(:,:,:,:,:,:,:,:,:,3);
wa02 = R(:,:,:,:,:,:,:,:,:,4);
wa01 = R(:,:,:,:,:,:,:,:,:,5);
wa14 = R(:,:,:,:,:,:,:,:,:,6);
wa35 = R(:,:,:,:,:,:,:,:,:,7);
wa13 = R(:,:,:,:,:,:,:,:,:,8);
wa23 = R(:,:,:,:,:,:,:,:,:,9);
we50 = R(:,:,:,:,:,:,:,:,:,10);
we30 = R(:,:,:,:,:,:,:,:,:,11);
we31 = R(:,:,:,:,:,:,:,:,:,12);
we10 = R(:,:,:,:,:,:,:,:,:,13);
we32 = R(:,:,:,:,:,:,:,:,:,14);

% Jacobian matrix
dN1dt_N1 = - A(2,1) - Gamma(1) - wa01 - wa13 - wa14 - we10 - 4*N1*k1012 - 4*N1*k1310 - 2*N2*k2101 - 2*N3*k3101 - N5*k5313;
dN1dt_N2 = A(3,2) + Gamma(2) - wa01 - 2*N3*k3101 + 2*k2101*(N0 - N2);
dN1dt_N3 = A(4,2) - wa01 + we31 - 2*N2*k2101 + N5*k5631 + 2*k3101*(N0 - N3);
dN1dt_N4 = A(5,2) - wa01 - 2*N2*k2101 - 2*N3*k3101;
dN1dt_N5 = A(6,2) - wa01 - 2*N2*k2101 - 2*N3*k3101 - N1*k5313 + N3*k5631 + 2*N5*k5751;
dN1dt_N6 = A(7,2) - wa01 - 2*N2*k2101 - 2*N3*k3101;
dN1dt_N7 = A(8,2) - wa01 - 2*N2*k2101 - 2*N3*k3101;
dN2dt_N1 = 2*N1*k1012 - wa02 + N2*k2101 - N5*k5203;
dN2dt_N2 = (N2 - N0)*k2101 - A(3,2) - Gamma(2) - wa02 - wa23 - A(3,1) - N5*k5203;
dN2dt_N3 = A(4,3) + Gamma(3) - wa02 + we32 + N2*k2101 - N5*k5203;
dN2dt_N4 = A(5,3) - wa02 + N2*k2101 - N5*k5203;
dN2dt_N5 = A(6,3) - wa02 + N2*k2101 + k5203*(N0 - N5);
dN2dt_N6 = A(7,3) - wa02 + N2*k2101 - N5*k5203;
dN2dt_N7 = A(8,3) - wa02 + N2*k2101 - N5*k5203;
dN3dt_N1 = wa13 - wa03 + 2*N1*k1310 + N3*k3101 - N5*k5203 + N5*k5313;
dN3dt_N2 = wa23 - wa03 + N3*k3101 - N5*k5203;
dN3dt_N3 = - A(4,2) - A(4,3) - Gamma(3) - wa03 - wa35 - we30 - we31 - we32 - A(4,1) + (N3 - N0)*k3101 - N5*k5203 - N5*k5631;
dN3dt_N4 = A(5,4) + Gamma(4) - wa03 + N3*k3101 - N5*k5203;
dN3dt_N5 = A(6,4) - wa03 + N3*k3101 + N1*k5313 - N3*k5631 + k5203*(N0 - N5);
dN3dt_N6 = A(7,4) - wa03 + N3*k3101 - N5*k5203;
dN3dt_N7 = A(8,4) - wa03 + N3*k3101 - N5*k5203;
dN4dt_N1 = wa14 - wa04 + N5*k5313;
dN4dt_N2 = - wa04;
dN4dt_N3 = - wa04;
dN4dt_N4 = - A(5,1) - A(5,2) - A(5,3) - A(5,4) - Gamma(4) - wa04;
dN4dt_N5 = A(6,5) + Gamma(5) - wa04 + N1*k5313;
dN4dt_N6 = A(7,5) - wa04;
dN4dt_N7 = A(8,5) - wa04;
dN5dt_N1 = N5*k5203 - wa05 - N5*k5313;
dN5dt_N2 = N5*k5203 - wa05;
dN5dt_N3 = wa35 - wa05 + N5*k5203 - N5*k5631;
dN5dt_N4 = N5*k5203 - wa05;
dN5dt_N5 = - A(6,2) - A(6,3) - A(6,4) - A(6,5) - Gamma(5) - wa05 - we50 - A(6,1) + (N5 - N0)*k5203 - N1*k5313 - N3*k5631 - 4*N5*k5751;
dN5dt_N6 = A(7,6) + Gamma(6) - wa05 + N5*k5203;
dN5dt_N7 = A(8,6) - wa05 + N5*k5203;
%dN6dt_N1 = 0;
%dN6dt_N2 = 0;
dN6dt_N3 = N5*k5631;
%dN6dt_N4 = 0;
dN6dt_N5 = N3*k5631;
dN6dt_N6 = - A(7,1) - A(7,2) - A(7,3) - A(7,4) - A(7,5) - A(7,6) - Gamma(6);
dN6dt_N7 = A(8,7) + Gamma(7);
%dN7dt_N1 = 0;
%dN7dt_N2 = 0;
%dN7dt_N3 = 0;
%dN7dt_N4 = 0;
dN7dt_N5 = 2*N5*k5751;
%dN7dt_N6 = 0;
dN7dt_N7 = - A(8,1) - A(8,2) - A(8,3) - A(8,4) - A(8,5) - A(8,6) - A(8,7) - Gamma(7);

% Some computations don't involve photon rates, so they have
% less-dimensional arrays, which makes array catenation into one Jacobian
% matrix fail. Extension of their dimensions is thus required.
num_x = size(R,3);
num_y = size(R,4);
M = size(R,8);
if size(N,8) ~= M % the number of parallelizations
    % Before extension, they're (1,1,num_x,num_y,1,1,1,1).
    dN6dt_N3 = repmat(dN6dt_N3,1,1,1,1,1,1,1,M);
    dN6dt_N5 = repmat(dN6dt_N5,1,1,1,1,1,1,1,M);
    dN7dt_N5 = repmat(dN7dt_N5,1,1,1,1,1,1,1,M);
end
% Before extension, they're scalars with the dimension (1,1,1,1,1,1,1,1).
dN6dt_N1 = zeros(1,1,num_x,num_y,1,1,1,M); % = repmat(dN6dt_N1,1,1,num_x,num_y,1,1,1,M);
dN6dt_N2 = zeros(1,1,num_x,num_y,1,1,1,M); % = repmat(dN6dt_N2,1,1,num_x,num_y,1,1,1,M);
dN6dt_N4 = zeros(1,1,num_x,num_y,1,1,1,M); % = repmat(dN6dt_N4,1,1,num_x,num_y,1,1,1,M);
dN6dt_N6 = repmat(dN6dt_N6,1,1,num_x,num_y,1,1,1,M);
dN6dt_N7 = repmat(dN6dt_N7,1,1,num_x,num_y,1,1,1,M);
dN7dt_N1 = zeros(1,1,num_x,num_y,1,1,1,M); % = repmat(dN7dt_N1,1,1,num_x,num_y,1,1,1,M);
dN7dt_N2 = zeros(1,1,num_x,num_y,1,1,1,M); % = repmat(dN7dt_N2,1,1,num_x,num_y,1,1,1,M);
dN7dt_N3 = zeros(1,1,num_x,num_y,1,1,1,M); % = repmat(dN7dt_N3,1,1,num_x,num_y,1,1,1,M);
dN7dt_N4 = zeros(1,1,num_x,num_y,1,1,1,M); % = repmat(dN7dt_N4,1,1,num_x,num_y,1,1,1,M);
dN7dt_N6 = zeros(1,1,num_x,num_y,1,1,1,M); % = repmat(dN7dt_N6,1,1,num_x,num_y,1,1,1,M);
dN7dt_N7 = repmat(dN7dt_N7,1,1,num_x,num_y,1,1,1,M);

Jacobian_Tm = cat(1,cat(2,dN1dt_N1,dN1dt_N2,dN1dt_N3,dN1dt_N4,dN1dt_N5,dN1dt_N6,dN1dt_N7),...
                    cat(2,dN2dt_N1,dN2dt_N2,dN2dt_N3,dN2dt_N4,dN2dt_N5,dN2dt_N6,dN2dt_N7),...
                    cat(2,dN3dt_N1,dN3dt_N2,dN3dt_N3,dN3dt_N4,dN3dt_N5,dN3dt_N6,dN3dt_N7),...
                    cat(2,dN4dt_N1,dN4dt_N2,dN4dt_N3,dN4dt_N4,dN4dt_N5,dN4dt_N6,dN4dt_N7),...
                    cat(2,dN5dt_N1,dN5dt_N2,dN5dt_N3,dN5dt_N4,dN5dt_N5,dN5dt_N6,dN5dt_N7),...
                    cat(2,dN6dt_N1,dN6dt_N2,dN6dt_N3,dN6dt_N4,dN6dt_N5,dN6dt_N6,dN6dt_N7),...
                    cat(2,dN7dt_N1,dN7dt_N2,dN7dt_N3,dN7dt_N4,dN7dt_N5,dN7dt_N6,dN7dt_N7));

end