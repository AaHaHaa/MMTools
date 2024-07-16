function func = FJ_Nd
%FJ_ND container for the coupled equations of Nd population evolutions and
%their Jacobian matrix
%
% 1. Boley et al., "E-band neodymium-doped fiber amplifier: model and
%    application," Appl. Opt. 58(9), 2320-2327 (2019)
% 2. Wu et al., "Modeling and Numerical Simulation of the Gain Spectrum of
%    1250nm-1350nm Neodymium-doped Broadband Fibre Amplifiers," Highlights
%    in Science, Engineering and Technology 46, 266-271 (2023)
% 3. Skrzypczak et al., "Comprehensive Rate Equation Analysis of 
%    Upconversion Luminescence Enhancement Due to BaCl2 Nanocrystals in
%    Neodymium-Doped Fluorozirconate-Based Glass Ceramics," J. Phys. Chem. 
%    C 118(24), 13087-13098 (2014)
% 4. Verlinden et al., "The Excited State Absorption Cross Section of
%    Neodymium-doped Silica Glass Fiber in the 1200-1500 nm Wavelength 
%    Range," Worcester Polytechnic Institute, PhD Thesis (2008)

func.F = @F;
func.J = @J;

end

%% Contained functions
function dNdt = F(N,N_total,A,Gamma,kijkl,R)
%F Population coupled equations among Nd's energy levels
%
% Input arguments:
%   N: population; (11,1,num_x,num_y,1,1,1,M)
%   N_total: total population; a (1,1,num_x,num_y) array
%   A: spontaneous transition rate/probabilities; a (12,12) matrix
%   Gamma: multiphonon decay rates; a (11,1) vector
%   kijkl: nonlinear coupling coefficients for upconversion or cross relaxation effects; a (1,5) vector
%   R: photon rate; (1,1,num_x,num_y,1,1,1,M,1,num_cross_sections)
%
% Output argument:
%   dNdt: variation of populations; a (11,1,num_x,num_y,1,1,1,M) array

% Populations of each energy level
N0  = N( 1,:,:,:,:,:,:,:); % ground state, 4I9/2
N1  = N( 2,:,:,:,:,:,:,:); % 4I11/2
N2  = N( 3,:,:,:,:,:,:,:); % 4I13/2
N3  = N( 4,:,:,:,:,:,:,:); % 4I15/2
N4  = N( 5,:,:,:,:,:,:,:); % 4F3/2
N5  = N( 6,:,:,:,:,:,:,:); % (2H9/2 +) 4F5/2
N6  = N( 7,:,:,:,:,:,:,:); % (4S3/2 +) 4F7/2
N7  = N( 8,:,:,:,:,:,:,:); % 4F9/2
N8  = N( 9,:,:,:,:,:,:,:); % (2G7/2 + 4G5/2 +) 2H11/2
N9  = N(10,:,:,:,:,:,:,:); % (4G9/2 + 4G7/2 +) 2K13/2
N10 = N(11,:,:,:,:,:,:,:); % (4G11/2 + 2(D,P)3/2 + 2G9/2 +) 2K15/2
N11 = N_total - sum(N,1); % (2D5/2 +) 2P1/2

% Nonlinear coupling coefficients
k4438   = kijkl(1);
k4429  = kijkl(2);
k440110 = kijkl(3);
k440011 = kijkl(4);
k4033   = kijkl(5);

% Stimulated photon rates (absorption and emission)
wa0011 = R(:,:,:,:,:,:,:,:,:, 1);
wa0010 = R(:,:,:,:,:,:,:,:,:, 2);
wa09   = R(:,:,:,:,:,:,:,:,:, 3);
wa08   = R(:,:,:,:,:,:,:,:,:, 4);
wa07   = R(:,:,:,:,:,:,:,:,:, 5);
wa06   = R(:,:,:,:,:,:,:,:,:, 6);
wa05   = R(:,:,:,:,:,:,:,:,:, 7);
wa04   = R(:,:,:,:,:,:,:,:,:, 8);
wa03   = R(:,:,:,:,:,:,:,:,:, 9);
wa49   = R(:,:,:,:,:,:,:,:,:,10);
we40   = R(:,:,:,:,:,:,:,:,:,11);
we41   = R(:,:,:,:,:,:,:,:,:,12);
we42   = R(:,:,:,:,:,:,:,:,:,13);

% Variation of populations
% dN0 + dN1 + dN2 + dN3 + dN4 + dN5 + dN6 + dN7 + dN8 + dN9 + dN10 + dN11 = dN_total = 0
dN0dt = A(2,1)*N1+A(3,1)*N2+A(4,1)*N3+A(5,1)*N4+A(6,1)*N5+A(7,1)*N6+A(8,1)*N7+A(9,1)*N8+A(10,1)*N9+A(11,1)*N10+A(12,1)*N11 ...
        + Gamma(1)*N1 ...
        + we40.*N4 ...
        - (wa03 + wa04 + wa05 + wa06 + wa07 + wa08 + wa09 + wa0010 + wa0011).*N0 ...
        + k440011*N4.^2 - k4033*N0.*N4;
dN1dt = A(3,2)*N2+A(4,2)*N3+A(5,2)*N4+A(6,2)*N5+A(7,2)*N6+A(8,2)*N7+A(9,2)*N8+A(10,2)*N9+A(11,2)*N10+A(12,2)*N11 ...
        + Gamma(2)*N2 ...
        - (A(2,1) + Gamma(1))*N1 ...
        + we41.*N4 ...
        + k440110*N4.^2;
dN2dt = A(4,3)*N3+A(5,3)*N4+A(6,3)*N5+A(7,3)*N6+A(8,3)*N7+A(9,3)*N8+A(10,3)*N9+A(11,3)*N10+A(12,3)*N11 ...
        + Gamma(3)*N3 ...
        - (A(3,1) + A(3,2) + Gamma(2))*N2 ...
        + we42.*N4 ...
        + k4429*N4.^2;
dN3dt = A(5,4)*N4+A(6,4)*N5+A(7,4)*N6+A(8,4)*N7+A(9,4)*N8+A(10,4)*N9+A(11,4)*N10+A(12,4)*N11 ...
        + Gamma(4)*N4 ...
        - (A(4,1) + A(4,2) + A(4,3) + Gamma(3))*N3 ...
        + wa03.*N0 ...
        + k4438*N4.^2 + 2*k4033*N0.*N4;
dN4dt = A(6,5)*N5+A(7,5)*N6+A(8,5)*N7+A(9,5)*N8+A(10,5)*N9+A(11,5)*N10+A(12,5)*N11 ...
        + Gamma(5)*N5 ...
        - (A(5,1) + A(5,2) + A(5,3) + A(5,4) + Gamma(4))*N4 ...
        - we40.*N4 - we41.*N4 - we42.*N4 ...
        + wa04.*N0 - wa49.*N4 ...
        - 2*(k4438+k4429+k440110+k440011)*N4.^2 - k4033*N0.*N4;
dN5dt = A(7,6)*N6+A(8,6)*N7+A(9,6)*N8+A(10,6)*N9+A(11,6)*N10+A(12,6)*N11 ...
        + Gamma(6)*N6 ...
        - (A(6,1) + A(6,2) + A(6,3) + A(6,4) + A(6,5) + Gamma(5))*N5 ...
        + wa05.*N0;
dN6dt = A(8,7)*N7+A(9,7)*N8+A(10,7)*N9+A(11,7)*N10+A(12,7)*N11 ...
        + Gamma(7)*N7 ...
        - (A(7,1) + A(7,2) + A(7,3) + A(7,4) + A(7,5) + A(7,6) + Gamma(6))*N6 ...
        + wa06.*N0;
dN7dt = A(9,8)*N8+A(10,8)*N9+A(11,8)*N10+A(12,8)*N11 ...
        + Gamma(8)*N8 ...
        - (A(8,1) + A(8,2) + A(8,3) + A(8,4) + A(8,5) + A(8,6) + A(8,7) + Gamma(7))*N7 ...
        + wa07.*N0;
dN8dt = A(10,9)*N9+A(11,9)*N10+A(12,9)*N11 ...
        + Gamma(9)*N9 ...
        - (A(9,1) + A(9,2) + A(9,3) + A(9,4) + A(9,5) + A(9,6) + A(9,7) + A(9,8) + Gamma(8))*N8 ...
        + wa08.*N0 ...
        + k4438*N4.^2;
dN9dt = A(11,10)*N10+A(12,10)*N11 ...
        + Gamma(10)*N10 ...
        - (A(10,1) + A(10,2) + A(10,3) + A(10,4) + A(10,5) + A(10,6) + A(10,7) + A(10,8) + A(10,9) + Gamma(9))*N9 ...
        + wa09.*N0 + wa49.*N4 ...
        + k4429*N4.^2;
dN10dt = A(12,11)*N11 ...
        + Gamma(11)*N11 ...
        - (A(11,1) + A(11,2) + A(11,3) + A(11,4) + A(11,5) + A(11,6) + A(11,7) + A(11,8) + A(11,9) + A(11,10) + Gamma(10))*N10 ...
        + wa0010.*N0 ...
        + k440110*N4.^2;
%dN11dt = - (A(12,1) + A(12,2) + A(12,3) + A(12,4) + A(12,5) + A(12,6) + A(12,7) + A(12,8) + A(12,9) + A(12,10) + A(12,11) + Gamma(11))*N11 ...
%        + wa0011.*N0 ...
%        + k440011*N4.^2;

dNdt = cat(1,dN0dt,dN1dt,dN2dt,dN3dt,dN4dt,dN5dt,dN6dt,dN7dt,dN8dt,dN9dt,dN10dt);

end

function Jacobian_Nd = J(N,A,Gamma,kijkl,R)
%J Jacobian matrix of the population coupled equations among Nd's energy
%levels
%
% Input arguments:
%   N: population, but it's a ratio here over N_total; (11,1,num_x,num_y,1,1,1,M)
%   A: spontaneous transition rate/probabilities; a (12,12) matrix
%   Gamma: multiphonon decay rates; a (11,1) vector
%   kijkl: nonlinear coupling coefficients for upconversion or cross relaxation effects; a (1,5) vector
%   R: photon rate; (1,1,num_x,num_y,1,1,1,M,1,num_cross_sections)
%
% Output argument:
%   Jacobian_Nd: Jacobian matrix; a (11,11,num_x,num_y,1,1,1,M) array

% Populations of each energy level
N0  = N( 1,:,:,:,:,:,:,:); % ground state, 4I9/2
%N1  = N( 2,:,:,:,:,:,:,:); % 4I11/2
%N2  = N( 3,:,:,:,:,:,:,:); % 4I13/2
%N3  = N( 4,:,:,:,:,:,:,:); % 4I15/2
N4  = N( 5,:,:,:,:,:,:,:); % 4F3/2
%N5  = N( 6,:,:,:,:,:,:,:); % 2H9/2 + 4F5/2
%N6  = N( 7,:,:,:,:,:,:,:); % 4S3/2 + 4F7/2
%N7  = N( 8,:,:,:,:,:,:,:); % 4F9/2
%N8  = N( 9,:,:,:,:,:,:,:); % 2G7/2 + 4G5/2
%N9  = N(10,:,:,:,:,:,:,:); % 4G9/2 + 4G7/2 + 2K13/2
%N10 = N(11,:,:,:,:,:,:,:); % 2G11/2 + 2(D,P)3/2 + 2G9/2 + 2K15/2
%N11 = N_total - sum(N,1); % 2D5/2 + 2P1/2

% Nonlinear coupling coefficients
k4438   = kijkl(1);
k4429  = kijkl(2);
k440110 = kijkl(3);
k440011 = kijkl(4);
k4033   = kijkl(5);

% Stimulated photon rates (absorption and emission)
wa0011 = R(:,:,:,:,:,:,:,:,:, 1);
wa0010 = R(:,:,:,:,:,:,:,:,:, 2);
wa09   = R(:,:,:,:,:,:,:,:,:, 3);
wa08   = R(:,:,:,:,:,:,:,:,:, 4);
wa07   = R(:,:,:,:,:,:,:,:,:, 5);
wa06   = R(:,:,:,:,:,:,:,:,:, 6);
wa05   = R(:,:,:,:,:,:,:,:,:, 7);
wa04   = R(:,:,:,:,:,:,:,:,:, 8);
wa03   = R(:,:,:,:,:,:,:,:,:, 9);
wa49   = R(:,:,:,:,:,:,:,:,:,10);
we40   = R(:,:,:,:,:,:,:,:,:,11);
we41   = R(:,:,:,:,:,:,:,:,:,12);
we42   = R(:,:,:,:,:,:,:,:,:,13);

% Jacobian matrix
dN0dt_N0 = -A(12,1) ...
           - (wa03 + wa04 + wa05 + wa06 + wa07 + wa08 + wa09 + wa0010 + wa0011) ...
           - k4033*N4;
dN0dt_N1 = A(2,1)-A(12,1) ...
           + Gamma(1);
dN0dt_N2 = A(3,1)-A(12,1);
dN0dt_N3 = A(4,1)-A(12,1);
dN0dt_N4 = A(5,1)-A(12,1) ...
           + we40 ...
           + 2*k440011*N4 - k4033*N0;
dN0dt_N5 = A(6,1)-A(12,1);
dN0dt_N6 = A(7,1)-A(12,1);
dN0dt_N7 = A(8,1)-A(12,1);
dN0dt_N8 = A(9,1)-A(12,1);
dN0dt_N9 = A(10,1)-A(12,1);
dN0dt_N10 = A(11,1)-A(12,1);
dN1dt_N0 = -A(12,2);
dN1dt_N1 = -A(12,2) ...
           - (A(2,1) + Gamma(1));
dN1dt_N2 = A(3,2)-A(12,2) ...
           + Gamma(2);
dN1dt_N3 = A(4,2)-A(12,2);
dN1dt_N4 = A(5,2)-A(12,2) ...
           + we41 ...
           + 2*k440110*N4;
dN1dt_N5 = A(6,2)-A(12,2);
dN1dt_N6 = A(7,2)-A(12,2);
dN1dt_N7 = A(8,2)-A(12,2);
dN1dt_N8 = A(9,2)-A(12,2);
dN1dt_N9 = A(10,2)-A(12,2);
dN1dt_N10 = A(11,2)-A(12,2);
dN2dt_N0 = -A(12,3);
dN2dt_N1 = -A(12,3);
dN2dt_N2 = -A(12,3) ...
           - (A(3,1) + A(3,2) + Gamma(2));
dN2dt_N3 = A(4,3)-A(12,3) ...
           + Gamma(3);
dN2dt_N4 = A(5,3)-A(12,3) ...
           + we42 ...
           + 2*k4429*N4;
dN2dt_N5 = A(6,3)-A(12,3);
dN2dt_N6 = A(7,3)-A(12,3);
dN2dt_N7 = A(8,3)-A(12,3);
dN2dt_N8 = A(9,3)-A(12,3);
dN2dt_N9 = A(10,3)-A(12,3);
dN2dt_N10 = A(11,3)-A(12,3);
dN3dt_N0 = -A(12,4) ...
           + wa03 ...
           + 2*k4033*N4;
dN3dt_N1 = -A(12,4);
dN3dt_N2 = -A(12,4);
dN3dt_N3 = -A(12,4) ...
           - (A(4,1) + A(4,2) + A(4,3) + Gamma(3));
dN3dt_N4 = A(5,4)-A(12,4) ...
           + Gamma(4) ...
           +2*k4033*N0 + 2*k4438*N4;
dN3dt_N5 = A(6,4)-A(12,4);
dN3dt_N6 = A(7,4)-A(12,4);
dN3dt_N7 = A(8,4)-A(12,4);
dN3dt_N8 = A(9,4)-A(12,4);
dN3dt_N9 = A(10,4)-A(12,4);
dN3dt_N10 = A(11,4)-A(12,4);
dN4dt_N0 = -A(12,5) ...
           + wa04 ...
           -k4033*N4;
dN4dt_N1 = -A(12,5);
dN4dt_N2 = -A(12,5);
dN4dt_N3 = -A(12,5);
dN4dt_N4 = -A(12,5) ...
           - (A(5,1) + A(5,2) + A(5,3) + A(5,4) + Gamma(4)) ...
           - wa49 ...
           - we40 - we41 - we42 ...
           - k4033*N0 - 2*(2*k4429 + 2*k4438 + 2*k440011 + 2*k440110)*N4;
dN4dt_N5 = A(6,5)-A(12,5) ...
           + Gamma(5);
dN4dt_N6 = A(7,5)-A(12,5);
dN4dt_N7 = A(8,5)-A(12,5);
dN4dt_N8 = A(9,5)-A(12,5);
dN4dt_N9 = A(10,5)-A(12,5);
dN4dt_N10 = A(11,5)-A(12,5);
dN5dt_N0 = -A(12,6) ...
           + wa05;
dN5dt_N1 = -A(12,6);
dN5dt_N2 = -A(12,6);
dN5dt_N3 = -A(12,6);
dN5dt_N4 = -A(12,6);
dN5dt_N5 = -A(12,6) ...
           - (A(6,1) + A(6,2) + A(6,3) + A(6,4) + A(6,5) + Gamma(5));
dN5dt_N6 = A(7,6)-A(12,6) ...
           + Gamma(6);
dN5dt_N7 = A(8,6)-A(12,6);
dN5dt_N8 = A(9,6)-A(12,6);
dN5dt_N9 = A(10,6)-A(12,6);
dN5dt_N10 = A(11,6)-A(12,6);
dN6dt_N0 = -A(12,7) ...
           + wa06;
dN6dt_N1 = -A(12,7);
dN6dt_N2 = -A(12,7);
dN6dt_N3 = -A(12,7);
dN6dt_N4 = -A(12,7);
dN6dt_N5 = -A(12,7);
dN6dt_N6 = -A(12,7) ...
           - (A(7,1) + A(7,2) + A(7,3) + A(7,4) + A(7,5) + A(7,6) + Gamma(6));
dN6dt_N7 = A(8,7)-A(12,7) ...
           + Gamma(7);
dN6dt_N8 = A(9,6)-A(12,6);
dN6dt_N9 = A(10,6)-A(12,6);
dN6dt_N10 = A(11,6)-A(12,6);
dN7dt_N0 = -A(12,8) ...
           + wa07;
dN7dt_N1 = -A(12,8);
dN7dt_N2 = -A(12,8);
dN7dt_N3 = -A(12,8);
dN7dt_N4 = -A(12,8);
dN7dt_N5 = -A(12,8);
dN7dt_N6 = -A(12,8);
dN7dt_N7 = -A(12,8) ...
           - (A(8,1) + A(8,2) + A(8,3) + A(8,4) + A(8,5) + A(8,6) + A(8,7) + Gamma(7));
dN7dt_N8 = A(9,8)-A(12,8) ...
           + Gamma(8);
dN7dt_N9 = A(10,8)-A(12,9);
dN7dt_N10 = A(11,8)-A(12,9);
dN8dt_N0 = -A(12,9) ...
           + wa08;
dN8dt_N1 = -A(12,9);
dN8dt_N2 = -A(12,9);
dN8dt_N3 = -A(12,9);
dN8dt_N4 = -A(12,9) ...
           + 2*k4438*N4;
dN8dt_N5 = -A(12,9);
dN8dt_N6 = -A(12,9);
dN8dt_N7 = -A(12,9);
dN8dt_N8 = -A(12,9) ...
           - (A(9,1) + A(9,2) + A(9,3) + A(9,4) + A(9,5) + A(9,6) + A(9,7) + A(9,8) + Gamma(8));
dN8dt_N9 = A(10,9)-A(12,9) ...
           +Gamma(9);
dN8dt_N10 = A(11,9)-A(12,9);
dN9dt_N0 = -A(12,10) ...
           + wa09;
dN9dt_N1 = -A(12,10);
dN9dt_N2 = -A(12,10);
dN9dt_N3 = -A(12,10);
dN9dt_N4 = -A(12,10) ...
           + wa49 ...
           + 2*k4429*N4;
dN9dt_N5 = -A(12,10);
dN9dt_N6 = -A(12,10);
dN9dt_N7 = -A(12,10);
dN9dt_N8 = -A(12,10);
dN9dt_N9 = -A(12,10) ...
           - (A(10,1) + A(10,2) + A(10,3) + A(10,4) + A(10,5) + A(10,6) + A(10,7) + A(10,8) + A(10,9) + Gamma(9));
dN9dt_N10 = A(11,10)-A(12,10) ...
            + Gamma(10);
dN10dt_N0 = -A(12,11) ...
            - Gamma(11) ...
            + wa0010;
dN10dt_N1 = -A(12,11) ...
            - Gamma(11);
dN10dt_N2 = -A(12,11) ...
            - Gamma(11);
dN10dt_N3 = -A(12,11) ...
            - Gamma(11);
dN10dt_N4 = -A(12,11) ...
            - Gamma(11) ...
           + 2*k440110*N4;
dN10dt_N5 = -A(12,11) ...
            - Gamma(11);
dN10dt_N6 = -A(12,11) ...
            - Gamma(11);
dN10dt_N7 = -A(12,11) ...
            - Gamma(11);
dN10dt_N8 = -A(12,11) ...
            - Gamma(11);
dN10dt_N9 = -A(12,11) ...
            - Gamma(11);
dN10dt_N10 = -A(12,11) ...
             - Gamma(11) ...
             - (A(11,1) + A(11,2) + A(11,3) + A(11,4) + A(11,5) + A(11,6) + A(11,7) + A(11,8) + A(11,9) + A(11,10) + Gamma(10));

% Some computations don't involve photon rates, so they have
% less-dimensional arrays, which makes array catenation into one Jacobian
% matrix fail. Extension of their dimensions is thus required.
num_x = size(R,3);
num_y = size(R,4);
M = size(R,8);
if size(N,8) ~= M % the number of parallelizations
    % Before extension, they're (1,1,num_x,num_y,1,1,1,1).
    dN3dt_N4 = repmat(dN3dt_N4,1,1,1,1,1,1,1,M);
    dN8dt_N4 = repmat(dN8dt_N4,1,1,1,1,1,1,1,M);
    dN10dt_N4 = repmat(dN10dt_N4,1,1,1,1,1,1,1,M);
end
% Before extension, they're scalars with the dimension (1,1,1,1,1,1,1,1).
dN0dt_N1 = repmat(dN0dt_N1,1,1,num_x,num_y,1,1,1,M);
dN0dt_N2 = repmat(dN0dt_N2,1,1,num_x,num_y,1,1,1,M);
dN0dt_N3 = repmat(dN0dt_N3,1,1,num_x,num_y,1,1,1,M);
dN0dt_N5 = repmat(dN0dt_N5,1,1,num_x,num_y,1,1,1,M);
dN0dt_N6 = repmat(dN0dt_N6,1,1,num_x,num_y,1,1,1,M);
dN0dt_N7 = repmat(dN0dt_N7,1,1,num_x,num_y,1,1,1,M);
dN0dt_N8 = repmat(dN0dt_N8,1,1,num_x,num_y,1,1,1,M);
dN0dt_N9 = repmat(dN0dt_N9,1,1,num_x,num_y,1,1,1,M);
dN0dt_N10 = repmat(dN0dt_N10,1,1,num_x,num_y,1,1,1,M);
dN1dt_N0 = repmat(dN1dt_N0,1,1,num_x,num_y,1,1,1,M);
dN1dt_N1 = repmat(dN1dt_N1,1,1,num_x,num_y,1,1,1,M);
dN1dt_N2 = repmat(dN1dt_N2,1,1,num_x,num_y,1,1,1,M);
dN1dt_N3 = repmat(dN1dt_N3,1,1,num_x,num_y,1,1,1,M);
dN1dt_N5 = repmat(dN1dt_N5,1,1,num_x,num_y,1,1,1,M);
dN1dt_N6 = repmat(dN1dt_N6,1,1,num_x,num_y,1,1,1,M);
dN1dt_N7 = repmat(dN1dt_N7,1,1,num_x,num_y,1,1,1,M);
dN1dt_N8 = repmat(dN1dt_N8,1,1,num_x,num_y,1,1,1,M);
dN1dt_N9 = repmat(dN1dt_N9,1,1,num_x,num_y,1,1,1,M);
dN1dt_N10 = repmat(dN1dt_N10,1,1,num_x,num_y,1,1,1,M);
dN2dt_N0 = repmat(dN2dt_N0,1,1,num_x,num_y,1,1,1,M);
dN2dt_N1 = repmat(dN2dt_N1,1,1,num_x,num_y,1,1,1,M);
dN2dt_N2 = repmat(dN2dt_N2,1,1,num_x,num_y,1,1,1,M);
dN2dt_N3 = repmat(dN2dt_N3,1,1,num_x,num_y,1,1,1,M);
dN2dt_N5 = repmat(dN2dt_N5,1,1,num_x,num_y,1,1,1,M);
dN2dt_N6 = repmat(dN2dt_N6,1,1,num_x,num_y,1,1,1,M);
dN2dt_N7 = repmat(dN2dt_N7,1,1,num_x,num_y,1,1,1,M);
dN2dt_N8 = repmat(dN2dt_N8,1,1,num_x,num_y,1,1,1,M);
dN2dt_N9 = repmat(dN2dt_N9,1,1,num_x,num_y,1,1,1,M);
dN2dt_N10 = repmat(dN2dt_N10,1,1,num_x,num_y,1,1,1,M);
dN3dt_N1 = repmat(dN3dt_N1,1,1,num_x,num_y,1,1,1,M);
dN3dt_N2 = repmat(dN3dt_N2,1,1,num_x,num_y,1,1,1,M);
dN3dt_N3 = repmat(dN3dt_N3,1,1,num_x,num_y,1,1,1,M);
dN3dt_N5 = repmat(dN3dt_N5,1,1,num_x,num_y,1,1,1,M);
dN3dt_N6 = repmat(dN3dt_N6,1,1,num_x,num_y,1,1,1,M);
dN3dt_N7 = repmat(dN3dt_N7,1,1,num_x,num_y,1,1,1,M);
dN3dt_N8 = repmat(dN3dt_N8,1,1,num_x,num_y,1,1,1,M);
dN3dt_N9 = repmat(dN3dt_N9,1,1,num_x,num_y,1,1,1,M);
dN3dt_N10 = repmat(dN3dt_N10,1,1,num_x,num_y,1,1,1,M);
dN4dt_N1 = repmat(dN4dt_N1,1,1,num_x,num_y,1,1,1,M);
dN4dt_N2 = repmat(dN4dt_N2,1,1,num_x,num_y,1,1,1,M);
dN4dt_N3 = repmat(dN4dt_N3,1,1,num_x,num_y,1,1,1,M);
dN4dt_N5 = repmat(dN4dt_N5,1,1,num_x,num_y,1,1,1,M);
dN4dt_N6 = repmat(dN4dt_N6,1,1,num_x,num_y,1,1,1,M);
dN4dt_N7 = repmat(dN4dt_N7,1,1,num_x,num_y,1,1,1,M);
dN4dt_N8 = repmat(dN4dt_N8,1,1,num_x,num_y,1,1,1,M);
dN4dt_N9 = repmat(dN4dt_N9,1,1,num_x,num_y,1,1,1,M);
dN4dt_N10 = repmat(dN4dt_N10,1,1,num_x,num_y,1,1,1,M);
dN5dt_N1 = repmat(dN5dt_N1,1,1,num_x,num_y,1,1,1,M);
dN5dt_N2 = repmat(dN5dt_N2,1,1,num_x,num_y,1,1,1,M);
dN5dt_N3 = repmat(dN5dt_N3,1,1,num_x,num_y,1,1,1,M);
dN5dt_N4 = repmat(dN5dt_N4,1,1,num_x,num_y,1,1,1,M);
dN5dt_N5 = repmat(dN5dt_N5,1,1,num_x,num_y,1,1,1,M);
dN5dt_N6 = repmat(dN5dt_N6,1,1,num_x,num_y,1,1,1,M);
dN5dt_N7 = repmat(dN5dt_N7,1,1,num_x,num_y,1,1,1,M);
dN5dt_N8 = repmat(dN5dt_N8,1,1,num_x,num_y,1,1,1,M);
dN5dt_N9 = repmat(dN5dt_N9,1,1,num_x,num_y,1,1,1,M);
dN5dt_N10 = repmat(dN5dt_N10,1,1,num_x,num_y,1,1,1,M);
dN6dt_N1 = repmat(dN6dt_N1,1,1,num_x,num_y,1,1,1,M);
dN6dt_N2 = repmat(dN6dt_N2,1,1,num_x,num_y,1,1,1,M);
dN6dt_N3 = repmat(dN6dt_N3,1,1,num_x,num_y,1,1,1,M);
dN6dt_N4 = repmat(dN6dt_N4,1,1,num_x,num_y,1,1,1,M);
dN6dt_N5 = repmat(dN6dt_N5,1,1,num_x,num_y,1,1,1,M);
dN6dt_N6 = repmat(dN6dt_N6,1,1,num_x,num_y,1,1,1,M);
dN6dt_N7 = repmat(dN6dt_N7,1,1,num_x,num_y,1,1,1,M);
dN6dt_N8 = repmat(dN6dt_N8,1,1,num_x,num_y,1,1,1,M);
dN6dt_N9 = repmat(dN6dt_N9,1,1,num_x,num_y,1,1,1,M);
dN6dt_N10 = repmat(dN6dt_N10,1,1,num_x,num_y,1,1,1,M);
dN7dt_N1 = repmat(dN7dt_N1,1,1,num_x,num_y,1,1,1,M);
dN7dt_N2 = repmat(dN7dt_N2,1,1,num_x,num_y,1,1,1,M);
dN7dt_N3 = repmat(dN7dt_N3,1,1,num_x,num_y,1,1,1,M);
dN7dt_N4 = repmat(dN7dt_N4,1,1,num_x,num_y,1,1,1,M);
dN7dt_N5 = repmat(dN7dt_N5,1,1,num_x,num_y,1,1,1,M);
dN7dt_N6 = repmat(dN7dt_N6,1,1,num_x,num_y,1,1,1,M);
dN7dt_N7 = repmat(dN7dt_N7,1,1,num_x,num_y,1,1,1,M);
dN7dt_N8 = repmat(dN7dt_N8,1,1,num_x,num_y,1,1,1,M);
dN7dt_N9 = repmat(dN7dt_N9,1,1,num_x,num_y,1,1,1,M);
dN7dt_N10 = repmat(dN7dt_N10,1,1,num_x,num_y,1,1,1,M);
dN8dt_N1 = repmat(dN8dt_N1,1,1,num_x,num_y,1,1,1,M);
dN8dt_N2 = repmat(dN8dt_N2,1,1,num_x,num_y,1,1,1,M);
dN8dt_N3 = repmat(dN8dt_N3,1,1,num_x,num_y,1,1,1,M);
dN8dt_N5 = repmat(dN8dt_N5,1,1,num_x,num_y,1,1,1,M);
dN8dt_N6 = repmat(dN8dt_N6,1,1,num_x,num_y,1,1,1,M);
dN8dt_N7 = repmat(dN8dt_N7,1,1,num_x,num_y,1,1,1,M);
dN8dt_N8 = repmat(dN8dt_N8,1,1,num_x,num_y,1,1,1,M);
dN8dt_N9 = repmat(dN8dt_N9,1,1,num_x,num_y,1,1,1,M);
dN8dt_N10 = repmat(dN8dt_N10,1,1,num_x,num_y,1,1,1,M);
dN9dt_N1 = repmat(dN9dt_N1,1,1,num_x,num_y,1,1,1,M);
dN9dt_N2 = repmat(dN9dt_N2,1,1,num_x,num_y,1,1,1,M);
dN9dt_N3 = repmat(dN9dt_N3,1,1,num_x,num_y,1,1,1,M);
dN9dt_N5 = repmat(dN9dt_N5,1,1,num_x,num_y,1,1,1,M);
dN9dt_N6 = repmat(dN9dt_N6,1,1,num_x,num_y,1,1,1,M);
dN9dt_N7 = repmat(dN9dt_N7,1,1,num_x,num_y,1,1,1,M);
dN9dt_N8 = repmat(dN9dt_N8,1,1,num_x,num_y,1,1,1,M);
dN9dt_N9 = repmat(dN9dt_N9,1,1,num_x,num_y,1,1,1,M);
dN9dt_N10 = repmat(dN9dt_N10,1,1,num_x,num_y,1,1,1,M);
dN10dt_N1 = repmat(dN10dt_N1,1,1,num_x,num_y,1,1,1,M);
dN10dt_N2 = repmat(dN10dt_N2,1,1,num_x,num_y,1,1,1,M);
dN10dt_N3 = repmat(dN10dt_N3,1,1,num_x,num_y,1,1,1,M);
dN10dt_N5 = repmat(dN10dt_N5,1,1,num_x,num_y,1,1,1,M);
dN10dt_N6 = repmat(dN10dt_N6,1,1,num_x,num_y,1,1,1,M);
dN10dt_N7 = repmat(dN10dt_N7,1,1,num_x,num_y,1,1,1,M);
dN10dt_N8 = repmat(dN10dt_N8,1,1,num_x,num_y,1,1,1,M);
dN10dt_N9 = repmat(dN10dt_N9,1,1,num_x,num_y,1,1,1,M);
dN10dt_N10 = repmat(dN10dt_N10,1,1,num_x,num_y,1,1,1,M);

Jacobian_Nd = cat(1,cat(2,dN0dt_N0,dN0dt_N1,dN0dt_N2,dN0dt_N3,dN0dt_N4,dN0dt_N5,dN0dt_N6,dN0dt_N7,dN0dt_N8,dN0dt_N9,dN0dt_N10),...
                    cat(2,dN1dt_N0,dN1dt_N1,dN1dt_N2,dN1dt_N3,dN1dt_N4,dN1dt_N5,dN1dt_N6,dN1dt_N7,dN1dt_N8,dN1dt_N9,dN1dt_N10),...
                    cat(2,dN2dt_N0,dN2dt_N1,dN2dt_N2,dN2dt_N3,dN2dt_N4,dN2dt_N5,dN2dt_N6,dN2dt_N7,dN2dt_N8,dN2dt_N9,dN2dt_N10),...
                    cat(2,dN3dt_N0,dN3dt_N1,dN3dt_N2,dN3dt_N3,dN3dt_N4,dN3dt_N5,dN3dt_N6,dN3dt_N7,dN3dt_N8,dN3dt_N9,dN3dt_N10),...
                    cat(2,dN4dt_N0,dN4dt_N1,dN4dt_N2,dN4dt_N3,dN4dt_N4,dN4dt_N5,dN4dt_N6,dN4dt_N7,dN4dt_N8,dN4dt_N9,dN4dt_N10),...
                    cat(2,dN5dt_N0,dN5dt_N1,dN5dt_N2,dN5dt_N3,dN5dt_N4,dN5dt_N5,dN5dt_N6,dN5dt_N7,dN5dt_N8,dN5dt_N9,dN5dt_N10),...
                    cat(2,dN6dt_N0,dN6dt_N1,dN6dt_N2,dN6dt_N3,dN6dt_N4,dN6dt_N5,dN6dt_N6,dN6dt_N7,dN6dt_N8,dN6dt_N9,dN6dt_N10),...
                    cat(2,dN7dt_N0,dN7dt_N1,dN7dt_N2,dN7dt_N3,dN7dt_N4,dN7dt_N5,dN7dt_N6,dN7dt_N7,dN7dt_N8,dN7dt_N9,dN7dt_N10),...
                    cat(2,dN8dt_N0,dN8dt_N1,dN8dt_N2,dN8dt_N3,dN8dt_N4,dN8dt_N5,dN8dt_N6,dN8dt_N7,dN8dt_N8,dN8dt_N9,dN8dt_N10),...
                    cat(2,dN9dt_N0,dN9dt_N1,dN9dt_N2,dN9dt_N3,dN9dt_N4,dN9dt_N5,dN9dt_N6,dN9dt_N7,dN9dt_N8,dN9dt_N9,dN9dt_N10),...
                    cat(2,dN10dt_N0,dN10dt_N1,dN10dt_N2,dN10dt_N3,dN10dt_N4,dN10dt_N5,dN10dt_N6,dN10dt_N7,dN10dt_N8,dN10dt_N9,dN10dt_N10));

end