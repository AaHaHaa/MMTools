function func = FJ_Er()
%FJ_ER container for the coupled equations of Er population evolutions and
%their Jacobian matrix
%
% 1. Morkel and Laming, "Theoretical modeling of erbium-doped fiber
%    amplifiers with excited-state absorption," Opt. Lett. 14(19), 1062-
%    1064 (1989)
% 2. Maciuc et al., "Rate equations for an erbium laser system: a numerical
%    approach," ROMOPTO 2000: Sixth Conference on Optics 4430, 136-146
%    (2001)
% 3. Henderson-Sapir et al., "New energy-transfer upconversion process in 
%    Er3+:ZBLAN mid-infrared fiber lasers," Opt. Express 24, 6869-6883 
%    (2016)
% 4. Jian-Feng Li and Stuart D. jackson, "Theoretical study and
%    optimization of a high power mid-infrared erbium-doped ZBLAN fibre 
%    laser," Chin. Phys. B 20(3), 034205 (2011)
% 5. Cai et al., "Numerical Analysis of a Dual-Wavelength-Clad-Pumped 3.5
%    μm Erbium-Doped Fluoride Fiber Laser," Appl. Sciences 12(15) (2022)
% 6. Pele et al., "Wavelength conversion in Er^3+ doped chalcogenide fibers
%    for optical gas sensors," Opt. Express 23(4), 4163-4172 (2015)

func.F = @F;
func.J = @J;

end

%% Contained functions
function dNdt = F(N,N_total,A,Gamma,kijkl,R)
%F Population coupled equations among Er's energy levels
%
% Input arguments:
%   N: population; (8,1,num_x,num_y,1,1,1,M)
%   N_total: total population; a (1,1,num_x,num_y) array
%   A: spontaneous transition rate/probabilities; a (9,9) matrix
%   Gamma: multiphonon decay rates; a (8,1) vector
%   kijkl: nonlinear coupling coefficients for upconversion or cross relaxation effects; a (1,4) vector
%   R: photon rate; (1,1,num_x,num_y,1,1,1,M,1,num_cross_sections)
%
% Output argument:
%   dNdt: variation of populations; a (8,1,num_x,num_y,1,1,1,M) array

% Populations of each energy level
N0 = N_total - sum(N,1); % ground state, 4I15/2
N1 = N(1,:,:,:,:,:,:,:); % 4I13/2
N2 = N(2,:,:,:,:,:,:,:); % 4I11/2
N3 = N(3,:,:,:,:,:,:,:); % 4I9/2
N4 = N(4,:,:,:,:,:,:,:); % 4F9/2
N5 = N(5,:,:,:,:,:,:,:); % 4S3/2 + 2H11/2
N6 = N(6,:,:,:,:,:,:,:); % 4F7/2
N7 = N(7,:,:,:,:,:,:,:); % 4F5/2
N8 = N(8,:,:,:,:,:,:,:); % 2H9/2

% Nonlinear coupling coefficients
k1103 = kijkl(1);
k2206 = kijkl(2);
k5031 = kijkl(3);
k4251 = kijkl(4);

% Stimulated photon rates (absorption and emission)
wa08 = R(:,:,:,:,:,:,:,:,:,1);
wa07 = R(:,:,:,:,:,:,:,:,:,2);
wa06 = R(:,:,:,:,:,:,:,:,:,3);
wa05 = R(:,:,:,:,:,:,:,:,:,4);
wa04 = R(:,:,:,:,:,:,:,:,:,5);
wa03 = R(:,:,:,:,:,:,:,:,:,6);
wa02 = R(:,:,:,:,:,:,:,:,:,7);
wa01 = R(:,:,:,:,:,:,:,:,:,8);
wa15 = R(:,:,:,:,:,:,:,:,:,9);
wa26 = R(:,:,:,:,:,:,:,:,:,10);
wa12 = R(:,:,:,:,:,:,:,:,:,11);
wa23 = R(:,:,:,:,:,:,:,:,:,12);
we10 = R(:,:,:,:,:,:,:,:,:,13);
we21 = R(:,:,:,:,:,:,:,:,:,14);
we32 = R(:,:,:,:,:,:,:,:,:,15);

% Variation of populations
% dN0 + dN1 + dN2 + dN3 + dN4 + dN5 + dN6 + dN7 + dN8 = dN_total = 0
%dN0dt = A(2,1)*N1+A(3,1)*N2+A(4,1)*N3+A(5,1)*N4+A(6,1)*N5+A(7,1)*N6+A(8,1)*N7+A(9,1)*N8 ...
%        + Gamma(1)*N1 ...
%        + we10.*N1 ...
%        - (wa01 + wa02 + wa03 + wa04 + wa05 + wa06 + wa07 + wa08).*N0 ...
%        + k1103*N1.^2 + k2206*N2.^2 - k5031*N0.*N5;
dN1dt = A(3,2)*N2+A(4,2)*N3+A(5,2)*N4+A(6,2)*N5+A(7,2)*N6+A(8,2)*N7+A(9,2)*N8 ...
        + Gamma(2)*N2 ...
        - (A(2,1) + Gamma(1))*N1 ...
        - we10.*N1 + we21.*N2 ...
        + wa01.*N0 - (wa12 + wa15).*N1 ...
        - 2*k1103*N1.^2 + k4251*N2.*N4 + k5031*N0.*N5;
dN2dt = A(4,3)*N3+A(5,3)*N4+A(6,3)*N5+A(7,3)*N6+A(8,3)*N7+A(9,3)*N8 ...
        + Gamma(3)*N3 ...
        - (A(3,1) + A(3,2) + Gamma(2))*N2 ...
        + we32.*N3 - we21.*N2 ...
        + wa02.*N0 + wa12.*N1 - wa23.*N2 - wa26.*N2 ...
        - 2*k2206*N2.^2 - k4251*N2.*N4;
dN3dt = A(5,4)*N4+A(6,4)*N5+A(7,4)*N6+A(8,4)*N7+A(9,4)*N8 ...
        + Gamma(4)*N4 ...
        - (A(4,1) + A(4,2) + A(4,3) + Gamma(3))*N3 ...
        - we32.*N3 ...
        + wa03.*N0 + wa23.*N2 ...
        + k1103*N1.^2 + k5031*N0.*N5;
dN4dt = A(6,5)*N5+A(7,5)*N6+A(8,5)*N7+A(9,5)*N8 ...
        + Gamma(5)*N5 ...
        - (A(5,1) + A(5,2) + A(5,3) + A(5,4) + Gamma(4))*N4 ...
        + wa04.*N0 ...
        - k4251*N2.*N4;
dN5dt = A(7,6)*N6+A(8,6)*N7+A(9,6)*N8 ...
        + Gamma(6)*N6 ...
        - (A(6,1) + A(6,2) + A(6,3) + A(6,4) + A(6,5) + Gamma(5))*N5 ...
        + wa05.*N0 + wa15.*N1 ...
        + k4251*N2.*N4- k5031*N0.*N5;
dN6dt = A(8,7)*N7+A(9,7)*N8 ...
        + Gamma(7)*N7 ...
        - (A(7,1) + A(7,2) + A(7,3) + A(7,4) + A(7,5) + A(7,6) + Gamma(6))*N6 ...
        + wa06.*N0 + wa26.*N2 ...
        + k2206*N2.^2;
dN7dt = A(9,8)*N8 ...
        + Gamma(8)*N8 ...
        - (A(8,1) + A(8,2) + A(8,3) + A(8,4) + A(8,5) + A(8,6) + A(8,7) + Gamma(7))*N7 ...
        + wa07.*N0;
dN8dt = - (A(9,1) + A(9,2) + A(9,3) + A(9,4) + A(9,5) + A(9,6) + A(9,7) + A(9,8) + Gamma(8))*N8 ...
        + wa08.*N0;

dNdt = cat(1,dN1dt,dN2dt,dN3dt,dN4dt,dN5dt,dN6dt,dN7dt,dN8dt);

end

function Jacobian_Er = J(N,N_total,A,Gamma,kijkl,R)
%J Jacobian matrix of the population coupled equations among Er's energy
%levels
%
% Input arguments:
%   N: population, but it's a ratio here over N_total; (8,1,num_x,num_y,1,1,1,M)
%   A: spontaneous transition rate/probabilities; a (9,9) matrix
%   Gamma: multiphonon decay rates; a (8,1) vector
%   kijkl: nonlinear coupling coefficients for upconversion or cross relaxation effects; a (1,4) vector
%   R: photon rate; (1,1,num_x,num_y,1,1,1,M,1,num_cross_sections)
%
% Output argument:
%   Jacobian_Er: Jacobian matrix; a (8,8,num_x,num_y,1,1,1,M) array

% Populations of each energy level
N0 = N_total - sum(N,1); % ground state, 4I15/2
N1 = N(1,:,:,:,:,:,:,:); % 4I13/2
N2 = N(2,:,:,:,:,:,:,:); % 4I11/2
%N3 = N(3,:,:,:,:,:,:,:); % 4I9/2
N4 = N(4,:,:,:,:,:,:,:); % 4F9/2
N5 = N(5,:,:,:,:,:,:,:); % 4S3/2 + 2H11/2
%N6 = N(6,:,:,:,:,:,:,:); % 4F7/2
%N7 = N(7,:,:,:,:,:,:,:); % 4F5/2
%N8 = N(8,:,:,:,:,:,:,:); % 2H9/2

% Nonlinear coupling coefficients
k1103 = kijkl(1);
k2206 = kijkl(2);
k5031 = kijkl(3);
k4251 = kijkl(4);

% Stimulated photon rates (absorption and emission)
wa08 = R(:,:,:,:,:,:,:,:,:,1);
wa07 = R(:,:,:,:,:,:,:,:,:,2);
wa06 = R(:,:,:,:,:,:,:,:,:,3);
wa05 = R(:,:,:,:,:,:,:,:,:,4);
wa04 = R(:,:,:,:,:,:,:,:,:,5);
wa03 = R(:,:,:,:,:,:,:,:,:,6);
wa02 = R(:,:,:,:,:,:,:,:,:,7);
wa01 = R(:,:,:,:,:,:,:,:,:,8);
wa15 = R(:,:,:,:,:,:,:,:,:,9);
wa26 = R(:,:,:,:,:,:,:,:,:,10);
wa12 = R(:,:,:,:,:,:,:,:,:,11);
wa23 = R(:,:,:,:,:,:,:,:,:,12);
we10 = R(:,:,:,:,:,:,:,:,:,13);
we21 = R(:,:,:,:,:,:,:,:,:,14);
we32 = R(:,:,:,:,:,:,:,:,:,15);

% Jacobian matrix
dN1dt_N1 = - A(2,1) - Gamma(1) - wa01 - wa12 - wa15 - we10 - 4*N1*k1103 - N5*k5031;
dN1dt_N2 = A(3,2) + Gamma(2) - wa01 + we21 + N4*k4251 - N5*k5031;
dN1dt_N3 = A(4,2) - wa01 - N5*k5031;
dN1dt_N4 = A(5,2) - wa01 + N2*k4251 - N5*k5031;
dN1dt_N5 = A(6,2) - wa01 + (N0 - N5)*k5031;
dN1dt_N6 = A(7,2) - wa01 - N5*k5031;
dN1dt_N7 = A(8,2) - wa01 - N5*k5031;
dN1dt_N8 = A(9,2) - wa01 - N5*k5031;
dN2dt_N1 = wa12 - wa02;
dN2dt_N2 = - A(3,1) - A(3,2) - Gamma(2) - wa02 - wa23 - wa26 - we21 - 4*N2*k2206 - N4*k4251;
dN2dt_N3 = A(4,3) + Gamma(3) - wa02 + we32;
dN2dt_N4 = A(5,3) - wa02 - N2*k4251;
dN2dt_N5 = A(6,3) - wa02;
dN2dt_N6 = A(7,3) - wa02;
dN2dt_N7 = A(8,3) - wa02;
dN2dt_N8 = A(9,3) - wa02;
dN3dt_N1 = 2*N1*k1103 - wa03 - N5*k5031;
dN3dt_N2 = wa23 - wa03 - N5*k5031;
dN3dt_N3 = - A(4,1) - A(4,2) - A(4,3) - Gamma(3) - wa03 - we32 - N5*k5031;
dN3dt_N4 = A(5,4) + Gamma(4) - wa03 - N5*k5031;
dN3dt_N5 = A(6,4) - wa03 + (N0 - N5)*k5031;
dN3dt_N6 = A(7,4) - wa03 - N5*k5031;
dN3dt_N7 = A(8,4) - wa03 - N5*k5031;
dN3dt_N8 = A(9,4) - wa03 - N5*k5031;
dN4dt_N1 = - wa04;
dN4dt_N2 = - wa04 - N4*k4251;
dN4dt_N3 = - wa04;
dN4dt_N4 = - A(5,1) - A(5,2) - A(5,3) - A(5,4) - Gamma(4) - wa04 - N2*k4251;
dN4dt_N5 = A(6,5) + Gamma(5) - wa04;
dN4dt_N6 = A(7,5) - wa04;
dN4dt_N7 = A(8,5) - wa04;
dN4dt_N8 = A(9,5) - wa04;
dN5dt_N1 = wa15 - wa05 + N5*k5031;
dN5dt_N2 = N4*k4251 - wa05 + N5*k5031;
dN5dt_N3 = N5*k5031 - wa05;
dN5dt_N4 = N2*k4251 - wa05 + N5*k5031;
dN5dt_N5 = (N5 -N0)*k5031 - A(6,2) - A(6,3) - A(6,4) - A(6,5) - Gamma(5) - wa05 - A(6,1);
dN5dt_N6 = A(7,6) + Gamma(6) - wa05 + N5*k5031;
dN5dt_N7 = A(8,6) - wa05 + N5*k5031;
dN5dt_N8 = A(9,6) - wa05 + N5*k5031;
dN6dt_N1 = - wa06;
dN6dt_N2 = wa26 - wa06 + 2*N2*k2206;
dN6dt_N3 = - wa06;
dN6dt_N4 = - wa06;
dN6dt_N5 = - wa06;
dN6dt_N6 = - A(7,1) - A(7,2) - A(7,3) - A(7,4) - A(7,5) - A(7,6) - Gamma(6) - wa06;
dN6dt_N7 = A(8,7) + Gamma(7) - wa06;
dN6dt_N8 = A(9,7) - wa06;
dN7dt_N1 = - wa07;
dN7dt_N2 = - wa07;
dN7dt_N3 = - wa07;
dN7dt_N4 = - wa07;
dN7dt_N5 = - wa07;
dN7dt_N6 = - wa07;
dN7dt_N7 = - A(8,1) - A(8,2) - A(8,3) - A(8,4) - A(8,5) - A(8,6) - A(8,7) - Gamma(7) - wa07;
dN7dt_N8 = A(9,8) + Gamma(8) - wa07;
dN8dt_N1 = - wa08;
dN8dt_N2 = - wa08;
dN8dt_N3 = - wa08;
dN8dt_N4 = - wa08;
dN8dt_N5 = - wa08;
dN8dt_N6 = - wa08;
dN8dt_N7 = - wa08;
dN8dt_N8 = - A(9,1) - A(9,2) - A(9,3) - A(9,4) - A(9,5) - A(9,6) - A(9,7) - A(9,8) - Gamma(8) - wa08;

Jacobian_Er = cat(1,cat(2,dN1dt_N1,dN1dt_N2,dN1dt_N3,dN1dt_N4,dN1dt_N5,dN1dt_N6,dN1dt_N7,dN1dt_N8),...
                    cat(2,dN2dt_N1,dN2dt_N2,dN2dt_N3,dN2dt_N4,dN2dt_N5,dN2dt_N6,dN2dt_N7,dN2dt_N8),...
                    cat(2,dN3dt_N1,dN3dt_N2,dN3dt_N3,dN3dt_N4,dN3dt_N5,dN3dt_N6,dN3dt_N7,dN3dt_N8),...
                    cat(2,dN4dt_N1,dN4dt_N2,dN4dt_N3,dN4dt_N4,dN4dt_N5,dN4dt_N6,dN4dt_N7,dN4dt_N8),...
                    cat(2,dN5dt_N1,dN5dt_N2,dN5dt_N3,dN5dt_N4,dN5dt_N5,dN5dt_N6,dN5dt_N7,dN5dt_N8),...
                    cat(2,dN6dt_N1,dN6dt_N2,dN6dt_N3,dN6dt_N4,dN6dt_N5,dN6dt_N6,dN6dt_N7,dN6dt_N8),...
                    cat(2,dN7dt_N1,dN7dt_N2,dN7dt_N3,dN7dt_N4,dN7dt_N5,dN7dt_N6,dN7dt_N7,dN7dt_N8),...
                    cat(2,dN8dt_N1,dN8dt_N2,dN8dt_N3,dN8dt_N4,dN8dt_N5,dN8dt_N6,dN8dt_N7,dN8dt_N8));

end