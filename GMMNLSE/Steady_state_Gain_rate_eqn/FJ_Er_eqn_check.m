% This code acts as a checker for symbolically calculating the rate
% equations of Er populations.

%% Define variables
syms N0 N1 N2 N3 N4 N5 N6 N7 N8; % populations
syms N_total; % total population = N0+N1+N2+N3+N4+N5+N6+N7+N8
syms A10 A20 A30 A40 A50 A60 A70 A80; % spontaneous emission terms
syms A11 A21 A31 A41 A51 A61 A71 A81;
syms A12 A22 A32 A42 A52 A62 A72 A82;
syms A13 A23 A33 A43 A53 A63 A73 A83;
syms A14 A24 A34 A44 A54 A64 A74 A84;
syms A15 A25 A35 A45 A55 A65 A75 A85;
syms A16 A26 A36 A46 A56 A66 A76 A86;
syms A17 A27 A37 A47 A57 A67 A77 A87;
syms A18 A28 A38 A48 A58 A68 A78 A88;
syms r1 r2 r3 r4 r5 r6 r7 r8; % multiphonon decay rates; here, "r" stands for "Gamma" symbolically which is more commonly used mathematically for decay
syms we10; % stimulated emission terms
syms we21 we32
syms wa01 wa02 wa03 wa04 wa05 wa06 wa07 wa08; % stimulated absorption terms
syms wa15 wa26 wa12 wa23;
syms k1103 k2206 k5031 k4251; % cross-interacting terms, such as cross relaxation and upconversion

N0 = N_total - N1-N2-N3-N4-N5-N6-N7-N8; % make N0 a dependent variable since N_total is given and fixed

%% Coupled equations, dNi/dt (i=0,1,...,8), for each population
dN0dt = A10*N1+A20*N2+A30*N3+A40*N4+A50*N5+A60*N6+A70*N7+A80*N8 ...
        + r1*N1 ...
        + we10*N1 ...
        - (wa01 + wa02 + wa03 + wa04 + wa05 + wa06 + wa07 + wa08)*N0 ...
        + k1103*N1^2 + k2206*N2^2 - k5031*N0*N5;
dN1dt = A21*N2+A31*N3+A41*N4+A51*N5+A61*N6+A71*N7+A81*N8 ...
        + r2*N2 ...
        - (A10 + r1)*N1 ...
        - we10*N1 + we21*N2 ...
        + wa01*N0 - (wa12 + wa15)*N1 ...
        - 2*k1103*N1^2 + k4251*N2*N4 + k5031*N0*N5;
dN2dt = A32*N3+A42*N4+A52*N5+A62*N6+A72*N7+A82*N8 ...
        + r3*N3 ...
        - (A20 + A21 + r2)*N2 ...
        + we32*N3 - we21*N2 ...
        + wa02*N0 + wa12*N1 - wa23*N2 - wa26*N2 ...
        - 2*k2206*N2^2 - k4251*N2*N4;
dN3dt = A43*N4+A53*N5+A63*N6+A73*N7+A83*N8 ...
        + r4*N4 ...
        - (A30 + A31 + A32 + r3)*N3 ...
        - we32*N3 ...
        + wa03*N0 + wa23*N2 ...
        + k1103*N1^2 + k5031*N0*N5;
dN4dt = A54*N5+A64*N6+A74*N7+A84*N8 ...
        + r5*N5 ...
        - (A40 + A41 + A42 + A43 + r4)*N4 ...
        + wa04*N0 ...
        - k4251*N2*N4;
dN5dt = A65*N6+A75*N7+A85*N8 ...
        + r6*N6 ...
        - (A50 + A51 + A52 + A53 + A54 + r5)*N5 ...
        + wa05*N0 + wa15*N1 ...
        + k4251*N2*N4- k5031*N0*N5;
dN6dt = A76*N7+A86*N8 ...
        + r7*N7 ...
        - (A60 + A61 + A62 + A63 + A64 + A65 + r6)*N6 ...
        + wa06*N0 + wa26*N2 ...
        + k2206*N2^2;
dN7dt = A87*N8 ...
        + r8*N8 ...
        - (A70 + A71 + A72 + A73 + A74 + A75 + A76 + r7)*N7 ...
        + wa07*N0;
dN8dt = - (A80 + A81 + A82 + A83 + A84 + A85 + A86 + A87 + r8)*N8 ...
        + wa08*N0;
  
%% Check
disp(simplify(dN0dt+dN1dt+dN2dt+dN3dt+dN4dt+dN5dt+dN6dt+dN7dt+dN8dt)); % This should be zero

% Vary the variables below for calculating the derivatives to find the
% Jacobian matrix
for i = 1:8
    for j = 1:8
        fprintf('dN%udt_N%u = %s\n',i,j,simplify(eval(sprintf('diff(dN%udt,N%u)',i,j)))); % vary it to diff(dN0,dN2), etc. for example
    end
end