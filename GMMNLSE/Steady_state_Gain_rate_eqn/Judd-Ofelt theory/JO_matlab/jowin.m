%  This program performs a Judd-Ofelt analysis on absorption spectra by performing a least-squares fit to the linestrength (in units cm^2).
%  A set of Judd-Ofelt parameters are found and used to predict branching ratios and transition probabilities for excited manifolds to all lower
%  manifolds. ED (Electric Dipole) and MD (Magnetic Dipole) transition probabilities are calculated for the full analysis. Radiative lifetimes
%  reflect both contributions (ED & MD) to the transition probability. 
%
% Program work in 3 different modes
% 
% Mode 0: when the mean wavelength and bandsum is available, we have to use this mode In this case we
% have to provide provide bandsum (bandsum),  bandsom energy levels(BEL) and mean wavelength(meanv),
% 
% Mode 1:When pure absolution cross section is available the code can compute and analyses the data. In
% this case. We have to provide absorption cross section boundaries(BSB) and bandsom energy levels(BEL), 
% select from excel files
% 
% Mode 2: when the Judd-Ofelt Parameters are available, we can analysis the life time , branching ratio and so
% on. We should provide Judd-Ofelt Parameters(omega)
% 
%For more information please refer to code documentation
%--------------------------------------------------------------------
%Writen by Siamak Dawazdah Emami
%Date:
%--------------------------------------------------------------------
clc
clear
close all
%--------------------Define variable here--------------------------


mode=1;     %0: provide bandsum (bandsum),  bandsom energy levels(BEL) and mean wavelength(meanv), 
                    %1: provide absorption cross section boundaries(BSB) and bandsom energy levels(BEL), select from  excel files  
                     %2: provide Judd-Ofelt Parameters(omega)
                 
rare=1;                   % rare earth element available data(1:thulium, 2:erbium,3:holmium, 4:terbium 5:dysprosium )
host= 'SIO2' ;                      % (YAG , LUAG, YSAG, YGG, LUGG, GGG, GSGG, YSGG, YSAG
                                           % BFAP, Q246, SIO2, SILICATE, ZBLAN, Y2O3, GEO2)
%mode 0 example
NBS2=6; %number of bandsum
BEL2=[6 7 9 10 11 12]
bandsum2=[2.865,2.815,10.707 ,19.692,34.82 ,86.911];        % bandsum matrix 
meanv2=[3590,470,682,789,1192,1708];               % mean value matrix 

%mode 1 example 
NBS=4; %number of bandsum
BSB=[0,507,1846,3866] ;  % Boundaries of bandsom
BEL=[8, 10, 11, 12]    ;       %bandsom energy levels, select from  excel files 

%mode2 example      
omega=[7.7479  2.5110   0.7599]

if mode ==2
Trp( omega,host,rare );
elseif mode==0
[omega2,slm]= JOFA( NBS2,mode,BSB,BEL2,bandsum2,meanv2,rare,host );
Trp( omega2,host,rare );
elseif mode==1
    bandsum=0
    meanv=0
[omega2,slm]= JOFA( NBS,mode,BSB,BEL,bandsum,meanv,rare,host );
Trp( omega2,host,rare );
end






