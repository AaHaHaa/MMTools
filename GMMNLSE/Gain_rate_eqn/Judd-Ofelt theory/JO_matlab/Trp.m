function Trp( omega,host,rare )
%T Summary of this function goes here
%   Detailed explanation goes here
if rare ==1
   load('thulium.mat')
   dfa='Thulium Dope Fiber';
   
  elseif rare ==2
    load('erbium.mat')
   dfa='Erbium Dope Fiber';
end

load('sellmeier.mat','material','sellmeier')
for n=1:16
    if strcmp(material{n},host)
        slm=sellmeier(n,:);
        break
    end
end
%--------------------JO TRANSITION PROBABILITY ----------------
fileID2 = fopen('transition.txt','w');
fprintf(fileID2,'%1s \r\n','                                            Judd Ofelt Analysis  ');
fprintf(fileID2,'%1s \r\n','       ');
fprintf(fileID2,'%12s %16s %12s %10s %19s %13s %15s \r\n','Transition', 'Wavelength',	'S(ED)',	'A(ED)',	'A(MD)',	'BETA', 'Life Time');
fprintf(fileID2,'%1s \r\n','_______________________________________________________________________________________________________');


s=0;
ni=length(energy)-1;
TSTRENGTH=zeros(1,50);
for i=1:ni
    
    for k=i+1:ni+1
        WT(k)=1e7./abs(energy(i)-energy(k)); % nm
         WM(k)=WT(k)/1000; % um
        N2=slm(1)+slm(2)*WM(k)^2/(WM(k)^2-slm(3))+slm(4)*WM(k)^2/(WM(k)^2-slm(5));
     
            N(k)=sqrt(slm(1)+slm(2).*WM(k).^2./(WM(k).^2-slm(3))+slm(4).*WM(k).^2./(WM(k).^2-slm(5)));
            NS(k)=N(k).*(3./(N(k).^2+2)).^2;
            NL(k)=N(k).*((N(k).^2+2)./3).^2;
            NDEX(k)=N(k);
        
        
  s=s+1;
        TSTRENGTH(k) =s1(s).*omega(1) +s2(s).*omega(2)+ s3(s).*omega(3) ;%LINESTRENGTH BY JUDD-OFELT EQUATION
        AED(k)=(7.235432e10*NL(k)) * (TSTRENGTH(k)*1e-20) / (num(i) * (WT(k)*.0000001) ^ 3);%ELECTRIC DIPOLE TRANSITION PROBABILITY
        AMD(k)=2.697348e10 * (NDEX(k))^3 / (num(i) * (WT(k)) ^ 3) * s4(s);
    end
    for KK=2: ni+1
    AEDSUM=0;
    AMDSUM=0;
        for JJ=1: KK
            AEDSUM=AED(JJ)+AEDSUM; %'SUM OF TRANSITION PROBABILITIES
            AMDSUM=AMD(JJ)+AMDSUM; %SUM OF MD TRANSITION PROBABILITIES
        end
        
    
    if AEDSUM~=0 
        LIF(KK)=(1./(AEDSUM+AMDSUM)).*1000; % RADIATIVE LIFETIME IN MILLISECONDS
    end
     

    
    
    
    end
    
    
    
for J=i+1:k    
BRATIO=(AED(J)+AMD(J))/(AEDSUM(1)+AMDSUM(1));%REM BRANCHING RATIO
    fprintf(fileID2,'%5s %1s %3s %15.3f %13.3f %15.3f  %13.3f %13.3f %13.3f  \r\n',layere{i},'-',layere{J},WT(J),TSTRENGTH(J),AED(J), AMD(J), BRATIO,LIF(J));

% SP=16-LEN(LEV$(I))-LEN(LEV$(J))-3
% PRINT#6,LEV$(I);" - ";LEV$(J);CHR$(9); 'SPACE$(SP);
% PRINT#6,USING("######.#",WT(J));CHR$(9);USING("##.####",TSTRENGTH(J));
% PRINT#6,CHR$(9);
% PRINT#6,USING("########.###",AED(J));CHR$(9);USING("#####.###",AMD(J));CHR$(9);USING("##.####",BRATIO);CHR$(9);USING("##########.####",LIF(J))
end
fprintf(fileID2,'%1s \r\n','  ');

 STRENGTH=zeros(1,k);
AED=zeros(1,k);
AMD=zeros(1,k);
end
fclose(fileID2)
end

