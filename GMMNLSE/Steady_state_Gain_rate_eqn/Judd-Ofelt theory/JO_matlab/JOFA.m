function [omega,slm]= JOFA( NBS,bandsum_mode,BSB,BEL,bandsum,meanv,rare,host )
%JOFA Summary of this function goes here
%   Detailed explanation goes here
 %-------------------- bandsum and mean wavelength calculation--------------------------    
 TAGS =2*NBS+1; %TOTAL ANGULAR MOMENTUM OF GROUND STATE
if bandsum_mode==1
      data=csvread('sigma_sa.csv'); %Read ecxel file related to absorption
    lamda=data(:,1)';
    sigma_sa=data(:,2)';

    for i=1:NBS-1
    lamdap=lamda(BSB(i)+1:BSB(i+1));
    sigma_sap=sigma_sa(BSB(i)+1:BSB(i+1));
    figure(i)
    area(lamdap,sigma_sap)
     xlabel( 'wavelength(nm) ')
     ylabel( 'Absorption cross section ')
    

    bandsum(i)=trapz(lamdap,sigma_sap)*1e24;
    meanv(i)=(sum(lamdap.*sigma_sap)./sum(sigma_sap))./1000;
    textd=['   Mean wavelength = ',num2str(meanv(i)*1000) ,' nm','       Bandsum =' ,num2str(bandsum(i)), '  10^-20 cm^2-nm'];
     text(min(lamdap),max(sigma_sap),textd)
    end
    lamdap=lamda(BSB(end)+1:end);
    sigma_sap=sigma_sa(BSB(end)+1:end);
    bandsum(NBS)=trapz(lamdap,sigma_sap)*1e24;
     figure(NBS)
     area(lamdap,sigma_sap)
     xlabel( 'wavelength(nm) ')
     ylabel( 'Absorption cross section ')
    meanv(NBS)=(sum(lamdap.*sigma_sap)./sum(sigma_sap))./1000;         
  
    textd=['   Mean wavelength = ',num2str(meanv(NBS)*1000) ,' nm','       Bandsum =' ,num2str(bandsum(NBS)),'  10^-20 cm^2-nm'];
     text(min(lamdap),max(sigma_sap),textd)
   
     if (NBS~=length(BSB)) || (NBS~=length(BEL))
        disp('error--------> number of bandsom should be the same as number or boundaries')
        return
    end

end 

 %-------------------- Load material and host data--------------------------                                                
if rare ==1
   load('thulium.mat');
   dfa='Thulium Dope Fiber';
elseif rare ==2
    load('erbium.mat');
   dfa='Erbium Dope Fiber';
  
end



for i=1:NBS
x(i,:)=raree(BEL(i),:);
end

 %-------------------- Sellmeier calculation--------------------------        
load('sellmeier.mat','material','sellmeier');
for n=1:16
    if strcmp(material{n},host)
        slm=sellmeier(n,:);
        break
    end
end

for i=1:NBS
N(i)=sqrt(slm(1)+slm(2)*meanv(i)^2/(meanv(i)^2-slm(3))+slm(4)*meanv(i)^2/(meanv(i)^2-slm(5)));
NS(i)=N(i)*(3/(N(i)^2+2))^2;
NL(i)=N(i)*((N(i)^2+2)/3)^2;
end

%-------------------- Omega and LINESTRENGTH(SE) calculation--------------------------       
SE = bandsum.*10.413487 .*NS.*TAGS./( meanv*1000); %(LINESTRENGTH: NOTE THAT 3hc/8pi^3e^2=10.413487)
omega=inv((x'*x))*(x'*SE'); %JUDD-OFELT LEAST SQUARES FIT FOR ABSORPTION LINE STRENGTHS
ST= (x*omega);

%--------------------FIND RADIATIVE LIFETIMES -------------------------    

EA = (7.2166e-10.*NL) .* ST' ./ (TAGS .* ((meanv*1000).*.0000001) .^ 3);
RL = 1000*(1 ./ EA);

%--------------------MEAN SQUARE ERROR -------------------------   

SSUM =(sum((SE - ST') .^ 2))/(TAGS-3);
RMS = sqrt(SSUM); %ROOT MEAN SQUARE ERROR
SERR=(sqrt(SSUM*diag(x)));


%------------------- OUTPUT FILE GENARATION-------

fileID = fopen('output.txt','w');
fprintf(fileID,'%1s \r\n','-----------------------------------------------------------------');
fprintf(fileID,'%1s \r\n','                    Output File ');
fprintf(fileID,'%1s \r\n','-----------------------------------------------------------------');
fprintf(fileID,'%1s \r\n',dfa);
fprintf(fileID,'%1s','Host is :');
fprintf(fileID,'%1s \r\n',host);
fprintf(fileID,'%1s \r\n','sellmeier coefficient');
fprintf(fileID,'%3s %5.2f %6s %5.2f %6s %5.2f  %6s %5.2f %6s %5.2f \r\n','A=', slm(1),'B=', slm(2),'C', slm(3),'D=', slm(4),'E=', slm(5));

fprintf(fileID,'%1s \r\n','  ');
fprintf(fileID,'%1s \r\n','                   Band sum & Mean value  ');
fprintf(fileID,'%1s \r\n','____________________________________________________________________');

fprintf(fileID,'%9s %15s %20s %18s \r\n','Transition', 'Band Sum' , 'Mean Wavelength','Refractive Index' );
for i=1:NBS
    fprintf(fileID,'%4s %18.8f %19.2f %15.2f  \r\n',layerM{BEL(i)}, bandsum(i), meanv(i)*1000,N(i));
end
fprintf(fileID,'%1s \r\n','____________________________________________________________________');

fprintf(fileID,'%1s  \r\n','AVERAGE WAVELENGTH OF TRANSITION IN nm' );
fprintf(fileID,'%1s  \r\n','INTEGRATED ABSORPTION OF MANIFOLD (BANDSUM) IN UNITS 10^-20 cm^2-nm' );
fprintf(fileID,'%1s \r\n','  ');

fprintf(fileID,'%1s \r\n','  ');
fprintf(fileID,'%1s \r\n','                                            JOF Analysis  ');
fprintf(fileID,'%1s \r\n','________________________________________________________________________________________');

fprintf(fileID,'%9s %10s %10s %10s %18s %15s %10s \r\n','Transition', 'U(2)',	'U(4)',	'U(6)',	'S(Experiment)',	'S(Theoritical)', 'Life time');

for i=1:NBS
    fprintf(fileID,'%4s %13.5f %10.5f %10.5f %14.5f %14.5f %14.5f \r\n',layerM{BEL(i)}, x(i,1), x(i,2),x(i,3),SE(i),ST(i),RL(i));
end
fprintf(fileID,'%1s \r\n','________________________________________________________________________________________');

fprintf(fileID,'%1s \r\n','  ');
fprintf(fileID,'%1s  \r\n','Judd-Ofelt Parameters' );
fprintf(fileID,'%12s %5.2f %12s %5.2f %12s %5.2f  \r\n','Omega(2)=', omega(1),'Omega(4)=', omega(2),'Omega(6)', omega(3));
fprintf(fileID,'%1s \r\n','  ');
fprintf(fileID,'%1s  \r\n','Standard Error In Judd-Ofelt Parameters' );
fprintf(fileID,'%12s %5.2f %12s %5.2f %12s %5.2f  \r\n','Omega(2)=', SERR(1),'Omega(4)=', SERR(2),'Omega(6)', SERR(3));

fprintf(fileID,'%1s \r\n','  ');
fprintf(fileID,'%1s  \r\n','Root Mean Square Error' );
fprintf(fileID,'%1s %5.5f   \r\n','RMS=', RMS);
fclose(fileID)

end

