clearvars; close all;

filename = 'Nd (Martin thesis 2006).txt';
%filename = 'Er (optiwave and Smirnov).txt';

%% Reading data from the specified file
delimiter = ',';
startRow = 1;
endRow = inf;

formatSpec = '%f%f%f%f%[^\n\r]';

fileID = fopen(filename,'r');

dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

fclose(fileID);

%% Defining the raw data (from the file)
lambda_raw = dataArray{:, 1};
absorption_raw = dataArray{:, 2};
emission_raw = dataArray{:, 3};
ESA_raw = dataArray{:, 4}; % excited-state absorption

% Smooth the data by smoothing the derivative, then integrate it back
% to make sure a C1 curve
[absorption_raw,emission_raw,ESA_raw] = smooth_cross_sections(lambda_raw,absorption_raw,emission_raw,ESA_raw);

%% Plot
figure;
h = plot(lambda_raw*1e9,[absorption_raw,emission_raw,ESA_raw]*1e24,'linewidth',2);
set(h(1),'Color','b'); set(h(2),'Color','r');
xlim([min(lambda_raw*1e9),max(lambda_raw*1e9)]);
set(gca,'fontsize',20);
xlabel('Wavelength (nm)');
ylabel('Cross section (pm^2)');
l = legend('Absorption','Emission','ESA');
set(l,'fontsize',12,'location','northeast');

%print(gcf,'Nd (Martin thesis 2006).jpg','-djpeg');
%print('Er (optiwave and Smirnov).jpg','-djpeg');