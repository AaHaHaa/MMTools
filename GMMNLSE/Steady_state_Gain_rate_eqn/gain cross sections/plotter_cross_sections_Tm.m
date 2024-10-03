clearvars; close all;

addpath('../../user_helpers/');

filename = 'Tm.txt';

%% Reading data from the specified file
delimiter = ',';
startRow = 1;
endRow = inf;

formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

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
 GSA05_raw = dataArray{:, 2}; % ground-state absorption (0 --> 5)
 GSA04_raw = dataArray{:, 3}; % ground-state absorption (0 --> 4)
 GSA03_raw = dataArray{:, 4}; % ground-state absorption (0 --> 3)
 GSA02_raw = dataArray{:, 5}; % ground-state absorption (0 --> 2)
 GSA01_raw = dataArray{:, 6}; % ground-state absorption (0 --> 1)
 ESA14_raw = dataArray{:, 7}; % excited-state absorption (1 --> 4)
 ESA35_raw = dataArray{:, 8}; % excited-state absorption (3 --> 5)
 ESA13_raw = dataArray{:, 9}; % excited-state absorption (1 --> 3)
 ESA23_raw = dataArray{:,10}; % excited-state absorption (2 --> 3)
emi50_raw  = dataArray{:,11}; % emission (5 --> 0)
emi30_raw  = dataArray{:,12}; % emission (3 --> 0)
emi31_raw  = dataArray{:,13}; % emission (3 --> 1)
emi10_raw  = dataArray{:,14}; % emission (1 --> 0)
emi32_raw  = dataArray{:,15}; % emission (3 --> 2)

% Smooth the data
[GSA05_raw,GSA04_raw,GSA03_raw,GSA02_raw,GSA01_raw,...
 ESA14_raw,ESA35_raw,ESA13_raw,ESA23_raw,...
 emi50_raw,emi30_raw,emi31_raw,emi10_raw,emi32_raw] = smooth_cross_sections(lambda_raw,...
                                                                            GSA05_raw,GSA04_raw,GSA03_raw,GSA02_raw,GSA01_raw,...
                                                                            ESA14_raw,ESA35_raw,ESA13_raw,ESA23_raw,...
                                                                            emi50_raw,emi30_raw,emi31_raw,emi10_raw,emi32_raw);

%% Plot
cross_sections = [GSA05_raw,GSA04_raw,GSA03_raw,GSA02_raw,GSA01_raw,...
                  ESA14_raw,ESA35_raw,ESA13_raw,ESA23_raw,...
                  emi50_raw,emi30_raw,emi31_raw,emi10_raw,emi32_raw];
figure;
ccc = distinguishable_colors(size(cross_sections,2));
h = plot(lambda_raw*1e9,cross_sections*1e24,'linewidth',2);
for i = 1:size(cross_sections,2)
    set(h(i),'Color',ccc(i,:));
end
xlim([min(lambda_raw*1e9),max(lambda_raw*1e9)]);
set(gca,'fontsize',20);
xlabel('Wavelength (nm)');
ylabel('Cross section (pm^2)');
l = legend('GSA05','GSA04','GSA03','GSA02','GSA01','ESA14','ESA35','ESA13','ESA23','emi50','emi30','emi31','emi10','emi32');
set(l,'fontsize',6,'location','northeast');

print(gcf,'Tm.jpg','-djpeg');