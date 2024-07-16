clearvars; close all;

addpath('../../user_helpers/');

filename = 'Ho.txt';

%% Reading data from the specified file
delimiter = ',';
startRow = 1;
endRow = inf;

formatSpec = '%f%f%f%f%f%f%f%[^\n\r]';

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
lambda_raw = dataArray{:,1};
 GSA04_raw = dataArray{:,2}; % ground-state absorption (0 --> 4)
 GSA03_raw = dataArray{:,3}; % ground-state absorption (0 --> 3)
 GSA02_raw = dataArray{:,4}; % ground-state absorption (0 --> 2)
 GSA01_raw = dataArray{:,5}; % ground-state absorption (0 --> 1)
 ESA14_raw = dataArray{:,6}; % excited-state absorption (1 --> 4)
emi10_raw  = dataArray{:,7}; % emission (1 --> 0)

% Smooth the data
[GSA04_raw,GSA03_raw,GSA02_raw,GSA01_raw,...
 ESA14_raw,...
 emi10_raw] = smooth_cross_sections(lambda_raw,...
                                    GSA04_raw,GSA03_raw,GSA02_raw,GSA01_raw,...
                                    ESA14_raw,...
                                    emi10_raw);

%% Plot
cross_sections = [GSA04_raw,GSA03_raw,GSA02_raw,GSA01_raw,...
                  ESA14_raw,...
                  emi10_raw];
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
l = legend('GSA04','GSA03','GSA02','GSA01','ESA14','emi10');
set(l,'fontsize',6,'location','northeast');

print(gcf,'Ho.jpg','-djpeg');