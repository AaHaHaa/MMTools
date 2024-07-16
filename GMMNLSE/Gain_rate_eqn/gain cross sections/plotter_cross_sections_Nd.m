clearvars; close all;

filename = 'Nd.txt';

%% Reading data from the specified file
delimiter = ',';
startRow = 1;
endRow = inf;

formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

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
GSA011_raw = dataArray{:, 2}; % ground-state absorption (0 --> 11)
GSA010_raw = dataArray{:, 3}; % ground-state absorption (0 --> 10)
 GSA09_raw = dataArray{:, 4}; % ground-state absorption (0 --> 9)
 GSA08_raw = dataArray{:, 5}; % ground-state absorption (0 --> 8)
 GSA07_raw = dataArray{:, 6}; % ground-state absorption (0 --> 7)
 GSA06_raw = dataArray{:, 7}; % ground-state absorption (0 --> 6)
 GSA05_raw = dataArray{:, 8}; % ground-state absorption (0 --> 5)
 GSA04_raw = dataArray{:, 9}; % ground-state absorption (0 --> 4)
 GSA03_raw = dataArray{:,10}; % ground-state absorption (0 --> 3)
 ESA49_raw = dataArray{:,11}; % excited-state absorption (4 --> 9)
emi40_raw  = dataArray{:,12}; % emission (4 --> 0)
emi41_raw  = dataArray{:,13}; % emission (4 --> 1)
emi42_raw  = dataArray{:,14}; % emission (4 --> 2)

% Smooth the data
[GSA011_raw,GSA010_raw,GSA09_raw,GSA08_raw,GSA07_raw,GSA06_raw,GSA05_raw,GSA04_raw,GSA03_raw,...
 ESA49_raw,...
 emi40_raw,emi41_raw,emi42_raw] = smooth_cross_sections(lambda_raw,...
                                                        GSA011_raw,GSA010_raw,GSA09_raw,GSA08_raw,GSA07_raw,GSA06_raw,GSA05_raw,GSA04_raw,GSA03_raw,...
                                                        ESA49_raw,...
                                                        emi40_raw,emi41_raw,emi42_raw);

%% Plot
cross_sections = [GSA011_raw,GSA010_raw,GSA09_raw,GSA08_raw,GSA07_raw,GSA06_raw,GSA05_raw,GSA04_raw,GSA03_raw,...
                  ESA49_raw,...
                  emi40_raw,emi41_raw,emi42_raw];
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
l = legend('GSA011','GSA010','GSA09','GSA08','GSA07','GSA06','GSA05','GSA04','GSA03','ESA49','emi40','emi41','emi42');
set(l,'fontsize',6,'location','northeast');

print(gcf,'Nd.jpg','-djpeg');