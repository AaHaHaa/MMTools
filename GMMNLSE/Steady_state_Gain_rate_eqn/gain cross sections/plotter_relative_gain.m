clearvars; close all;

%filename = 'optiwave Er.txt';
%filename = 'Liekki Yb_AV_20160530.txt';
filename = 'Yb_Gen_VIII_Cross_Section (Nufern).txt';

%% Reading data from the specified file
delimiter = ',';
startRow = 1;
endRow = inf;

formatSpec = '%f%f%f%[^\n\r]';

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

% Smooth the data by smoothing the derivative, then integrate it back
% to make sure a C1 curve
[absorption_raw,emission_raw] = smooth_cross_sections(lambda_raw,absorption_raw,emission_raw);

%% Calculate the gain/loss
N2 = linspace(0,0.5,10);
N1 = 1 - N2;
gain = emission_raw.*N2 - absorption_raw.*N1;

figure;
plot(lambda_raw*1e9,gain*1e24);