function [GSA04,GSA03,GSA02,GSA01,...
          ESA14,...
          emi10] = read_cross_sections_Ho(filename, lambda)
%READ_CROSS_SECTIONS_HO This reads the absorption and emission cross 
%sections of Ho
%
% The data comes from
%
% 1. Hemming et al., "A review of recent progress in holmium-doped silica
%    fibre sources," Opt. Fiber Technol. 20(6), 621-630 (2014)
% 2. Simakov et al., "High gain holmium-doped fibre amplifiers," Opt.
%    Express 24(13), 13946-13956 (2016)
% 3. Alharbi et al., "Performance Optimization of Holmium Doped Fiber
%    Amplifiers for Optical Communication Applications in 2–2.15 um
%    Wavelength Range," Photonics 9(4) (2022)

% "lambda" needs to be a column vector.
if size(lambda,1)==1
    lambda = lambda.';
end
% "lambda" needs to be in the order of small to large because of extrapolation below
% If not, sort it here and the result will be reverse back to the unsorted original one after computations.
if length(lambda)>1
    lambda_sort = true;
    [lambda,sort_idx] = sort(lambda);
    [~,reverse_sort_idx] = sort(sort_idx);
else
    lambda_sort = false;
end

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
lambda_raw = dataArray{:, 1};
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

%% Interpolating data for our wavelength grid
lambda_raw_min = lambda_raw(1);
lambda_raw_max = lambda_raw(end);
lambda_inside_lambda_raw = lambda(lambda<=lambda_raw_max & lambda>=lambda_raw_min); % interpolation
lambda_outside_lambda_raw_left  = lambda(lambda<lambda_raw_min); % extrapolation
lambda_outside_lambda_raw_right = lambda(lambda>lambda_raw_max); % extrapolation

% interpolation
GSA04 = [zeros(size(lambda_outside_lambda_raw_left)); ...
           interp1(lambda_raw,GSA04_raw, lambda_inside_lambda_raw, 'pchip'); ...
         zeros(size(lambda_outside_lambda_raw_right))];
GSA03 = [zeros(size(lambda_outside_lambda_raw_left)); ...
           interp1(lambda_raw,GSA03_raw, lambda_inside_lambda_raw, 'pchip'); ...
         zeros(size(lambda_outside_lambda_raw_right))];
GSA02 = [zeros(size(lambda_outside_lambda_raw_left)); ...
           interp1(lambda_raw,GSA02_raw, lambda_inside_lambda_raw, 'pchip'); ...
         zeros(size(lambda_outside_lambda_raw_right))];
GSA01 = [zeros(size(lambda_outside_lambda_raw_left)); ...
           interp1(lambda_raw,GSA01_raw, lambda_inside_lambda_raw, 'pchip'); ...
         zeros(size(lambda_outside_lambda_raw_right))];
ESA14 = [zeros(size(lambda_outside_lambda_raw_left)); ...
           interp1(lambda_raw,ESA14_raw, lambda_inside_lambda_raw, 'pchip'); ...
         zeros(size(lambda_outside_lambda_raw_right))];
emi10 = [zeros(size(lambda_outside_lambda_raw_left)); ...
           interp1(lambda_raw,emi10_raw, lambda_inside_lambda_raw, 'pchip'); ...
         zeros(size(lambda_outside_lambda_raw_right))];

% extrapolation with an exponential decay
% The left part
if ~isempty(lambda_outside_lambda_raw_left)
    % B.C.: Slope at the edge
    GSA04_raw_slope_left = (GSA04_raw(2)-GSA04_raw(1))/(lambda_raw(2)-lambda_raw(1));
    GSA03_raw_slope_left = (GSA03_raw(2)-GSA03_raw(1))/(lambda_raw(2)-lambda_raw(1));
    GSA02_raw_slope_left = (GSA02_raw(2)-GSA02_raw(1))/(lambda_raw(2)-lambda_raw(1));
    GSA01_raw_slope_left = (GSA01_raw(2)-GSA01_raw(1))/(lambda_raw(2)-lambda_raw(1));
    ESA14_raw_slope_left = (ESA14_raw(2)-ESA14_raw(1))/(lambda_raw(2)-lambda_raw(1));
    emi10_raw_slope_left = (emi10_raw(2)-emi10_raw(1))/(lambda_raw(2)-lambda_raw(1));
    % Compute extrapolations
    GSA04(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),GSA04_raw(1),max(0,GSA04_raw_slope_left),lambda_outside_lambda_raw_left);
    GSA03(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),GSA03_raw(1),max(0,GSA03_raw_slope_left),lambda_outside_lambda_raw_left);
    GSA02(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),GSA02_raw(1),max(0,GSA02_raw_slope_left),lambda_outside_lambda_raw_left);
    GSA01(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),GSA01_raw(1),max(0,GSA01_raw_slope_left),lambda_outside_lambda_raw_left);
    ESA14(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),ESA14_raw(1),max(0,ESA14_raw_slope_left),lambda_outside_lambda_raw_left);
    emi10(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),emi10_raw(1),max(0,emi10_raw_slope_left),lambda_outside_lambda_raw_left);
end
% The right part
if ~isempty(lambda_outside_lambda_raw_right)
    % B.C.: Slope at the edge
    GSA04_raw_slope_right = (GSA04_raw(end)-GSA04_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    GSA03_raw_slope_right = (GSA03_raw(end)-GSA03_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    GSA02_raw_slope_right = (GSA02_raw(end)-GSA02_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    GSA01_raw_slope_right = (GSA01_raw(end)-GSA01_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    ESA14_raw_slope_right = (ESA14_raw(end)-ESA14_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    emi10_raw_slope_right = (emi10_raw(end)-emi10_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    % Compute extrapolations
    GSA04(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),GSA04_raw(end),min(0,GSA04_raw_slope_right),lambda_outside_lambda_raw_right);
    GSA03(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),GSA03_raw(end),min(0,GSA03_raw_slope_right),lambda_outside_lambda_raw_right);
    GSA02(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),GSA02_raw(end),min(0,GSA02_raw_slope_right),lambda_outside_lambda_raw_right);
    GSA01(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),GSA01_raw(end),min(0,GSA01_raw_slope_right),lambda_outside_lambda_raw_right);
    ESA14(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),ESA14_raw(end),min(0,ESA14_raw_slope_right),lambda_outside_lambda_raw_right);
    emi10(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),emi10_raw(end),min(0,emi10_raw_slope_right),lambda_outside_lambda_raw_right);
end

% Reverse back the original order of "lambda"
if lambda_sort
    GSA04 = GSA04(reverse_sort_idx);
    GSA03 = GSA03(reverse_sort_idx);
    GSA02 = GSA02(reverse_sort_idx);
    GSA01 = GSA01(reverse_sort_idx);
    ESA14 = ESA14(reverse_sort_idx);
    emi10 = emi10(reverse_sort_idx);
end

end

function y = exponential_decay(x_edge,A,B,x)
%EXPONENTIAL_DECAY
%   It extrapolates with the function exp((+-)a*(x-b)).
%
%   B.C.:
%       continuity of points: exp((+-)a*(x_edge-b))       = A
%       continuity of slope:  (+-)a*exp((+-)a*(x_edge-b)) = B
%   
%   solution:
%       (+-)a = B/A
%       b = x_edge - ((+-)log(A))/a = x_edge - A*log(A)/B

if A == 0
    y = zeros(size(x));
elseif B == 0
    % choose somewhere near the center of spectrum as "b" and find "a"
    b = 1e-6; % wavelength: 1um
    pma = log(A)/(x_edge-b);
    
    % extrapolation
    y = exp(pma*(x-b));
else
    % coefficients of exponential decay function
    pma = B/A;
    b = x_edge - A*log(A)/B;

    % extrapolation
    y = exp(pma*(x-b));
end

y = y.*exp(-abs(x-x_edge)/5e-8); % make it decay within 50 nm

end