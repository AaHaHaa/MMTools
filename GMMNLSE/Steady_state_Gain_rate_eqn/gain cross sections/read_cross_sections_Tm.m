function [GSA05,GSA04,GSA03,GSA02,GSA01,...
          ESA14,ESA35,ESA13,ESA23,...
          emi50,emi30,emi31,emi10,emi32] = read_cross_sections_Tm(filename, lambda)
%READ_CROSS_SECTIONS_TM This reads the absorption and emission cross 
%sections of Tm
%
% The data comes from
%
% 1. Walsh and Barnes, "Comparison of Tm:ZBLAN and Tm:silica fiber lasers;
%    Spectroscopy and tunable pulsed laser operation around 1.9 μm," Appl.
%    Phys. B 78(3), 325-333 (2004)
% 2. Peterka et al., "Theoretical modeling of fiber laser at 810 nm based
%    on thulium-doped silica fibers with enhanced ^3H_4 level lifetime,"
%    Opt. Express 19(3), 2773-2781 (2011)
% 3. Muravyev et al., "Dual-band Tm^3+-doped tellurite fiber amplifier and
%    laser at 1.9 μm and 2.3 μm," Sci. Rep. 8(1), 16164 (2018)
% 4. Jackson and King, "Theoretical modeling of Tm-doped silica fiber
%    lasers," J. Light. Technol. 17(5), 948-956 (1999)

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

%% Interpolating data for our wavelength grid
lambda_raw_min = lambda_raw(1);
lambda_raw_max = lambda_raw(end);
lambda_inside_lambda_raw = lambda(lambda<=lambda_raw_max & lambda>=lambda_raw_min); % interpolation
lambda_outside_lambda_raw_left  = lambda(lambda<lambda_raw_min); % extrapolation
lambda_outside_lambda_raw_right = lambda(lambda>lambda_raw_max); % extrapolation

% interpolation
GSA05 = [zeros(size(lambda_outside_lambda_raw_left)); ...
           interp1(lambda_raw,GSA05_raw, lambda_inside_lambda_raw, 'pchip'); ...
         zeros(size(lambda_outside_lambda_raw_right))];
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
ESA35 = [zeros(size(lambda_outside_lambda_raw_left)); ...
           interp1(lambda_raw,ESA35_raw, lambda_inside_lambda_raw, 'pchip'); ...
         zeros(size(lambda_outside_lambda_raw_right))];
ESA13 = [zeros(size(lambda_outside_lambda_raw_left)); ...
           interp1(lambda_raw,ESA13_raw, lambda_inside_lambda_raw, 'pchip'); ...
         zeros(size(lambda_outside_lambda_raw_right))];
ESA23 = [zeros(size(lambda_outside_lambda_raw_left)); ...
           interp1(lambda_raw,ESA23_raw, lambda_inside_lambda_raw, 'pchip'); ...
         zeros(size(lambda_outside_lambda_raw_right))];
emi50 = [zeros(size(lambda_outside_lambda_raw_left)); ...
           interp1(lambda_raw,emi50_raw, lambda_inside_lambda_raw, 'pchip'); ...
         zeros(size(lambda_outside_lambda_raw_right))];
emi30 = [zeros(size(lambda_outside_lambda_raw_left)); ...
           interp1(lambda_raw,emi30_raw, lambda_inside_lambda_raw, 'pchip'); ...
         zeros(size(lambda_outside_lambda_raw_right))];
emi31 = [zeros(size(lambda_outside_lambda_raw_left)); ...
           interp1(lambda_raw,emi31_raw, lambda_inside_lambda_raw, 'pchip'); ...
         zeros(size(lambda_outside_lambda_raw_right))];
emi10 = [zeros(size(lambda_outside_lambda_raw_left)); ...
           interp1(lambda_raw,emi10_raw, lambda_inside_lambda_raw, 'pchip'); ...
         zeros(size(lambda_outside_lambda_raw_right))];
emi32 = [zeros(size(lambda_outside_lambda_raw_left)); ...
           interp1(lambda_raw,emi32_raw, lambda_inside_lambda_raw, 'pchip'); ...
         zeros(size(lambda_outside_lambda_raw_right))];

% extrapolation with an exponential decay
% The left part
if ~isempty(lambda_outside_lambda_raw_left)
    % B.C.: Slope at the edge
    GSA05_raw_slope_left = (GSA05_raw(2)-GSA05_raw(1))/(lambda_raw(2)-lambda_raw(1));
    GSA04_raw_slope_left = (GSA04_raw(2)-GSA04_raw(1))/(lambda_raw(2)-lambda_raw(1));
    GSA03_raw_slope_left = (GSA03_raw(2)-GSA03_raw(1))/(lambda_raw(2)-lambda_raw(1));
    GSA02_raw_slope_left = (GSA02_raw(2)-GSA02_raw(1))/(lambda_raw(2)-lambda_raw(1));
    GSA01_raw_slope_left = (GSA01_raw(2)-GSA01_raw(1))/(lambda_raw(2)-lambda_raw(1));
    ESA14_raw_slope_left = (ESA14_raw(2)-ESA14_raw(1))/(lambda_raw(2)-lambda_raw(1));
    ESA35_raw_slope_left = (ESA35_raw(2)-ESA35_raw(1))/(lambda_raw(2)-lambda_raw(1));
    ESA13_raw_slope_left = (ESA13_raw(2)-ESA13_raw(1))/(lambda_raw(2)-lambda_raw(1));
    ESA23_raw_slope_left = (ESA23_raw(2)-ESA23_raw(1))/(lambda_raw(2)-lambda_raw(1));
    emi50_raw_slope_left = (emi50_raw(2)-emi50_raw(1))/(lambda_raw(2)-lambda_raw(1));
    emi30_raw_slope_left = (emi30_raw(2)-emi30_raw(1))/(lambda_raw(2)-lambda_raw(1));
    emi31_raw_slope_left = (emi31_raw(2)-emi31_raw(1))/(lambda_raw(2)-lambda_raw(1));
    emi10_raw_slope_left = (emi10_raw(2)-emi10_raw(1))/(lambda_raw(2)-lambda_raw(1));
    emi32_raw_slope_left = (emi32_raw(2)-emi32_raw(1))/(lambda_raw(2)-lambda_raw(1));
    % Compute extrapolations
    GSA05(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),GSA05_raw(1),max(0,GSA05_raw_slope_left),lambda_outside_lambda_raw_left);
    GSA04(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),GSA04_raw(1),max(0,GSA04_raw_slope_left),lambda_outside_lambda_raw_left);
    GSA03(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),GSA03_raw(1),max(0,GSA03_raw_slope_left),lambda_outside_lambda_raw_left);
    GSA02(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),GSA02_raw(1),max(0,GSA02_raw_slope_left),lambda_outside_lambda_raw_left);
    GSA01(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),GSA01_raw(1),max(0,GSA01_raw_slope_left),lambda_outside_lambda_raw_left);
    ESA14(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),ESA14_raw(1),max(0,ESA14_raw_slope_left),lambda_outside_lambda_raw_left);
    ESA35(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),ESA35_raw(1),max(0,ESA35_raw_slope_left),lambda_outside_lambda_raw_left);
    ESA13(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),ESA13_raw(1),max(0,ESA13_raw_slope_left),lambda_outside_lambda_raw_left);
    ESA23(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),ESA23_raw(1),max(0,ESA23_raw_slope_left),lambda_outside_lambda_raw_left);
    emi50(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),emi50_raw(1),max(0,emi50_raw_slope_left),lambda_outside_lambda_raw_left);
    emi30(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),emi30_raw(1),max(0,emi30_raw_slope_left),lambda_outside_lambda_raw_left);
    emi31(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),emi31_raw(1),max(0,emi31_raw_slope_left),lambda_outside_lambda_raw_left);
    emi10(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),emi10_raw(1),max(0,emi10_raw_slope_left),lambda_outside_lambda_raw_left);
    emi32(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),emi32_raw(1),max(0,emi32_raw_slope_left),lambda_outside_lambda_raw_left);
end
% The right part
if ~isempty(lambda_outside_lambda_raw_right)
    % B.C.: Slope at the edge
    GSA05_raw_slope_right = (GSA05_raw(end)-GSA05_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    GSA04_raw_slope_right = (GSA04_raw(end)-GSA04_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    GSA03_raw_slope_right = (GSA03_raw(end)-GSA03_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    GSA02_raw_slope_right = (GSA02_raw(end)-GSA02_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    GSA01_raw_slope_right = (GSA01_raw(end)-GSA01_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    ESA14_raw_slope_right = (ESA14_raw(end)-ESA14_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    ESA35_raw_slope_right = (ESA35_raw(end)-ESA35_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    ESA13_raw_slope_right = (ESA13_raw(end)-ESA13_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    ESA23_raw_slope_right = (ESA23_raw(end)-ESA23_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    emi50_raw_slope_right = (emi50_raw(end)-emi50_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    emi30_raw_slope_right = (emi30_raw(end)-emi30_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    emi31_raw_slope_right = (emi31_raw(end)-emi31_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    emi10_raw_slope_right = (emi10_raw(end)-emi10_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    emi32_raw_slope_right = (emi32_raw(end)-emi32_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    % Compute extrapolations
    GSA05(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),GSA05_raw(end),min(0,GSA05_raw_slope_right),lambda_outside_lambda_raw_right);
    GSA04(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),GSA04_raw(end),min(0,GSA04_raw_slope_right),lambda_outside_lambda_raw_right);
    GSA03(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),GSA03_raw(end),min(0,GSA03_raw_slope_right),lambda_outside_lambda_raw_right);
    GSA02(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),GSA02_raw(end),min(0,GSA02_raw_slope_right),lambda_outside_lambda_raw_right);
    GSA01(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),GSA01_raw(end),min(0,GSA01_raw_slope_right),lambda_outside_lambda_raw_right);
    ESA14(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),ESA14_raw(end),min(0,ESA14_raw_slope_right),lambda_outside_lambda_raw_right);
    ESA35(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),ESA35_raw(end),min(0,ESA35_raw_slope_right),lambda_outside_lambda_raw_right);
    ESA13(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),ESA13_raw(end),min(0,ESA13_raw_slope_right),lambda_outside_lambda_raw_right);
    ESA23(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),ESA23_raw(end),min(0,ESA23_raw_slope_right),lambda_outside_lambda_raw_right);
    emi50(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),emi50_raw(end),min(0,emi50_raw_slope_right),lambda_outside_lambda_raw_right);
    emi30(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),emi30_raw(end),min(0,emi30_raw_slope_right),lambda_outside_lambda_raw_right);
    emi31(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),emi31_raw(end),min(0,emi31_raw_slope_right),lambda_outside_lambda_raw_right);
    emi10(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),emi10_raw(end),min(0,emi10_raw_slope_right),lambda_outside_lambda_raw_right);
    emi32(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),emi32_raw(end),min(0,emi32_raw_slope_right),lambda_outside_lambda_raw_right);
end

% Reverse back the original order of "lambda"
if lambda_sort
    GSA05 = GSA05(reverse_sort_idx);
    GSA04 = GSA04(reverse_sort_idx);
    GSA03 = GSA03(reverse_sort_idx);
    GSA02 = GSA02(reverse_sort_idx);
    GSA01 = GSA01(reverse_sort_idx);
    ESA14 = ESA14(reverse_sort_idx);
    ESA35 = ESA35(reverse_sort_idx);
    ESA13 = ESA13(reverse_sort_idx);
    ESA23 = ESA23(reverse_sort_idx);
    emi50 = emi50(reverse_sort_idx);
    emi30 = emi30(reverse_sort_idx);
    emi31 = emi31(reverse_sort_idx);
    emi10 = emi10(reverse_sort_idx);
    emi32 = emi32(reverse_sort_idx);
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