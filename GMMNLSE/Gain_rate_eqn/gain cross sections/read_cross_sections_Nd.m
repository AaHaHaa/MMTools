function [GSA011,GSA010,GSA09,GSA08,GSA07,GSA06,GSA05,GSA04,GSA03,...
          ESA49,...
          emi40,emi41,emi42] = read_cross_sections_Nd(filename, lambda)
%READ_CROSS_SECTIONS_Nd This reads the absorption and emission cross 
%sections of Nd
%
% The data comes from
%
% 1. Verlinden et al., "The Excited State Absorption Cross Section of
%    Neodymium-doped Silica Glass Fiber in the 1200-1500 nm Wavelength 
%    Range," Worcester Polytechnic Institute, PhD Thesis (2008)
% 2. George Robert Geddes, "Amplified Spontaneous Emission in Nd-Doped
%    Fiber Amplifiers," Worcester Polytechnic Institute (2011)

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

%% Interpolating data for our wavelength grid
lambda_raw_min = lambda_raw(1);
lambda_raw_max = lambda_raw(end);
lambda_inside_lambda_raw = lambda(lambda<=lambda_raw_max & lambda>=lambda_raw_min); % interpolation
lambda_outside_lambda_raw_left  = lambda(lambda<lambda_raw_min); % extrapolation
lambda_outside_lambda_raw_right = lambda(lambda>lambda_raw_max); % extrapolation

% interpolation
GSA011 = [zeros(size(lambda_outside_lambda_raw_left)); ...
           interp1(lambda_raw,GSA011_raw, lambda_inside_lambda_raw, 'pchip'); ...
         zeros(size(lambda_outside_lambda_raw_right))];
GSA010 = [zeros(size(lambda_outside_lambda_raw_left)); ...
           interp1(lambda_raw,GSA010_raw, lambda_inside_lambda_raw, 'pchip'); ...
         zeros(size(lambda_outside_lambda_raw_right))];
GSA09 = [zeros(size(lambda_outside_lambda_raw_left)); ...
           interp1(lambda_raw,GSA09_raw, lambda_inside_lambda_raw, 'pchip'); ...
         zeros(size(lambda_outside_lambda_raw_right))];
GSA08 = [zeros(size(lambda_outside_lambda_raw_left)); ...
           interp1(lambda_raw,GSA08_raw, lambda_inside_lambda_raw, 'pchip'); ...
         zeros(size(lambda_outside_lambda_raw_right))];
GSA07 = [zeros(size(lambda_outside_lambda_raw_left)); ...
           interp1(lambda_raw,GSA07_raw, lambda_inside_lambda_raw, 'pchip'); ...
         zeros(size(lambda_outside_lambda_raw_right))];
GSA06 = [zeros(size(lambda_outside_lambda_raw_left)); ...
           interp1(lambda_raw,GSA06_raw, lambda_inside_lambda_raw, 'pchip'); ...
         zeros(size(lambda_outside_lambda_raw_right))];
GSA05 = [zeros(size(lambda_outside_lambda_raw_left)); ...
           interp1(lambda_raw,GSA05_raw, lambda_inside_lambda_raw, 'pchip'); ...
         zeros(size(lambda_outside_lambda_raw_right))];
GSA04 = [zeros(size(lambda_outside_lambda_raw_left)); ...
           interp1(lambda_raw,GSA04_raw, lambda_inside_lambda_raw, 'pchip'); ...
         zeros(size(lambda_outside_lambda_raw_right))];
GSA03 = [zeros(size(lambda_outside_lambda_raw_left)); ...
           interp1(lambda_raw,GSA03_raw, lambda_inside_lambda_raw, 'pchip'); ...
         zeros(size(lambda_outside_lambda_raw_right))];
ESA49 = [zeros(size(lambda_outside_lambda_raw_left)); ...
           interp1(lambda_raw,ESA49_raw, lambda_inside_lambda_raw, 'pchip'); ...
         zeros(size(lambda_outside_lambda_raw_right))];
emi40 = [zeros(size(lambda_outside_lambda_raw_left)); ...
           interp1(lambda_raw,emi40_raw, lambda_inside_lambda_raw, 'pchip'); ...
         zeros(size(lambda_outside_lambda_raw_right))];
emi41 = [zeros(size(lambda_outside_lambda_raw_left)); ...
           interp1(lambda_raw,emi41_raw, lambda_inside_lambda_raw, 'pchip'); ...
         zeros(size(lambda_outside_lambda_raw_right))];
emi42 = [zeros(size(lambda_outside_lambda_raw_left)); ...
           interp1(lambda_raw,emi42_raw, lambda_inside_lambda_raw, 'pchip'); ...
         zeros(size(lambda_outside_lambda_raw_right))];

% extrapolation with an exponential decay
% The left part
if ~isempty(lambda_outside_lambda_raw_left)
    % B.C.: Slope at the edge
    GSA011_raw_slope_left = (GSA011_raw(2)-GSA011_raw(1))/(lambda_raw(2)-lambda_raw(1));
    GSA010_raw_slope_left = (GSA010_raw(2)-GSA010_raw(1))/(lambda_raw(2)-lambda_raw(1));
    GSA09_raw_slope_left = (GSA09_raw(2)-GSA09_raw(1))/(lambda_raw(2)-lambda_raw(1));
    GSA08_raw_slope_left = (GSA08_raw(2)-GSA08_raw(1))/(lambda_raw(2)-lambda_raw(1));
    GSA07_raw_slope_left = (GSA07_raw(2)-GSA07_raw(1))/(lambda_raw(2)-lambda_raw(1));
    GSA06_raw_slope_left = (GSA06_raw(2)-GSA06_raw(1))/(lambda_raw(2)-lambda_raw(1));
    GSA05_raw_slope_left = (GSA05_raw(2)-GSA05_raw(1))/(lambda_raw(2)-lambda_raw(1));
    GSA04_raw_slope_left = (GSA04_raw(2)-GSA04_raw(1))/(lambda_raw(2)-lambda_raw(1));
    GSA03_raw_slope_left = (GSA03_raw(2)-GSA03_raw(1))/(lambda_raw(2)-lambda_raw(1));
    ESA49_raw_slope_left = (ESA49_raw(2)-ESA49_raw(1))/(lambda_raw(2)-lambda_raw(1));
    emi40_raw_slope_left = (emi40_raw(2)-emi40_raw(1))/(lambda_raw(2)-lambda_raw(1));
    emi41_raw_slope_left = (emi41_raw(2)-emi41_raw(1))/(lambda_raw(2)-lambda_raw(1));
    emi42_raw_slope_left = (emi42_raw(2)-emi42_raw(1))/(lambda_raw(2)-lambda_raw(1));
    % Compute extrapolations
    GSA011(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),GSA011_raw(1),max(0,GSA011_raw_slope_left),lambda_outside_lambda_raw_left);
    GSA010(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),GSA010_raw(1),max(0,GSA010_raw_slope_left),lambda_outside_lambda_raw_left);
    GSA09(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),GSA09_raw(1),max(0,GSA09_raw_slope_left),lambda_outside_lambda_raw_left);
    GSA08(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),GSA08_raw(1),max(0,GSA08_raw_slope_left),lambda_outside_lambda_raw_left);
    GSA07(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),GSA07_raw(1),max(0,GSA07_raw_slope_left),lambda_outside_lambda_raw_left);
    GSA06(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),GSA06_raw(1),max(0,GSA06_raw_slope_left),lambda_outside_lambda_raw_left);
    GSA05(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),GSA05_raw(1),max(0,GSA05_raw_slope_left),lambda_outside_lambda_raw_left);
    GSA04(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),GSA04_raw(1),max(0,GSA04_raw_slope_left),lambda_outside_lambda_raw_left);
    GSA03(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),GSA03_raw(1),max(0,GSA03_raw_slope_left),lambda_outside_lambda_raw_left);
    ESA49(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),ESA49_raw(1),max(0,ESA49_raw_slope_left),lambda_outside_lambda_raw_left);
    emi40(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),emi40_raw(1),max(0,emi40_raw_slope_left),lambda_outside_lambda_raw_left);
    emi41(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),emi41_raw(1),max(0,emi41_raw_slope_left),lambda_outside_lambda_raw_left);
    emi42(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),emi42_raw(1),max(0,emi42_raw_slope_left),lambda_outside_lambda_raw_left);
end
% The right part
if ~isempty(lambda_outside_lambda_raw_right)
    % B.C.: Slope at the edge
    GSA011_raw_slope_right = (GSA011_raw(end)-GSA011_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    GSA010_raw_slope_right = (GSA010_raw(end)-GSA010_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    GSA09_raw_slope_right = (GSA09_raw(end)-GSA09_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    GSA08_raw_slope_right = (GSA08_raw(end)-GSA08_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    GSA07_raw_slope_right = (GSA07_raw(end)-GSA07_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    GSA06_raw_slope_right = (GSA06_raw(end)-GSA06_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    GSA05_raw_slope_right = (GSA05_raw(end)-GSA05_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    GSA04_raw_slope_right = (GSA04_raw(end)-GSA04_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    GSA03_raw_slope_right = (GSA03_raw(end)-GSA03_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    ESA49_raw_slope_right = (ESA49_raw(end)-ESA49_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    emi40_raw_slope_right = (emi40_raw(end)-emi40_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    emi41_raw_slope_right = (emi41_raw(end)-emi41_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    emi42_raw_slope_right = (emi42_raw(end)-emi42_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    % Compute extrapolations
    GSA011(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),GSA011_raw(end),min(0,GSA011_raw_slope_right),lambda_outside_lambda_raw_right);
    GSA010(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),GSA010_raw(end),min(0,GSA010_raw_slope_right),lambda_outside_lambda_raw_right);
    GSA09(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),GSA09_raw(end),min(0,GSA09_raw_slope_right),lambda_outside_lambda_raw_right);
    GSA08(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),GSA08_raw(end),min(0,GSA08_raw_slope_right),lambda_outside_lambda_raw_right);
    GSA07(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),GSA07_raw(end),min(0,GSA07_raw_slope_right),lambda_outside_lambda_raw_right);
    GSA06(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),GSA06_raw(end),min(0,GSA06_raw_slope_right),lambda_outside_lambda_raw_right);
    GSA05(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),GSA05_raw(end),min(0,GSA05_raw_slope_right),lambda_outside_lambda_raw_right);
    GSA04(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),GSA04_raw(end),min(0,GSA04_raw_slope_right),lambda_outside_lambda_raw_right);
    GSA03(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),GSA03_raw(end),min(0,GSA03_raw_slope_right),lambda_outside_lambda_raw_right);
    ESA49(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),ESA49_raw(end),min(0,ESA49_raw_slope_right),lambda_outside_lambda_raw_right);
    emi40(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),emi40_raw(end),min(0,emi40_raw_slope_right),lambda_outside_lambda_raw_right);
    emi41(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),emi41_raw(end),min(0,emi41_raw_slope_right),lambda_outside_lambda_raw_right);
    emi42(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),emi42_raw(end),min(0,emi42_raw_slope_right),lambda_outside_lambda_raw_right);
end

% Reverse back the original order of "lambda"
if lambda_sort
    GSA011 = GSA011(reverse_sort_idx);
    GSA010 = GSA010(reverse_sort_idx);
    GSA09 = GSA09(reverse_sort_idx);
    GSA08 = GSA08(reverse_sort_idx);
    GSA07 = GSA07(reverse_sort_idx);
    GSA06 = GSA06(reverse_sort_idx);
    GSA05 = GSA05(reverse_sort_idx);
    GSA04 = GSA04(reverse_sort_idx);
    GSA03 = GSA03(reverse_sort_idx);
    ESA49 = ESA49(reverse_sort_idx);
    emi40 = emi40(reverse_sort_idx);
    emi41 = emi41(reverse_sort_idx);
    emi42 = emi42(reverse_sort_idx);
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