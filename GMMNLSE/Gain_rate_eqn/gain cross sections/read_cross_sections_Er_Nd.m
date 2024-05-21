function [absorption,emission,ESA] = read_cross_sections_Er(filename, lambda)
%READ_CROSS_SECTIONS_ER_Nd This reads the absorption and emission cross 
%sections of ER and Nd

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

%% Interpolating data for our wavelength grid
lambda_raw_min = lambda_raw(1);
lambda_raw_max = lambda_raw(end);
lambda_inside_lambda_raw = lambda(lambda<=lambda_raw_max & lambda>=lambda_raw_min); % interpolation
lambda_outside_lambda_raw_left  = lambda(lambda<lambda_raw_min); % extrapolation
lambda_outside_lambda_raw_right = lambda(lambda>lambda_raw_max); % extrapolation

% interpolation
absorption = [zeros(size(lambda_outside_lambda_raw_left)); ...
                interp1(lambda_raw,absorption_raw, lambda_inside_lambda_raw, 'pchip'); ...
              zeros(size(lambda_outside_lambda_raw_right))];
emission   = [zeros(size(lambda_outside_lambda_raw_left)); ...
                interp1(lambda_raw,emission_raw  , lambda_inside_lambda_raw, 'pchip'); ...
              zeros(size(lambda_outside_lambda_raw_right))];
ESA   = [zeros(size(lambda_outside_lambda_raw_left)); ...
           interp1(lambda_raw,ESA_raw  , lambda_inside_lambda_raw, 'pchip'); ...
         zeros(size(lambda_outside_lambda_raw_right))];

% extrapolation with an exponential decay
% The left part
if ~isempty(lambda_outside_lambda_raw_left)
    % B.C.: Slope at the edge
    absorption_raw_slope_left = (absorption_raw(2)-absorption_raw(1))/(lambda_raw(2)-lambda_raw(1));
    emission_raw_slope_left = (emission_raw(2)-emission_raw(1))/(lambda_raw(2)-lambda_raw(1));
    ESA_raw_slope_left = (ESA_raw(2)-ESA_raw(1))/(lambda_raw(2)-lambda_raw(1));
    % Compute extrapolations
    absorption(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),absorption_raw(1),max(0,absorption_raw_slope_left),lambda_outside_lambda_raw_left);
    emission  (1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),  emission_raw(1),  max(0,emission_raw_slope_left),lambda_outside_lambda_raw_left);
    ESA       (1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),       ESA_raw(1),       max(0,ESA_raw_slope_left),lambda_outside_lambda_raw_left);
end
% The right part
if ~isempty(lambda_outside_lambda_raw_right)
    % B.C.: Slope at the edge
    absorption_raw_slope_right = (absorption_raw(end)-absorption_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    emission_raw_slope_right = (emission_raw(end)-emission_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    ESA_raw_slope_right = (ESA_raw(end)-ESA_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    % Compute extrapolations
    absorption(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),absorption_raw(end),min(0,absorption_raw_slope_right),lambda_outside_lambda_raw_right);
    emission  (end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),  emission_raw(end),  min(0,emission_raw_slope_right),lambda_outside_lambda_raw_right);
    ESA       (end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),       ESA_raw(end),       min(0,ESA_raw_slope_right),lambda_outside_lambda_raw_right);
end

% Reverse back the original order of "lambda"
if lambda_sort
    absorption = absorption(reverse_sort_idx);
    emission = emission(reverse_sort_idx);
    ESA = ESA(reverse_sort_idx);
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