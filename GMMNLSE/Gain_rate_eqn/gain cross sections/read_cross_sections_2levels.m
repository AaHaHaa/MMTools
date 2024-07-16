%%%% Function to read Yb cross sections from a file %%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Author: Michel Olivier                                                %
%   Affiliation: COPL, Universit� Laval, Qu�bec, Canada                   %
%                                                                         %
%   Date of creation: 2017-07-05                                          %
%   Date of last modification: 2017-11-21                                 %
%                                                                         %
%   Rem.:This function reads the Yb absorption and emission cross         %
%   sections from a file and it interpolates to define the cross          %
%   sections according to the simulation frequency grid.                  %
%                                                                         %
%   Output:                                                               %
%   absorption: Yb absorption cross section array (float)                 %
%   emission: Yb emission cross section array (float)                     %
%                                                                         %
%   Input:                                                                %
%   filename: name of the file to be read (string)                        %
%   lambda: array of wavelengths used in simulation (float)               %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Yb's Nufern data source:
% https://www.coherent.com/resources/application-note/components-and-accessories/specialty-optical-fibers/yb-absorption-emission.pdf
%
% Modified by Yi-Hao Chen
% date       | modification
%  6/28/2018 | Added exponential-decay extrapolation instead of just making them zero.
% 10/27/2018 | Change the function name from "read_Yb_cross_sections" to "read_cross_sections"
%              because I added Er cross section data
% 6/27/2021  | Smooth the measured data due to weird sharp change of derivative at around 1050 nm for the Yb absorption cross section in "Liekki Yb_AV_20160530.txt".
% 5/20/2024  | Restrict this function only to Yb and Nd. 
%              A separate function is created to read Er data due to the cross section of its excited-state absorption.

function [absorption,emission] = read_cross_sections_2levels(filename, lambda)
%READ_CROSS_SECTIONS_2LEVELS This reads the absorption and emission cross 
%sections of a two-level system, such as Yb.

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

% extrapolation with an exponential decay
% The left part
if ~isempty(lambda_outside_lambda_raw_left)
    % B.C.: Slope at the edge
    absorption_raw_slope_left = (absorption_raw(2)-absorption_raw(1))/(lambda_raw(2)-lambda_raw(1));
    emission_raw_slope_left = (emission_raw(2)-emission_raw(1))/(lambda_raw(2)-lambda_raw(1));
    % Compute extrapolations
    absorption(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),absorption_raw(1),max(0,absorption_raw_slope_left),lambda_outside_lambda_raw_left);
    emission  (1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),  emission_raw(1),  max(0,emission_raw_slope_left),lambda_outside_lambda_raw_left);
end
% The right part
if ~isempty(lambda_outside_lambda_raw_right)
    % B.C.: Slope at the edge
    absorption_raw_slope_right = (absorption_raw(end)-absorption_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    emission_raw_slope_right = (emission_raw(end)-emission_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    % Compute extrapolations
    absorption(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),absorption_raw(end),min(0,absorption_raw_slope_right),lambda_outside_lambda_raw_right);
    emission  (end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),  emission_raw(end),  min(0,emission_raw_slope_right),lambda_outside_lambda_raw_right);
end

% Reverse back the original order of "lambda"
if lambda_sort
    absorption = absorption(reverse_sort_idx);
    emission = emission(reverse_sort_idx);
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