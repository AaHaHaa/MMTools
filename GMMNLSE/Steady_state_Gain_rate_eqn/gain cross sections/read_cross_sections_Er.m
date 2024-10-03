function [GSA08,GSA07,GSA06,GSA05,GSA04,GSA03,GSA02,GSA01,...
          ESA15,ESA26,ESA12,ESA23,...
          emi10,emi21,emi32] = read_cross_sections_Er(filename, lambda)
%READ_CROSS_SECTIONS_ER This reads the absorption and emission cross 
%sections of Er
%
% The data comes from
%
% 1. Verlinden et al., "The Excited State Absorption Cross Section of
%    Neodymium-doped Silica Glass Fiber in the 1200-1500 nm Wavelength 
%    Range," Worcester Polytechnic Institute, PhD Thesis (2008)
% 2. Zemon et al., "Excited-state absorption cross sections in the 800-nm
%    band for Er-doped, Al/P-silica fibers: Measurements and amplifier 
%    modeling ," IEEE Photonics Technology Letters 3(7), 621-624 (1991)
% 3. Quimby et al., "Excited-state absorption at 980 nm in erbium-doped
%    glass," Fiber Laser Sources and Amplifiers III 1581, 72-79 (1992)
% 4. Lin et al., "Enhanced mid-infrared emissions of Er^3+ at 2.7 μm via
%    Nd^3+ sensitization in chalcohalide glass," Opt. Lett. 36(10), 1815-
%    1817 (2011)
% 5. Smirnov et al., "Thermal Switching of Lasing Regimes in Heavily Doped
%    Er^3+ Fiber Lasers," ACS Photonics 5(12), 5038-5046 (2018)
% 6. Yamasaki et al., "Optical properties of Er^3+ heavily doped silica
%    glass fabricated by zeolite method," J. Non-Cryst. Solids 543, 120149
%    (2020)
% 7. Pele et al., "Wavelength conversion in Er^3+ doped chalcogenide fibers
%    for optical gas sensors," Opt. Express 23(4), 4163-4172 (2015)

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

formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

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
 GSA08_raw = dataArray{:, 2}; % ground-state absorption (0 --> 8)
 GSA07_raw = dataArray{:, 3}; % ground-state absorption (0 --> 7)
 GSA06_raw = dataArray{:, 4}; % ground-state absorption (0 --> 6)
 GSA05_raw = dataArray{:, 5}; % ground-state absorption (0 --> 5)
 GSA04_raw = dataArray{:, 6}; % ground-state absorption (0 --> 4)
 GSA03_raw = dataArray{:, 7}; % ground-state absorption (0 --> 3)
 GSA02_raw = dataArray{:, 8}; % ground-state absorption (0 --> 2)
 GSA01_raw = dataArray{:, 9}; % ground-state absorption (0 --> 1)
 ESA15_raw = dataArray{:,10}; % excited-state absorption (1 --> 5)
 ESA26_raw = dataArray{:,11}; % excited-state absorption (2 --> 6)
 ESA12_raw = dataArray{:,12}; % excited-state absorption (1 --> 2)
 ESA23_raw = dataArray{:,13}; % excited-state absorption (2 --> 3)
emi10_raw  = dataArray{:,14}; % emission (1 --> 0)
emi21_raw  = dataArray{:,15}; % emission (2 --> 1)
emi32_raw  = dataArray{:,16}; % emission (3 --> 2)

% Smooth the data
[GSA08_raw,GSA07_raw,GSA06_raw,GSA05_raw,GSA04_raw,GSA03_raw,GSA02_raw,GSA01_raw,...
 ESA15_raw,ESA26_raw,ESA12_raw,ESA23_raw,...
 emi10_raw,emi21_raw,emi32_raw] = smooth_cross_sections(lambda_raw,...
                                                        GSA08_raw,GSA07_raw,GSA06_raw,GSA05_raw,GSA04_raw,GSA03_raw,GSA02_raw,GSA01_raw,...
                                                        ESA15_raw,ESA26_raw,ESA12_raw,ESA23_raw,...
                                                        emi10_raw,emi21_raw,emi32_raw);

%% Interpolating data for our wavelength grid
lambda_raw_min = lambda_raw(1);
lambda_raw_max = lambda_raw(end);
lambda_inside_lambda_raw = lambda(lambda<=lambda_raw_max & lambda>=lambda_raw_min); % interpolation
lambda_outside_lambda_raw_left  = lambda(lambda<lambda_raw_min); % extrapolation
lambda_outside_lambda_raw_right = lambda(lambda>lambda_raw_max); % extrapolation

% interpolation
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
GSA02 = [zeros(size(lambda_outside_lambda_raw_left)); ...
           interp1(lambda_raw,GSA02_raw, lambda_inside_lambda_raw, 'pchip'); ...
         zeros(size(lambda_outside_lambda_raw_right))];
GSA01 = [zeros(size(lambda_outside_lambda_raw_left)); ...
           interp1(lambda_raw,GSA01_raw, lambda_inside_lambda_raw, 'pchip'); ...
         zeros(size(lambda_outside_lambda_raw_right))];
ESA15 = [zeros(size(lambda_outside_lambda_raw_left)); ...
           interp1(lambda_raw,ESA15_raw, lambda_inside_lambda_raw, 'pchip'); ...
         zeros(size(lambda_outside_lambda_raw_right))];
ESA26 = [zeros(size(lambda_outside_lambda_raw_left)); ...
           interp1(lambda_raw,ESA26_raw, lambda_inside_lambda_raw, 'pchip'); ...
         zeros(size(lambda_outside_lambda_raw_right))];
ESA12 = [zeros(size(lambda_outside_lambda_raw_left)); ...
           interp1(lambda_raw,ESA12_raw, lambda_inside_lambda_raw, 'pchip'); ...
         zeros(size(lambda_outside_lambda_raw_right))];
ESA23 = [zeros(size(lambda_outside_lambda_raw_left)); ...
           interp1(lambda_raw,ESA23_raw, lambda_inside_lambda_raw, 'pchip'); ...
         zeros(size(lambda_outside_lambda_raw_right))];
emi10 = [zeros(size(lambda_outside_lambda_raw_left)); ...
           interp1(lambda_raw,emi10_raw, lambda_inside_lambda_raw, 'pchip'); ...
         zeros(size(lambda_outside_lambda_raw_right))];
emi21 = [zeros(size(lambda_outside_lambda_raw_left)); ...
           interp1(lambda_raw,emi21_raw, lambda_inside_lambda_raw, 'pchip'); ...
         zeros(size(lambda_outside_lambda_raw_right))];
emi32 = [zeros(size(lambda_outside_lambda_raw_left)); ...
           interp1(lambda_raw,emi32_raw, lambda_inside_lambda_raw, 'pchip'); ...
         zeros(size(lambda_outside_lambda_raw_right))];

% extrapolation with an exponential decay
% The left part
if ~isempty(lambda_outside_lambda_raw_left)
    % B.C.: Slope at the edge
    GSA08_raw_slope_left = (GSA08_raw(2)-GSA08_raw(1))/(lambda_raw(2)-lambda_raw(1));
    GSA07_raw_slope_left = (GSA07_raw(2)-GSA07_raw(1))/(lambda_raw(2)-lambda_raw(1));
    GSA06_raw_slope_left = (GSA06_raw(2)-GSA06_raw(1))/(lambda_raw(2)-lambda_raw(1));
    GSA05_raw_slope_left = (GSA05_raw(2)-GSA05_raw(1))/(lambda_raw(2)-lambda_raw(1));
    GSA04_raw_slope_left = (GSA04_raw(2)-GSA04_raw(1))/(lambda_raw(2)-lambda_raw(1));
    GSA03_raw_slope_left = (GSA03_raw(2)-GSA03_raw(1))/(lambda_raw(2)-lambda_raw(1));
    GSA02_raw_slope_left = (GSA02_raw(2)-GSA02_raw(1))/(lambda_raw(2)-lambda_raw(1));
    GSA01_raw_slope_left = (GSA01_raw(2)-GSA01_raw(1))/(lambda_raw(2)-lambda_raw(1));
    ESA15_raw_slope_left = (ESA15_raw(2)-ESA15_raw(1))/(lambda_raw(2)-lambda_raw(1));
    ESA26_raw_slope_left = (ESA26_raw(2)-ESA26_raw(1))/(lambda_raw(2)-lambda_raw(1));
    ESA12_raw_slope_left = (ESA12_raw(2)-ESA12_raw(1))/(lambda_raw(2)-lambda_raw(1));
    ESA23_raw_slope_left = (ESA23_raw(2)-ESA23_raw(1))/(lambda_raw(2)-lambda_raw(1));
    emi10_raw_slope_left = (emi10_raw(2)-emi10_raw(1))/(lambda_raw(2)-lambda_raw(1));
    emi21_raw_slope_left = (emi21_raw(2)-emi21_raw(1))/(lambda_raw(2)-lambda_raw(1));
    emi32_raw_slope_left = (emi32_raw(2)-emi32_raw(1))/(lambda_raw(2)-lambda_raw(1));
    % Compute extrapolations
    GSA08(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),GSA08_raw(1),max(0,GSA08_raw_slope_left),lambda_outside_lambda_raw_left);
    GSA07(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),GSA07_raw(1),max(0,GSA07_raw_slope_left),lambda_outside_lambda_raw_left);
    GSA06(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),GSA06_raw(1),max(0,GSA06_raw_slope_left),lambda_outside_lambda_raw_left);
    GSA05(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),GSA05_raw(1),max(0,GSA05_raw_slope_left),lambda_outside_lambda_raw_left);
    GSA04(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),GSA04_raw(1),max(0,GSA04_raw_slope_left),lambda_outside_lambda_raw_left);
    GSA03(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),GSA03_raw(1),max(0,GSA03_raw_slope_left),lambda_outside_lambda_raw_left);
    GSA02(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),GSA02_raw(1),max(0,GSA02_raw_slope_left),lambda_outside_lambda_raw_left);
    GSA01(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),GSA01_raw(1),max(0,GSA01_raw_slope_left),lambda_outside_lambda_raw_left);
    ESA15(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),ESA15_raw(1),max(0,ESA15_raw_slope_left),lambda_outside_lambda_raw_left);
    ESA26(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),ESA26_raw(1),max(0,ESA26_raw_slope_left),lambda_outside_lambda_raw_left);
    ESA12(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),ESA12_raw(1),max(0,ESA12_raw_slope_left),lambda_outside_lambda_raw_left);
    ESA23(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),ESA23_raw(1),max(0,ESA23_raw_slope_left),lambda_outside_lambda_raw_left);
    emi10(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),emi10_raw(1),max(0,emi10_raw_slope_left),lambda_outside_lambda_raw_left);
    emi21(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),emi21_raw(1),max(0,emi21_raw_slope_left),lambda_outside_lambda_raw_left);
    emi32(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),emi32_raw(1),max(0,emi32_raw_slope_left),lambda_outside_lambda_raw_left);
end
% The right part
if ~isempty(lambda_outside_lambda_raw_right)
    % B.C.: Slope at the edge
    GSA08_raw_slope_right = (GSA08_raw(end)-GSA08_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    GSA07_raw_slope_right = (GSA07_raw(end)-GSA07_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    GSA06_raw_slope_right = (GSA06_raw(end)-GSA06_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    GSA05_raw_slope_right = (GSA05_raw(end)-GSA05_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    GSA04_raw_slope_right = (GSA04_raw(end)-GSA04_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    GSA03_raw_slope_right = (GSA03_raw(end)-GSA03_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    GSA02_raw_slope_right = (GSA02_raw(end)-GSA02_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    GSA01_raw_slope_right = (GSA01_raw(end)-GSA01_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    ESA15_raw_slope_right = (ESA15_raw(end)-ESA15_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    ESA26_raw_slope_right = (ESA26_raw(end)-ESA26_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    ESA12_raw_slope_right = (ESA12_raw(end)-ESA12_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    ESA23_raw_slope_right = (ESA23_raw(end)-ESA23_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    emi10_raw_slope_right = (emi10_raw(end)-emi10_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    emi21_raw_slope_right = (emi21_raw(end)-emi21_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    emi32_raw_slope_right = (emi32_raw(end)-emi32_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    % Compute extrapolations
    GSA08(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),GSA08_raw(end),min(0,GSA08_raw_slope_right),lambda_outside_lambda_raw_right);
    GSA07(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),GSA07_raw(end),min(0,GSA07_raw_slope_right),lambda_outside_lambda_raw_right);
    GSA06(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),GSA06_raw(end),min(0,GSA06_raw_slope_right),lambda_outside_lambda_raw_right);
    GSA05(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),GSA05_raw(end),min(0,GSA05_raw_slope_right),lambda_outside_lambda_raw_right);
    GSA04(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),GSA04_raw(end),min(0,GSA04_raw_slope_right),lambda_outside_lambda_raw_right);
    GSA03(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),GSA03_raw(end),min(0,GSA03_raw_slope_right),lambda_outside_lambda_raw_right);
    GSA02(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),GSA02_raw(end),min(0,GSA02_raw_slope_right),lambda_outside_lambda_raw_right);
    GSA01(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),GSA01_raw(end),min(0,GSA01_raw_slope_right),lambda_outside_lambda_raw_right);
    ESA15(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),ESA15_raw(end),min(0,ESA15_raw_slope_right),lambda_outside_lambda_raw_right);
    ESA26(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),ESA26_raw(end),min(0,ESA26_raw_slope_right),lambda_outside_lambda_raw_right);
    ESA12(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),ESA12_raw(end),min(0,ESA12_raw_slope_right),lambda_outside_lambda_raw_right);
    ESA23(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),ESA23_raw(end),min(0,ESA23_raw_slope_right),lambda_outside_lambda_raw_right);
    emi10(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),emi10_raw(end),min(0,emi10_raw_slope_right),lambda_outside_lambda_raw_right);
    emi21(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),emi21_raw(end),min(0,emi21_raw_slope_right),lambda_outside_lambda_raw_right);
    emi32(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),emi32_raw(end),min(0,emi32_raw_slope_right),lambda_outside_lambda_raw_right);
end

% Reverse back the original order of "lambda"
if lambda_sort
    GSA08 = GSA08(reverse_sort_idx);
    GSA07 = GSA07(reverse_sort_idx);
    GSA06 = GSA06(reverse_sort_idx);
    GSA05 = GSA05(reverse_sort_idx);
    GSA04 = GSA04(reverse_sort_idx);
    GSA03 = GSA03(reverse_sort_idx);
    GSA02 = GSA02(reverse_sort_idx);
    GSA01 = GSA01(reverse_sort_idx);
    ESA15 = ESA15(reverse_sort_idx);
    ESA26 = ESA26(reverse_sort_idx);
    ESA12 = ESA12(reverse_sort_idx);
    ESA23 = ESA23(reverse_sort_idx);
    emi10 = emi10(reverse_sort_idx);
    emi21 = emi21(reverse_sort_idx);
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