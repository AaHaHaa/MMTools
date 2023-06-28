% ----------------------------------------------------------------------- %
% FUNCTION "whitejet": defines a new colormap with the same colors        %
% that "jet", but it also replaces the green tones with white ones. This  %
% useful when a signed metric is depicted, and its null values are useless.
% The color structure is the following:                                   %
%                                                                         %
%           DR  R       Y       W       C       B   DB                    %
%           |---|-------|-------|?------|?------|---|                     %
%           0  0.1     0.3     0.5     0.7     0.9  1                     %
% where:                                                                  %
%       - DR:   Deep Red    (RGB: 0.5 0 0)                                %
%       - R:    Red         (RGB: 1 0 0)                                  %
%       - Y:    Yellow      (RGB: 1 1 0)                                  %
%       - W:    White       (RGB: 1 1 1)                                  %
%       - C:    Cyan        (RGB: 0 1 1)                                  %
%       - B:    Blue        (RGB: 0 0 1)                                  %
%       - DB:   Deep Blue   (RGB: 0 0 0.5)                                %
%                                                                         %
%   Input parameters:                                                     %
%       - m:    Number of points (recommended: m > 64, min value: m = 7). %
%                                                                         %
%   Output variables:                                                     %
%       - J:    Colormap in RGB values (dimensions [mx3]).                %
% ----------------------------------------------------------------------- %
%   Example of use:                                                       %
%       C = 2.*rand(5,100)-1;                                             %
%       imagesc(C);                                                       %
%       colormap(whitejet);                                               %
%       colorbar;                                                         %
% ----------------------------------------------------------------------- %
%       - Author:   Víctor Martínez-Cagigal                               %
%       - Date:     16/05/2018                                            %
%       - Version:  1.0                                                   %
%       - E-mail:   victor.martinez (at) gib (dot) tel (dot) uva (dot) es %
%                                                                         %
%       Biomedical Engineering Group (University of Valladolid), Spain    %
% ----------------------------------------------------------------------- %
function J = whitejet_centered(m)

if nargin < 1
   f = get(groot,'CurrentFigure');
   if isempty(f)
      m = size(get(groot,'DefaultFigureColormap'),1);
   else
      m = size(f.Colormap,1);
   end
end

% Colors
color_palette = [1/2 0 0;   % Deep red
                 1 0 0;     % Red
                 1 1 0;     % Yellow
                 1 1 1;     % White
                 0 1 1;     % Cyan
                 0 0 1;     % Blue
                 0 0 1/2];  % Deep blue
             
% Compute distributions along the samples
color_dist = cumsum([0 1/10 1/5 1/5 1/5 1/5 1/10]);
color_samples = round((m-1)*color_dist)+1;

% Make the gradients
J = zeros(m,3);
J(color_samples,:) = color_palette(1:7,:);
diff_samples = diff(color_samples)-1;
for d = 1:1:length(diff_samples)
    if diff_samples(d)~=0
        color1 = color_palette(d,:);
        color2 = color_palette(d+1,:);
        G = zeros(diff_samples(d),3);
        for idx_rgb = 1:1:3
            g = linspace(color1(idx_rgb), color2(idx_rgb), diff_samples(d)+2);
            g([1, length(g)]) = [];
            G(:,idx_rgb) = g';
        end
        J(color_samples(d)+1:color_samples(d+1)-1,:) = G;
    end
end
J = flipud(J);




