function varargout = confplot(varargin)
%CONFPLOT Linear plot with continuous confidence/error boundaries.
%   (Modified by Yi-Hao Chen, 4/13/2022)
%
%   CONFPLOT(X,Y,L,U) plots the graph of vector X vs. vector Y with
%   'continuous' confidence/error boundaries specified by the vectors
%   L and U.  L and U contain the lower and upper error ranges for each
%   point in Y. The vectors X,Y,L and U must all be the same length.  
%
%   CONFPLOT(X,Y,E) or CONFPLOT(Y,E) plots Y with error bars [Y-E Y+E].
%   CONFPLOT(...,'LineSpec') uses the color and linestyle specified by
%   the string 'LineSpec'.  See PLOT for possibilities.
%
%   H = CONFPLOT(...) returns a vector of line handles.
%
%   For example,
%      x = 1:0.1:10;
%      y = sin(x);
%      e = std(y)*ones(size(x));
%      confplot(x,y,e)
%   draws symmetric continuous confidence/error boundaries of unit standard deviation.
%
%   See also ERRORBAR, SEMILOGX, SEMILOGY, LOGLOG, PLOTYY, GRID, CLF, CLC, TITLE,
%   XLABEL, YLABEL, AXIS, AXES, HOLD, COLORDEF, LEGEND, SUBPLOT, STEM.
%
%     © 2002 - Michele Giugliano, PhD (http://www.giugliano.info) (Bern, Monday Nov 4th, 2002 - 19:02)
%    (bug-reports to michele@giugliano.info)
%   $Revision: 1.0 $  $Date: 2002/11/11 14:36:08 $
%
% -------------------------------------------------------------------------
% Yi-Hao's modification:
%    Instead of plotting two areas with gray and white to generate the
%    shaded region, I use one area(x,y) with y being a matrix, [z2,z1-z2].
%    Michele's original version has the problem of covering those in the
%    white region. This becomes a problem if confplot() is overlapped on
%    top of an nonempty figure, such as after "hold on" or within yyaxis
%    left or right.
%    The matrix version of MATLAB's area() has the capability of setting
%    the white region as transparency but maintaining the gray color within
%    only the narrow standard-deviation region.

spec_argin_idx = 0;
for i = 1:nargin
    if ischar(varargin{i}) || isstring(varargin{i})
        spec_argin_idx = i;
        nargin_vector = spec_argin_idx - 1;
        break;
    end
end

if nargin_vector < 2
	disp('ERROR: not enough input arguments!');
	return;
end

x = [];  y = [];  z1 = [];  z2 = [];

switch nargin_vector
	case 2
        y  = varargin{1};
        z1 = y + varargin{2};
        z2 = y - varargin{2};
        x  = 1:length(y);
	case 3
        x  = varargin{1};
        y  = varargin{2};
        z1 = y + varargin{3};
        z2 = y - varargin{3};
	case 4
        x  = varargin{1};
        y  = varargin{2};
        z1 = y + varargin{4};
        z2 = y - varargin{3};
end
% column vector for all inputs
sx = size(x); if sx(2) > sx(1), x = x.'; end
sy = size(y); if sy(2) > sy(1), y = y.'; end
sz1 = size(z1); if sz1(2) > sz1(1), z1 = z1.'; end
sz2 = size(z2); if sz2(2) > sz2(1), z2 = z2.'; end

% draw the shaded area
p = plot(x,y,x,z1,x,z2);    YLIM = get(gca,'YLim');    delete(p); % get ylim; otherwise, area() starts from 0 in y values
axes_area = area(x,[z2,z1-z2]);
set(axes_area,'LineStyle','none','ShowBaseLine','off');
%set(axes_area(2),'FaceColor',[0.9 0.9 0.9]); % preset to gray color
set(axes_area(1),'FaceColor','none');
% draw the line of mean, y
hold on;
if spec_argin_idx ~= 0
    spec = sprintf('p = plot(x,y,varargin{spec_argin_idx}');
    for i = spec_argin_idx+1:nargin
        spec = sprintf('%s,varargin{%d}',spec,i);
    end
    spec = sprintf('%s);',spec);
    eval(spec);
else
    p = plot(x,y);
end
hold off;
LineColor = get(p,'Color'); LineColor(LineColor==0) = 0.01;
set(axes_area(2),'FaceColor',(LineColor).^(1/10)); % set the area color to the light-dark of the line color

%set(gca,'Layer','top','XGrid','on','YGrid','on');             
set(gca,'Layer','top');
ylim(YLIM);

% Output
H = {p, axes_area};
varargout = H;