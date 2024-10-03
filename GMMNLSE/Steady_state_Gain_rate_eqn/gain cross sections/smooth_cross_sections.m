function varargout = smooth_cross_sections(data_wl,varargin)
%SMOOTH_CROSS_SECTIONS It smoothens the cross sections
%
%   data_wl: the measured wavelength
%   varargout: container cell for the measured cross sections

%% Smooth

varargout = cell(1,length(varargin));
for i = 1:length(varargin)
    varargout{i} = smooth_varargin(data_wl,varargin{i},5,1);
    
    % remove negative values
    varargout{i}(varargout{i}<0) = 0;
end

end

%% Helper function
function output = smooth_varargin(wl,cross_sections,smooth_range,smooth_rep)

for i = 1:smooth_rep
    output = smooth(wl,cross_sections,smooth_range,'sgolay',3);
end

end