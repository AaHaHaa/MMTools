function [data_abs,data_emi,varargout] = smooth_cross_sections(data_wl,data_abs,data_emi,varargin)
%SMOOTH_CROSS_SECTIONS It smoothens the cross sections by smoothening their
%derivatives and integrating them back
%
%   data_wl: the measured wavelength
%   data_abs: the measured absorption cross section
%   data-emi: the measured emission cross section

%% Smooth

[data_abs,data_emi] = smooth_nk(data_wl,data_abs,data_emi,5,2);

% remove negative values
data_abs(data_abs<0) = 0;
data_emi(data_emi<0) = 0;

% If there are more input arguments
varargout = cell(1,length(varargin));
for i = 1:length(varargin)
    varargout{i} = smooth_varargin(data_wl,varargin{i},5,2);
    
    % remove negative values
    varargout{i}(varargout{i}<0) = 0;
end

end

%% Helper functions
function [abs,emi] = smooth_nk(wl,abs,emi,smooth_range,smooth_rep)

for i = 1:smooth_rep
    abs = smooth(wl,abs,smooth_range,'sgolay',3);
    emi = smooth(wl,emi,smooth_range,'sgolay',3);
end

end

function output = smooth_varargin(wl,more_argin,smooth_range,smooth_rep)

for i = 1:smooth_rep
    output = smooth(wl,more_argin,smooth_range,'sgolay',3);
end

end