function [data_abs,data_emi] = smooth_cross_sections(data_wl,data_abs,data_emi)
%SMOOTH_CROSS_SECTIONS It smoothens the cross sections by smoothening their
%derivatives and integrating them back
%
%   data_wl: the measured wavelength
%   data_abs: the measured absorption cross section
%   data-emi: the measured emission cross section

%% Smooth

[Ddata_abs,Ddata_emi] = smooth_nk(data_wl(2:end),diff(data_abs)./diff(data_wl),diff(data_emi)./diff(data_wl),5,2);
data_abs = cumtrapz(data_wl,[0;Ddata_abs]) + data_abs(1);
data_emi = cumtrapz(data_wl,[0;Ddata_emi]) + data_emi(1);

% remove negative values
data_abs(data_abs<0) = 0;
data_emi(data_emi<0) = 0;

end

function [abs,emi] = smooth_nk(wl,abs,emi,smooth_range,smooth_rep)

for i = 1:smooth_rep
    abs = smooth(wl,abs,smooth_range,'sgolay',3);
    emi = smooth(wl,emi,smooth_range,'sgolay',3);
end

end