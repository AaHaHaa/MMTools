function Aeff = get_Aeff( SR )
%GET_AEFF Gets Aeff for each mode from SR values
%
%   SR: (num_modes,num_modes,num_modes,num_modes,num_fiber)
%   Aeff: (1,num_modes,num_fiber)

num_fiber = size(SR,5);
num_modes = size(SR,1);

Aeff = zeros(1,num_modes,num_fiber);
for i = 1:num_fiber
    for m = 1:num_modes
        Aeff(1,m,i) = 1/SR(m,m,m,m,i);
    end
end

end

