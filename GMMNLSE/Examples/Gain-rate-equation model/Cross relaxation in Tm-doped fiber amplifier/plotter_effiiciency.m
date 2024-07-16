clearvars; close all;

pump_power = 10:5:30; % W

all_absorbed_pump_power = zeros(length(pump_power),1);
all_output_power = zeros(length(pump_power),1);
for i = 1:length(pump_power)
    filename = sprintf('Tm_CW_%uW.mat',pump_power(i));
    load(filename);
    
    all_absorbed_pump_power(i) = absorbed_pump_power;
    all_output_power(i) = output_power;
end

[ linefitcoeff,corr ] = line_fitting( all_absorbed_pump_power,all_output_power );