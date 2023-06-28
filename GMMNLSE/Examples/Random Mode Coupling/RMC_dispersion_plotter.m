close all;  clearvars;

addpath('../../user_helpers/');

load('RMC.mat');
    
each_energy = permute(trapz(abs(output_field.fields).^2)*output_field.dt/1e3,[3,2,1]);
each_ratio = each_energy./sum(each_energy,2);

figure;
plot(output_field.z,each_ratio);
xlabel('Propagation length (m)');
ylabel('Energy ratio');