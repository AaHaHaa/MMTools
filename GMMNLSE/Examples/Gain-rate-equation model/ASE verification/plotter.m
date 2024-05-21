clearvars; close all;

%% Load paper data
data = readmatrix('ASE data.csv');
exp_pump_power = data(:,1);
exp_backward_ASE_power = data(:,2);
exp_forward_ASE_power = data(:,3);

%%
pump_power = [20,30,60,80,100,120];

fitting_ratio = 1.58; % after trial and error

forward_ASE_power = zeros(length(pump_power),1);
backward_ASE_power = zeros(length(pump_power),1);
for i = 1:length(pump_power)
    load(sprintf('ASE_evolutions_Pump%3uW.mat',pump_power(i)));
    
    f = fftshift(f);
    
    forward_ASE_power(i) = trapz(f,prop_output.Power.ASE.forward(:,:,end));
    backward_ASE_power(i) = trapz(f,prop_output.Power.ASE.backward(:,:,1));
end

figure;
h = plot(pump_power*fitting_ratio,[forward_ASE_power,backward_ASE_power],'linewidth',2);
set(gca,'fontsize',20);
set(h(1),'Color','b'); set(h(2),'Color','r');
xlabel('Pump power (W)');
ylabel('ASE power (W)');

hold on;
h = plot(exp_pump_power,[exp_forward_ASE_power,exp_backward_ASE_power],'linewidth',2,'linestyle','--');
set(h(1),'Color','b'); set(h(2),'Color','r');
hold off;