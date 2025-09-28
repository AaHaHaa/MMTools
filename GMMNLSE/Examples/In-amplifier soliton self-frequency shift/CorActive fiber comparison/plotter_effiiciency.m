clearvars; close all;

pump_power = 1:10; % W

all_output_power = zeros(length(pump_power),1);
for i = 1:length(pump_power)
    filename = sprintf('Tm_CW_%uW.mat',pump_power(i));
    load(filename);
    
    all_output_power(i) = output_power;
end

figure;
data = readmatrix('CorActive data.csv');
h = plot(data(:,1),data(:,2),'LineWidth',2,'Color','b');
hold on;
h2 = plot(pump_power,all_output_power,'LineWidth',2,'Color','r');
hold off;
set(gca,'fontsize',20);
xlabel('Injected pump power (W)');
ylabel('Signal power (W)');
legend('data','sim','Location','southeast');