% Note that the final pump power and the inversion, N2, are the same at
% each point of the gain fiber for co- and counter-pumping because the gain
% can't respond to the pulse in time and sees only the average effect
% including both the forward and backward propagating pulses.
close all; clearvars;

filename = 'linear_Mamyshev_oscillator_noASE.mat';

load(filename);

addpath('../../../user_helpers/');

analyze_field(t,f,output_field(:,1,end),'Treacy-t',pi/6,1e-6,true,false);

func = analyze_sim;
func.analyze_fields(t,f,field{end},saved_z,splice_z);

pump_plot.forward = cat(3,zeros(1,1,length(saved_z)/2),pump{end}(1,1,length(saved_z)/2+1:end));
pump_plot.backward = cat(3,pump{end}(1,1,1:length(saved_z)/2),zeros(1,1,length(saved_z)/2));
func.analyze_gain(saved_z,splice_z,pump_plot,N2{end});