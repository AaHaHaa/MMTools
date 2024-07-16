% Note that the final pump power and the inversion, N2, are the same at
% each point of the gain fiber for co- and counter-pumping because the gain
% can't respond to the pulse in time and sees only the average effect
% including both the forward and backward propagating pulses.
close all; clearvars;

filename = 'linear_Mamyshev_oscillator_withASE.mat';

load(filename);

addpath('../../../user_helpers/');

ASE_out = struct('spectrum',ASE.forward{end},'t_rep',gain_rate_eqn{1}.t_rep);
analyze_field(t,f,output_field(:,1,end),'Treacy-t',pi/6,1e-6,true,false,ASE_out);

func = analyze_sim;
func.analyze_fields(t,f,field{end},saved_z,splice_z);

pump_plot.forward = cat(3,zeros(1,1,length(saved_z)/2),pump{end}(1,1,length(saved_z)/2+1:end));
pump_plot.backward = cat(3,pump{end}(1,1,1:length(saved_z)/2),zeros(1,1,length(saved_z)/2));
func.analyze_gain(saved_z,splice_z,pump_plot,1-N2{end});

ASE_plot.forward = ASE.forward{end}; ASE_plot.backward = ASE.backward{end};
func.analyze_ASE(f,ASE_plot,saved_z,splice_z);