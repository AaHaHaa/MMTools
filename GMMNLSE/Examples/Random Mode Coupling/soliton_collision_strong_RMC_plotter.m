% It plots the field for each simulation. Besides, it also calculates the
% total energy to check if it's conserved. Because I found out that if the
% running precision for linear mode coupling isn't set to "double", the
% energy may not be conserved.

close all;  clearvars;

addpath('../../user_helpers/');

filename = {'soliton_collision_noRMC.mat',...
            'soliton_collision_strong_RMC_1.mat','soliton_collision_strong_RMC_2.mat'};

propagation_length = 3790e3; % approximately 80 soliton periods, numbers taken from "soliton_collision" code
num_save = 100;

for i = 1:length(filename)

    load(filename{i});

    for n = 1:length(output_field)
        %{
        % Total energy
        figure;
        energy = permute(sum(trapz(abs(output_field(n).fields).^2),2),[3 1 2])*output_field(n).dt/10^3; % energy in nJ
        plot(energy);
        ax = gca;
        ylim([0.9 1.1]*max(ax.YLim));
        xlabel('Propagation distance (km)');
        ylabel('Energy (nJ)');
        title('Check the energy conservation');
        %}
        % Propagation
        figure;
        abs_field = zeros(size(output_field(n).fields,1),size(output_field(n).fields,3));
        linecolor = distinguishable_colors(2);
        for zidx = 1:size(output_field(n).fields,3)
            for k = 1:size(output_field(n).fields,2)
                abs_field(:,zidx) = abs_field(:,zidx) + abs(output_field(n).fields(:,k,zidx)).^2;
            end
            %plot3(zidx*ones(size(abs_field,1),1),t,abs_field(:,zidx));
            %plot3(j*ones(size(abs_field,1),1),t,abs(output_field(n).fields(:,1,j)),'color','r');
            %hold on; axis tight;
            %plot3(j*ones(size(abs_field,1),1),t,abs(output_field(n).fields(:,3,j)),'color','b');
        end
        %hold off;

        pcolor(t,linspace(0,propagation_length,num_save+1),abs_field.'); colorbar; grid off; shading interp;
        xlabel('time (ps)');
        ylabel('Propagation distance (km)');
        %title(sprintf('Pulse Power (W)\n(%s%s)',title_string_model,title_string_RMC));

        % Save the figure
        sep = strfind(filename{i},'.');
        %print(gcf,filename{i}(1:sep-1), '-depsc');
    end
end