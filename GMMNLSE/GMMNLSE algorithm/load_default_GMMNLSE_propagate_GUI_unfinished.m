function [fiber,sim] = load_default_GMMNLSE_propagate_GUI()
%LOAD_DEFAULT_GMMNLSE_PROPAGATE_GUI Summary of this function goes here
%

% Load GUI global variables
app = GUI_variables();

% Create UIFigure and components
createComponents();

% Execute the startup function
startupFcn();

% Wait until people finish/close it.
uiwait(app.UIFigure);

    %% Wrap up variables into fiber and sim
    function wrapUp()
        if app.singlemodeButton.Value % scalar computation
            app.betas = app.betas(1:end,1);
        end
        
        % fiber
        fiber = struct('betas',app.betas,...
                       'n2',app.Nonlinearrefractiveindexn2m2WEditField.Value,...
                       'SR',app.SR,...
                       'L0',app.FiberlengthmEditField.Value);
            % Gain model
            if ismember(app.GainmodelDropDown.Value,1:3)
                fiber.dB_gain = app.GaindBmEditField.Value;
                fiber.gain_coeff = app.gaincoeff1mEditField.Value;
                fiber.gain_fwhm = app.FWHMnmEditField.Value*1e-9; % m
                switch app.GainmodelDropDown.Value
                    case 1
                        fiber.saturation_energy = app.SaturationenergynJEditField.Value;
                    case {2,3}
                        fiber.saturation_intensity = app.SaturationintensityJm2EditField.Value;
                end
            end
        
        % sim
        sim = struct('betas',app.sim_betas,...
                     'f0',app.CenterfrequencyTHzEditField.Value,...
                     'sw',app.IncludeshocktermButton.Value,...
                     'save_period',app.SavingstepsizemEditField.Value,...
                     'scalar',app.scalarButton.Value,...
                     'ellipticity',app.EllipticityEditField.Value,...
                     'single_yes',app.singleButton.Value,...
                     'gpu_yes',app.UseGPUButton.Value,...
                     'Raman_model',app.RamanmodelDropDown.Value,...
                     'step_method',app.SteppingalgorithmDropDown.Value,...
                     'gain_model',app.GainmodelDropDown.Value,...
                     'rmc',struct('model',app.RmcmodelDropDown.Value),...
                     'pulse_centering',app.CenterthepulsewrtthetimewindowButton.Value,...
                     'num_photon_noise_per_band',app.PhotonnoiseperfreqbandEditField.Value,...
                     'progress_bar',app.ShowprogressbarButton.Value);
            % MPA
            if isequal(app.SteppingalgorithmDropDown.Value,'MPA')
                sim.MPA = struct('M',app.NumberofparallelizationsEditField.Value,...
                                 'n_tot_max',app.MaxiterationsEditField.Value,...
                                 'n_tot_min',app.MiniterationsEditField.Value,...
                                 'tol',app.ToleranceEditField.Value);
            end
            % Single-mode/Multimode
            switch app.SimulationtypeButtonGroup.SelectedObject.Text
                case 'single mode'
                    sim.midx = 1;
                case 'multimode'
                    sim.midx = unique(uint8(eval(['[',app.modeindexEditField.Value,']'])));
            end
            % Adaptive algorithm
            sim.adaptive_deltaZ.model = app.adaptivestepsizeButton.Value;
            if app.adaptivestepsizeButton.Value
                sim.adaptive_deltaZ.threshold = app.ThresholdEditField.Value;
            else % fixed step size
                sim.deltaZ = app.StepsizemEditField.Value;
            end
            % Raman model
            if ismember(app.RamanmodelDropDown.Value,1:2)
                sim.Raman_sponRS = app.IncludespontaneousRamanButton.Value;
            end
            % GPU
            if app.UseGPUButton.Value
                sim.gpuDevice = struct('Index',app.GPUindextouseEditField.Value,...
                                       'Device',gpuDevice(app.GPUindextouseEditField.Value));
                sim.cuda_dir_path = fullfile(app.upper_folder,app.CudafolderEditField.Value);
            end
            % Progress bar
            if app.ShowprogressbarButton.Value
                sim.progress_bar_name = app.ProgressbarnameEditField.Value;
            end
    end

    %% Create a callback creation function used in MATLAB app designer
    function callback_handle = createCallbackFcn(app,function_handle,dummy) %#ok
        callback_handle = function_handle;
    end

    %%
    % =========================================================================
    % Below are copied from load_GMMNLSE.mlapp
    % =========================================================================

    %% GUI global variables
    function app = GUI_variables()
        app.c = 299792458; % m/s; speed of light
        app.upper_folder = [];
        app.betas = []; % Taylor-expansion coefficients of the propagation constant
        app.Aeff = []; % m^2; effective area of the SMF
        app.SR = []; % S tensors

        app.MM_folder = '';
        app.midx = 1; % mode index

        app.sim_betas = [];
    end

    %% Create UIFigure and components
    function createComponents()

        % Create UIFigure and hide until all components are created
        app.UIFigure = uifigure('Visible', 'off');
        app.UIFigure.Position = [100 100 751 636];
        app.UIFigure.Name = 'MATLAB App';
        app.UIFigure.CloseRequestFcn = @CloseReqFcn;

        % Create GMMNLSE_propagateinputparametersLabel
        app.GMMNLSE_propagateinputparametersLabel = uilabel(app.UIFigure);
        app.GMMNLSE_propagateinputparametersLabel.FontSize = 18;
        app.GMMNLSE_propagateinputparametersLabel.FontWeight = 'bold';
        app.GMMNLSE_propagateinputparametersLabel.Position = [118 600 364 23];
        app.GMMNLSE_propagateinputparametersLabel.Text = 'GMMNLSE_propagate() input parameters';

        % Create CenterfrequencyTHzEditFieldLabel
        app.CenterfrequencyTHzEditFieldLabel = uilabel(app.UIFigure);
        app.CenterfrequencyTHzEditFieldLabel.HorizontalAlignment = 'right';
        app.CenterfrequencyTHzEditFieldLabel.Tooltip = {''};
        app.CenterfrequencyTHzEditFieldLabel.Position = [13 567 131 22];
        app.CenterfrequencyTHzEditFieldLabel.Text = 'Center frequency (THz)';

        % Create CenterfrequencyTHzEditField
        app.CenterfrequencyTHzEditField = uieditfield(app.UIFigure, 'numeric');
        app.CenterfrequencyTHzEditField.Limits = [0 Inf];
        app.CenterfrequencyTHzEditField.ValueChangedFcn = createCallbackFcn(app, @CenterfrequencyTHzEditFieldValueChanged, true);
        app.CenterfrequencyTHzEditField.Tooltip = {'sim.f0'};
        app.CenterfrequencyTHzEditField.Position = [152 567 50 22];
        app.CenterfrequencyTHzEditField.Value = 291.06;

        % Create CenterwavelengthnmEditFieldLabel
        app.CenterwavelengthnmEditFieldLabel = uilabel(app.UIFigure);
        app.CenterwavelengthnmEditFieldLabel.HorizontalAlignment = 'right';
        app.CenterwavelengthnmEditFieldLabel.Tooltip = {''};
        app.CenterwavelengthnmEditFieldLabel.Position = [210 567 134 22];
        app.CenterwavelengthnmEditFieldLabel.Text = 'Center wavelength (nm)';

        % Create CenterwavelengthnmEditField
        app.CenterwavelengthnmEditField = uieditfield(app.UIFigure, 'numeric');
        app.CenterwavelengthnmEditField.ValueChangedFcn = createCallbackFcn(app, @CenterwavelengthnmEditFieldValueChanged, true);
        app.CenterwavelengthnmEditField.Tooltip = {'sim.lambda0'};
        app.CenterwavelengthnmEditField.Position = [354 567 55 22];
        app.CenterwavelengthnmEditField.Value = 1030;

        % Create TabGroup2
        app.TabGroup2 = uitabgroup(app.UIFigure);
        app.TabGroup2.Position = [13 6 489 240];

        % Create GainmodelTab
        app.GainmodelTab = uitab(app.TabGroup2);
        app.GainmodelTab.Title = 'Gain model';

        % Create GaindBmEditFieldLabel
        app.GaindBmEditFieldLabel = uilabel(app.GainmodelTab);
        app.GaindBmEditFieldLabel.HorizontalAlignment = 'right';
        app.GaindBmEditFieldLabel.Enable = 'off';
        app.GaindBmEditFieldLabel.Position = [11 176 71 22];
        app.GaindBmEditFieldLabel.Text = 'Gain (dB/m)';

        % Create GaindBmEditField
        app.GaindBmEditField = uieditfield(app.GainmodelTab, 'numeric');
        app.GaindBmEditField.ValueChangedFcn = createCallbackFcn(app, @GaindBmEditFieldValueChanged, true);
        app.GaindBmEditField.Enable = 'off';
        app.GaindBmEditField.Tooltip = {'fiber.dB_gain'};
        app.GaindBmEditField.Position = [94 176 32 22];
        app.GaindBmEditField.Value = 30;

        % Create gaincoeff1mEditFieldLabel
        app.gaincoeff1mEditFieldLabel = uilabel(app.GainmodelTab);
        app.gaincoeff1mEditFieldLabel.HorizontalAlignment = 'right';
        app.gaincoeff1mEditFieldLabel.Enable = 'off';
        app.gaincoeff1mEditFieldLabel.Position = [140 176 88 22];
        app.gaincoeff1mEditFieldLabel.Text = 'gain coeff (1/m)';

        % Create gaincoeff1mEditField
        app.gaincoeff1mEditField = uieditfield(app.GainmodelTab, 'numeric');
        app.gaincoeff1mEditField.Limits = [0 Inf];
        app.gaincoeff1mEditField.ValueChangedFcn = createCallbackFcn(app, @gaincoeff1mEditFieldValueChanged, true);
        app.gaincoeff1mEditField.Enable = 'off';
        app.gaincoeff1mEditField.Tooltip = {'fiber.gain_coeff=fiber.dB_gain*log(10)/(10*fiber.L0)'};
        app.gaincoeff1mEditField.Position = [241 176 61 22];
        app.gaincoeff1mEditField.Value = 6.90775527898214;

        % Create FWHMnmEditFieldLabel
        app.FWHMnmEditFieldLabel = uilabel(app.GainmodelTab);
        app.FWHMnmEditFieldLabel.HorizontalAlignment = 'right';
        app.FWHMnmEditFieldLabel.Enable = 'off';
        app.FWHMnmEditFieldLabel.Position = [12 144 71 22];
        app.FWHMnmEditFieldLabel.Text = 'FWHM (nm)';

        % Create FWHMnmEditField
        app.FWHMnmEditField = uieditfield(app.GainmodelTab, 'numeric');
        app.FWHMnmEditField.Limits = [0 Inf];
        app.FWHMnmEditField.Enable = 'off';
        app.FWHMnmEditField.Tooltip = {'fiber.gain_fwhm'};
        app.FWHMnmEditField.Position = [93 144 32 22];
        app.FWHMnmEditField.Value = 40;

        % Create TabGroup3
        app.TabGroup3 = uitabgroup(app.GainmodelTab);
        app.TabGroup3.Position = [18 6 325 134];

        % Create EnergyTab
        app.EnergyTab = uitab(app.TabGroup3);
        app.EnergyTab.Title = 'Energy';

        % Create SaturationenergynJEditFieldLabel
        app.SaturationenergynJEditFieldLabel = uilabel(app.EnergyTab);
        app.SaturationenergynJEditFieldLabel.HorizontalAlignment = 'right';
        app.SaturationenergynJEditFieldLabel.Enable = 'off';
        app.SaturationenergynJEditFieldLabel.Position = [12 70 124 22];
        app.SaturationenergynJEditFieldLabel.Text = 'Saturation energy (nJ)';

        % Create SaturationenergynJEditField
        app.SaturationenergynJEditField = uieditfield(app.EnergyTab, 'numeric');
        app.SaturationenergynJEditField.Limits = [0 Inf];
        app.SaturationenergynJEditField.Enable = 'off';
        app.SaturationenergynJEditField.Tooltip = {'fiber.saturation_energy'};
        app.SaturationenergynJEditField.Position = [147 70 159 22];
        app.SaturationenergynJEditField.Value = 0.335683395724872;

        % Create IntensityTab
        app.IntensityTab = uitab(app.TabGroup3);
        app.IntensityTab.Title = 'Intensity';

        % Create SaturationintensityJm2EditFieldLabel
        app.SaturationintensityJm2EditFieldLabel = uilabel(app.IntensityTab);
        app.SaturationintensityJm2EditFieldLabel.HorizontalAlignment = 'right';
        app.SaturationintensityJm2EditFieldLabel.Enable = 'off';
        app.SaturationintensityJm2EditFieldLabel.Position = [17 70 151 22];
        app.SaturationintensityJm2EditFieldLabel.Text = 'Saturation intensity (J/m^2)';

        % Create SaturationintensityJm2EditField
        app.SaturationintensityJm2EditField = uieditfield(app.IntensityTab, 'numeric');
        app.SaturationintensityJm2EditField.Limits = [0 Inf];
        app.SaturationintensityJm2EditField.Enable = 'off';
        app.SaturationintensityJm2EditField.Tooltip = {'fiber.saturation_intensity'};
        app.SaturationintensityJm2EditField.Position = [183 70 123 22];
        app.SaturationintensityJm2EditField.Value = 11.1187662317349;

        % Create DerivationfromscratchTab
        app.DerivationfromscratchTab = uitab(app.TabGroup3);
        app.DerivationfromscratchTab.Title = 'Derivation from scratch';

        % Create UpperstatelifetimesEditFieldLabel
        app.UpperstatelifetimesEditFieldLabel = uilabel(app.DerivationfromscratchTab);
        app.UpperstatelifetimesEditFieldLabel.HorizontalAlignment = 'right';
        app.UpperstatelifetimesEditFieldLabel.Enable = 'off';
        app.UpperstatelifetimesEditFieldLabel.Position = [11 79 127 22];
        app.UpperstatelifetimesEditFieldLabel.Text = 'Upper-state lifetime (s)';

        % Create UpperstatelifetimesEditField
        app.UpperstatelifetimesEditField = uieditfield(app.DerivationfromscratchTab, 'numeric');
        app.UpperstatelifetimesEditField.Limits = [0 Inf];
        app.UpperstatelifetimesEditField.ValueChangedFcn = createCallbackFcn(app, @UpperstatelifetimesEditFieldValueChanged, true);
        app.UpperstatelifetimesEditField.Enable = 'off';
        app.UpperstatelifetimesEditField.Tooltip = {'fiber.gain_tau'};
        app.UpperstatelifetimesEditField.Position = [153 79 70 22];
        app.UpperstatelifetimesEditField.Value = 0.00084;

        % Create LaserrepetitionrateHzEditFieldLabel
        app.LaserrepetitionrateHzEditFieldLabel = uilabel(app.DerivationfromscratchTab);
        app.LaserrepetitionrateHzEditFieldLabel.HorizontalAlignment = 'right';
        app.LaserrepetitionrateHzEditFieldLabel.Enable = 'off';
        app.LaserrepetitionrateHzEditFieldLabel.Position = [11 52 138 22];
        app.LaserrepetitionrateHzEditFieldLabel.Text = 'Laser repetition rate (Hz)';

        % Create LaserrepetitionrateHzEditField
        app.LaserrepetitionrateHzEditField = uieditfield(app.DerivationfromscratchTab, 'numeric');
        app.LaserrepetitionrateHzEditField.Limits = [0 Inf];
        app.LaserrepetitionrateHzEditField.ValueChangedFcn = createCallbackFcn(app, @LaserrepetitionrateHzEditFieldValueChanged, true);
        app.LaserrepetitionrateHzEditField.Enable = 'off';
        app.LaserrepetitionrateHzEditField.Tooltip = {'fiber.t_rep'};
        app.LaserrepetitionrateHzEditField.Position = [154 52 70 22];
        app.LaserrepetitionrateHzEditField.Value = 30000000;

        % Create Totalcrosssectionm2emissionabsorptionEditFieldLabel
        app.Totalcrosssectionm2emissionabsorptionEditFieldLabel = uilabel(app.DerivationfromscratchTab);
        app.Totalcrosssectionm2emissionabsorptionEditFieldLabel.HorizontalAlignment = 'right';
        app.Totalcrosssectionm2emissionabsorptionEditFieldLabel.Enable = 'off';
        app.Totalcrosssectionm2emissionabsorptionEditFieldLabel.Tooltip = {'fiber.gain_cross_section'};
        app.Totalcrosssectionm2emissionabsorptionEditFieldLabel.Position = [4 14 138 28];
        app.Totalcrosssectionm2emissionabsorptionEditFieldLabel.Text = {'Total cross section (m^2)'; '(emission + absorption)'};

        % Create Totalcrosssectionm2emissionabsorptionEditField
        app.Totalcrosssectionm2emissionabsorptionEditField = uieditfield(app.DerivationfromscratchTab, 'numeric');
        app.Totalcrosssectionm2emissionabsorptionEditField.Limits = [0 Inf];
        app.Totalcrosssectionm2emissionabsorptionEditField.ValueChangedFcn = createCallbackFcn(app, @Totalcrosssectionm2emissionabsorptionEditFieldValueChanged, true);
        app.Totalcrosssectionm2emissionabsorptionEditField.Enable = 'off';
        app.Totalcrosssectionm2emissionabsorptionEditField.Tooltip = {'fiber.gain_cross_section'};
        app.Totalcrosssectionm2emissionabsorptionEditField.Position = [157 20 70 22];
        app.Totalcrosssectionm2emissionabsorptionEditField.Value = 6.883e-25;

        % Create MPATab
        app.MPATab = uitab(app.TabGroup2);
        app.MPATab.Tooltip = {'Parallelization algorithm for multimode simulations'};
        app.MPATab.Title = 'MPA';

        % Create NumberofparallelizationsEditFieldLabel
        app.NumberofparallelizationsEditFieldLabel = uilabel(app.MPATab);
        app.NumberofparallelizationsEditFieldLabel.HorizontalAlignment = 'right';
        app.NumberofparallelizationsEditFieldLabel.Enable = 'off';
        app.NumberofparallelizationsEditFieldLabel.Tooltip = {''};
        app.NumberofparallelizationsEditFieldLabel.Position = [12 176 144 22];
        app.NumberofparallelizationsEditFieldLabel.Text = 'Number of parallelizations';

        % Create NumberofparallelizationsEditField
        app.NumberofparallelizationsEditField = uieditfield(app.MPATab, 'numeric');
        app.NumberofparallelizationsEditField.Limits = [1 20];
        app.NumberofparallelizationsEditField.ValueChangedFcn = createCallbackFcn(app, @NumberofparallelizationsEditFieldValueChanged, true);
        app.NumberofparallelizationsEditField.Enable = 'off';
        app.NumberofparallelizationsEditField.Tooltip = {'sim.MPA.M'};
        app.NumberofparallelizationsEditField.Position = [161 176 32 22];
        app.NumberofparallelizationsEditField.Value = 10;

        % Create MaxiterationsEditFieldLabel
        app.MaxiterationsEditFieldLabel = uilabel(app.MPATab);
        app.MaxiterationsEditFieldLabel.HorizontalAlignment = 'right';
        app.MaxiterationsEditFieldLabel.Enable = 'off';
        app.MaxiterationsEditFieldLabel.Tooltip = {''};
        app.MaxiterationsEditFieldLabel.Position = [12 144 80 22];
        app.MaxiterationsEditFieldLabel.Text = 'Max iterations';

        % Create MaxiterationsEditField
        app.MaxiterationsEditField = uieditfield(app.MPATab, 'numeric');
        app.MaxiterationsEditField.Limits = [1 Inf];
        app.MaxiterationsEditField.Enable = 'off';
        app.MaxiterationsEditField.Tooltip = {'sim.MPA.n_tot_max'};
        app.MaxiterationsEditField.Position = [98 143 34 22];
        app.MaxiterationsEditField.Value = 20;

        % Create MiniterationsEditFieldLabel
        app.MiniterationsEditFieldLabel = uilabel(app.MPATab);
        app.MiniterationsEditFieldLabel.HorizontalAlignment = 'right';
        app.MiniterationsEditFieldLabel.Enable = 'off';
        app.MiniterationsEditFieldLabel.Tooltip = {''};
        app.MiniterationsEditFieldLabel.Position = [12 111 77 22];
        app.MiniterationsEditFieldLabel.Text = 'Min iterations';

        % Create MiniterationsEditField
        app.MiniterationsEditField = uieditfield(app.MPATab, 'numeric');
        app.MiniterationsEditField.Limits = [1 Inf];
        app.MiniterationsEditField.Enable = 'off';
        app.MiniterationsEditField.Tooltip = {'sim.MPA.n_tot_min'};
        app.MiniterationsEditField.Position = [98 108 34 22];
        app.MiniterationsEditField.Value = 2;

        % Create ToleranceEditFieldLabel
        app.ToleranceEditFieldLabel = uilabel(app.MPATab);
        app.ToleranceEditFieldLabel.HorizontalAlignment = 'right';
        app.ToleranceEditFieldLabel.Enable = 'off';
        app.ToleranceEditFieldLabel.Tooltip = {''};
        app.ToleranceEditFieldLabel.Position = [13 77 58 22];
        app.ToleranceEditFieldLabel.Text = 'Tolerance';

        % Create ToleranceEditField
        app.ToleranceEditField = uieditfield(app.MPATab, 'numeric');
        app.ToleranceEditField.Limits = [0 Inf];
        app.ToleranceEditField.Enable = 'off';
        app.ToleranceEditField.Tooltip = {'sim.MPA.tol'};
        app.ToleranceEditField.Position = [78 77 56 22];
        app.ToleranceEditField.Value = 0.0001;

        % Create AdaptivealgorithmTab
        app.AdaptivealgorithmTab = uitab(app.TabGroup2);
        app.AdaptivealgorithmTab.Title = 'Adaptive algorithm';

        % Create ThresholdEditFieldLabel
        app.ThresholdEditFieldLabel = uilabel(app.AdaptivealgorithmTab);
        app.ThresholdEditFieldLabel.HorizontalAlignment = 'right';
        app.ThresholdEditFieldLabel.Tooltip = {''};
        app.ThresholdEditFieldLabel.Position = [13 176 59 22];
        app.ThresholdEditFieldLabel.Text = 'Threshold';

        % Create ThresholdEditField
        app.ThresholdEditField = uieditfield(app.AdaptivealgorithmTab, 'numeric');
        app.ThresholdEditField.Limits = [0 1];
        app.ThresholdEditField.Tooltip = {'sim.adaptive_deltaZ.threshold'; 'Recommended value is less than 1e-5. A value larger than 1e-3 is too large.'};
        app.ThresholdEditField.Position = [79 176 55 22];
        app.ThresholdEditField.Value = 1e-06;

        % Create maxdeltaZEditFieldLabel
        app.maxdeltaZEditFieldLabel = uilabel(app.AdaptivealgorithmTab);
        app.maxdeltaZEditFieldLabel.HorizontalAlignment = 'right';
        app.maxdeltaZEditFieldLabel.Position = [10 142 65 22];
        app.maxdeltaZEditFieldLabel.Text = 'max deltaZ';

        % Create maxdeltaZEditField
        app.maxdeltaZEditField = uieditfield(app.AdaptivealgorithmTab, 'numeric');
        app.maxdeltaZEditField.Limits = [0 Inf];
        app.maxdeltaZEditField.ValueChangedFcn = createCallbackFcn(app, @maxdeltaZEditFieldValueChanged, true);
        app.maxdeltaZEditField.Position = [79 143 57 22];

        % Create RamanTab
        app.RamanTab = uitab(app.TabGroup2);
        app.RamanTab.Title = 'Raman';

        % Create IncludespontaneousRamanButton
        app.IncludespontaneousRamanButton = uibutton(app.RamanTab, 'state');
        app.IncludespontaneousRamanButton.Tooltip = {'sim.Raman_sponRS'};
        app.IncludespontaneousRamanButton.Text = 'Include spontaneous Raman';
        app.IncludespontaneousRamanButton.Position = [18 176 168 22];
        app.IncludespontaneousRamanButton.Value = true;

        % Create GPUTab
        app.GPUTab = uitab(app.TabGroup2);
        app.GPUTab.Title = 'GPU';

        % Create GPUindextouseEditFieldLabel
        app.GPUindextouseEditFieldLabel = uilabel(app.GPUTab);
        app.GPUindextouseEditFieldLabel.HorizontalAlignment = 'right';
        app.GPUindextouseEditFieldLabel.Tooltip = {''};
        app.GPUindextouseEditFieldLabel.Position = [14 176 100 22];
        app.GPUindextouseEditFieldLabel.Text = 'GPU index to use';

        % Create GPUindextouseEditField
        app.GPUindextouseEditField = uieditfield(app.GPUTab, 'numeric');
        app.GPUindextouseEditField.Limits = [1 Inf];
        app.GPUindextouseEditField.ValueChangedFcn = createCallbackFcn(app, @GPUindextouseEditFieldValueChanged, true);
        app.GPUindextouseEditField.Tooltip = {'sim.gpuDevice.Index'};
        app.GPUindextouseEditField.Position = [129 176 37 22];
        app.GPUindextouseEditField.Value = 1;

        % Create CudafolderEditFieldLabel
        app.CudafolderEditFieldLabel = uilabel(app.GPUTab);
        app.CudafolderEditFieldLabel.HorizontalAlignment = 'right';
        app.CudafolderEditFieldLabel.Position = [14 139 68 22];
        app.CudafolderEditFieldLabel.Text = 'Cuda folder';

        % Create CudafolderEditField
        app.CudafolderEditField = uieditfield(app.GPUTab, 'text');
        app.CudafolderEditField.Position = [97 139 117 22];
        app.CudafolderEditField.Value = 'GMMNLSE/cuda';

        % Create ProgressbarTab
        app.ProgressbarTab = uitab(app.TabGroup2);
        app.ProgressbarTab.Title = 'Progress bar';

        % Create ProgressbarnameEditFieldLabel
        app.ProgressbarnameEditFieldLabel = uilabel(app.ProgressbarTab);
        app.ProgressbarnameEditFieldLabel.HorizontalAlignment = 'right';
        app.ProgressbarnameEditFieldLabel.Position = [12 176 108 22];
        app.ProgressbarnameEditFieldLabel.Text = 'Progress-bar name';

        % Create ProgressbarnameEditField
        app.ProgressbarnameEditField = uieditfield(app.ProgressbarTab, 'text');
        app.ProgressbarnameEditField.Tooltip = {'sim.progress_bar_name'};
        app.ProgressbarnameEditField.Position = [135 176 73 22];

        % Create SimulationtypeButtonGroup
        app.SimulationtypeButtonGroup = uibuttongroup(app.UIFigure);
        app.SimulationtypeButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @SimulationtypeButtonGroupSelectionChanged, true);
        app.SimulationtypeButtonGroup.Title = 'Simulation type';
        app.SimulationtypeButtonGroup.Position = [13 256 123 75];

        % Create singlemodeButton
        app.singlemodeButton = uiradiobutton(app.SimulationtypeButtonGroup);
        app.singlemodeButton.Text = 'single mode';
        app.singlemodeButton.Position = [11 29 87 22];
        app.singlemodeButton.Value = true;

        % Create multimodeButton
        app.multimodeButton = uiradiobutton(app.SimulationtypeButtonGroup);
        app.multimodeButton.Text = 'multimode';
        app.multimodeButton.Position = [11 7 78 22];

        % Create IncludeshocktermButton
        app.IncludeshocktermButton = uibutton(app.UIFigure, 'state');
        app.IncludeshocktermButton.Tooltip = {'sim.sw'};
        app.IncludeshocktermButton.Text = 'Include shock term';
        app.IncludeshocktermButton.Position = [522 546 116 22];
        app.IncludeshocktermButton.Value = true;

        % Create StepsizemEditFieldLabel
        app.StepsizemEditFieldLabel = uilabel(app.UIFigure);
        app.StepsizemEditFieldLabel.HorizontalAlignment = 'right';
        app.StepsizemEditFieldLabel.Enable = 'off';
        app.StepsizemEditFieldLabel.Tooltip = {''};
        app.StepsizemEditFieldLabel.Position = [583 224 76 22];
        app.StepsizemEditFieldLabel.Text = 'Step size (m)';

        % Create StepsizemEditField
        app.StepsizemEditField = uieditfield(app.UIFigure, 'numeric');
        app.StepsizemEditField.Limits = [0 Inf];
        app.StepsizemEditField.ValueChangedFcn = createCallbackFcn(app, @StepsizemEditFieldValueChanged, true);
        app.StepsizemEditField.Enable = 'off';
        app.StepsizemEditField.Tooltip = {'sim.deltaZ'; 'This is used only with a fixed step size.'};
        app.StepsizemEditField.Position = [675 224 69 22];
        app.StepsizemEditField.Value = 0.001;

        % Create SavingstepsizemEditFieldLabel
        app.SavingstepsizemEditFieldLabel = uilabel(app.UIFigure);
        app.SavingstepsizemEditFieldLabel.HorizontalAlignment = 'right';
        app.SavingstepsizemEditFieldLabel.Tooltip = {''};
        app.SavingstepsizemEditFieldLabel.Position = [546 193 114 22];
        app.SavingstepsizemEditFieldLabel.Text = 'Saving step size (m)';

        % Create SavingstepsizemEditField
        app.SavingstepsizemEditField = uieditfield(app.UIFigure, 'numeric');
        app.SavingstepsizemEditField.Limits = [0 Inf];
        app.SavingstepsizemEditField.ValueChangedFcn = createCallbackFcn(app, @SavingstepsizemEditFieldValueChanged, true);
        app.SavingstepsizemEditField.Tooltip = {'sim.save_period'; '0 = only save input and output (saving step = fiber length)'};
        app.SavingstepsizemEditField.Position = [675 193 70 22];

        % Create ScalarPolarizedsimulationsButtonGroup
        app.ScalarPolarizedsimulationsButtonGroup = uibuttongroup(app.UIFigure);
        app.ScalarPolarizedsimulationsButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @ScalarPolarizedsimulationsButtonGroupSelectionChanged, true);
        app.ScalarPolarizedsimulationsButtonGroup.Tooltip = {'sim.scalar'; 'Make sure your fields, betas, and SR have the correct dimensions.'};
        app.ScalarPolarizedsimulationsButtonGroup.Title = 'Scalar/Polarized simulations';
        app.ScalarPolarizedsimulationsButtonGroup.Position = [150 256 164 75];

        % Create scalarButton
        app.scalarButton = uiradiobutton(app.ScalarPolarizedsimulationsButtonGroup);
        app.scalarButton.Text = 'scalar';
        app.scalarButton.Position = [11 29 54 22];
        app.scalarButton.Value = true;

        % Create polarizedButton
        app.polarizedButton = uiradiobutton(app.ScalarPolarizedsimulationsButtonGroup);
        app.polarizedButton.Text = 'polarized';
        app.polarizedButton.Position = [11 7 71 22];

        % Create EllipticityEditFieldLabel
        app.EllipticityEditFieldLabel = uilabel(app.UIFigure);
        app.EllipticityEditFieldLabel.HorizontalAlignment = 'right';
        app.EllipticityEditFieldLabel.Tooltip = {''};
        app.EllipticityEditFieldLabel.Position = [614 162 52 22];
        app.EllipticityEditFieldLabel.Text = 'Ellipticity';

        % Create EllipticityEditField
        app.EllipticityEditField = uieditfield(app.UIFigure, 'numeric');
        app.EllipticityEditField.Limits = [0 1];
        app.EllipticityEditField.ValueChangedFcn = createCallbackFcn(app, @EllipticityEditFieldValueChanged, true);
        app.EllipticityEditField.Tooltip = {'sim.ellipticity'; '1 = linear polarization; 0 = circular polarization'};
        app.EllipticityEditField.Position = [675 162 41 22];

        % Create AdaptivealgorithmButtonGroup
        app.AdaptivealgorithmButtonGroup = uibuttongroup(app.UIFigure);
        app.AdaptivealgorithmButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @AdaptivealgorithmButtonGroupSelectionChanged, true);
        app.AdaptivealgorithmButtonGroup.Tooltip = {'sim.adaptive_deltaZ.model'};
        app.AdaptivealgorithmButtonGroup.Title = 'Adaptive algorithm';
        app.AdaptivealgorithmButtonGroup.Position = [325 256 164 75];

        % Create adaptivestepsizeButton
        app.adaptivestepsizeButton = uiradiobutton(app.AdaptivealgorithmButtonGroup);
        app.adaptivestepsizeButton.Text = 'adaptive step size';
        app.adaptivestepsizeButton.Position = [11 29 118 22];
        app.adaptivestepsizeButton.Value = true;

        % Create fixedstepsizeButton
        app.fixedstepsizeButton = uiradiobutton(app.AdaptivealgorithmButtonGroup);
        app.fixedstepsizeButton.Text = 'fixed step size';
        app.fixedstepsizeButton.Position = [11 7 98 22];

        % Create NumericalprecisionButtonGroup
        app.NumericalprecisionButtonGroup = uibuttongroup(app.UIFigure);
        app.NumericalprecisionButtonGroup.Tooltip = {'sim.single_yes'; 'Single precision doesn''t help boost the computation speed a lot.'};
        app.NumericalprecisionButtonGroup.Title = 'Numerical precision';
        app.NumericalprecisionButtonGroup.Position = [499 256 123 75];

        % Create singleButton
        app.singleButton = uiradiobutton(app.NumericalprecisionButtonGroup);
        app.singleButton.Text = 'single';
        app.singleButton.Position = [11 29 54 22];

        % Create doubleButton
        app.doubleButton = uiradiobutton(app.NumericalprecisionButtonGroup);
        app.doubleButton.Text = 'double';
        app.doubleButton.Position = [11 7 58 22];
        app.doubleButton.Value = true;

        % Create RamanmodelDropDownLabel
        app.RamanmodelDropDownLabel = uilabel(app.UIFigure);
        app.RamanmodelDropDownLabel.HorizontalAlignment = 'right';
        app.RamanmodelDropDownLabel.Position = [522 399 80 22];
        app.RamanmodelDropDownLabel.Text = 'Raman model';

        % Create RamanmodelDropDown
        app.RamanmodelDropDown = uidropdown(app.UIFigure);
        app.RamanmodelDropDown.Items = {'No Raman', 'isotropic Raman', 'anisotropic Raman'};
        app.RamanmodelDropDown.ItemsData = [0 1 2];
        app.RamanmodelDropDown.Tooltip = {'sim.Raman_model'};
        app.RamanmodelDropDown.Position = [617 399 126 22];
        app.RamanmodelDropDown.Value = 1;

        % Create SteppingalgorithmDropDownLabel
        app.SteppingalgorithmDropDownLabel = uilabel(app.UIFigure);
        app.SteppingalgorithmDropDownLabel.HorizontalAlignment = 'right';
        app.SteppingalgorithmDropDownLabel.Position = [522 426 106 22];
        app.SteppingalgorithmDropDownLabel.Text = 'Stepping algorithm';

        % Create SteppingalgorithmDropDown
        app.SteppingalgorithmDropDown = uidropdown(app.UIFigure);
        app.SteppingalgorithmDropDown.Items = {'split-step', 'RK4IP', 'MPA'};
        app.SteppingalgorithmDropDown.ValueChangedFcn = createCallbackFcn(app, @SteppingalgorithmDropDownValueChanged, true);
        app.SteppingalgorithmDropDown.Tooltip = {'sim.step_method'};
        app.SteppingalgorithmDropDown.Position = [643 426 100 22];
        app.SteppingalgorithmDropDown.Value = 'RK4IP';

        % Create UseGPUButton
        app.UseGPUButton = uibutton(app.UIFigure, 'state');
        app.UseGPUButton.ValueChangedFcn = createCallbackFcn(app, @UseGPUButtonValueChanged, true);
        app.UseGPUButton.Tooltip = {'sim.gpu_yes'};
        app.UseGPUButton.Text = 'Use GPU';
        app.UseGPUButton.Position = [522 515 100 22];
        app.UseGPUButton.Value = true;

        % Create CenterthepulsewrtthetimewindowButton
        app.CenterthepulsewrtthetimewindowButton = uibutton(app.UIFigure, 'state');
        app.CenterthepulsewrtthetimewindowButton.Tooltip = {'sim.pulse_centering'; 'The offset time delay will be saved in the output for recovery.'};
        app.CenterthepulsewrtthetimewindowButton.Text = 'Center the pulse w.r.t. the time window';
        app.CenterthepulsewrtthetimewindowButton.Position = [522 456 222 22];
        app.CenterthepulsewrtthetimewindowButton.Value = true;

        % Create PhotonnoiseperfreqbandLabel
        app.PhotonnoiseperfreqbandLabel = uilabel(app.UIFigure);
        app.PhotonnoiseperfreqbandLabel.HorizontalAlignment = 'right';
        app.PhotonnoiseperfreqbandLabel.Tooltip = {''};
        app.PhotonnoiseperfreqbandLabel.Position = [502 133 157 22];
        app.PhotonnoiseperfreqbandLabel.Text = '#Photon noise per freq band';

        % Create PhotonnoiseperfreqbandEditField
        app.PhotonnoiseperfreqbandEditField = uieditfield(app.UIFigure, 'numeric');
        app.PhotonnoiseperfreqbandEditField.Limits = [0 Inf];
        app.PhotonnoiseperfreqbandEditField.Tooltip = {''};
        app.PhotonnoiseperfreqbandEditField.Position = [674 133 42 22];

        % Create GainmodelDropDownLabel
        app.GainmodelDropDownLabel = uilabel(app.UIFigure);
        app.GainmodelDropDownLabel.HorizontalAlignment = 'right';
        app.GainmodelDropDownLabel.Position = [522 372 67 22];
        app.GainmodelDropDownLabel.Text = 'Gain model';

        % Create GainmodelDropDown
        app.GainmodelDropDown = uidropdown(app.UIFigure);
        app.GainmodelDropDown.Items = {'No gain', 'SM gain', 'New gain', 'Taylor gain', 'Rate-equation gain'};
        app.GainmodelDropDown.ItemsData = [0 1 2 3 4];
        app.GainmodelDropDown.ValueChangedFcn = createCallbackFcn(app, @GainmodelDropDownValueChanged, true);
        app.GainmodelDropDown.Tooltip = {'sim.gain_model'};
        app.GainmodelDropDown.Position = [604 372 139 22];
        app.GainmodelDropDown.Value = 0;

        % Create FiberPanel
        app.FiberPanel = uipanel(app.UIFigure);
        app.FiberPanel.Title = 'Fiber';
        app.FiberPanel.Position = [16 341 495 220];

        % Create FiberlengthmEditFieldLabel
        app.FiberlengthmEditFieldLabel = uilabel(app.FiberPanel);
        app.FiberlengthmEditFieldLabel.HorizontalAlignment = 'right';
        app.FiberlengthmEditFieldLabel.Tooltip = {''};
        app.FiberlengthmEditFieldLabel.Position = [8 131 90 22];
        app.FiberlengthmEditFieldLabel.Text = 'Fiber length (m)';

        % Create FiberlengthmEditField
        app.FiberlengthmEditField = uieditfield(app.FiberPanel, 'numeric');
        app.FiberlengthmEditField.Limits = [0 Inf];
        app.FiberlengthmEditField.ValueChangedFcn = createCallbackFcn(app, @FiberlengthmEditFieldValueChanged, true);
        app.FiberlengthmEditField.Tooltip = {'fiber.L0'};
        app.FiberlengthmEditField.Position = [104 131 48 22];
        app.FiberlengthmEditField.Value = 1;

        % Create Nonlinearrefractiveindexn2m2WEditFieldLabel
        app.Nonlinearrefractiveindexn2m2WEditFieldLabel = uilabel(app.FiberPanel);
        app.Nonlinearrefractiveindexn2m2WEditFieldLabel.HorizontalAlignment = 'right';
        app.Nonlinearrefractiveindexn2m2WEditFieldLabel.Tooltip = {''};
        app.Nonlinearrefractiveindexn2m2WEditFieldLabel.Position = [8 73 208 22];
        app.Nonlinearrefractiveindexn2m2WEditFieldLabel.Text = 'Nonlinear refractive index n2 (m^2/W)';

        % Create Nonlinearrefractiveindexn2m2WEditField
        app.Nonlinearrefractiveindexn2m2WEditField = uieditfield(app.FiberPanel, 'numeric');
        app.Nonlinearrefractiveindexn2m2WEditField.Tooltip = {'fiber.n2'};
        app.Nonlinearrefractiveindexn2m2WEditField.Position = [130 52 85 22];
        app.Nonlinearrefractiveindexn2m2WEditField.Value = 2.3e-20;

        % Create FibermaterialDropDownLabel
        app.FibermaterialDropDownLabel = uilabel(app.FiberPanel);
        app.FibermaterialDropDownLabel.HorizontalAlignment = 'right';
        app.FibermaterialDropDownLabel.Position = [8 165 79 22];
        app.FibermaterialDropDownLabel.Text = 'Fiber material';

        % Create FibermaterialDropDown
        app.FibermaterialDropDown = uidropdown(app.FiberPanel);
        app.FibermaterialDropDown.Items = {'silica', 'calcogenide', 'ZBLAN'};
        app.FibermaterialDropDown.ValueChangedFcn = createCallbackFcn(app, @FibermaterialDropDownValueChanged, true);
        app.FibermaterialDropDown.Tooltip = {'fiber.material'; 'This controls which Raman model to use.'};
        app.FibermaterialDropDown.Position = [102 165 100 22];
        app.FibermaterialDropDown.Value = 'silica';

        % Create TabGroup
        app.TabGroup = uitabgroup(app.FiberPanel);
        app.TabGroup.Position = [263 4 224 191];

        % Create SinglemodeTab
        app.SinglemodeTab = uitab(app.TabGroup);
        app.SinglemodeTab.Title = 'Single mode';

        % Create MFDumEditField_2Label
        app.MFDumEditField_2Label = uilabel(app.SinglemodeTab);
        app.MFDumEditField_2Label.HorizontalAlignment = 'right';
        app.MFDumEditField_2Label.Position = [6 136 59 22];
        app.MFDumEditField_2Label.Text = 'MFD (um)';

        % Create MFDumEditField_2
        app.MFDumEditField_2 = uieditfield(app.SinglemodeTab, 'numeric');
        app.MFDumEditField_2.Limits = [0 Inf];
        app.MFDumEditField_2.ValueChangedFcn = createCallbackFcn(app, @MFDumEditField_2ValueChanged, true);
        app.MFDumEditField_2.Tooltip = {'fiber.MFD'; 'Mode-field diameter'};
        app.MFDumEditField_2.Position = [76 136 41 22];
        app.MFDumEditField_2.Value = 5.95;

        % Create BetasButtonGroup
        app.BetasButtonGroup = uibuttongroup(app.SinglemodeTab);
        app.BetasButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @BetasButtonGroupSelectionChanged, true);
        app.BetasButtonGroup.Tooltip = {'fiber.betas'; 'Taylor-expansion coefficients of the propagation constant over the center angular frequency'};
        app.BetasButtonGroup.Title = 'Betas';
        app.BetasButtonGroup.Position = [6 7 215 126];

        % Create Beta0Button
        app.Beta0Button = uitogglebutton(app.BetasButtonGroup);
        app.Beta0Button.Text = 'Beta0';
        app.Beta0Button.Position = [2 71 39 22];
        app.Beta0Button.Value = true;

        % Create Beta1Button
        app.Beta1Button = uitogglebutton(app.BetasButtonGroup);
        app.Beta1Button.Text = 'Beta1';
        app.Beta1Button.Position = [41 71 39 22];

        % Create Beta2Button
        app.Beta2Button = uitogglebutton(app.BetasButtonGroup);
        app.Beta2Button.Text = 'Beta2';
        app.Beta2Button.Position = [80 71 39 22];

        % Create Beta3Button
        app.Beta3Button = uitogglebutton(app.BetasButtonGroup);
        app.Beta3Button.Text = 'Beta3';
        app.Beta3Button.Position = [118 71 39 22];

        % Create Beta4Button
        app.Beta4Button = uitogglebutton(app.BetasButtonGroup);
        app.Beta4Button.Text = 'Beta4';
        app.Beta4Button.Position = [156 71 47 22];

        % Create mLabel
        app.mLabel = uilabel(app.BetasButtonGroup);
        app.mLabel.Position = [103 41 41 22];
        app.mLabel.Text = '1/m';

        % Create mLabel_3
        app.mLabel_3 = uilabel(app.BetasButtonGroup);
        app.mLabel_3.Enable = 'off';
        app.mLabel_3.Position = [103 8 41 22];
        app.mLabel_3.Text = '1/m';

        % Create DefaultButton
        app.DefaultButton = uibutton(app.BetasButtonGroup, 'push');
        app.DefaultButton.ButtonPushedFcn = createCallbackFcn(app, @DefaultButtonPushed, true);
        app.DefaultButton.Tooltip = {'This code saves two betas based on normal or anomalous dispersion.'};
        app.DefaultButton.Position = [152 21 55 22];
        app.DefaultButton.Text = 'Default';

        % Create XEditFieldLabel
        app.XEditFieldLabel = uilabel(app.BetasButtonGroup);
        app.XEditFieldLabel.HorizontalAlignment = 'right';
        app.XEditFieldLabel.Position = [7 41 12 22];
        app.XEditFieldLabel.Text = 'X';

        % Create XEditField
        app.XEditField = uieditfield(app.BetasButtonGroup, 'numeric');
        app.XEditField.ValueChangedFcn = createCallbackFcn(app, @XEditFieldValueChanged, true);
        app.XEditField.Position = [27 41 73 22];

        % Create YEditFieldLabel
        app.YEditFieldLabel = uilabel(app.BetasButtonGroup);
        app.YEditFieldLabel.HorizontalAlignment = 'right';
        app.YEditFieldLabel.Enable = 'off';
        app.YEditFieldLabel.Position = [-6 8 25 22];
        app.YEditFieldLabel.Text = 'Y';

        % Create YEditField
        app.YEditField = uieditfield(app.BetasButtonGroup, 'numeric');
        app.YEditField.ValueChangedFcn = createCallbackFcn(app, @YEditFieldValueChanged, true);
        app.YEditField.Enable = 'off';
        app.YEditField.Position = [27 8 73 22];

        % Create MultimodeTab
        app.MultimodeTab = uitab(app.TabGroup);
        app.MultimodeTab.Title = 'Multimode';

        % Create modeindexEditFieldLabel
        app.modeindexEditFieldLabel = uilabel(app.MultimodeTab);
        app.modeindexEditFieldLabel.HorizontalAlignment = 'right';
        app.modeindexEditFieldLabel.Enable = 'off';
        app.modeindexEditFieldLabel.Position = [9 136 68 22];
        app.modeindexEditFieldLabel.Text = 'mode index';

        % Create modeindexEditField
        app.modeindexEditField = uieditfield(app.MultimodeTab, 'text');
        app.modeindexEditField.ValueChangedFcn = createCallbackFcn(app, @modeindexEditFieldValueChanged, true);
        app.modeindexEditField.Enable = 'off';
        app.modeindexEditField.Tooltip = {'sim.midx'};
        app.modeindexEditField.Position = [92 136 73 22];
        app.modeindexEditField.Value = '1';

        % Create MM_folderButton
        app.MM_folderButton = uibutton(app.MultimodeTab, 'push');
        app.MM_folderButton.ButtonPushedFcn = createCallbackFcn(app, @MM_folderButtonPushed, true);
        app.MM_folderButton.Enable = 'off';
        app.MM_folderButton.Tooltip = {'fiber.MM_folder'};
        app.MM_folderButton.Position = [9 103 97 22];
        app.MM_folderButton.Text = 'MM_folder';

        % Create BetasfilenameEditFieldLabel
        app.BetasfilenameEditFieldLabel = uilabel(app.MultimodeTab);
        app.BetasfilenameEditFieldLabel.HorizontalAlignment = 'right';
        app.BetasfilenameEditFieldLabel.Enable = 'off';
        app.BetasfilenameEditFieldLabel.Tooltip = {''};
        app.BetasfilenameEditFieldLabel.Position = [9 69 85 22];
        app.BetasfilenameEditFieldLabel.Text = 'Betas filename';

        % Create BetasfilenameEditField
        app.BetasfilenameEditField = uieditfield(app.MultimodeTab, 'text');
        app.BetasfilenameEditField.ValueChangedFcn = createCallbackFcn(app, @BetasfilenameEditFieldValueChanged, true);
        app.BetasfilenameEditField.Enable = 'off';
        app.BetasfilenameEditField.Tooltip = {'fiber.betas_filename'};
        app.BetasfilenameEditField.Position = [124 69 73 22];
        app.BetasfilenameEditField.Value = 'betas.mat';

        % Create StensorsfilenameEditFieldLabel
        app.StensorsfilenameEditFieldLabel = uilabel(app.MultimodeTab);
        app.StensorsfilenameEditFieldLabel.HorizontalAlignment = 'right';
        app.StensorsfilenameEditFieldLabel.Enable = 'off';
        app.StensorsfilenameEditFieldLabel.Position = [9 41 106 22];
        app.StensorsfilenameEditFieldLabel.Text = 'S-tensors filename';

        % Create StensorsfilenameEditField
        app.StensorsfilenameEditField = uieditfield(app.MultimodeTab, 'text');
        app.StensorsfilenameEditField.ValueChangedFcn = createCallbackFcn(app, @StensorsfilenameEditFieldValueChanged, true);
        app.StensorsfilenameEditField.Enable = 'off';
        app.StensorsfilenameEditField.Tooltip = {'fiber.S_tensors_filename'};
        app.StensorsfilenameEditField.Position = [124 41 89 22];
        app.StensorsfilenameEditField.Value = 'S_tensors_6modes.mat';

        % Create ShowprogressbarButton
        app.ShowprogressbarButton = uibutton(app.UIFigure, 'state');
        app.ShowprogressbarButton.ValueChangedFcn = createCallbackFcn(app, @ShowprogressbarButtonValueChanged, true);
        app.ShowprogressbarButton.Tooltip = {'sim.progress_bar'};
        app.ShowprogressbarButton.Text = 'Show progress bar';
        app.ShowprogressbarButton.Position = [522 485 116 22];
        app.ShowprogressbarButton.Value = true;

        % Create TimewindowBetasButtonGroup
        app.TimewindowBetasButtonGroup = uibuttongroup(app.UIFigure);
        app.TimewindowBetasButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @TimewindowBetasButtonGroupSelectionChanged, true);
        app.TimewindowBetasButtonGroup.Tooltip = {'sim.betas'; 'This defaults to the beta0 and beta1 of the first mode.'};
        app.TimewindowBetasButtonGroup.Title = 'Time-window Betas';
        app.TimewindowBetasButtonGroup.Position = [546 29 175 92];

        % Create Beta0Button_2
        app.Beta0Button_2 = uitogglebutton(app.TimewindowBetasButtonGroup);
        app.Beta0Button_2.Text = 'Beta0';
        app.Beta0Button_2.Position = [2 37 39 22];
        app.Beta0Button_2.Value = true;

        % Create Beta1Button_2
        app.Beta1Button_2 = uitogglebutton(app.TimewindowBetasButtonGroup);
        app.Beta1Button_2.Text = 'Beta1';
        app.Beta1Button_2.Position = [41 37 39 22];

        % Create mLabel_2
        app.mLabel_2 = uilabel(app.TimewindowBetasButtonGroup);
        app.mLabel_2.Position = [147 7 26 22];
        app.mLabel_2.Text = '1/m';

        % Create ValueEditFieldLabel
        app.ValueEditFieldLabel = uilabel(app.TimewindowBetasButtonGroup);
        app.ValueEditFieldLabel.HorizontalAlignment = 'right';
        app.ValueEditFieldLabel.Position = [9 7 35 22];
        app.ValueEditFieldLabel.Text = 'Value';

        % Create ValueEditField
        app.ValueEditField = uieditfield(app.TimewindowBetasButtonGroup, 'numeric');
        app.ValueEditField.ValueChangedFcn = createCallbackFcn(app, @ValueEditFieldValueChanged, true);
        app.ValueEditField.Position = [52 7 86 22];

        % Create RmcmodelDropDownLabel
        app.RmcmodelDropDownLabel = uilabel(app.UIFigure);
        app.RmcmodelDropDownLabel.HorizontalAlignment = 'right';
        app.RmcmodelDropDownLabel.Position = [522 341 66 22];
        app.RmcmodelDropDownLabel.Text = 'Rmc model';

        % Create RmcmodelDropDown
        app.RmcmodelDropDown = uidropdown(app.UIFigure);
        app.RmcmodelDropDown.Items = {'No random mode coupling', 'unitary matrix R', 'iQA addition', 'Manakov'};
        app.RmcmodelDropDown.ItemsData = [0 1 2 3];
        app.RmcmodelDropDown.Position = [603 341 140 22];
        app.RmcmodelDropDown.Value = 0;

        % Create ClosethisforexportingtheparametersLabel
        app.ClosethisforexportingtheparametersLabel = uilabel(app.UIFigure);
        app.ClosethisforexportingtheparametersLabel.Position = [491 588 213 22];
        app.ClosethisforexportingtheparametersLabel.Text = 'Close this for exporting the parameters';

        % Show the figure after all components are created
        app.UIFigure.Visible = 'on';
    end

    %% Code that executes after component creation
    function startupFcn()
        if ispc
            sep = '\';
        else % unix
            sep = '/';
        end
        current_path = mfilename('fullpath');
        sep_pos = strfind(current_path,sep);
        app.upper_folder = current_path(1:sep_pos(end));

        read_betas();
        read_SR();

        UseGPUButtonValueChanged(app.UseGPUButton);
        maxdeltaZEditFieldValueChanged(app.maxdeltaZEditField);
    end

    function CloseReqFcn(dummy, event) %#ok
        % Wrap up variables into fiber and sim
        wrapUp();
        
        % Delete UIFigure when app is deleted
        delete(app.UIFigure);
    end

    % Value changed function: CenterfrequencyTHzEditField
    function CenterfrequencyTHzEditFieldValueChanged(new_value, event) %#ok
        app.CenterfrequencyTHzEditField = new_value;
        
        app.CenterwavelengthnmEditField.Value = app.c/app.CenterfrequencyTHzEditField.Value*1e-3; % nm
    end

    % Value changed function: CenterwavelengthnmEditField
    function CenterwavelengthnmEditFieldValueChanged(new_value, event) %#ok
        app.CenterwavelengthnmEditField = new_value;
        
        app.CenterfrequencyTHzEditField.Value = app.c/app.CenterwavelengthnmEditField.Value*1e-3; % THz
    end

    % Selection changed function: SimulationtypeButtonGroup
    function SimulationtypeButtonGroupSelectionChanged(new_value, event) %#ok
        app.SimulationtypeButtonGroup = new_value;
        
        switch app.SimulationtypeButtonGroup.SelectedObject.Text
            case 'single mode'
                % Enable
                app.MFDumEditField_2Label.Enable = 'on';
                app.MFDumEditField_2.Enable = 'on';
                matlab_year = version('-release'); % ButtonGroup has the "Enable" property only after version 2019.
                if str2double(matlab_year(1:4)) > 2019
                    app.BetasButtonGroup.Enable = 'on';
                else
                    app.XEditFieldLabel.Enable = 'on';
                    app.XEditField.Enable = 'on';
                    app.mLabel.Enable = 'on';
                    app.DefaultButton.Enable = 'on';
                end
                % Disable
                app.modeindexEditFieldLabel.Enable = 'off';
                app.modeindexEditField.Enable = 'off';
                app.MM_folderButton.Enable = 'off';
                app.BetasfilenameEditFieldLabel.Enable = 'off';
                app.BetasfilenameEditField.Enable = 'off';
                app.StensorsfilenameEditFieldLabel.Enable = 'off';
                app.StensorsfilenameEditField.Enable = 'off';
            case 'multimode'
                % Enable
                app.modeindexEditFieldLabel.Enable = 'on';
                app.modeindexEditField.Enable = 'on';
                app.MM_folderButton.Enable = 'on';
                app.BetasfilenameEditFieldLabel.Enable = 'on';
                app.BetasfilenameEditField.Enable = 'on';
                app.StensorsfilenameEditFieldLabel.Enable = 'on';
                app.StensorsfilenameEditField.Enable = 'on';
                % Disable
                app.MFDumEditField_2Label.Enable = 'off';
                app.MFDumEditField_2.Enable = 'off';
                matlab_year = version('-release'); % ButtonGroup has the "Enable" property only after version 2019.
                if str2double(matlab_year(1:4)) > 2019
                    app.BetasButtonGroup.Enable = 'off';
                else
                    app.XEditFieldLabel.Enable = 'off';
                    app.XEditField.Enable = 'off';
                    app.mLabel.Enable = 'off';
                    app.DefaultButton.Enable = 'off';
                end
        end
    end

    % Selection changed function: BetasButtonGroup
    function BetasButtonGroupSelectionChanged(new_value, event) %#ok
        app.BetasButtonGroup = new_value;
        
        switch app.BetasButtonGroup.SelectedObject.Text
            case 'Beta0'
                app.XEditField.Value = app.betas(1,1);
                app.YEditField.Value = app.betas(1,2);
                app.mLabel.Text = '1/m';
                app.mLabel_3.Text = '1/m';
            case 'Beta1'
                app.XEditField.Value = app.betas(2,1);
                app.YEditField.Value = app.betas(2,2);
                app.mLabel.Text = 'ps/m';
                app.mLabel_3.Text = 'ps/m';
            case 'Beta2'
                app.XEditField.Value = app.betas(3,1);
                app.YEditField.Value = app.betas(3,2);
                app.mLabel.Text = 'ps^2/m';
                app.mLabel_3.Text = 'ps^2/m';
            case 'Beta3'
                app.XEditField.Value = app.betas(4,1);
                app.YEditField.Value = app.betas(4,2);
                app.mLabel.Text = 'ps^3/m';
                app.mLabel_3.Text = 'ps^3/m';
            case 'Beta4'
                app.XEditField.Value = app.betas(5,1);
                app.YEditField.Value = app.betas(5,2);
                app.mLabel.Text = 'ps^4/m';
                app.mLabel_3.Text = 'ps^4/m';
        end
    end

    % Button pushed function: DefaultButton
    function DefaultButtonPushed(new_value, event) %#ok
        
        if app.CenterwavelengthnmEditField.Value < 1300 % lambda(zero dispersion)=1.3um; assume 1030nm for positive dispersion
            app.betas = repmat([8.8268e6; 4.8821e3; 0.0209; 32.9e-6; -26.7e-9],1,2);
            app.MFDumEditField_2.Value = 5.95; % um; 1030nm from Thorlabs 1060XP
        else % assume 1550nm for anomalous dispersion
            app.betas = repmat([5.8339e6; 4.8775e3; -0.0123; 0.1049e-6; -378.3e-9],1,2);
            app.MFDumEditField_2.Value = 8.09; % um; 1550nm from Thorlabs 1060XP
        end

        % Update
        switch app.BetasButtonGroup.SelectedObject.Text
            case 'Beta0'
                app.XEditField.Value = app.betas(1,1);
                app.YEditField.Value = app.betas(1,2);
            case 'Beta1'
                app.XEditField.Value = app.betas(2,1);
                app.YEditField.Value = app.betas(2,2);
            case 'Beta2'
                app.XEditField.Value = app.betas(3,1);
                app.YEditField.Value = app.betas(3,2);
            case 'Beta3'
                app.XEditField.Value = app.betas(4,1);
                app.YEditField.Value = app.betas(4,2);
            case 'Beta4'
                app.XEditField.Value = app.betas(5,1);
                app.YEditField.Value = app.betas(5,2);
        end

        % Update time-window betas
        app.sim_betas = app.betas([1,2],1);
        switch app.TimewindowBetasButtonGroup.SelectedObject.Text
            case 'Beta0'
                app.ValueEditField.Value = app.sim_betas(1);
            case 'Beta1'
                app.ValueEditField.Value = app.sim_betas(2);
        end
    end

    % Value changed function: MFDumEditField_2
    function MFDumEditField_2ValueChanged(new_value, event) %#ok
        app.MFDumEditField_2 = new_value;
        app.Aeff = pi*(app.MFDumEditField_2.Value/2)^2*1e-12; % m^2; effective area of the SMF
        app.SR = 1/app.Aeff;
    end

    % Button pushed function: MM_folderButton
    function MM_folderButtonPushed(new_value, event) %#ok
        MM_folder_string = uigetdir(app.upper_folder,'Multimode parameter folder');
        figure(app.UIFigure);
        
        if MM_folder_string == 0
            errordlg('Please specify the folder where betas and SR are saved.','MM_folder Error','replace');
        else
            app.MM_folder = MM_folder_string;

            % Load betas and SR from files
            read_betas();
            read_SR();
        end
    end

    % Value changed function: modeindexEditField
    function modeindexEditFieldValueChanged(new_value, event) %#ok
        app.modeindexEditField = new_value;
        
        try
            midx_string = unique(uint8(eval(['[',app.modeindexEditField.Value,']'])));
            if size(midx_string,1) ~= 1
                ME = MException('modeindexEditFieldValueChanged:midxWrongSize','');
                throw(ME);
            end
            if any(midx_string == 0)
                ME = MException('modeindexEditFieldValueChanged:midxWrongValue','');
                throw(ME);
            end

            app.midx = midx_string;

        catch
            errordlg('mode index needs to be a 1-by-n integer array.','midx Error','replace');
        end

        % If there's no error,
        % I do this again in case the input contains decimals for which uint8() removes the decimal parts.
        % If there's an error,
        % put the original midx back
        app.modeindexEditField.Value = num2str(app.midx);

        % Load betas and SR from files
        read_betas();
        read_SR();
    end

    % Value changed function: XEditField
    function XEditFieldValueChanged(new_value, event) %#ok
        app.XEditField = new_value;
        
        switch app.BetasButtonGroup.SelectedObject.Text
            case 'Beta0'
                app.betas(1,1) = app.XEditField.Value;
            case 'Beta1'
                app.betas(2,1) = app.XEditField.Value;
            case 'Beta2'
                app.betas(3,1) = app.XEditField.Value;
            case 'Beta3'
                app.betas(4,1) = app.XEditField.Value;
            case 'Beta4'
                app.betas(5,1) = app.XEditField.Value;
        end
    end

    % Value changed function: YEditField
    function YEditFieldValueChanged(new_value, event) %#ok
        app.YEditField = new_value;
        
        switch app.BetasButtonGroup.SelectedObject.Text
            case 'Beta0'
                app.betas(1,2) = app.YEditField.Value;
            case 'Beta1'
                app.betas(2,2) = app.YEditField.Value;
            case 'Beta2'
                app.betas(3,2) = app.YEditField.Value;
            case 'Beta3'
                app.betas(4,2) = app.YEditField.Value;
            case 'Beta4'
                app.betas(5,2) = app.YEditField.Value;
        end
    end

    % Value changed function: BetasfilenameEditField
    function BetasfilenameEditFieldValueChanged(new_value, event) %#ok
        app.BetasfilenameEditField = new_value;
        
        % Load betas from the specified file
        read_betas();
    end

    % Value changed function: StensorsfilenameEditField
    function StensorsfilenameEditFieldValueChanged(new_value, event) %#ok
        app.StensorsfilenameEditField = new_value;
        
        % Load SR from the specified file
        read_SR();
    end

    % Value changed function: GaindBmEditField
    function GaindBmEditFieldValueChanged(new_value, event) %#ok
        app.GaindBmEditField = new_value;
        
        app.gaincoeff1mEditField.Value = app.GaindBmEditField.Value*log(10)/(10*app.FiberlengthmEditField.Value); % 1/m
    end

    % Value changed function: gaincoeff1mEditField
    function gaincoeff1mEditFieldValueChanged(new_value, event) %#ok
        app.gaincoeff1mEditField = new_value;
        
        app.GaindBmEditField.Value = app.gaincoeff1mEditField.Value*(10*app.FiberlengthmEditField.Value)/log(10); % dB/m
    end

    % Value changed function: UpperstatelifetimesEditField
    function UpperstatelifetimesEditFieldValueChanged(new_value, event) %#ok
        app.UpperstatelifetimesEditField = new_value;
        
        h = 6.626e-34; % Planck constant

        app.SaturationintensityJm2EditField.Value = h*app.c/(app.CenterwavelengthnmEditField.Value*1e-9)/app.Totalcrosssectionm2emissionabsorptionEditField.Value/app.UpperstatelifetimesEditField.Value/app.LaserrepetitionrateHzEditField.Value; % J/m^2
        app.SaturationenergynJEditField.Value = app.SaturationintensityJm2EditField.Value*app.Aeff*1e9; % nJ
    end

    % Value changed function: LaserrepetitionrateHzEditField
    function LaserrepetitionrateHzEditFieldValueChanged(new_value, event) %#ok
        app.LaserrepetitionrateHzEditField = new_value;
        
        h = 6.626e-34; % Planck constant

        app.SaturationintensityJm2EditField.Value = h*app.c/(app.CenterwavelengthnmEditField.Value*1e-9)/app.Totalcrosssectionm2emissionabsorptionEditField.Value/app.UpperstatelifetimesEditField.Value/app.LaserrepetitionrateHzEditField.Value; % J/m^2
        app.SaturationenergynJEditField.Value = app.SaturationintensityJm2EditField.Value*app.Aeff*1e9; % nJ
    end

    % Value changed function: 
    % Totalcrosssectionm2emissionabsorptionEditField
    function Totalcrosssectionm2emissionabsorptionEditFieldValueChanged(new_value, event) %#ok
        app.Totalcrosssectionm2emissionabsorptionEditField = new_value;
        
        h = 6.626e-34; % Planck constant

        app.SaturationintensityJm2EditField.Value = h*app.c/(app.CenterwavelengthnmEditField.Value*1e-9)/app.Totalcrosssectionm2emissionabsorptionEditField.Value/app.UpperstatelifetimesEditField.Value/app.LaserrepetitionrateHzEditField.Value; % J/m^2
        app.SaturationenergynJEditField.Value = app.SaturationintensityJm2EditField.Value*app.Aeff*1e9; % nJ
    end

    % Value changed function: SavingstepsizemEditField
    function SavingstepsizemEditFieldValueChanged(new_value, event) %#ok
        app.SavingstepsizemEditField = new_value;
        
        if app.SavingstepsizemEditField.Value ~= 0
            reset_saving_period = false;
            num_saves_total = app.FiberlengthmEditField.Value/app.SavingstepsizemEditField.Value;
            if rem(num_saves_total,1) && rem(num_saves_total+eps(num_saves_total),1) && rem(num_saves_total-eps(num_saves_total),1)
                errordlg(sprintf('The save period is %4.2f m and the fiber length is %4.2f m, which are not commensurate', app.SavingstepsizemEditField.Value,app.FiberlengthmEditField.Value),'Saving step size Incommemsurate Error','replace');
                reset_saving_period = true; % reset
            end

            if isequal(app.AdaptivealgorithmButtonGroup.SelectedObject.Text,'fixed step size')
                if isequal(app.SteppingalgorithmDropDown.Value,'MPA')
                    large_step = app.StepsizemEditField.Value*app.NumberofparallelizationsEditField.Value;
                else
                    large_step = app.StepsizemEditField.Value;
                end
                num_steps_persave = app.SavingstepsizemEditField.Value/large_step;
                if rem(num_steps_persave,1) && rem(num_steps_persave+eps(num_steps_persave),1) && rem(num_steps_persave-eps(num_steps_persave),1)
                    errordlg(sprintf('The step size is %4.2f m and the save period is %4.2f m, which are not commensurate', large_step,app.SavingstepsizemEditField.Value),'Saving step size Incommemsurate Error','replace');
                    reset_saving_period = true;
                end
            end

            if reset_saving_period
                app.SavingstepsizemEditField.Value = 0; % reset to 0
            end
        end
        
        maxdeltaZEditFieldValueChanged(app.maxdeltaZEditField); % check the max adaptive step size again
    end

    % Value changed function: EllipticityEditField
    function EllipticityEditFieldValueChanged(new_value, event) %#ok
        app.EllipticityEditField = new_value;
        
        if isequal(app.ScalarPolarizedsimulationsButtonGroup.SelectedObject.Text,'scalar')
            if app.EllipticityEditField.Value ~= 0 && app.EllipticityEditField.Value ~= 1
                errordlg('Ellipticity can only be 0 (linearly polarized) or 1 (circularly polarized) in scalar computations.','Ellipticity Error','replace');
                app.EllipticityEditField.Value = 0; % reset to linear polarization
            end
        end
    end

    % Selection changed function: TimewindowBetasButtonGroup
    function TimewindowBetasButtonGroupSelectionChanged(new_value, event) %#ok
        app.TimewindowBetasButtonGroup = new_value;
        
        switch app.TimewindowBetasButtonGroup.SelectedObject.Text
            case 'Beta0'
                app.ValueEditField.Value = app.sim_betas(1);
                app.mLabel_2.Text = '1/m';
            case 'Beta1'
                app.ValueEditField.Value = app.sim_betas(2);
                app.mLabel_2.Text = 'ps/m';
        end
    end

    % Selection changed function: 
    % ScalarPolarizedsimulationsButtonGroup
    function ScalarPolarizedsimulationsButtonGroupSelectionChanged(new_value, event) %#ok
        app.ScalarPolarizedsimulationsButtonGroup = new_value;
        
        if app.singlemodeButton.Value
            switch app.ScalarPolarizedsimulationsButtonGroup.SelectedObject.Text
                case 'scalar'
                    app.YEditFieldLabel.Enable = 'off';
                    app.YEditField.Enable = 'off';

                    if app.EllipticityEditField.Value ~= 0 && app.EllipticityEditField.Value ~= 1
                        errordlg('Ellipticity can only be 0 (linearly polarized) or 1 (circularly polarized) in scalar computations.','Ellipticity Error','replace');
                        app.EllipticityEditField.Value = 0; % reset to linear polarization
                    end
                case 'polarized'
                    app.YEditFieldLabel.Enable = 'on';
                    app.YEditField.Enable = 'on';
            end
        end
    end

    % Selection changed function: AdaptivealgorithmButtonGroup
    function AdaptivealgorithmButtonGroupSelectionChanged(new_value, event) %#ok
        app.AdaptivealgorithmButtonGroup = new_value;
        
        switch app.AdaptivealgorithmButtonGroup.SelectedObject.Text
            case 'adaptive step size'
                app.StepsizemEditFieldLabel.Enable = 'off';
                app.StepsizemEditField.Enable = 'off';
                app.ThresholdEditFieldLabel.Enable = 'on';
                app.ThresholdEditField.Enable = 'on';
                app.maxdeltaZEditField.Enable = 'on';
                app.maxdeltaZEditFieldLabel.Enable = 'on';

                app.SteppingalgorithmDropDown.Value = 'RK4IP'; % Adaptive algorithm supports only RK4IP.
                app.NumberofparallelizationsEditFieldLabel.Enable = 'off';
                app.NumberofparallelizationsEditField.Enable = 'off';
                app.MaxiterationsEditFieldLabel.Enable = 'off';
                app.MaxiterationsEditField.Enable = 'off';
                app.MiniterationsEditFieldLabel.Enable = 'off';
                app.MiniterationsEditField.Enable = 'off';
                app.ToleranceEditFieldLabel.Enable = 'off';
                app.ToleranceEditField.Enable = 'off';
            case 'fixed step size'
                app.StepsizemEditFieldLabel.Enable = 'on';
                app.StepsizemEditField.Enable = 'on';
                app.ThresholdEditFieldLabel.Enable = 'off';
                app.ThresholdEditField.Enable = 'off';
                app.maxdeltaZEditField.Enable = 'off';
                app.maxdeltaZEditFieldLabel.Enable = 'off';
        end
    end

    % Value changed function: UseGPUButton
    function UseGPUButtonValueChanged(new_value, event) %#ok
        app.UseGPUButton = new_value;
        
        if app.UseGPUButton.Value
            try
                gpuDevice(); % check if it has an available GPU

                app.GPUindextouseEditFieldLabel.Enable = 'on';
                app.GPUindextouseEditField.Enable = 'on';
                app.CudafolderEditFieldLabel.Enable = 'on';
                app.CudafolderEditField.Enable = 'on';
            catch
                app.UseGPUButton.Value = false;

                app.GPUindextouseEditFieldLabel.Enable = 'off';
                app.GPUindextouseEditField.Enable = 'off';
                app.CudafolderEditFieldLabel.Enable = 'off';
                app.CudafolderEditField.Enable = 'off';
            end
        else
            app.GPUindextouseEditFieldLabel.Enable = 'off';
            app.GPUindextouseEditField.Enable = 'off';
            app.CudafolderEditFieldLabel.Enable = 'off';
            app.CudafolderEditField.Enable = 'off';
        end
    end

    % Value changed function: ShowprogressbarButton
    function ShowprogressbarButtonValueChanged(new_value, event) %#ok
        app.ShowprogressbarButton = new_value;
        
        if app.ShowprogressbarButton.Value
            app.ProgressbarnameEditFieldLabel.Enable = 'on';
            app.ProgressbarnameEditField.Enable = 'on';
        else
            app.ProgressbarnameEditFieldLabel.Enable = 'off';
            app.ProgressbarnameEditField.Enable = 'off';
        end
    end

    % Value changed function: SteppingalgorithmDropDown
    function SteppingalgorithmDropDownValueChanged(new_value, event) %#ok
        app.SteppingalgorithmDropDown = new_value;
        
        switch app.AdaptivealgorithmButtonGroup.SelectedObject.Text
            case 'adaptive step size'
                app.SteppingalgorithmDropDown.Value = 'RK4IP';

                warndlg('Adaptive-step algorithm supports only RK4IP.','Stepping algorithm','replace');
            case 'fixed step size'
                if isequal(app.SteppingalgorithmDropDown.Value,'MPA')
                    app.NumberofparallelizationsEditFieldLabel.Enable = 'on';
                    app.NumberofparallelizationsEditField.Enable = 'on';
                    app.MaxiterationsEditFieldLabel.Enable = 'on';
                    app.MaxiterationsEditField.Enable = 'on';
                    app.MiniterationsEditFieldLabel.Enable = 'on';
                    app.MiniterationsEditField.Enable = 'on';
                    app.ToleranceEditFieldLabel.Enable = 'on';
                    app.ToleranceEditField.Enable = 'on';
                end
        end
    end

    % Value changed function: GainmodelDropDown
    function GainmodelDropDownValueChanged(new_value, event) %#ok
        app.GainmodelDropDown = new_value;
        
        if ismember(app.GainmodelDropDown.Value,1:3)
            app.GaindBmEditFieldLabel.Enable = 'on';
            app.GaindBmEditField.Enable = 'on';
            app.gaincoeff1mEditFieldLabel.Enable = 'on';
            app.gaincoeff1mEditField.Enable = 'on';
            app.FWHMnmEditFieldLabel.Enable = 'on';
            app.FWHMnmEditField.Enable = 'on';
            switch app.GainmodelDropDown.Value
                case 1
                    app.SaturationenergynJEditFieldLabel.Enable = 'on';
                    app.SaturationenergynJEditField.Enable = 'on';

                    app.SaturationintensityJm2EditFieldLabel.Enable = 'off';
                    app.SaturationintensityJm2EditField.Enable = 'off';
                case {2,3}
                    app.SaturationintensityJm2EditFieldLabel.Enable = 'on';
                    app.SaturationintensityJm2EditField.Enable = 'on';

                    app.SaturationenergynJEditFieldLabel.Enable = 'off';
                    app.SaturationenergynJEditField.Enable = 'off';
            end
            app.UpperstatelifetimesEditFieldLabel.Enable = 'on';
            app.UpperstatelifetimesEditField.Enable = 'on';
            app.LaserrepetitionrateHzEditFieldLabel.Enable = 'on';
            app.LaserrepetitionrateHzEditField.Enable = 'on';
            app.Totalcrosssectionm2emissionabsorptionEditFieldLabel.Enable = 'on';
            app.Totalcrosssectionm2emissionabsorptionEditField.Enable = 'on';
        else
            app.GaindBmEditFieldLabel.Enable = 'off';
            app.GaindBmEditField.Enable = 'off';
            app.gaincoeff1mEditFieldLabel.Enable = 'off';
            app.gaincoeff1mEditField.Enable = 'off';
            app.FWHMnmEditFieldLabel.Enable = 'off';
            app.FWHMnmEditField.Enable = 'off';
            app.SaturationenergynJEditFieldLabel.Enable = 'off';
            app.SaturationenergynJEditField.Enable = 'off';
            app.SaturationintensityJm2EditFieldLabel.Enable = 'off';
            app.SaturationintensityJm2EditField.Enable = 'off';
            app.UpperstatelifetimesEditFieldLabel.Enable = 'off';
            app.UpperstatelifetimesEditField.Enable = 'off';
            app.LaserrepetitionrateHzEditFieldLabel.Enable = 'off';
            app.LaserrepetitionrateHzEditField.Enable = 'off';
            app.Totalcrosssectionm2emissionabsorptionEditFieldLabel.Enable = 'off';
            app.Totalcrosssectionm2emissionabsorptionEditField.Enable = 'off';
        end
    end

    % Value changed function: GPUindextouseEditField
    function GPUindextouseEditFieldValueChanged(new_value, event) %#ok
        app.GPUindextouseEditField = new_value;
        
        num_GPU = gpuDeviceCount();

        if app.GPUindextouseEditField.Value == 0 && num_GPU ~= 0
            errordlg('It starts from 1.','GPU index Error','replace');
            app.GPUindextouseEditField.Value = 1;
        end
        if app.GPUindextouseEditField.Value > num_GPU
            errordlg(sprintf('It needs to be smaller than or equal to %u.',num_GPU),'GPU index Error','replace');
            app.GPUindextouseEditField.Value = num_GPU;
        end
    end

    % Value changed function: ValueEditField
    function ValueEditFieldValueChanged(new_value, event) %#ok
        app.ValueEditField = new_value;
        
        switch app.BetasButtonGroup.SelectedObject.Text
            case 'Beta0'
                app.sim_betas(1) = app.ValueEditField.Value;
            case 'Beta1'
                app.sim_betas(2) = app.ValueEditField.Value;
        end
    end

    % Value changed function: FibermaterialDropDown
    function FibermaterialDropDownValueChanged(new_value, event) %#ok
        app.FibermaterialDropDown = new_value;
        
        switch app.FibermaterialDropDown.Value
            case 'silica'
                app.Nonlinearrefractiveindexn2m2WEditField.Value = 2.3e-20; % m^2/W
            case 'calcogenide'
                app.Nonlinearrefractiveindexn2m2WEditField.Value = 1.7e-18; % m^2/W
            case 'ZBLAN'
                app.Nonlinearrefractiveindexn2m2WEditField.Value = 2.1e-20; % m^2/W
        end
    end

    % Value changed function: StepsizemEditField
    function StepsizemEditFieldValueChanged(new_value, event) %#ok
        app.StepsizemEditField = new_value;
        
        if isequal(app.SteppingalgorithmDropDown.Value,'MPA')
            large_step = app.NumberofparallelizationsEditField.Value*app.StepsizemEditField.Value;
        else
            large_step = app.StepsizemEditField.Value;
        end
        num_steps_total = app.FiberlengthmEditField.Value/large_step;
        if rem(num_steps_total,1) && rem(num_steps_total+eps(num_steps_total),1) && rem(num_steps_total-eps(num_steps_total),1)
            errordlg(sprintf('The step size is %4.2f m and the fiber length is %4.2f m, which are not commensurate', large_step,app.FiberlengthmEditField.Value),'Step size Incommemsurate Error','replace');
            app.StepsizemEditField.Value = app.FiberlengthmEditField.Value/round(num_steps_total);
            if isequal(app.SteppingalgorithmDropDown.Value,'MPA')
                app.StepsizemEditField.Value = app.StepsizemEditField.Value/app.NumberofparallelizationsEditField.Value;
            end
        end
    end

    % Value changed function: NumberofparallelizationsEditField
    function NumberofparallelizationsEditFieldValueChanged(new_value, event) %#ok
        app.NumberofparallelizationsEditField = new_value;
        
        StepsizemEditFieldValueChanged(app.StepsizemEditField);
        SavingstepsizemEditFieldValueChanged(app.SavingstepsizemEditField);
    end

    % Value changed function: FiberlengthmEditField
    function FiberlengthmEditFieldValueChanged(new_value, event) %#ok
        app.FiberlengthmEditField = new_value;
        
        SavingstepsizemEditFieldValueChanged(app.SavingstepsizemEditField);
        StepsizemEditFieldValueChanged(app.StepsizemEditField);
        maxdeltaZEditFieldValueChanged(app.maxdeltaZEditField); % check the max adaptive step size again

        GaindBmEditFieldValueChanged(app.GaindBmEditField);
        gaincoeff1mEditFieldValueChanged(app.gaincoeff1mEditField);
    end

    % Value changed function: maxdeltaZEditField
    function maxdeltaZEditFieldValueChanged(new_value, event) %#ok
        app.maxdeltaZEditField = new_value;
        
        if app.SavingstepsizemEditField.Value == 0
            save_step_size = app.FiberlengthmEditField.Value;
        else
            save_step_size = app.SavingstepsizemEditField.Value;
        end

        % The saved step size can't be too small; otherwise, it'll skip
        % saving the data and send an error.
        max_maxDeltaZ = save_step_size/10;

        if app.maxdeltaZEditField.Value > max_maxDeltaZ || app.maxdeltaZEditField.Value == 0
            app.maxdeltaZEditField.Value = max_maxDeltaZ;
        end
    end

    %% Global functions
    function read_betas()
        % Find the default betas

        switch app.SimulationtypeButtonGroup.SelectedObject.Text
            case 'single mode'
                DefaultButtonPushed();
            case 'multimode'
                try
                    % Load betas from the specified file
                    loaded_variables = load(fullfile(app.MM_folder,app.BetasfilenameEditField.Value),'betas'); % in fs^n/mm
                    if max(app.midx) > size(loaded_variables.betas,2)
                        ME = MException('read_betas:midxWrongValue',sprintf('Mode index can''t be larger than the mode number in "%s".',app.BetasfilenameEditField.Value));
                        throw(ME);
                    end
                    app.betas = loaded_variables.betas;

                    unit_conversion = 0.001.^(-1:size(app.betas, 1)-2)'; % The imported values are in fs^n/mm, but the simulation uses ps^n/m
                    app.betas = app.betas(:,app.midx).*unit_conversion;
                catch ME
                    switch ME.identifier
                        case 'read_betas:midxWrongValue'
                            errordlg(ME.message,'Load Betas Error','replace');
                        case 'MATLAB:load:couldNotReadFile'
                            errordlg(sprintf('Please check if MM_folder contains "betas" variable in the "%s".',app.BetasfilenameEditField.Value),'Load Betas Error','replace');
                    end
                end
        end

        % Update time-window betas
        app.sim_betas = app.betas([1,2],1);
        switch app.TimewindowBetasButtonGroup.SelectedObject.Text
            case 'Beta0'
                app.ValueEditField.Value = app.sim_betas(1);
            case 'Beta1'
                app.ValueEditField.Value = app.sim_betas(2);
        end
    end
    function read_SR()
        switch app.SimulationtypeButtonGroup.SelectedObject.Text
            case 'single mode'
                MFDumEditField_2ValueChanged(app.MFDumEditField_2);
            case 'multimode'
                try
                    % Load SR from the specified file
                    loaded_variables = load(fullfile(app.MM_folder,app.StensorsfilenameEditField.Value),'SR'); % in 1/m^2
                    if max(app.midx) > max(size(loaded_variables.SR))
                        ME = MException('read_SR:midxWrongValue',sprintf('Mode index can''t be larger than the mode number in "%s".',app.StensorsfilenameEditField.Value));
                        throw(ME);
                    end
                    app.SR = loaded_variables.SR(app.midx,app.midx,app.midx,app.midx);

                    app.Aeff = 1/app.SR(1,1,1,1); % m^2
                    
                    % Update gain values
                    UpperstatelifetimesEditFieldValueChanged(app.UpperstatelifetimesEditField);
                    LaserrepetitionrateHzEditFieldValueChanged(app.LaserrepetitionrateHzEditField);
                    Totalcrosssectionm2emissionabsorptionEditFieldValueChanged(app.Totalcrosssectionm2emissionabsorptionEditField);
                catch ME
                    switch ME.identifier
                        case 'read_SR:midxWrongValue'
                            errordlg(ME.message,'Load SR Error','replace');
                        case 'MATLAB:load:couldNotReadFile'
                            errordlg(sprintf('Please check if MM_folder contains "SR" variable in the "%s".',app.StensorsfilenameEditField.Value),'Load SR Error','replace');
                    end
                end
        end
    end

end