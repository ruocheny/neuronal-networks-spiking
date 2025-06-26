clear all; close all; clc;

rpt_num = 1; % number of trials per setting
T = 500000; input_type = "lognormal";
input_strength_list = [10000];
noise_scalar_list = [0.6];
amp_EE_list = [0.07];
indir = './spikes_multi_input_exponential';

%% params for network
Ne_xy_dim = 30;           Ni_xy_dim = 15;
Ne = Ne_xy_dim^2;         Ni = Ni_xy_dim^2;
isWeighted = 1; weightScalar = 32;
sigma_EE = 10; sigma_IE = 10; sigma_EI = 10; sigma_II = 10;
amp_EI = 0.27; amp_II = 1.08;

%% params for firing model
input_strength_scalar_E = 5; input_strength_scalar_I = 2;
input_loc_x_start = 6; input_loc_x_end = 25;
input_loc_y_start = 6; input_loc_y_end = 25;


%% iterate over all settings
idx_rpt = 1;

for idx_input_strength = 1:length(input_strength_list)
    for idx_amp_EE = 1:length(amp_EE_list)
        for idx_noise_scalar = 1:length(noise_scalar_list)
        
            
            %% dataname
            input_strength = input_strength_list(idx_input_strength);
            noise_scalar = noise_scalar_list(idx_noise_scalar);
            amp_EE = amp_EE_list(idx_amp_EE);
            amp_IE = 4*amp_EE;

            network_filename = sprintf('net_isW_%d_%d_amp_EE_%d_EI_%d_IE_%d_II_%d', ...
                isWeighted, weightScalar, amp_EE*1000, amp_EI*1000, amp_IE*1000, amp_II*1000);
            
            input_noise_scalar_E = 5*noise_scalar;
            input_noise_scalar_I = 2*noise_scalar;
            
            spike_filename = sprintf('exp_spike_T_%d_signal_%d_%d_%d_X_%d_%d_Y_%d_%d_noise_%d_%d', ...
                T, input_strength, ...
                input_strength_scalar_E, input_strength_scalar_I,...
                input_loc_x_start, input_loc_x_end, input_loc_y_start, input_loc_y_end, ...
                input_noise_scalar_E*10, input_noise_scalar_I*10);
            spike_filename_alter = sprintf('exp_spike_T_%d_signal_%d_%d_%d_X_%d_%d_Y_%d_%d_noise_%.6e_%.6e', ...
                T, input_strength, ...
                input_strength_scalar_E, input_strength_scalar_I,...
                input_loc_x_start, input_loc_x_end, input_loc_y_start, input_loc_y_end, ...
                input_noise_scalar_E*10, input_noise_scalar_I*10);
            
            
            %% load ISI
            
            load(sprintf('%s/%s_%s_rpt_%d.mat',indir, spike_filename,network_filename,idx_rpt),...
                    'firings', 'spiking_time_length', 'y');
            
            
            % convert firings to time series
            % spiking_time{i} is the spiking time series for the i-th neuron
            spiking_time = cell(Ne+Ni,1);
            for row = 1:size(firings,1)
                spiking_time{firings(row,2)} = [spiking_time{firings(row,2)}, firings(row,1)];
            end
            
            ISI = cell(Ne+Ni,1);
            for u = 1:Ne+Ni
                ISI{u,1} = diff(spiking_time{u,1});
            end
            
            %% MFDFA
            for u=1:Ne+Ni
                [MFDFA_Fq{u}, ...
                    MFDFA_Hq{u}, ...
                    MFDFA_qRegLine{u}, ...
                    MFDFA_scale{u}, ...
                    MFDFA_tq{u}, ...
                    MFDFA_hq{u}, ...
                    MFDFA_Dq{u}]...
                    = MFDFA(ISI{u,1});
            end
            
            
            save(sprintf('MFDFA_multi_input_exponential/%s_%s.mat', spike_filename, network_filename), ...
                'MFDFA_Fq', 'MFDFA_Hq', 'MFDFA_qRegLine', 'MFDFA_scale', ...
                'MFDFA_tq', 'MFDFA_hq', 'MFDFA_Dq', 'spiking_time_length');
        end
    end
end




