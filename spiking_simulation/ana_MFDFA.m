clear all; close all; clc;

rpt_num = 1;
T = 500000; input_type = "lognormal";
input_signal_start_t = 5000; input_signal_interval_t = 50000;
input_strength_list = 5000*[2];
noise_scalar_list = [0.6];
amp_EE_list = [0.07];
indir = './MFDFA_multi_input_exponential';

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

edges = 0:0.01:1;

idx_noise_scalar = 1;
noise_scalar = noise_scalar_list(idx_noise_scalar);
input_noise_scalar_E = 5*noise_scalar;
input_noise_scalar_I = 2*noise_scalar;

q=-5:0.25:5;


%% plot all in one figure

ls = [":", ":", "--", "--", "-", "-"];
cmap = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250]];

figure(1); hold on;
set(gcf,'position',[0,0,600,450]);
set(gca,'FontSize',20)
set(gca,'TickDir','out');

figure(2); hold on;
set(gcf,'position',[0,0,600,450]);
set(gca,'FontSize',20)
set(gca,'TickDir','out');


for idx_amp_EE = 1:length(amp_EE_list)
    for idx_input_strength = 1:length(input_strength_list)
        
        input_strength = input_strength_list(idx_input_strength);
        
        amp_EE = amp_EE_list(idx_amp_EE);
        amp_IE = 4*amp_EE;

        network_filename = sprintf('net_isW_%d_%d_amp_EE_%d_EI_%d_IE_%d_II_%d', ...
            isWeighted, weightScalar, amp_EE*1000, amp_EI*1000, amp_IE*1000, amp_II*1000);
        spike_filename = sprintf('exp_spike_T_%d_signal_%d_%d_%d_X_%d_%d_Y_%d_%d_noise_%d_%d', ...
            T, input_strength, ...
            input_strength_scalar_E, input_strength_scalar_I,...
            input_loc_x_start, input_loc_x_end, input_loc_y_start, input_loc_y_end, ...
            input_noise_scalar_E*10, input_noise_scalar_I*10);
        load(sprintf('%s/%s_%s.mat', indir, spike_filename, network_filename), 'MFDFA_hq','MFDFA_Dq', 'MFDFA_tq','MFDFA_Hq','spiking_time_length');
        
        alpha = [];
        f_alpha = [];
        tau = []; % q-order mass exponent
        H = []; % q-order Hurst exponent
        D = []; % generalized fractal dimension
        for u=1:Ne
            if (~any(isnan(MFDFA_hq{u}))) && (~any(isnan(MFDFA_Dq{u})))
                alpha = [alpha; MFDFA_hq{u}];
                f_alpha = [f_alpha; MFDFA_Dq{u}];
                tau = [tau; MFDFA_tq{u}];
                H = [H;MFDFA_Hq{u}];
                D = [D; MFDFA_tq{u} ./ (q-1)];
            end
        end
        
        
        %% spectrum
        figure(1); 
        x = mean(alpha, 1);
        y = mean(f_alpha, 1);
        
        plot(x(21:end), y(21:end), 'LineStyle', ls(idx_input_strength), 'LineWidth', 3, ...
            'Color', cmap(idx_amp_EE,:), 'DisplayName', ...
            sprintf('\\alpha=%.2f, A=%dK', amp_EE, input_strength/1000));
        

        %% q-order Hurst exponent
        figure(2); 
        x = q;
        y = mean(H, 1);
        ystd = std(H, 1);
        
       plot(x(21:end), y(21:end), 'LineStyle', ls(idx_input_strength), 'LineWidth', 3, ...
                'Color', cmap(idx_amp_EE,:), 'DisplayName', ...
                sprintf('\\alpha=%.2f, A=%dK', amp_EE, input_strength/1000));
       
    end
end

figure(1);
ylim([0.15,1]);
xlabel('q-order singularity exponent');
ylabel('q-order singularity dimension');
legend();

figure(2);
ylim([0.55,1.4]);
xlabel('q');
ylabel('q-order Hurst exponent');
legend();

