clear all; close all; clc;

indir = './networks';
outdir = './spikes_multi_input_exponential';

num_rpt = 1;

%% params for network
% make sure here Ne_xy_dim%Ni_xy_dim = 0 for later computation
% Excitatory neurons      %Inhibitatory neurons
Ne_xy_dim = 30;           Ni_xy_dim = 15;
Ne = Ne_xy_dim^2;         Ni = Ni_xy_dim^2;

isWeighted = 1; weightScalar = 32;

% standard deviation are chosed to be the same (similar to Reyes code)
sigma_EE = 10; sigma_IE = 10; sigma_EI = 10; sigma_II = 10;

amp_EI = 0.27; amp_II = 1.08;
amp_EE_list = [0.07];


%% params for firing model

T = 500000; input_type = "lognormal";
input_strength_list = [10000];
noise_scalar_list = [0.6];
input_strength_scalar_E = 5; input_strength_scalar_I = 2;
input_loc_x_start = 6; input_loc_x_end = 25;
input_loc_y_start = 6; input_loc_y_end = 25;

% % interval of stimulus is generated from exponential distribution with mean of mu
% % the interval is further bounded within a certain range (40000, 70000)
num_input = 10; exp_mu = 5;
input_signal_t_interval = [exp_mu];
while length(input_signal_t_interval) < num_input
    tmp = exprnd(exp_mu);
    if (4 < tmp) && (tmp < 7)
        input_signal_t_interval = [input_signal_t_interval, tmp];
    end
end
input_signal_t = cumsum(input_signal_t_interval);
input_signal_t = floor(input_signal_t * T/num_input/exp_mu);


%% izhikevich model parameters
% Excitatory neurons      %Inhibitatory neurons
re = rand(Ne,1);          ri=rand(Ni,1);
a = [0.02*ones(Ne,1);     0.02+0.08*ri];
b = [0.2*ones(Ne,1);      0.25-0.05*ri];
c = [-65+15*re.^2;        -65*ones(Ni,1)];
d = [8-6*re.^2;           2*ones(Ni,1)];


%% start simulation

for amp_EE = amp_EE_list
    amp_IE = 4*amp_EE;
    
    network_filename = sprintf('net_isW_%d_%d_amp_EE_%d_EI_%d_IE_%d_II_%d', ...
        isWeighted, weightScalar, amp_EE*1000, amp_EI*1000, amp_IE*1000, amp_II*1000);
    
    load(sprintf('%s/%s.mat', indir, network_filename), 'S');
    
    for input_strength = input_strength_list
        for noise_scalar = noise_scalar_list
            
            input_noise_scalar_E = 5*noise_scalar;
            input_noise_scalar_I = 2*noise_scalar;
            
            spike_filename = sprintf('exp_spike_T_%d_signal_%d_%d_%d_X_%d_%d_Y_%d_%d_noise_%d_%d', ...
                T, input_strength, ...
                input_strength_scalar_E, input_strength_scalar_I,...
                input_loc_x_start, input_loc_x_end, input_loc_y_start, input_loc_y_end, ...
                input_noise_scalar_E*10, input_noise_scalar_I*10);
            
            generateSpikes_input(outdir, network_filename, spike_filename, ...
                num_rpt, T, S, Ne, Ni, Ne_xy_dim, Ni_xy_dim, ...
                input_type, input_strength, input_strength_scalar_E, input_strength_scalar_I, ...
                input_noise_scalar_E, input_noise_scalar_I,...
                input_loc_x_start, input_loc_x_end, input_loc_y_start, input_loc_y_end, ...
                input_signal_t, a, b, c, d)

        end
    end
    
end



