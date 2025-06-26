clear all; close all; clc;

outdir = './networks';

%% parameters
% make sure here Ne_xy_dim%Ni_xy_dim = 0 for later computation
% Excitatory neurons      %Inhibitatory neurons
Ne_xy_dim = 30;           Ni_xy_dim = 15;
Ne = Ne_xy_dim^2;         Ni = Ni_xy_dim^2;


%% weighted
isWeighted = 1; weightScalar = 32;
sigma_EE = 10; sigma_IE = 10; sigma_EI = 10; sigma_II = 10;

amp_EI = 0.27; amp_II = 1.08;

amp_EE_list = [0.07];

for amp_EE = amp_EE_list
    for amp_IE = [4*amp_EE]
        
        network_filename = sprintf('net_isW_%d_%d_amp_EE_%d_EI_%d_IE_%d_II_%d', ...
            isWeighted, weightScalar, amp_EE*1000, amp_EI*1000, amp_IE*1000, amp_II*1000);
        
        %% generate network of connection (binary)
        
        S = generateNetwork(Ne, Ni, Ne_xy_dim, Ni_xy_dim, ...
            amp_EE, amp_EI, amp_IE, amp_II, ...
            sigma_EE, sigma_EI, sigma_IE, sigma_II,isWeighted, weightScalar);
        
        % save new network
        save(sprintf('%s/%s.mat', outdir, network_filename),...
            'S', 'Ne', 'Ni',...
            'amp_EE', 'amp_EI', 'amp_IE', 'amp_II',...
            'sigma_EE','sigma_EE','sigma_EE','sigma_EE');
        
    end
end



