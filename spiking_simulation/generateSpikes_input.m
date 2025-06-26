function generateSpikes_input(outdir, network_filename, spike_filename, ...
    num_rpt, T, S, Ne, Ni, Ne_xy_dim, Ni_xy_dim, ...
    input_type, input_strength, input_strength_scalar_E, input_strength_scalar_I, ...
    input_noise_scalar_E, input_noise_scalar_I,...
    input_loc_x_start, input_loc_x_end, input_loc_y_start, input_loc_y_end, ...
    input_signal_t, a, b, c, d)

for idx_rpt = 1:num_rpt
        
    %% input mask
    [mask_E, mask_I] = generateInputMask(Ne, Ni, Ne_xy_dim, Ni_xy_dim, ...
        input_loc_x_start, input_loc_x_end, input_loc_y_start, input_loc_y_end);
    
    %% lognormal input
    pd = makedist('Lognormal','mu',7.5,'sigma',1);

    y = zeros(T,1);
    for i_t = 1:length(input_signal_t)
        t_start = input_signal_t(i_t);
        x_0 = (1:1:T-t_start)';
        y(t_start+1:T) = y(t_start+1:T) + pdf(pd,x_0)*input_strength;
    end
    
    
    %% start simulation
    % code source of Izhikevich model: https://www.izhikevich.org/publications/net.m
    
    v=-65*ones(Ne+Ni,1);
    u=b.*v;
    firings=[];
    
    % firing = [t,fired neuron idx;...]
    % I_t: thalamic input
    
    for t=1:T
        
        fired=find(v>=30); % indices of spikes
        firings=[firings; t+0*fired,fired];
        v(fired)=c(fired);
        u(fired)=u(fired)+d(fired);
        
        if input_type == "random"
            % random thalamic input
            I=[5*randn(Ne,1);2*randn(Ni,1)];
            
        elseif input_type == "lognormal"
            
            % lognormal input
            input_signal_to_E = input_strength_scalar_E * y(t);
            input_signal_to_E = input_signal_to_E .* mask_E;
            input_signal_to_I = input_strength_scalar_I * y(t);
            input_signal_to_I = input_signal_to_I .* mask_I;
            I_signal = [input_signal_to_E;input_signal_to_I];
            
            % noise input
            input_noise_to_E = [input_noise_scalar_E*randn(Ne,1)]; % standard Gaussian noise
            input_noise_to_I = [input_noise_scalar_I*randn(Ni,1)];
            I_noise = [input_noise_to_E;input_noise_to_I];
            
            I = I_signal + I_noise;
        end
        
        I_t(:,t) = I;
        
        I=I+sum(S(:,fired),2); % add the influence from spiked neuron
        
        v=v+0.5*(0.04*v.^2+5*v+140-u+I); % step 0.5 ms
        v=v+0.5*(0.04*v.^2+5*v+140-u+I); % for numerical
        u=u+a.*(b.*v-u); % stability
    end
    
    %% analysis
    
    % convert firings to time series
    % spiking_time{i} is the spiking time series for the i-th neuron
    spiking_time = cell(Ne+Ni,1);
    for row = 1:size(firings,1)
        spiking_time{firings(row,2)} = [spiking_time{firings(row,2)}, firings(row,1)];
    end
    
    spiking_time_length = zeros(Ne+Ni,1);
    for i = 1:Ne+Ni
        spiking_time_length(i) = length(spiking_time{i});
    end
    spiking_time_length_aver_E = mean(spiking_time_length(1:Ne));
    spiking_time_length_aver_I = mean(spiking_time_length(Ne+1:Ne+Ni));
    
    fprintf('E: %f, I: %f \n', spiking_time_length_aver_E, spiking_time_length_aver_I);
    
    save(sprintf('%s/%s_%s_rpt_%d.mat',outdir, spike_filename,network_filename,idx_rpt),...
        'firings', 'spiking_time_length',...
        'spiking_time_length_aver_E', 'spiking_time_length_aver_I', 'y');
    
end
