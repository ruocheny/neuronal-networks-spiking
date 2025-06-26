function S = generateNetwork(Ne, Ni, Ne_xy_dim, Ni_xy_dim, ...
    amp_EE, amp_EI, amp_IE, amp_II, ...
    sigma_EE, sigma_EI, sigma_IE, sigma_II, isWeighted, weightScalar)


%% neuron connection S, S(i,j): directed connection from neuron j to neuron i

% define coordinate of E cells and I cells in the sheet
E_coor = zeros(Ne,2);
cnt = 0;
for i=1:Ne_xy_dim
    for j=1:Ne_xy_dim
        cnt = cnt + 1;
        E_coor(cnt,:) = [i,j];
    end
end

I_loc_scalar = Ne_xy_dim/Ni_xy_dim;
I_coor = zeros(Ni,2);
cnt = 0;
for i=1:Ni_xy_dim
    for j=1:Ni_xy_dim
        cnt = cnt + 1;
        I_coor(cnt,:) = [i*I_loc_scalar,j*I_loc_scalar];
    end
end


S = zeros(Ne+Ni,Ne+Ni); % first Ne rows/columns are excitatory; last Ni are inhibitatory

% E to E (positive; top left block)
for i=1:Ne          % E cells
    for j=1:Ne      % E cells
        distance = pdist2(E_coor(i,:), E_coor(j,:),'euclidean');
        prob = amp_EE * exp(-(distance^2)/(2*sigma_EE^2));
        S(i,j) = generateStrength(prob, distance, isWeighted, weightScalar);
    end
end

% E to I (postive; bottom left block)
for i = Ne+1:Ne+Ni  % I cells
    for j = 1:Ne    % E cells
        distance = pdist2(I_coor(i-Ne,:), E_coor(j,:),'euclidean');
        prob = amp_EI * exp(-(distance^2)/(2*sigma_EI^2));
        S(i,j) = generateStrength(prob, distance, isWeighted, weightScalar);
    end
end


% I to E (the corresponding element in S should be negative;; top right block)
for i = 1:Ne                % E cells
    for j = Ne+1:Ne+Ni      % I cells
        distance = pdist2(E_coor(i,:), I_coor(j-Ne,:),'euclidean');
        prob = amp_IE * exp(-(distance^2)/(2*sigma_IE^2));
        S(i,j) = (-1)*generateStrength(prob, distance, isWeighted, weightScalar);
    end
end

% I to I (what should the sign of the corresponding element be? also negative?; bottom right block)
for i = Ne+1:Ne+Ni          % I cells
    for j = Ne+1:Ne+Ni      % I cells
        distance = pdist2(I_coor(i-Ne,:), I_coor(j-Ne,:),'euclidean');
        prob = amp_II * exp(-(distance^2)/(2*sigma_II^2));
        S(i,j) = (-1)*generateStrength(prob, distance, isWeighted, weightScalar);
    end
end



function strength = generateStrength(prob, distance, isWeighted, weightScalar)

tmp = rand;

if tmp<=prob
    if isWeighted
        strength = 1/(1+distance)*weightScalar;
    else
        strength = 1; % binary case
    end
else
    strength = 0;
end

