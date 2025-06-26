function [mask_E, mask_I] = generateInputMask(Ne, Ni, Ne_xy_dim, Ni_xy_dim, ...
            input_loc_x_start, input_loc_x_end, input_loc_y_start, input_loc_y_end)

mask_E = zeros(Ne,1);
mask_I = zeros(Ni,1);

for k = 1:Ne
    x = ceil(k/Ne_xy_dim);
    y = mod(k-1, Ne_xy_dim) + 1;
    if (x>=input_loc_x_start)&&(x<=input_loc_x_end) && (y>=input_loc_y_start)&&(y<=input_loc_y_end)
        mask_E(k) = 1;
    end
end

I_loc_scalar = Ne_xy_dim/Ni_xy_dim;
for k = 1:Ni
    x = ceil(k/Ni_xy_dim);
    y = mod(k-1, Ni_xy_dim) + 1;
    x = x * I_loc_scalar;
    y = y * I_loc_scalar;
    if (x>=input_loc_x_start)&&(x<=input_loc_x_end) && (y>=input_loc_y_start)&&(y<=input_loc_y_end)
        mask_I(k) = 1;
    end
end
