function [kgrid_AS, medium_AS, source_AS] = convert_2d_to_axisymmetric(kgrid, medium, source)

    swap_axes = 0;

    swapXY = @(A) permute(A, [2, 1]);

    % Swap medium and mask if needed so longer dimension is rows
    if swap_axes
        medium.sound_speed = swapXY(medium.sound_speed);
        medium.density = swapXY(medium.density);
        medium.alpha_coeff = swapXY(medium.alpha_coeff);
    end

    medium_AS.sound_speed = medium.sound_speed(:,floor(size(medium.sound_speed,2)/2)+1:end);
    medium_AS.density = medium.density(:,floor(size(medium.density,2)/2)+1:end);
    medium_AS.alpha_coeff = medium.alpha_coeff(:,floor(size(medium.alpha_coeff,2)/2)+1:end);

    %% retain other inputs
    kgrid_AS = kgrid;
    source_AS = source;

% % Convert the simulation grid, medium, and source to axisymmetric setup by halving the shorter dimension.
% % If y-dimension is longer than x-dimension, swap axes first so longer dimension is rows.
% % Then halve the shorter dimension (columns).
% % Finally, swap back to original orientation.
% %
% % Input:
% %   kgrid  - original kWaveGrid object (2D)
% %   medium - structure with 2D fields: sound_speed, density, alpha_coeff
% %   source - structure with p_mask (2D logical) and p (N_points x Nt)
% %
% % Output:
% %   kgrid_AS  - adjusted kWaveGrid with halved shorter dimension
% %   medium_AS - symmetrized medium fields
% %   source_AS - adjusted source with halved p_mask and corresponding p
% 
%     Nx = kgrid.Nx;
%     Ny = kgrid.Ny;
%     dx = kgrid.dx;
%     dy = kgrid.dy;
% 
%     swap_axes = Ny > Nx; % true if y longer than x
% 
%     swapXY = @(A) permute(A, [2, 1]);
% 
%     % Swap medium and mask if needed so longer dimension is rows
%     if swap_axes
%         Nx_tmp = Nx; Nx = Ny; Ny = Nx_tmp;
%         dx_tmp = dx; dx = dy; dy = dx_tmp;
% 
%         medium.sound_speed = swapXY(medium.sound_speed);
%         medium.density = swapXY(medium.density);
%         medium.alpha_coeff = swapXY(medium.alpha_coeff);
% 
%         source.p_mask = swapXY(source.p_mask);
%     end
% 
%     % --- Halve the shorter dimension (columns) ---
% 
%     % Define coordinate vectors for current orientation
%     x_vec = linspace(-dx*(Nx/2), dx*(Nx/2-1), Nx); % longer dim (rows)
%     y_vec = linspace(0, dy*(Ny-1), Ny);           % shorter dim (cols)
% 
%     % Find midpoint index of shorter dimension (columns)
%     half_idx = floor(Ny/2);
% 
%     % Halve the shorter dimension by taking first half columns
%     medium_AS.sound_speed = medium.sound_speed(:, 1:half_idx);
%     medium_AS.density = medium.density(:, 1:half_idx);
%     medium_AS.alpha_coeff = medium.alpha_coeff(:, 1:half_idx);
% 
% 
% 
% %     source_AS.p_mask = source.p_mask(:, 1:half_idx);
% % 
% %     % --- Adjust source.p accordingly ---
% % 
% %     % Find linear indices of 1s in original mask and halved mask
% %     full_mask_ind = find(source.p_mask);
% %     halved_mask_ind = find(source_AS.p_mask);
% % 
% %     % Map halved mask indices to rows in source.p
% %     [~, loc] = ismember(halved_mask_ind, full_mask_ind);
% % 
% %     if any(loc == 0)
% %         error('Some indices in halved mask not found in original mask!');
% %     end
% % 
% %     source_AS.p = source.p(loc, :);
% % 
% %     % --- Create new kWaveGrid with halved shorter dimension ---
% % 
% %     kgrid_AS = kWaveGrid(Nx, dx, half_idx, dy);
% %     kgrid_AS.dt = kgrid.dt;
% %     kgrid_AS.Nt = kgrid.Nt;
% %     kgrid_AS.t_array = kgrid.t_array;
% 
%     % --- Swap back if axes were swapped ---
% 
% %     if swap_axes
% %         medium_AS.sound_speed = swapXY(medium_AS.sound_speed);
% %         medium_AS.density = swapXY(medium_AS.density);
% %         medium_AS.alpha_coeff = swapXY(medium_AS.alpha_coeff);
% % 
% %         source_AS.p_mask = swapXY(source_AS.p_mask);
% %         % source_AS.p remains unchanged (signal matrix)
% %     end
% end
% 
% % function [kgrid_AS, medium_AS, source_AS] = convert_2d_to_axisymmetric(kgrid, medium, source)
% % 
% % % Convert the simulation grid, medium, and source to a setup symmetric about x=0.
% % % The x-dimension is halved and symmetrized, y-dimension remains the same.
% % 
% %     % Original 2D grid parameters
% %     Nx = kgrid.Nx;
% %     Ny = kgrid.Ny;
% %     dx = kgrid.dx;
% %     dy = kgrid.dy;
% %     
% %     % Original x and y vectors (assuming centered around zero)
% %     x_vec = linspace(-dx*(Nx/2), dx*(Nx/2-1), Nx);
% %     y_vec = linspace(0, dy*(Ny-1), Ny);  % y assumed positive or unchanged
% %     
% %     % --- Step 1: Define new grid symmetric about x=0 ---
% %     
% %     % Find index of x=0 (or closest)
% %     [~, zero_idx] = min(abs(x_vec));
% %     
% %     % Extract positive half of x_vec including zero
% %     x_vec_AS = x_vec(zero_idx:end);
% %     Nx_AS = length(x_vec_AS);
% %     
% %     % y dimension remains the same
% %     y_vec_AS = y_vec;
% %     Ny_AS = length(y_vec_AS);
% %     
% %     % Create new kWaveGrid with halved x and full y
% %     kgrid_AS = kWaveGrid(Nx_AS, dx, Ny_AS, dy);
% %     kgrid_AS.dt = kgrid.dt;
% %     kgrid_AS.Nt = kgrid.Nt;
% %     kgrid_AS.t_array = kgrid.t_array;
% %     
% %     % --- Step 2: Symmetrize medium and source data about x=0 ---
% %     
% %     fields_to_symmetrize = {'sound_speed', 'density', 'alpha_coeff'};
% %     
% %     for i = 1:length(fields_to_symmetrize)
% %         field = fields_to_symmetrize{i};
% %         
% %         % Extract negative and positive halves of the field along x (rows)
% %         neg_half = medium.(field)(1:zero_idx-1, :); % x < 0
% %         pos_half = medium.(field)(zero_idx:end, :); % x >= 0
% %         
% %         % Flip negative half upside down to mirror about x=0
% %         neg_half_flipped = flipud(neg_half);
% %         
% %         % Average the two halves to enforce symmetry about x=0
% %         medium_AS.(field) = (pos_half + neg_half_flipped) / 2;
% %     end
% %     
% %     % --- Step 3: Adjust source similarly ---
% %     
% %     % Extract the top half of the source mask
% %     neg_half_mask = source.p_mask(1:zero_idx-1, :);
% %     source_AS.p_mask = neg_half_mask;
% %     
% %     % identify corrspeonding input time series for remaining data
% %     [row_sub, col_sub] = find(neg_half_mask);
% %     orig_row_sub = row_sub;  % because neg_half_mask = source.p_mask(1:zero_idx-1, :)
% %     orig_col_sub = col_sub;
% %     orig_mask_size = size(source.p_mask);
% %     orig_linear_ind = sub2ind(orig_mask_size, orig_row_sub, orig_col_sub);
% %     full_mask_ind = find(source.p_mask);  % linear indices of all 1s in original mask
% %     
% %     [~, loc] = ismember(orig_linear_ind, full_mask_ind);
% %     source_AS.p = source.p(loc, :);
% %     source_AS.p_mask = neg_half_mask;
% %     
% % end