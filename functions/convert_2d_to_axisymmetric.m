function [kgrid_AS, medium_AS, source_AS] = convert_2d_to_axisymmetric(kgrid, medium, source)

% Convert the simulation grid, medium, and source to a setup symmetric about x=0.
% The x-dimension is halved and symmetrized, y-dimension remains the same.

    % Original 2D grid parameters
    Nx = kgrid.Nx;
    Ny = kgrid.Ny;
    dx = kgrid.dx;
    dy = kgrid.dy;
    
    % Original x and y vectors (assuming centered around zero)
    x_vec = linspace(-dx*(Nx/2), dx*(Nx/2-1), Nx);
    y_vec = linspace(0, dy*(Ny-1), Ny);  % y assumed positive or unchanged
    
    % --- Step 1: Define new grid symmetric about x=0 ---
    
    % Find index of x=0 (or closest)
    [~, zero_idx] = min(abs(x_vec));
    
    % Extract positive half of x_vec including zero
    x_vec_AS = x_vec(zero_idx:end);
    Nx_AS = length(x_vec_AS);
    
    % y dimension remains the same
    y_vec_AS = y_vec;
    Ny_AS = length(y_vec_AS);
    
    % Create new kWaveGrid with halved x and full y
    kgrid_AS = kWaveGrid(Nx_AS, dx, Ny_AS, dy);
    kgrid_AS.dt = kgrid.dt;
    kgrid_AS.Nt = kgrid.Nt;
    kgrid_AS.t_array = kgrid.t_array;
    
    % --- Step 2: Symmetrize medium and source data about x=0 ---
    
    fields_to_symmetrize = {'sound_speed', 'density', 'alpha_coeff'};
    
    for i = 1:length(fields_to_symmetrize)
        field = fields_to_symmetrize{i};
        
        % Extract negative and positive halves of the field along x (rows)
        neg_half = medium.(field)(1:zero_idx-1, :); % x < 0
        pos_half = medium.(field)(zero_idx:end, :); % x >= 0
        
        % Flip negative half upside down to mirror about x=0
        neg_half_flipped = flipud(neg_half);
        
        % Average the two halves to enforce symmetry about x=0
        medium_AS.(field) = (pos_half + neg_half_flipped) / 2;
    end
    
    % --- Step 3: Adjust source similarly ---
    
    % For source.p (pressure field), symmetrize along x (rows) (if asymmetric, this may increase transducer by 1 voxel)
    zero_idx_p = ceil(size(source.p,1)/2);
    neg_half_p = source.p(1:zero_idx_p, :);
    source_AS.p = neg_half_p;
    
    % For source.p_mask (logical mask), do not combine halves (in case they are not symmetric)
    neg_half_mask = source.p_mask(1:zero_idx-1, :);
    source_AS.p_mask = neg_half_mask;
    
end