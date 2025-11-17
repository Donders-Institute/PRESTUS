clear all; close all;    
% Create empty kWaveArray object
    karray = kWaveArray('BLITolerance', 0.1, 'UpsamplingRate', 10, 'BLIType', 'sinc');
    element_height = 6;
    element_width = 6;

    tran_info = readtable('TR64pcd1d65r75_Stanford_fromED.xlsx');
    
    phys_positions = table2array(tran_info(1:64, 2:4));
    theta = table2array(tran_info(1:64, 6));
    phi = table2array(tran_info(1:64, 7));
    rotations = zeros(size(phys_positions, 1), 3);
    extracted_rot = zeros(size(phys_positions, 1), 3);
    for i = 1:size(phys_positions, 1)
        phys_pos = phys_positions(i, :);
        
        % Convert to k-Wave rotations
        tx = 0;
        ty = 90 - theta(i, :);
        tz = phi(i, :);

        rotations(i, :) = calculateTiltTowardsOrigin(phys_pos);
        extracted_rot(i, :) = [tx, ty, tz];

        %rotations(i, :) = [0, 0, 0];

        % Add the rectangular element with position and orientation
        karray.addRectElement(phys_pos, element_height, element_width,  rotations(i, :));
    end

    visualizeArrayWithNormals(phys_positions, rotations, element_height, element_width);
%end

function [tx, ty, tz] = calculateTiltTowardsOrigin(xyz)
    % CALCULATETILTTOWARDSORIGIN Calculates rotation angles (in degrees) needed
    % to make a square at position (x, y, z) have its normal pointing towards origin (0,0,0)
    %
    % Inputs:
    %   x, y, z - coordinates of the square's position
    %
    % Outputs:
    %   tx - rotation around x-axis (in degrees)
    %   ty - rotation around y-axis (in degrees)
    %   tz - rotation around z-axis (in degrees)
    %
    % Assumptions:
    %   - The square initially has its normal pointing along the positive z-axis [0,0,1]
    %   - Rotations are applied in the order: Z, then Y, then X (ZYX order)
    x = xyz(1);
    y = xyz(2);
    z = xyz(3);

    % Calculate distance from position to origin
    distance = sqrt(x^2 + y^2 + z^2);
    
    % Check if position is too close to origin
    if distance < 1e-10
        error('Position is too close to origin. Cannot calculate tilt.');
    end
    
    % Calculate the normalized vector from the position to the origin
    % This is the target normal direction
    nx = -x / distance;
    ny = -y / distance;
    nz = -z / distance;
    
    % Create rotation matrix that transforms [0,0,1] to [nx,ny,nz]
    % This is done using the formula for rotation matrix from one vector to another
    source = [0; 0; 1];  % Initial normal direction
    target = [nx; ny; nz];  % Target normal direction
    
    % If vectors are nearly parallel or opposite, handle specially
    dot_product = dot(source, target);
    if abs(dot_product - 1) < 1e-10
        % Vectors are nearly identical, no rotation needed
        tx = 0;
        ty = 0;
        tz = 0;
        return;
    elseif abs(dot_product + 1) < 1e-10
        % Vectors are nearly opposite, 180-degree rotation around x-axis
        tx = 180;
        ty = 0;
        tz = 0;
        return;
    end
    
    % For all other cases, calculate the rotation matrix
    % Cross product gives rotation axis
    v = cross(source, target);
    s = sqrt(sum(v.^2));  % sin of angle
    
    % Calculate rotation matrix using Rodrigues' rotation formula
    c = dot_product;  % cos of angle
    C = 1 - c;
    
    x = v(1);
    y = v(2);
    z = v(3);
    
    % Normalize rotation axis
    if s > 0
        x = x / s;
        y = y / s;
        z = z / s;
    end
    
    % Rodrigues' rotation formula
    R = [
        x^2*C+c,    x*y*C-z*s,  x*z*C+y*s;
        y*x*C+z*s,  y^2*C+c,    y*z*C-x*s;
        z*x*C-y*s,  z*y*C+x*s,  z^2*C+c
    ];
    
    % Extract Euler angles from the rotation matrix (ZYX convention)
    % This extracts the angles in the order Z, Y, X
    % Check for gimbal lock
    if abs(R(3,1)) >= 1.0
        % Gimbal lock case
        tz = 0;  % Set one angle to zero
        if R(3,1) < 0
            ty = pi/2;
            tx = atan2(R(1,2), R(1,3));
        else
            ty = -pi/2;
            tx = atan2(-R(1,2), -R(1,3));
        end
    else
        ty = -asin(R(3,1));
        tx = atan2(R(3,2)/cos(ty), R(3,3)/cos(ty));
        tz = atan2(R(2,1)/cos(ty), R(1,1)/cos(ty));
    end
    tx = tx * 180/pi;
    ty = ty * 180/pi;
    tz = tz * 180/pi;
end

function visualizeArrayWithNormals(positions, rotations, Lx, Ly)
    % positions: Nx3 array of element positions
    % rotations: Nx3 array of rotation angles [tx, ty, tz]
    % Lx, Ly: dimensions of each element
    
    figure;
    hold on;
    
    % Plot each element
    for i = 1:size(positions, 1)
        pos = positions(i, :);
        rot = rotations(i, :);
        
        % Plot element
        plotElement(pos, rot, Lx, Ly);
        
        % Calculate and plot normal vector
        plotNormal(pos, rot);
    end
    
    % Finalize plot
    axis equal;
    grid on;
    xlabel('x [m]');
    ylabel('y [m]');
    zlabel('z [m]');
    title('Array with Element Normals');
    view(3);
end

function plotElement(pos, rot, Lx, Ly)
    % Convert rotations to radians
    tx = rot(1) * pi/180;
    ty = rot(2) * pi/180;
    tz = rot(3) * pi/180;
    
    % Define corners of unrotated element
    corners = [
        -Lx/2, -Ly/2, 0;
        Lx/2, -Ly/2, 0;
        Lx/2, Ly/2, 0;
        -Lx/2, Ly/2, 0;
        -Lx/2, -Ly/2, 0  % Close the shape
    ];
    
    % Apply rotations
    R_x = [1 0 0; 0 cos(tx) -sin(tx); 0 sin(tx) cos(tx)];
    R_y = [cos(ty) 0 sin(ty); 0 1 0; -sin(ty) 0 cos(ty)];
    R_z = [cos(tz) -sin(tz) 0; sin(tz) cos(tz) 0; 0 0 1];
    
    % Combined rotation matrix (Z-Y-X order)
    R = R_x * R_y * R_z;
    
    % Apply rotation and translation
    rotated_corners = zeros(size(corners));
    for j = 1:size(corners, 1)
        rotated_corners(j, :) = (R * corners(j, :)')' + pos;
    end
    
    % Plot element
    patch(rotated_corners(:, 1), rotated_corners(:, 2), rotated_corners(:, 3), 'b', 'FaceAlpha', 0.3);
end

function plotNormal(pos, rot)
    normal_length = 1;
    % Convert rotations to radians
    tx = rot(1) * pi/180;
    ty = rot(2) * pi/180;
    tz = rot(3) * pi/180;
    
    % Initial normal vector points along Z-axis [0, 0, 1]
    normal = [0; 0; 1];
    
    % Create rotation matrices
    % Note: kWave applies rotations in Z-Y-X order (extrinsic)
    % but we need to apply them in X-Y-Z order when multiplying matrices
    R_x = [1 0 0; 0 cos(tx) -sin(tx); 0 sin(tx) cos(tx)];
    R_y = [cos(ty) 0 sin(ty); 0 1 0; -sin(ty) 0 cos(ty)];
    R_z = [cos(tz) -sin(tz) 0; sin(tz) cos(tz) 0; 0 0 1];
    
    % Combined rotation matrix (Z-Y-X order for extrinsic rotations)
    R = R_x * R_y * R_z;
    
    % Apply rotation to normal
    rotated_normal = R * normal;
    
    % Plot normal vector
    quiver3(pos(1), pos(2), pos(3), ...
        rotated_normal(1) * normal_length, ...
        rotated_normal(2) * normal_length, ...
        rotated_normal(3) * normal_length, ...
        'r', 'LineWidth', 2);
end

%function [transducer_pars, karray] = setup_stanford_transducer_design(transducer_pars, kgrid, trans_pos, focus_pos, karray, water_c0, n_sim_dims, display_tran)
    
