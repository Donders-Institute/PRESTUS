function positions = fibonacci_spherical_cap(N, R, D)
    if nargin == 0
        % Example/demo code
        N = 64;
        R = 75;
        D = 55;
        positions = fibonacci_spherical_cap(N, R, D);

        % Plotting
        figure;
        scatter3(positions(:,1), positions(:,2), positions(:,3), 'filled');
        xlabel('x'); ylabel('y'); zlabel('z');
        title('Fibonacci Spiral on Spherical Cap');
        axis equal; grid on;
        return
    end
    % N  = number of elements
    % R  = radius of curvature of the bowl (mm)
    % D  = aperture diameter (mm)

    % Aperture radius
    a = D / 2;

    % Golden angle in radians
    ga = pi * (3 - sqrt(5));

    % Preallocate
    positions = zeros(N, 3);

    for i = 1:N
        r = a * sqrt((i - 0.5) / N);  % Uniform density in aperture
        theta = ga * (i - 1);

        x = r * cos(theta);
        y = r * sin(theta);

        if (x^2 + y^2) > a^2
            continue  % Skip points outside aperture
        end

        % Project onto spherical cap
        z = R - sqrt(R^2 - x^2 - y^2);

        positions(i, :) = [x, y, z];
    end

    % Remove unused (zeroed) rows if any were skipped
    positions = positions(any(positions, 2), :);

    N = size(positions, 1);
    avg_dist = 0;
    for i = 1:N
        dists = vecnorm(positions - positions(i,:), 2, 2);  % Euclidean distances
        dists(i) = inf;  % Ignore self
        avg_dist = avg_dist + min(dists);  % Nearest neighbor distance
    end
    avg_pitch = avg_dist / N;  % Average center-to-center spacing
    
    % Define kerf and compute element size
    kerf = 0.5;  % in mm
    element_size = avg_pitch - kerf;
    
    fprintf('Estimated avg pitch: %.3f mm\n', avg_pitch);
    fprintf('Recommended element size (with kerf %.3f mm): %.3f mm\n', kerf, element_size);
end

