function show_binary_mask_transducer(karray, kgrid, parameters)
% SHOW_BINARY_MASK_TRANSDUCER Visualizes the binary mask of a transducer in the simulation grid.
%
% The resulting figure is saved to the debug directory defined in the
% simulation parameters for inspection and debugging purposes.
%
% Input:
%   karray     - kWaveArray object defining the transducer geometry.
%   kgrid      - kWaveGrid object describing the simulation grid.
%   parameters - Struct containing simulation parameters.
%
% Output:
%   (none) The generated visualization is saved to disk.

    % Generate binary mask of transducer elements on the simulation grid
    source_mask = karray.getArrayBinaryMask(kgrid);

    % Construct base filename for saved figures
    base_filename = sprintf('sub-%03d_%s_binary_transducer_mask%s', ...
        parameters.subject_id, parameters.simulation_medium, ...
        parameters.results_filename_affix);

    if parameters.n_sim_dims == 3

        % -----------------------------------------------------------------
        % 3D visualization using isosurface
        % -----------------------------------------------------------------

        h = figure('Name', '3D Transducer Binary Mask');
        hold on;

        % Extract isosurface from binary mask
        p = patch(isosurface(source_mask, 0.5));  % 0.5 is the threshold for binary mask
        isonormals(source_mask, p);  % Add surface normals for better visualization

        set(p, 'FaceColor', 'r', 'EdgeColor', 'none');

        camlight; 
        lighting gouraud;

        axis equal;
        xlabel('x-grid'); ylabel('y-grid'); zlabel('z-grid');
        title('3D Visualization of Transducer Binary Mask');

        % Save figure
        saveas(h, fullfile(parameters.debug_dir, [base_filename '.fig']));
        saveas(h, fullfile(parameters.debug_dir, [base_filename '.png']));

        close(h);
    else
        % -----------------------------------------------------------------
        % 2D visualization
        % -----------------------------------------------------------------

        h = figure('Name', 'Transducer in Simulation Grid');

        imagesc(source_mask);
        axis equal tight;
        colormap(hot);
        colorbar;
        title('Transducer Binary Mask in Simulation Grid');

        % Save figure
        saveas(h, fullfile(parameters.debug_dir, [base_filename '.png']));

        close(h);
    end
end