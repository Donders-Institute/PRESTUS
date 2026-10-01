function tf = should_save_output(io, specific_field)
% SHOULD_SAVE_OUTPUT  Decide whether a particular pipeline output should be saved
%
% Precedence (highest to lowest):
%   1. io.<specific_field>  — per-step flag, when set (non-empty)
%   2. io.save_matrices     — global flag, when set (non-empty); applies only
%                             to intermediate matrix/cache outputs (see below)
%   3. built-in default     — per-output default listed below
%
% Leaving a per-step flag empty ([] in YAML) means "inherit": setting
% io.save_matrices = 1 then saves every intermediate matrix, and 0 saves
% none, while any per-step flag that is set explicitly still wins.
%
% Intermediate outputs controlled by io.save_matrices (built-in default):
%   save_source_matrices   (1)  kwave_source.mat
%   save_acoustic_matrices (0)  acoustic results / cache
%   save_thermal_matrices  (1)  heating_res.mat
%   save_grid_cache        (0)  grid cache
%   save_medium_cache      (0)  medium cache
%   save_matrices          (0)  the global flag itself (e.g. uncertainty cleanup)
%
% Other outputs, not affected by io.save_matrices (built-in default):
%   save_MNI               (1)  MNI-space NIfTIs
%   save_property_maps     (0)  per-property NIfTIs
%   any other field        (1)
%
% Use as:
%   tf = should_save_output(parameters.io, 'save_acoustic_matrices')
%
% Input:
%   io             - parameters.io struct
%   specific_field - name of the per-step flag field
%
% Output:
%   tf - logical scalar; true = save, false = skip
%
% See also: SIMULATION_NIFTI, PRESTUS_PIPELINE

arguments
    io             (1,1) struct
    specific_field (1,:) char
end

    matrix_defaults = struct( ...
        'save_source_matrices',   true, ...
        'save_acoustic_matrices', false, ...
        'save_thermal_matrices',  true, ...
        'save_grid_cache',        false, ...
        'save_medium_cache',      false, ...
        'save_matrices',          false);
    other_defaults = struct( ...
        'save_MNI',               true, ...
        'save_property_maps',     false);

    % 1. explicitly set per-step flag
    if is_set(io, specific_field)
        tf = logical(io.(specific_field));
        return;
    end

    % 2. global flag (intermediate matrices only)
    is_matrix = isfield(matrix_defaults, specific_field);
    if is_matrix && is_set(io, 'save_matrices')
        tf = logical(io.save_matrices);
        return;
    end

    % 3. built-in default
    if is_matrix
        tf = matrix_defaults.(specific_field);
    elseif isfield(other_defaults, specific_field)
        tf = other_defaults.(specific_field);
    else
        tf = true;
    end
end

function tf = is_set(io, field)
    tf = isfield(io, field) && ~isempty(io.(field));
end
