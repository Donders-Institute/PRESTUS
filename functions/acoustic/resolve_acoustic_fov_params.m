function [parameters, has_acoustic_fov] = resolve_acoustic_fov_params(parameters, trans_pos, focus_pos)
% RESOLVE_ACOUSTIC_FOV_PARAMS  Resolve 'auto' FOV settings and check readiness.
%
% Converts string sentinel values in parameters.grid.acoustic_fov_length_mm
% and parameters.grid.acoustic_fov_diameter_mm to numeric mm values, then
% returns whether both fields are present and numeric (has_acoustic_fov).
%
%   'auto' for acoustic_fov_length_mm  ->  focal_depth + grid.z_beyond_focus_mm
%                                          (default z_beyond = 40 mm)
%   'auto' for acoustic_fov_diameter_mm -> max(transducer(1).elem_od_mm)
%
% Use as:
%   [parameters, has_acoustic_fov] = resolve_acoustic_fov_params(parameters, trans_pos, focus_pos)

arguments
    parameters (1,1) struct
    trans_pos  (1,:) {mustBeNumeric}
    focus_pos  (1,:) {mustBeNumeric}
end

is_auto = @(key) isfield(parameters.grid, key) && ischar(parameters.grid.(key)) && ...
    strcmpi(parameters.grid.(key), 'auto');
length_auto   = is_auto('acoustic_fov_length_mm');
diameter_auto = is_auto('acoustic_fov_diameter_mm');

if isfield(parameters.grid, 'acoustic_fov_length_mm') && ...
        ischar(parameters.grid.acoustic_fov_length_mm) && ...
        strcmpi(parameters.grid.acoustic_fov_length_mm, 'auto') && ...
        ~isempty(trans_pos) && ~isempty(focus_pos)
    focal_depth_mm = norm(focus_pos - trans_pos) * parameters.grid.resolution_mm;
    z_beyond_mm    = 40.0;
    if isfield(parameters.grid, 'z_beyond_focus_mm') && ~isempty(parameters.grid.z_beyond_focus_mm)
        z_beyond_mm = parameters.grid.z_beyond_focus_mm;
    end
    parameters.grid.acoustic_fov_length_mm = focal_depth_mm + z_beyond_mm;
    fprintf('[pipeline] acoustic_fov_length_mm auto -> %.1f mm (focal depth %.1f + z_beyond %.1f)\n', ...
        parameters.grid.acoustic_fov_length_mm, focal_depth_mm, z_beyond_mm);
end

if isfield(parameters.grid, 'acoustic_fov_diameter_mm') && ...
        ischar(parameters.grid.acoustic_fov_diameter_mm) && ...
        strcmpi(parameters.grid.acoustic_fov_diameter_mm, 'auto')
    if isfield(parameters.transducer(1), 'elem_od_mm') && ~isempty(parameters.transducer(1).elem_od_mm)
        parameters.grid.acoustic_fov_diameter_mm = max(parameters.transducer(1).elem_od_mm);
    elseif isfield(parameters.transducer(1), 'annular') && isfield(parameters.transducer(1).annular, 'elem_od_mm') && ~isempty(parameters.transducer(1).annular.elem_od_mm)
        parameters.grid.acoustic_fov_diameter_mm = max(parameters.transducer(1).annular.elem_od_mm);
    else
        error('pipeline:autoFovDiameter', ...
            ['acoustic_fov_diameter_mm=''auto'' requires transducer(1).elem_od_mm. ' ...
             'Set acoustic_fov_diameter_mm explicitly instead.']);
    end
    fprintf('[pipeline] acoustic_fov_diameter_mm auto -> %.1f mm (transducer outer diameter)\n', ...
        parameters.grid.acoustic_fov_diameter_mm);
end

% Several transducers: resolve 'auto' for each one from its own geometry, so
% acoustic_grid_fov can size each transducer's part of the FOV box.
n_tx = numel(parameters.transducer);
if n_tx > 1 && (length_auto || diameter_auto) && ...
        isnumeric(parameters.grid.acoustic_fov_length_mm) && ...
        isnumeric(parameters.grid.acoustic_fov_diameter_mm)
    per.length_mm   = repmat(parameters.grid.acoustic_fov_length_mm,   1, n_tx);
    per.diameter_mm = repmat(parameters.grid.acoustic_fov_diameter_mm, 1, n_tx);
    z_beyond_mm = 40.0;
    if isfield(parameters.grid, 'z_beyond_focus_mm') && ~isempty(parameters.grid.z_beyond_focus_mm)
        z_beyond_mm = parameters.grid.z_beyond_focus_mm;
    end
    for ti = 2:n_tx
        tr = parameters.transducer(ti);
        if length_auto && isfield(tr, 'trans_pos') && isfield(tr, 'focus_pos') && ...
                ~isempty(tr.trans_pos) && ~isempty(tr.focus_pos)
            per.length_mm(ti) = norm(tr.focus_pos(:)' - tr.trans_pos(:)') * ...
                parameters.grid.resolution_mm + z_beyond_mm;
        end
        if diameter_auto
            od = [];
            if isfield(tr, 'elem_od_mm') && ~isempty(tr.elem_od_mm)
                od = tr.elem_od_mm;
            elseif isfield(tr, 'annular') && isfield(tr.annular, 'elem_od_mm')
                od = tr.annular.elem_od_mm;
            end
            if ~isempty(od)
                per.diameter_mm(ti) = max(od);
            end
        end
    end
    parameters.state.acoustic_fov_per_tx = per;
    fprintf('[pipeline] acoustic FOV per transducer: length [%s] mm, diameter [%s] mm\n', ...
        num2str(per.length_mm, '%.1f '), num2str(per.diameter_mm, '%.1f '));
end

has_acoustic_fov = isfield(parameters.grid, 'acoustic_fov_diameter_mm') && ...
    ~isempty(parameters.grid.acoustic_fov_diameter_mm) && ...
    isnumeric(parameters.grid.acoustic_fov_diameter_mm) && ...
    isfield(parameters.grid, 'acoustic_fov_length_mm') && ...
    ~isempty(parameters.grid.acoustic_fov_length_mm) && ...
    isnumeric(parameters.grid.acoustic_fov_length_mm);

end
