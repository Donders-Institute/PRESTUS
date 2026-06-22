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
    else
        error('pipeline:autoFovDiameter', ...
            ['acoustic_fov_diameter_mm=''auto'' requires transducer(1).elem_od_mm. ' ...
             'Set acoustic_fov_diameter_mm explicitly instead.']);
    end
    fprintf('[pipeline] acoustic_fov_diameter_mm auto -> %.1f mm (transducer outer diameter)\n', ...
        parameters.grid.acoustic_fov_diameter_mm);
end

has_acoustic_fov = isfield(parameters.grid, 'acoustic_fov_diameter_mm') && ...
    ~isempty(parameters.grid.acoustic_fov_diameter_mm) && ...
    isnumeric(parameters.grid.acoustic_fov_diameter_mm) && ...
    isfield(parameters.grid, 'acoustic_fov_length_mm') && ...
    ~isempty(parameters.grid.acoustic_fov_length_mm) && ...
    isnumeric(parameters.grid.acoustic_fov_length_mm);

end
