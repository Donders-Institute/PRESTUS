function isppa = babelbrain_insitu_isppa(parameters, sensor_data, kwave_medium, medium_masks)
% BABELBRAIN_INSITU_ISPPA  In-brain peak ISPPA of an acoustic field, as BabelBrain's thermal Isppa input
%
% BabelBrain's BHTE rescales the acoustic field so that max p^2/(2*rho*c)
% inside the brain equals the Isppa it is given (CalculateTemperatureEffects.py,
% AnalyzeLosses). Passing the free-water calibration ISPPA there drops the
% transcranial loss; this function returns the value that makes BabelBrain
% thermal deposit the same field amplitude k-Wave thermal consumes.
%
% Brain-mask precedence matches acoustic_analysis: BabelBrain's own material
% map (sensor_data.babelbrain_material_map) first, then medium_masks. rho*c is
% taken per voxel from kwave_medium when it matches the field grid, otherwise
% from medium_properties.brain.
%
% Use as:
%   isppa = babelbrain_insitu_isppa(parameters, sensor_data, kwave_medium, medium_masks)
%
% Input:
%   parameters   - PRESTUS config (medium_properties used for labels / fallback rho*c)
%   sensor_data  - struct with p_max_all [Pa] (optionally babelbrain_material_map)
%   kwave_medium - struct with sound_speed, density (scalar or grid-sized); may be struct()
%   medium_masks - PRESTUS label volume on the same grid, or []
%
% Output:
%   isppa        - peak in-brain intensity [W/cm²]
%
% See also: BABELBRAIN_THERMAL_WRAPPER, ACOUSTIC_ANALYSIS

arguments
    parameters   (1,1) struct
    sensor_data  (1,1) struct
    kwave_medium (1,1) struct = struct()
    medium_masks = []
end

if ~isfield(sensor_data, 'p_max_all') || isempty(sensor_data.p_max_all)
    error('babelbrain_insitu_isppa:noPressure', ...
        ['sensor_data.p_max_all is missing; cannot derive the in-situ ISPPA for BabelBrain. ' ...
         'Set babelbrain.isppa_w_per_cm2 explicitly.']);
end
p = double(sensor_data.p_max_all);
dims = size(p);

brain = [];
if isfield(sensor_data, 'babelbrain_material_map') && ...
        isequal(size(sensor_data.babelbrain_material_map), dims)
    m = bb_tissuemask_binary(sensor_data.babelbrain_material_map);
    brain = m.brain;
end
if ~any(brain(:)) && ~isempty(medium_masks) && isequal(size(medium_masks), dims)
    m = tissuemask_binary(parameters, medium_masks);
    brain = m.brain;
end
if ~any(brain(:))
    error('babelbrain_insitu_isppa:noBrainMask', ...
        ['No brain voxels found on the acoustic grid (no matching babelbrain_material_map ' ...
         'or medium_masks). BabelBrain normalises heating to the in-brain peak ISPPA, so ' ...
         'set babelbrain.isppa_w_per_cm2 explicitly for this medium.']);
end

rho_c = [];
if isfield(kwave_medium, 'sound_speed') && isfield(kwave_medium, 'density')
    rho_c = double(kwave_medium.sound_speed) .* double(kwave_medium.density);
    if ~isscalar(rho_c) && ~isequal(size(rho_c), dims)
        rho_c = [];
    end
end
if isempty(rho_c)
    mp = parameters.medium_properties.brain;
    rho_c = mp.sound_speed * mp.density;
end

I = p.^2 ./ (2 * rho_c) * 1e-4;   % [W/cm²]
isppa = max(I(brain));

end
