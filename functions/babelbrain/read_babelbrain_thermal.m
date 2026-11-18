function results_heating = read_babelbrain_thermal(thermal_h5, parameters, planimg, focus_pos, focal_length_mm, trans_pos)
% READ_BABELBRAIN_THERMAL  Read BabelBrain thermal HDF5 into PRESTUS results_heating
%
% Reads the ThermalField HDF5 written by BabelBrain's CalculateTemperatureEffects,
% resamples the 3D temperature/CEM43/pressure volumes from BabelBrain's simulation
% grid to the PRESTUS acoustic grid space, and returns a results_heating struct
% compatible with PRESTUS thermal_analysis and nifti_thermal.
%
% BabelBrain thermal output fields used (from CalculateTemperatureEffects.py):
%   TempEndFUS   - [Nx x Ny x Nz] temperature at end of sonication [°C]
%   FinalTemp    - [Nx x Ny x Nz] final temperature after cooling [°C]
%   DoseEndFUS   - [Nx x Ny x Nz] raw CEM43 accumulation at end of sonication [CEM43·s]
%   FinalDose    - [Nx x Ny x Nz] raw CEM43 accumulation after cooling [CEM43·s]
%                  (divided by 60 on read to convert to [min] — BabelBrain's BHTE
%                   accumulates with dt in seconds; scalar summaries CEMBrain/CEMSkin
%                   are already in minutes but 3D arrays are not)
%   p_map        - [Nx x Ny x Nz] peak pressure map [Pa]
%   x_vec / y_vec / z_vec - spatial coordinates [mm]
%   TempProfileTarget / TimeProfileTarget - temperature time course at focus
%   MaxIsppa / MaxIspta / Isppa / Ispta - intensity metrics [W/cm²]
%   affine       - 4x4 affine from BabelBrain grid to RAS space
%
% The coordinate vectors (x_vec, y_vec, z_vec) are used to construct the
% BabelBrain→PRESTUS-grid resampling transform. If the BabelBrain grid matches
% the PRESTUS grid resolution, volumes are placed directly; otherwise linear
% resampling is applied via affine_resample_3d.
%
% Use as:
%   results_heating = read_babelbrain_thermal(thermal_h5, parameters, planimg)
%   results_heating = read_babelbrain_thermal(..., focus_pos, focal_length_mm, trans_pos)
%
% Input:
%   thermal_h5   - path to BabelBrain *ThermalField*.h5
%   parameters   - PRESTUS config (used for grid.dims and grid.resolution_mm)
%   planimg      - planimg struct with transf (T1→grid affine) and t1_header
%   focus_pos    - (optional) [1x3] focus voxel; anchors the depth resample on
%                  /TargetLocation (as READ_BABELBRAIN_OUTPUT)
%   focal_length_mm - (optional) preset FocalLength [mm]; anchor fallback
%   trans_pos    - (optional) [1x3] transducer voxel; sets the beam direction
%
% Output:
%   results_heating - struct compatible with thermal_simulation output: p_amp,
%                     maxT, heating_endT, CEM43, CEM43_end (ISO CEM43 fields
%                     NaN, focal-plane fields empty), timeseries, MaxIsppa,
%                     MaxIspta, Isppa, Ispta, dx_step, babelbrain_thermal_source
%
% See also: BABELBRAIN_THERMAL_WRAPPER, THERMAL_ANALYSIS, NIFTI_THERMAL

arguments
    thermal_h5      (1,:) char
    parameters      (1,1) struct
    planimg         (1,1) struct
    focus_pos       (1,3) double  = [NaN NaN NaN]
    focal_length_mm (1,1) double  = NaN
    trans_pos             double  = []   % [1x3] transducer voxel (same grid as focus_pos); sets beam-axis direction
end

if ~isfile(thermal_h5)
    error('read_babelbrain_thermal: file not found: %s', thermal_h5);
end

%% Read from HDF5
% BabelBrain writes via SaveToH5py which produces one dataset per key.
read_field = @(key) double(h5read(thermal_h5, ['/' key]));

% h5read reverses Python axis order: MATLAB dim 1 = BB z (depth), dim 2 = BB y, dim 3 = BB x.
% PRESTUS grid: dim 1 = lateral X (BB x), dim 2 = lateral Y (BB y), dim 3 = axial (BB z).
% Permute all 3D volumes immediately to PRESTUS axis order [N_x, N_y, N_z].
perm321 = @(v) permute(v, [3, 2, 1]);
T_end_fus  = perm321(read_field('TempEndFUS'));   % temperature at end of sonication
T_final    = perm321(read_field('FinalTemp'));    % final temperature after cooling
% BabelBrain accumulates CEM43 with the BHTE time step in SECONDS; the 3D
% DoseEndFUS/FinalDose arrays are therefore in CEM43·s, not CEM43·min.
cem_end    = perm321(read_field('DoseEndFUS')) / 60;   % CEM43 at end of sonication [min]
cem_final  = perm321(read_field('FinalDose'))  / 60;   % final CEM43 [min]

p_map      = perm321(read_field('p_map'));        % peak pressure [Pa]

% Coordinate vectors [mm] for the thermal output domain.
% Unlike the acoustic DataForSim.h5 (whose arrays are stored depth-flipped
% relative to z_vec, see read_babelbrain_output), the thermal arrays are in
% z_vec order: CalculateTemperatureEffects np.flip()s p_amp / MaterialMap on
% load (axis 2) and saves its results without flipping back (thermal
% MaterialMap == flip(acoustic MaterialMap, 3)). So z_vec is used as is;
% reversing it would mirror the heating about the focal marker. After permute([3,2,1]):
%   dim 1 = BB x (lateral X), 0 at focal spot
%   dim 2 = BB y (lateral Y), 0 at focal spot
%   dim 3 = BB z (depth), FocalLength at focal spot
x_vec_raw = read_field('x_vec');   % BB x [mm]
y_vec_raw = read_field('y_vec');   % BB y [mm]
z_vec_raw = read_field('z_vec');   % BB z [mm]
x_mm = x_vec_raw;                 % PRESTUS dim 1: lateral X, 0 at focal spot [mm]
y_mm = y_vec_raw;                 % PRESTUS dim 2: lateral Y, 0 at focal spot [mm]
z_mm = z_vec_raw;                 % PRESTUS dim 3: depth, FocalLength at focal [mm]

% Inventory root-level datasets once; some fields may only exist nested under
% AllData/item_0 depending on the BabelBrain version.
info_        = h5info(thermal_h5);
root_names_  = {};
if ~isempty(info_.Datasets)
    root_names_ = {info_.Datasets.Name};
end
nested_prefix_ = '/AllData/item_0/';
nested_names_  = {};
try
    nested_info_ = h5info(thermal_h5, '/AllData/item_0');
    if ~isempty(nested_info_.Datasets)
        nested_names_ = {nested_info_.Datasets.Name};
    end
catch
end

read_any_ = @(key) read_field_any(thermal_h5, key, root_names_, nested_prefix_, nested_names_);

% TempProfileTarget/TimeProfileTarget: root or AllData/item_0
temp_profile = read_any_('TempProfileTarget');   % [1 x Nt] °C at focal target
time_profile = read_any_('TimeProfileTarget');   % [1 x Nt] s

% TemperaturePoints: [Npts x Nt] — rows are skin(1), brain(2), skull(3), target(4).
% Row 4 may be absent when the target coincides with the brain maximum.
% May be absent in some BabelBrain versions; fall back to empty gracefully.
if any(strcmp(root_names_, 'TemperaturePoints'))
    temp_points_raw = read_field('TemperaturePoints');
elseif any(strcmp(nested_names_, 'TemperaturePoints'))
    temp_points_raw = double(h5read(thermal_h5, [nested_prefix_ 'TemperaturePoints']));
else
    warning('read_babelbrain_thermal:missingTemperaturePoints', ...
        'TemperaturePoints not found in %s; per-tissue time courses will be empty.', thermal_h5);
    temp_points_raw = zeros(0, 0);
end

% Intensity scalars: root or AllData/item_0; NaN when absent.
MaxIsppa = read_scalar_any_(thermal_h5, 'MaxIsppa', root_names_, nested_prefix_, nested_names_);
MaxIspta = read_scalar_any_(thermal_h5, 'MaxIspta', root_names_, nested_prefix_, nested_names_);
Isppa    = read_scalar_any_(thermal_h5, 'Isppa',    root_names_, nested_prefix_, nested_names_);
Ispta    = read_scalar_any_(thermal_h5, 'Ispta',    root_names_, nested_prefix_, nested_names_);

fprintf('[read_babelbrain_thermal] Grid (PRESTUS order): [%d %d %d]  dx_lat=%.3f dz_depth=%.3f mm\n', ...
    size(T_end_fus,1), size(T_end_fus,2), size(T_end_fus,3), ...
    mean(diff(x_mm(:)')), abs(mean(diff(z_mm(:)'))));
fprintf('[read_babelbrain_thermal] Peak T rise: %.2f °C  MaxIsppa: %.2f W/cm²\n', ...
    max(T_end_fus(:)) - min(T_end_fus(:)), MaxIsppa);

%% Resample BabelBrain volumes → PRESTUS grid
% After permute([3,2,1]) and coord reassignment (see above):
%   x_mm (PRESTUS dim 1): lateral X, 0 at focal spot [mm]
%   y_mm (PRESTUS dim 2): lateral Y, 0 at focal spot [mm]
%   z_mm (PRESTUS dim 3): depth (flip-corrected), FocalLength at focal spot [mm]
%
% Same tissue-marker anchor as read_babelbrain_output: /TargetLocation (the
% label-5 voxel written at focus_pos) is mapped back onto focus_pos, so the
% thermal volumes align with medium_masks exactly as the acoustic field does.
% focal_length_mm is only a FALLBACK for files without /TargetLocation.
% /TargetLocation indexes z_vec_raw (= z_mm here).
focal_length_mm_effective = focal_length_mm;
if any(strcmp(root_names_, 'TargetLocation')) || any(strcmp(nested_names_, 'TargetLocation'))
    try
        target_loc_bb = double(read_any_('TargetLocation')) + 1;   % 0-based (numpy) -> 1-based (MATLAB)
        z_vec_raw_col = z_vec_raw(:)';
        if numel(z_vec_raw_col) >= target_loc_bb(3)
            z_mm_at_target = z_vec_raw_col(target_loc_bb(3));   % marker z
            focal_length_mm_effective = z_mm_at_target;
            if isnan(focal_length_mm) || abs(z_mm_at_target - focal_length_mm) > 1e-6
                fprintf(['[read_babelbrain_thermal] Anchor: /TargetLocation tissue marker at %.2f mm ' ...
                    '-> focus_pos (passed fallback anchor was %.2f mm).\n'], ...
                    z_mm_at_target, focal_length_mm);
            else
                fprintf('[read_babelbrain_thermal] Anchor: /TargetLocation tissue marker at %.2f mm -> focus_pos.\n', ...
                    z_mm_at_target);
            end
        end
    catch ME
        warning('read_babelbrain_thermal:targetLocationReadFailed', ...
            'Could not read /TargetLocation for anchor: %s', ME.message);
    end
end
if focal_length_mm_effective == focal_length_mm && ~isnan(focal_length_mm) && ...
        ~(any(strcmp(root_names_, 'TargetLocation')) || any(strcmp(nested_names_, 'TargetLocation')))
    fprintf(['[read_babelbrain_thermal] /TargetLocation unavailable — falling back to the passed ' ...
        'anchor %.2f mm. Tissue alignment is only guaranteed when this equals the marker z ' ...
        '(preset FocalLength).\n'], focal_length_mm);
end

if any(isnan(focus_pos)) || isnan(focal_length_mm_effective)
    error(['read_babelbrain_thermal: focus_pos and focal_length_mm are required for ' ...
        'BabelBrain→PRESTUS resampling but were not provided (and no /TargetLocation ' ...
        'was available in the file to derive an anchor automatically).\n' ...
        'Pass them as the 4th and 5th arguments.']);
end

prestus_dims = parameters.grid.dims;
res_mm       = parameters.grid.resolution_mm;

% Beam-axis direction (see read_babelbrain_output): −z beams need the depth
% resample mirrored. Default +1 when trans_pos is unavailable.
if nargin >= 6 && ~isempty(trans_pos) && numel(trans_pos) >= 3 && ~any(isnan(focus_pos))
    beam_sign = sign(focus_pos(3) - trans_pos(3));
    if beam_sign == 0, beam_sign = 1; end
else
    beam_sign = 1;
end

resample_vol = @(vol) resample_bb_to_prestus( ...
    vol, x_mm(:)', y_mm(:)', z_mm(:)', planimg, prestus_dims, res_mm, focus_pos, focal_length_mm_effective, beam_sign);

% Temperatures are resampled as rises over the baseline, so voxels outside
% BabelBrain's thermal domain (zero-filled by the resampler) read the baseline
% temperature, not 0 degC (which appeared as a -37 degC "heatrise").
T_base = 37;
if any(strcmp(root_names_, 'BaselineTemperature')) || any(strcmp(nested_names_, 'BaselineTemperature'))
    T_base = read_any_('BaselineTemperature');
    T_base = T_base(1);
end
T_end_fus_g = resample_vol(T_end_fus - T_base) + T_base;
T_final_g   = resample_vol(T_final   - T_base) + T_base;
cem_end_g   = resample_vol(cem_end);
cem_final_g = resample_vol(cem_final);
p_map_g     = resample_vol(p_map);

% Hard check: every resampled volume must match the declared PRESTUS grid dims
exp = prestus_dims(:)';
for chk_ = {{'T_end_fus', T_end_fus_g}, {'T_final', T_final_g}, ...
             {'cem_end', cem_end_g}, {'cem_final', cem_final_g}, {'p_map', p_map_g}}
    act = size(chk_{1}{2});
    if ~isequal(act, exp)
        error(['read_babelbrain_thermal: resampled %s has size [%d %d %d] but expected ' ...
            'PRESTUS grid [%d %d %d].\nCheck planimg affine transforms and BabelBrain coordinate vectors.'], ...
            chk_{1}{1}, act(1), act(2), act(3), exp(1), exp(2), exp(3));
    end
end
clear chk_ act exp;

%% Assemble results struct
% Field names match thermal_simulation output exactly so that
% thermal_analysis and nifti_thermal work without modification.
% Required fields (as thermal_simulation returns them):
%   maxT, focal_planeT, heating_endT, CEM43, focal_planeCEM43,
%   CEM43_end, timeseries, CEM43_iso, focal_planeCEM43_iso, CEM43_iso_end

% Store absolute temperature (matching thermal_simulation convention).
% thermal_analysis subtracts temp_0 itself when computing riseT_* fields.
results_heating.p_amp                = p_map_g;
results_heating.maxT                 = T_end_fus_g;
results_heating.heating_endT         = T_final_g;
results_heating.focal_planeT         = [];    % BabelBrain has no focal-plane slice; filled as empty
results_heating.CEM43                = cem_end_g;
results_heating.CEM43_end            = cem_final_g;
results_heating.focal_planeCEM43     = [];
% BabelBrain does not compute ISO-formula CEM43 (R-value method); set to NaN
% so downstream code does not silently treat it as ISO-compliant.
results_heating.CEM43_iso            = NaN(size(cem_end_g));
results_heating.CEM43_iso_end        = NaN(size(cem_final_g));
results_heating.focal_planeCEM43_iso = [];
% Build timeseries struct.  BabelBrain stores per-tissue time courses in
% TemperaturePoints [Npts x Nt]: row 1=skin, row 2=brain, row 3=skull,
% row 4=target (only present when target ≠ brain maximum).
t_ax = time_profile(:)';
results_heating.timeseries           = struct('T', struct(), 'Tdiff', struct(), ...
    'CEM43', struct(), 'CEM43_iso', struct());
results_heating.timeseries.TempProfileTarget = temp_profile(:)';
results_heating.timeseries.TimeProfileTarget = t_ax;
tissue_labels = {'skin', 'brain', 'skull', 'target'};
for ti = 1:min(size(temp_points_raw, 1), numel(tissue_labels))
    lbl = tissue_labels{ti};
    results_heating.timeseries.T.(lbl)      = temp_points_raw(ti, :);
    results_heating.timeseries.Tdiff.(lbl)  = temp_points_raw(ti, :) - temp_points_raw(ti, 1);
end
results_heating.dx_step              = res_mm * 1e-3;   % [m]
results_heating.MaxIsppa             = MaxIsppa;
results_heating.MaxIspta             = MaxIspta;
results_heating.Isppa                = Isppa;
results_heating.Ispta                = Ispta;
results_heating.babelbrain_thermal_source = thermal_h5;

end

%% Local helpers

function vol_out = resample_bb_to_prestus(vol_bb, x_mm, y_mm, z_mm, planimg, prestus_dims, res_mm, focus_pos, focal_length_mm, beam_sign)
    if nargin < 10 || isempty(beam_sign)
        beam_sign = 1;
    end
    vol_out = resample_babelbrain_vol(vol_bb, x_mm, y_mm, z_mm, planimg, prestus_dims, res_mm, focus_pos, focal_length_mm, 'linear', beam_sign);
end

function val = read_field_any(h5, key, root_names, nested_prefix, nested_names)
% Read a scalar or array field from root or nested AllData/item_0 group.
    if any(strcmp(root_names, key))
        val = double(h5read(h5, ['/' key]));
    elseif any(strcmp(nested_names, key))
        val = double(h5read(h5, [nested_prefix key]));
    else
        error('read_babelbrain_thermal:missingField', ...
            'Required field ''%s'' not found at root or under %s in:\n  %s', ...
            key, nested_prefix, h5);
    end
end

function val = read_scalar_any_(h5, key, root_names, nested_prefix, nested_names)
% Read an optional scalar from root or nested group; return NaN if absent.
    if any(strcmp(root_names, key))
        val = double(h5read(h5, ['/' key]));
    elseif any(strcmp(nested_names, key))
        val = double(h5read(h5, [nested_prefix key]));
    else
        warning('read_babelbrain_thermal:missingScalar', ...
            'Optional field ''%s'' not found in thermal HDF5; returning NaN.', key);
        val = NaN;
    end
end
