function results_heating = read_babelbrain_thermal(thermal_h5, parameters, planimg, focus_pos, focal_length_mm)
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
%
% Input:
%   thermal_h5   - path to BabelBrain *ThermalField*.h5
%   parameters   - PRESTUS config (used for grid.dims and grid.resolution_mm)
%   planimg      - planimg struct with transf (T1→grid affine) and t1_header
%
% Output:
%   results_heating - struct compatible with thermal_simulation output
%
% See also: BABELBRAIN_THERMAL_WRAPPER, THERMAL_ANALYSIS, NIFTI_THERMAL

arguments
    thermal_h5      (1,:) char
    parameters      (1,1) struct
    planimg         (1,1) struct
    focus_pos       (1,3) double  = [NaN NaN NaN]
    focal_length_mm (1,1) double  = NaN
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
% BabelBrain thermal also applies np.flip on load (axis=2); z_vec is NOT flipped
% → flip z_vec here.  After permute([3,2,1]):
%   dim 1 = BB x (lateral X), 0 at focal spot
%   dim 2 = BB y (lateral Y), 0 at focal spot
%   dim 3 = BB z (depth, flip-corrected), FocalLength at focal spot
x_vec_raw = read_field('x_vec');   % BB x [mm]
y_vec_raw = read_field('y_vec');   % BB y [mm]
z_vec_raw = read_field('z_vec');   % BB z [mm], un-flipped
x_mm = x_vec_raw;                 % PRESTUS dim 1: lateral X, 0 at focal spot [mm]
y_mm = y_vec_raw;                 % PRESTUS dim 2: lateral Y, 0 at focal spot [mm]
z_mm = z_vec_raw(end:-1:1);      % PRESTUS dim 3: depth (flip-corrected), FocalLength at focal [mm]

temp_profile      = read_field('TempProfileTarget');   % [1 x Nt] °C at focal target
time_profile      = read_field('TimeProfileTarget');   % [1 x Nt] s
% TemperaturePoints: [4 x Nt] — rows are skin(1), brain(2), skull(3), target(4).
% Row 4 may be absent when the target coincides with the brain maximum
% (in that case BabelBrain uses IndTarget=2 and only 3 monitor points exist).
temp_points_raw   = read_field('TemperaturePoints');   % [Npts x Nt]

MaxIsppa = double(h5read(thermal_h5, '/MaxIsppa'));
MaxIspta = double(h5read(thermal_h5, '/MaxIspta'));
Isppa    = double(h5read(thermal_h5, '/Isppa'));
Ispta    = double(h5read(thermal_h5, '/Ispta'));

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
% Resampling is anchored at focus_pos (PRESTUS grid voxels) using focal_length_mm.
if any(isnan(focus_pos)) || isnan(focal_length_mm)
    error(['read_babelbrain_thermal: focus_pos and focal_length_mm are required for ' ...
        'BabelBrain→PRESTUS resampling but were not provided.\n' ...
        'Pass them as the 4th and 5th arguments.']);
end

prestus_dims = parameters.grid.dims;
res_mm       = parameters.grid.resolution_mm;

resample_vol = @(vol) resample_bb_to_prestus( ...
    vol, x_mm(:)', y_mm(:)', z_mm(:)', planimg, prestus_dims, res_mm, focus_pos, focal_length_mm);

T_end_fus_g = resample_vol(T_end_fus);
T_final_g   = resample_vol(T_final);
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
% Required fields (from thermal_simulation.m line 554–563):
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
    results_heating.timeseries.Tdiff.(lbl)  = [0, diff(temp_points_raw(ti, :))];
end
results_heating.dx_step              = res_mm * 1e-3;   % [m]
results_heating.MaxIsppa             = MaxIsppa;
results_heating.MaxIspta             = MaxIspta;
results_heating.Isppa                = Isppa;
results_heating.Ispta                = Ispta;
results_heating.babelbrain_thermal_source = thermal_h5;

end

%% Local helpers

function vol_out = resample_bb_to_prestus(vol_bb, x_mm, y_mm, z_mm, planimg, prestus_dims, res_mm, focus_pos, focal_length_mm)
    vol_out = resample_babelbrain_vol(vol_bb, x_mm, y_mm, z_mm, planimg, prestus_dims, res_mm, focus_pos, focal_length_mm);
end
