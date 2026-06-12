function results_heating = read_babelbrain_thermal(thermal_h5, parameters, planimg)
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
%   DoseEndFUS   - [Nx x Ny x Nz] CEM43 at end of sonication [min]
%   FinalDose    - [Nx x Ny x Nz] final CEM43 after cooling [min]
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
    thermal_h5   (1,:) char
    parameters   (1,1) struct
    planimg      (1,1) struct
end

if ~isfile(thermal_h5)
    error('read_babelbrain_thermal: file not found: %s', thermal_h5);
end

%% Read from HDF5
% BabelBrain writes via SaveToH5py which produces one dataset per key.
read_field = @(key) double(h5read(thermal_h5, ['/' key]));

T_end_fus  = read_field('TempEndFUS');   % temperature at end of sonication
T_final    = read_field('FinalTemp');    % final temperature after cooling
cem_end    = read_field('DoseEndFUS');   % CEM43 at end of sonication
cem_final  = read_field('FinalDose');    % final CEM43

p_map      = read_field('p_map');        % peak pressure [Pa]

x_mm = read_field('x_vec');   % [1 x N1] mm
y_mm = read_field('y_vec');   % [1 x N2] mm
z_mm = read_field('z_vec');   % [1 x N3] mm

temp_profile  = read_field('TempProfileTarget');  % [1 x Nt] °C
time_profile  = read_field('TimeProfileTarget');  % [1 x Nt] s

MaxIsppa = double(h5read(thermal_h5, '/MaxIsppa'));
MaxIspta = double(h5read(thermal_h5, '/MaxIspta'));
Isppa    = double(h5read(thermal_h5, '/Isppa'));
Ispta    = double(h5read(thermal_h5, '/Ispta'));

fprintf('[read_babelbrain_thermal] Grid: [%d %d %d]  dx=%.3f mm\n', ...
    size(T_end_fus,1), size(T_end_fus,2), size(T_end_fus,3), ...
    mean(diff(x_mm(:)')));
fprintf('[read_babelbrain_thermal] Peak T rise: %.2f °C  MaxIsppa: %.2f W/cm²\n', ...
    max(T_end_fus(:)) - min(T_end_fus(:)), MaxIsppa);

%% Resample BabelBrain volumes → PRESTUS grid
% BabelBrain coordinates (x_vec/y_vec/z_vec) are in RAS mm.
% planimg.transf maps T1 voxels → PRESTUS acoustic grid voxels.
% Build affine from BabelBrain RAS mm → PRESTUS grid voxels and resample.

prestus_dims = parameters.grid.dims;
res_mm       = parameters.grid.resolution_mm;

resample_vol = @(vol) resample_bb_to_prestus( ...
    vol, x_mm(:)', y_mm(:)', z_mm(:)', planimg, prestus_dims, res_mm);

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

T_rise_g     = T_end_fus_g - parameters.thermal.temp_0.brain;
T_final_rise = T_final_g   - parameters.thermal.temp_0.brain;

results_heating.maxT                 = T_rise_g;
results_heating.heating_endT         = T_final_rise;
results_heating.focal_planeT         = [];    % BabelBrain has no focal-plane slice; filled as empty
results_heating.CEM43                = cem_end_g;
results_heating.CEM43_end            = cem_final_g;
results_heating.focal_planeCEM43     = [];
% BabelBrain does not compute iso-contour CEM43; use the same map (conservative)
results_heating.CEM43_iso            = cem_end_g;
results_heating.CEM43_iso_end        = cem_final_g;
results_heating.focal_planeCEM43_iso = [];
% Build a minimal timeseries struct compatible with thermal_plot_sim.
% BabelBrain provides only a single time-course at the focal target, not
% per-layer timeseries. Initialise an empty struct via thermal_update_timeseries
% so field names are correct, then attach the focal time course separately.
results_heating.timeseries           = struct('T', struct(), 'Tdiff', struct(), ...
    'CEM43', struct(), 'CEM43_iso', struct());
results_heating.timeseries.TempProfileTarget = temp_profile(:)';
results_heating.timeseries.TimeProfileTarget = time_profile(:)';
results_heating.dx_step              = res_mm * 1e-3;   % [m]
results_heating.MaxIsppa             = MaxIsppa;
results_heating.MaxIspta             = MaxIspta;
results_heating.Isppa                = Isppa;
results_heating.Ispta                = Ispta;
results_heating.babelbrain_thermal_source = thermal_h5;

end

%% Local helpers

function vol_out = resample_bb_to_prestus(vol_bb, x_mm, y_mm, z_mm, planimg, prestus_dims, res_mm)
% Resample a BabelBrain 3D volume (in RAS mm coordinates) to the PRESTUS
% acoustic grid via affine_resample_3d.
%
% BabelBrain grid origin in RAS mm: [x_mm(1), y_mm(1), z_mm(1)]
% BabelBrain voxel spacing: mean(diff(x_mm)) mm (assumed isotropic)
%
% planimg.inv_transf maps PRESTUS-grid voxels → T1 voxels (row-vector convention).
% We need BabelBrain voxels → PRESTUS grid voxels.

    dx_bb = mean(diff(x_mm));   % mm per voxel in BabelBrain grid

    % Build RAS-mm → BabelBrain-voxel affine (4x4, column-major for MATLAB)
    % BabelBrain voxel index (1-based) = (ras_mm - origin) / dx + 1
    % Written as affine: P_vox = M * P_ras
    origin = [x_mm(1), y_mm(1), z_mm(1)];
    M_ras_to_bb = diag([1/dx_bb, 1/dx_bb, 1/dx_bb, 1]);
    M_ras_to_bb(1:3, 4) = (-origin ./ dx_bb + 1)';  % 1-based voxel offset

    % planimg.transf: T1-voxels → acoustic-grid-voxels (row-vector convention)
    % planimg.t1_header.Transform: MATLAB affine3d (T1 RAS mm → T1 voxels)
    if ~isempty(planimg) && isfield(planimg, 't1_header') && ~isempty(planimg.t1_header)
        % T1 RAS mm → T1 voxels via inverse of t1_header.Transform
        T_ras_to_t1vox = inv(planimg.t1_header.Transform.T');  % 4x4
        % T1 voxels → acoustic grid voxels
        T_t1vox_to_acgrid = planimg.transf;  % 4x4
        % Combined: BabelBrain RAS mm → acoustic grid voxels
        M_bb_to_acgrid = T_t1vox_to_acgrid * T_ras_to_t1vox * inv(M_ras_to_bb);
    else
        % Non-layered / phantom: diagonal affine assuming same RAS origin
        % and converting mm to voxel indices via resolution_mm
        M_bb_to_acgrid = diag([dx_bb/res_mm, dx_bb/res_mm, dx_bb/res_mm, 1]);
    end

    vol_out = affine_resample_3d(single(vol_bb), M_bb_to_acgrid, prestus_dims, 'linear', 0);
    vol_out = double(vol_out);
end
