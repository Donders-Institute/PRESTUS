function sensor_data = read_babelbrain_output(h5_path, use_refocus, parameters, planimg, focus_pos, focal_length_mm, trans_pos)
% READ_BABELBRAIN_OUTPUT  Read BabelBrain DataForSim.h5 into PRESTUS sensor_data
%
% Reads the peak pressure field from a BabelBrain *DataForSim.h5 output file,
% resamples it from BabelBrain's internal simulation grid to the PRESTUS acoustic
% grid when the two resolutions differ, and packages the result into the
% sensor_data struct expected by acoustic_analysis and thermal_simulation.
%
% BabelBrain runs at a grid step derived from its PPW setting:
%   dx_bb = c / (freq x ppw)
% This may differ from parameters.grid.resolution_mm (the PRESTUS acoustic grid
% step). When they differ, p_amp is resampled onto the PRESTUS grid using the
% same affine_resample_3d approach used by read_babelbrain_thermal. When
% parameters / planimg are not supplied, no resampling is performed.
%
% BabelBrain HDF5 output keys (babel_integration_base.py Step 9/10):
%   /p_amp          - [Nx x Ny x Nz] float32, peak positive pressure [Pa]
%   /p_amp_refocus  - same, after refocusing back-propagation (optional)
%   /p_complex      - [Nx x Ny x Nz] complex64, steady-state Fourier pressure [Pa]
%                     stored by h5py as a compound type with fields 'r' and 'i'
%   /p_complex_refocus - same, after refocusing (optional)
%   /x_vec, /y_vec, /z_vec - voxel centre coordinates [m] (cropped simulation grid)
%
% Use as:
%   sensor_data = read_babelbrain_output(h5_path)
%   sensor_data = read_babelbrain_output(h5_path, true)
%   sensor_data = read_babelbrain_output(h5_path, use_refocus, parameters, planimg)
%   sensor_data = read_babelbrain_output(..., focus_pos, focal_length_mm, trans_pos)
%
% Input:
%   h5_path      - path to *DataForSim.h5 produced by RUN_SIM().RunCases()
%   use_refocus  - (optional, default false) use /p_amp_refocus if available
%   parameters   - (optional) PRESTUS config struct; required for resampling
%                  (uses grid.resolution_mm and grid.dims)
%   planimg      - (optional) planimg struct from grid_tissue_setup; required
%                  for affine-correct resampling in layered (subject) mode
%   focus_pos    - (optional) [1x3] focus voxel; anchors the depth resample on
%                  BabelBrain's /TargetLocation marker
%   focal_length_mm - (optional) preset FocalLength [mm]; anchor fallback when
%                  /TargetLocation is absent
%   trans_pos    - (optional) [1x3] transducer voxel; sets the beam direction
%
% Output:
%   sensor_data  - struct with fields:
%                    .p_max_all   [Nx x Ny x Nz] double [Pa], peak pressure
%                    .p_complex   [Nx x Ny x Nz] complex double [Pa], Fourier
%                                 pressure (amplitude + phase); present when
%                                 /p_complex exists in the HDF5 file
%                    .babelbrain_material_map  BabelBrain labels on the PRESTUS
%                                 grid (when resampled and /MaterialMap exists)
%                    .babelbrain_x_vec/_y_vec/_z_vec  BabelBrain grid coordinates
%                    .babelbrain_source  h5_path
%                  Volumes are at parameters.grid.resolution_mm resolution when
%                  resampled.
%
% See also: BABELBRAIN_WRAPPER, ACOUSTIC_WRAPPER, READ_BABELBRAIN_THERMAL

arguments
    h5_path          (1,:) char
    use_refocus      (1,1) logical = false
    parameters       (1,1) struct  = struct()
    planimg                        = []
    focus_pos        (1,3) double  = [NaN NaN NaN]
    focal_length_mm  (1,1) double  = NaN
    trans_pos              double  = []   % [1x3] transducer voxel (same grid as focus_pos); sets beam-axis direction for the depth resample
end

% Beam-axis direction along dim 3: +1 when physical depth increases with dim 3
% (focus deeper in index than the transducer), −1 when it decreases. Must match
% the depth orientation write_babelbrain_input used (same focus_pos vs trans_pos
% test) so the resampled field is not mirrored. Defaults to +1 when trans_pos is
% unavailable (legacy beam-+z assumption).
if ~isempty(trans_pos) && numel(trans_pos) >= 3 && ~any(isnan(focus_pos))
    beam_sign = sign(focus_pos(3) - trans_pos(3));
    if beam_sign == 0, beam_sign = 1; end
else
    beam_sign = 1;
end

ensure_hdf5_plugin_path();

if ~isfile(h5_path)
    error('read_babelbrain_output: file not found: %s', h5_path);
end

% Inventory datasets in the file (needed for refocus fallback and p_complex check)
info     = h5info(h5_path);
ds_names = {info.Datasets.Name};

% Select amplitude pressure key
if ~any(strcmp(ds_names, 'p_amp'))
    error('read_babelbrain_output:missingDataset', ...
        '/p_amp dataset not found in BabelBrain output:\n  %s', h5_path);
end
p_key = '/p_amp';
if use_refocus
    if any(strcmp(ds_names, 'p_amp_refocus'))
        p_key = '/p_amp_refocus';
        fprintf('[read_babelbrain_output] Using refocused pressure map.\n');
    else
        warning('read_babelbrain_output: /p_amp_refocus not found; falling back to /p_amp.');
    end
end

% h5read reverses Python dimension order: MATLAB dim 1 = BabelBrain z (depth, Python axis 2),
% MATLAB dim 2 = BabelBrain y (Python axis 1), MATLAB dim 3 = BabelBrain x (Python axis 0).
% PRESTUS grid: dim 1 = lateral X (BB x), dim 2 = lateral Y (BB y), dim 3 = axial (BB z).
% Permute immediately to PRESTUS axis order [N_x, N_y, N_z].
p_amp_raw = double(h5read(h5_path, p_key));   % [N_z x N_y x N_x] after h5read reversal
p_amp     = permute(p_amp_raw, [3, 2, 1]);    % [N_x x N_y x N_z], PRESTUS axis order

% Read spatial coordinate vectors
% If parameters are supplied (resampling may be needed), fail hard on missing coords
% rather than silently skipping the resample step with a corrupt output.
has_params = isfield(parameters, 'grid') ...
             && isfield(parameters.grid, 'resolution_mm') ...
             && isfield(parameters.grid, 'dims');

xdim = []; ydim = []; zdim = [];
coord_read_ok = false;
try
    % Coordinate vectors [m].  BabelBrain names: x=lateral-X, y=lateral-Y, z=depth.
    % After permute([3,2,1]) above: PRESTUS dim 1 = BB x, dim 2 = BB y, dim 3 = BB z.
    % BabelBrain applied np.flip(p_amp, axis=2) on the depth axis before saving;
    % z_vec was NOT correspondingly flipped — flip it here so coord ↔ data align.
    x_vec_bb = h5read(h5_path, '/x_vec');   % BB x: lateral X, 0 at focal spot [m]
    y_vec_bb = h5read(h5_path, '/y_vec');   % BB y: lateral Y, 0 at focal spot [m]
    z_vec_bb = h5read(h5_path, '/z_vec');   % BB z: depth, FocalLength at focal spot [m]
    xdim = x_vec_bb;               % PRESTUS dim 1 = BB x (lateral X) [m]
    ydim = y_vec_bb;               % PRESTUS dim 2 = BB y (lateral Y) [m]
    zdim = z_vec_bb(end:-1:1);    % PRESTUS dim 3 = BB z (depth, flip-corrected) [m]
    coord_read_ok = true;
    dx_bb_mm = mean(diff(xdim)) * 1e3;
    dy_bb_mm = mean(diff(ydim)) * 1e3;
    dz_bb_mm = mean(diff(zdim)) * 1e3;   % negative (depth decreases along dim 3 after flip)
    fprintf('[read_babelbrain_output] BabelBrain grid: [%d %d %d]  dx_lat=%.4f dy_lat=%.4f dz_depth=%.4f mm\n', ...
        size(p_amp,1), size(p_amp,2), size(p_amp,3), dx_bb_mm, dy_bb_mm, dz_bb_mm);
    sensor_data.babelbrain_x_vec = x_vec_bb;
    sensor_data.babelbrain_y_vec = y_vec_bb;
    sensor_data.babelbrain_z_vec = z_vec_bb;
catch ME
    dx_bb_mm = NaN;
    p_dims = size(p_amp);
    if has_params
        expected = parameters.grid.dims(:)';
        if isequal(p_dims, expected)
            warning(['read_babelbrain_output: coordinate vectors /x_vec /y_vec /z_vec not found in:\n  %s\n' ...
                'BabelBrain output dims [%d %d %d] match the PRESTUS grid — skipping resampling.\n' ...
                'Delete the cached .h5 file and re-run to get coordinate outputs.'], ...
                h5_path, p_dims(1), p_dims(2), p_dims(3));
        else
            error(['read_babelbrain_output: coordinate vectors /x_vec /y_vec /z_vec not found in:\n  %s\n' ...
                'Error: %s\n' ...
                'BabelBrain output dims [%d %d %d] do not match PRESTUS grid [%d %d %d].\n' ...
                'Delete the cached .h5 and re-run to enable resampling.'], ...
                h5_path, ME.message, p_dims(1), p_dims(2), p_dims(3), expected(1), expected(2), expected(3));
        end
    end
    fprintf('[read_babelbrain_output] BabelBrain grid: [%d %d %d]  (no coordinate vectors; no resampling)\n', ...
        p_dims(1), p_dims(2), p_dims(3));
end

% Resample z-axis anchor — TISSUE-MARKER anchored (PRESTUS-space alignment).
%
% The resample must invert the forward mask transform (write_babelbrain_input)
% so that every returned volume, field and MaterialMap alike, lands back on the
% PRESTUS grid. The one voxel known exactly on both sides is the label-5 focal
% marker: PRESTUS wrote it at focus_pos, and BabelBrain reports it as
% /TargetLocation. Mapping z_vec_bb(TargetLocation) onto focus_pos therefore
% aligns tissue by construction. The field is aligned too, because
% run_babelbrain_case.py neutralizes ZDomainStart so the source sits at the
% requested TxMechanicalAdjustmentZ; the remaining source snap (<= 1 native
% voxel) is reported in the sidecar diagnostic below.
%
% focal_length_mm (preset FocalLength) is only a FALLBACK when /TargetLocation
% is absent. /TargetLocation indexes the RAW, pre-flip z_vec_bb, not the
% flip-corrected z axis used for the data arrays.
% See doc_babelbrain.md (Output alignment).
target_loc_bb = [];
focal_length_mm_effective = focal_length_mm;
if coord_read_ok && any(strcmp(ds_names, 'TargetLocation'))
    try
        target_loc_bb = double(h5read(h5_path, '/TargetLocation')) + 1;   % 0-based (numpy) -> 1-based (MATLAB)
        if numel(z_vec_bb) >= target_loc_bb(3)
            z_mm_at_target = z_vec_bb(target_loc_bb(3)) * 1e3;   % raw z_vec_bb — marker z
            focal_length_mm_effective = z_mm_at_target;
            if isnan(focal_length_mm) || abs(z_mm_at_target - focal_length_mm) > 1e-6
                fprintf(['[read_babelbrain_output] Anchor: /TargetLocation tissue marker at ' ...
                    '%.2f mm -> focus_pos (passed fallback anchor was %.2f mm).\n'], ...
                    z_mm_at_target, focal_length_mm);
            else
                fprintf('[read_babelbrain_output] Anchor: /TargetLocation tissue marker at %.2f mm -> focus_pos.\n', ...
                    z_mm_at_target);
            end
        end
    catch ME
        warning('read_babelbrain_output:targetLocationReadFailed', ...
            'Could not read /TargetLocation for anchor: %s', ME.message);
    end
end
if isempty(target_loc_bb) && ~isnan(focal_length_mm)
    fprintf(['[read_babelbrain_output] /TargetLocation unavailable — falling back to the passed ' ...
        'anchor %.2f mm. Tissue alignment is only guaranteed when this equals the marker z ' ...
        '(preset FocalLength).\n'], focal_length_mm);
end

% Source-placement DIAGNOSTIC (no anchor correction applied). The sidecar
% written by run_babelbrain_case.py records the realized-minus-requested source
% position (residual after the ZDomainStart neutralization; ideally <= 1 native
% voxel of whole-voxel snap-back). A large residual means the field was solved
% with the source off the PRESTUS-intended position — the tissue resample is
% still exact, but source-to-skull distances in the field differ from the
% intent, so surface it prominently.
[h5_dir, h5_name] = fileparts(h5_path);
tgt = regexprep(h5_name, '_\d+kHz.*$', '');   % '<target>_<freq>kHz_...' -> '<target>'
anchor_json = fullfile(h5_dir, [tgt '_prestus_anchor.json']);
if isfile(anchor_json)
    try
        aj = jsondecode(fileread(anchor_json));
        if isfield(aj, 'source_residual_mm') && isfinite(aj.source_residual_mm)
            res_mm_native = NaN;
            if coord_read_ok
                res_mm_native = abs(mean(diff(z_vec_bb))) * 1e3;   % native voxel [mm]
            end
            fprintf(['[read_babelbrain_output] Source placement residual (realized - requested): ' ...
                '%+.3f mm (native voxel %.3f mm).\n'], aj.source_residual_mm, res_mm_native);
            if isfinite(res_mm_native) && abs(aj.source_residual_mm) > 1.5 * res_mm_native
                warning('read_babelbrain_output:sourcePlacementResidual', ...
                    ['BabelBrain placed the source %+.2f mm off the requested position ' ...
                     '(> 1.5 native voxels). The tissue-anchored resample is unaffected, but ' ...
                     'the FIELD was solved with the transducer at this offset (focus shifted ' ...
                     'accordingly) — check the run log for "Need to apply correction to ' ...
                     'reposition Tx" (bowl inside the domain: write_babelbrain_input rim trim) ' ...
                     'and ZDomainStart lines (is the neutralization patch active?).'], ...
                    aj.source_residual_mm);
            end
        elseif isfield(aj, 'anchor_offset_mm') && isfinite(aj.anchor_offset_mm)
            % Legacy sidecar (pre tissue-anchor rework): offset = CoC - marker.
            % No longer applied to the anchor; report only.
            fprintf(['[read_babelbrain_output] Legacy anchor sidecar found (CoC - marker = %+.3f mm); ' ...
                'not applied — the resample is tissue-marker anchored. Re-run the simulation to ' ...
                'refresh the sidecar diagnostics.\n'], aj.anchor_offset_mm);
        end
    catch ME
        warning('read_babelbrain_output:anchorSidecarFailed', ...
            'Could not read anchor sidecar %s: %s', anchor_json, ME.message);
    end
end

% Resample onto PRESTUS acoustic grid when resolutions differ.
% Use abs() because dz_bb_mm is negative (depth flip makes it decreasing).
needs_resample = has_params && coord_read_ok && ~isnan(dx_bb_mm) && ~isempty(xdim) ...
                 && (abs(abs(dx_bb_mm) - parameters.grid.resolution_mm) > 1e-3 ...
                  || abs(abs(dy_bb_mm) - parameters.grid.resolution_mm) > 1e-3 ...
                  || abs(abs(dz_bb_mm) - parameters.grid.resolution_mm) > 1e-3);

% Validate required resampling parameters (focus_pos and a usable focal-length
% anchor). These anchor BabelBrain's internal coordinate system to the PRESTUS
% grid. focal_length_mm_effective prefers BabelBrain's own /TargetLocation
% (ground truth) over the passed-in focal_length_mm estimate; see above.
if needs_resample
    if any(isnan(focus_pos)) || isnan(focal_length_mm_effective)
        error(['read_babelbrain_output: focus_pos and focal_length_mm are required for ' ...
            'BabelBrain→PRESTUS resampling but were not provided (and no /TargetLocation ' ...
            'was available in the file to derive an anchor automatically).\n' ...
            'Pass them as the 5th and 6th arguments, e.g.:\n' ...
            '  read_babelbrain_output(h5, refocus, params, planimg, focus_pos, focal_length_mm)']);
    end
end

bb_debug = debug_on(parameters, 'babelbrain');

if needs_resample
    res_mm = parameters.grid.resolution_mm;
    fprintf('[read_babelbrain_output] Resampling %.4f mm -> %.4f mm  (target dims: [%d %d %d])\n', ...
        abs(dx_bb_mm), res_mm, parameters.grid.dims(1), parameters.grid.dims(2), parameters.grid.dims(3));
    % After permute([3,2,1]) and coord reassignment:
    %   xdim (= x_vec): lateral X, 0 at focal spot [m]  → PRESTUS dim 1
    %   ydim (= y_vec): lateral Y, 0 at focal spot [m]  → PRESTUS dim 2
    %   zdim (= flipud(z_vec)): depth, focal_length_mm_effective at focal spot [m] → PRESTUS dim 3
    x_mm = xdim(:)' * 1e3;   % lateral X [mm], 0 at focal spot
    y_mm = ydim(:)' * 1e3;   % lateral Y [mm], 0 at focal spot
    z_mm = zdim(:)' * 1e3;   % depth     [mm], focal_length_mm_effective at focal spot
    p_amp_bb_native = p_amp;   % pre-resample, BabelBrain-native grid (kept for debug FOV plot)
    p_amp = resample_bb_to_prestus(p_amp, x_mm, y_mm, z_mm, planimg, parameters.grid.dims, res_mm, focus_pos, focal_length_mm_effective, 'linear', beam_sign);
    fprintf('[read_babelbrain_output] Resampled grid: [%d %d %d]\n', ...
        size(p_amp,1), size(p_amp,2), size(p_amp,3));

    % Always resample BabelBrain's MaterialMap onto the PRESTUS grid (not only in
    % debug): acoustic_analysis uses it to mask the field with BabelBrain's OWN
    % segmentation, the medium the field was solved on, so voxels at the
    % skull-brain interface are attributed as BabelBrain labelled them (a
    % reflection hotspot in BabelBrain's cortical bone does not count as
    % brain). Mode
    % (majority) resampling — integer label data on a finer BB grid. BB labels:
    % 0=water 1=skin 2=cortical 3=trabecular 4=brain.
    if any(strcmp(ds_names, 'MaterialMap'))
        try
            mat_raw_am = double(h5read(h5_path, '/MaterialMap'));
            mat_bb_am  = permute(mat_raw_am, [3, 2, 1]);
            sensor_data.babelbrain_material_map = resample_bb_to_prestus(mat_bb_am, ...
                x_mm, y_mm, z_mm, planimg, parameters.grid.dims, res_mm, ...
                focus_pos, focal_length_mm_effective, 'mode', beam_sign);
        catch ME
            warning('read_babelbrain_output:materialMapResampleFailed', ...
                'Could not resample /MaterialMap for tissue masking: %s', ME.message);
        end
    end
    % Hard check: resampled dims must match the declared PRESTUS grid
    act = size(p_amp); exp = parameters.grid.dims(:)';
    if ~isequal(act, exp)
        error(['read_babelbrain_output: resampled p_amp size [%d %d %d] does not match ' ...
            'expected PRESTUS grid dims [%d %d %d].\n' ...
            'Check focus_pos, focal_length_mm, and BabelBrain coordinate vectors.'], ...
            act(1), act(2), act(3), exp(1), exp(2), exp(3));
    end
    if bb_debug
        % MaterialMap ground truth (when present): BabelBrain's own tissue-index
        % volume, on the SAME grid/indexing as p_amp (both flipped identically
        % before saving; see babel_integration_base.py Step10_GetResults). Resample
        % with MODE (majority vote) — this is integer label data, and the
        % BabelBrain-native grid is typically FINER than the PRESTUS grid, so
        % point-sampled 'nearest' would drop thin structures (e.g. sub-voxel
        % cortical bone) that miss the single sample point but occupy real
        % volume in the coarse cell. Mode preserves them; linear would smear
        % tissue indices into meaningless fractional values.
        material_bb_native = [];
        material_resampled = [];
        if any(strcmp(ds_names, 'MaterialMap'))
            try
                mat_raw = double(h5read(h5_path, '/MaterialMap'));
                material_bb_native = permute(mat_raw, [3, 2, 1]);
                material_resampled = resample_bb_to_prestus(material_bb_native, x_mm, y_mm, z_mm, ...
                    planimg, parameters.grid.dims, res_mm, focus_pos, focal_length_mm_effective, 'mode', beam_sign);
            catch ME
                warning('read_babelbrain_output:materialMapReadFailed', ...
                    'Could not read/resample /MaterialMap: %s', ME.message);
            end
        end

        % anchor_mismatch_mm: how far the passed fallback anchor differed from
        % the /TargetLocation marker actually used (0 in the normal case —
        % see the tissue-marker anchor block above). Uses the same raw
        % z_vec_bb lookup as that block, NOT z_mm/zdim.
        if ~isnan(focal_length_mm) && ~isempty(target_loc_bb) && numel(z_vec_bb) >= target_loc_bb(3)
            anchor_mismatch_mm = z_vec_bb(target_loc_bb(3)) * 1e3 - focal_length_mm;
        else
            anchor_mismatch_mm = NaN;
        end

        % FOV-crop offset for back-transforming the (possibly cropped) acoustic
        % grid to T1 voxel space (see the MaterialMap NIfTI block in
        % plot_babelbrain_fov_debug).
        if isfield(planimg, 'fov_offset_ac') && ~isempty(planimg.fov_offset_ac)
            fov_offset = planimg.fov_offset_ac;
        else
            fov_offset = [1 1 1];
        end

        try
            plot_babelbrain_fov_debug(p_amp_bb_native, p_amp, x_mm, y_mm, z_mm, ...
                res_mm, focus_pos, focal_length_mm_effective, parameters, ...
                material_bb_native, material_resampled, target_loc_bb, anchor_mismatch_mm, planimg, fov_offset, beam_sign);
        catch ME
            warning('read_babelbrain_output:fovDebugFailed', ...
                'plot_babelbrain_fov_debug failed: %s', ME.message);
        end

        % T1w-space debug copy of the BabelBrain OUTPUT pressure field, so it can
        % be overlaid directly on the T1 / medium_masks (the INPUT material map is
        % written the same way inside plot_babelbrain_fov_debug). Both are
        % back-transformed from the acoustic grid into T1 voxel space via the same
        % path nifti_to_t1w uses for medium_masks, guaranteeing overlay.
        try
            bb_dbg_dir = parameters.io.dir_img;
            if isfield(parameters.io, 'dir_debug_babelbrain') && isfolder(parameters.io.dir_debug_babelbrain)
                bb_dbg_dir = parameters.io.dir_debug_babelbrain;
            end
            p_out = fullfile(bb_dbg_dir, sprintf('sub-%03d_babelbrain_pressure%s', ...
                parameters.subject_id, parameters.io.output_affix));
            if debug_on(parameters, 'babelbrain', 'volumes')
                write_bb_t1w_debug(p_amp, p_out, planimg, fov_offset, 'single', 32, 'linear');
                debug_register(parameters, 'babelbrain', [p_out '.nii.gz'], ...
                    'BabelBrain output pressure amplitude in T1w space', 'volume');
            end
        catch ME
            warning('read_babelbrain_output:pressureT1wDebugFailed', ...
                'Could not write T1w pressure debug NIfTI: %s', ME.message);
        end
    end
elseif has_params && coord_read_ok && ~isnan(dx_bb_mm)
    fprintf('[read_babelbrain_output] Grid resolutions match (%.4f mm); no resampling needed.\n', dx_bb_mm);
    % Grids already aligned — carry the MaterialMap through unresampled for
    % BabelBrain-frame tissue masking in acoustic_analysis (see the resample
    % branch above for the rationale).
    if any(strcmp(ds_names, 'MaterialMap'))
        try
            mat_raw_nm = double(h5read(h5_path, '/MaterialMap'));
            sensor_data.babelbrain_material_map = permute(mat_raw_nm, [3, 2, 1]);
        catch ME
            warning('read_babelbrain_output:materialMapReadFailed', ...
                'Could not read /MaterialMap for tissue masking: %s', ME.message);
        end
    end
end

sensor_data.p_max_all         = p_amp;
sensor_data.babelbrain_source = h5_path;

p_max_pa = max(p_amp(:));
fprintf('[read_babelbrain_output] p_max = %.1f Pa (%.2f MPa)\n', p_max_pa, p_max_pa/1e6);

% Read complex Fourier pressure field when present.
% h5py stores numpy complex64 as a compound HDF5 type with fields 'r' (real)
% and 'i' (imaginary). MATLAB's h5read returns a struct with .r and .i.
cx_key = '/p_complex';
if use_refocus && any(strcmp(ds_names, 'p_complex_refocus'))
    cx_key = '/p_complex_refocus';
end
if any(strcmp(ds_names, cx_key(2:end)))
    raw = h5read(h5_path, cx_key);
    if isstruct(raw) && isfield(raw, 'r') && isfield(raw, 'i')
        p_cx = double(raw.r) + 1i * double(raw.i);
    else
        p_cx = double(raw);  % already numeric complex (uncommon)
    end
    % Apply the same permute as p_amp: h5read [N_z,N_y,N_x] → PRESTUS [N_x,N_y,N_z]
    p_cx = permute(p_cx, [3, 2, 1]);
    if needs_resample
        p_cx = resample_bb_to_prestus(real(p_cx), x_mm, y_mm, z_mm, planimg, parameters.grid.dims, res_mm, focus_pos, focal_length_mm_effective, 'linear', beam_sign) + ...
            1i * resample_bb_to_prestus(imag(p_cx), x_mm, y_mm, z_mm, planimg, parameters.grid.dims, res_mm, focus_pos, focal_length_mm_effective, 'linear', beam_sign);
    end
    sensor_data.p_complex = p_cx;
end

end

%% Local helpers

function vol_out = resample_bb_to_prestus(vol_bb, x_mm, y_mm, z_mm, planimg, prestus_dims, res_mm, focus_pos, focal_length_mm, method, beam_sign)
    if nargin < 10 || isempty(method)
        method = 'linear';
    end
    if nargin < 11 || isempty(beam_sign)
        beam_sign = 1;
    end
    vol_out = resample_babelbrain_vol(vol_bb, x_mm, y_mm, z_mm, planimg, prestus_dims, res_mm, focus_pos, focal_length_mm, method, beam_sign);
end

function plot_babelbrain_fov_debug(p_amp_bb_native, p_amp_resampled, x_mm, y_mm, z_mm, ...
        res_mm, focus_pos, focal_length_mm, parameters, material_bb_native, material_resampled, ...
        target_loc_bb, anchor_mismatch_mm, planimg, fov_offset, beam_sign)
% PLOT_BABELBRAIN_FOV_DEBUG  Diagnose BabelBrain -> PRESTUS domain-coverage limits
%   and report the tissue-marker anchor and fallback-anchor mismatch.
%
% Two independent failure modes can look identical in the final pressure
% profile (a band of exact-zero pressure near the transducer/skull), so this
% figure separates them:
%
%   1. DOMAIN TRUNCATION — BabelBrain's returned volume genuinely doesn't cover
%      part of the PRESTUS grid (its own domain sizing / narrow-beam cropping
%      fell short). resample_babelbrain_vol zero-fills those out-of-bounds
%      voxels. Fix: extend BabelBrain's domain / disable narrow-beam cropping.
%      Shown as the shaded region in panel 3.
%
%   2. FOCAL-ANCHOR ERROR — the resample anchors the /TargetLocation tissue
%      marker (the label-5 voxel PRESTUS wrote at focus_pos) back onto
%      focus_pos, making the output tissue-aligned with medium_masks by
%      construction. anchor_mismatch_mm reports how far the passed FALLBACK
%      anchor (preset FocalLength from babelbrain_wrapper) sat from the marker
%      actually used — normally 0; a nonzero value only matters for runs where
%      /TargetLocation was missing and the fallback anchored the resample.
%
% Panel 1: BabelBrain's own pressure profile in its native grid, with z=0 at
% the marker (focal spot). Panel 2: the same axis for the MaterialMap (tissue
% index), so the skull's actual position in BabelBrain's own output is
% directly visible. Panel 3: the resampled PRESTUS-grid pressure profile with
% the zero-filled (out-of-BB-domain) region shaded.

    if nargin < 16 || isempty(beam_sign)
        beam_sign = 1;
    end

    % BabelBrain-native central-axis profile: index nearest lateral (0,0) [focal spot].
    [~, ix_bb] = min(abs(x_mm));
    [~, iy_bb] = min(abs(y_mm));
    z_rel_focus_bb = z_mm - focal_length_mm;              % depth rel. focal spot [mm], BB native
    profile_bb     = squeeze(p_amp_bb_native(ix_bb, iy_bb, :));

    has_material = ~isempty(material_bb_native) && ~isempty(material_resampled);
    if has_material
        material_line_bb = squeeze(material_bb_native(ix_bb, iy_bb, :));
    end
    has_target = nargin >= 12 && ~isempty(target_loc_bb) && numel(z_mm) >= target_loc_bb(3);
    if nargin < 13
        anchor_mismatch_mm = NaN;
    end

    % PRESTUS resampled central-axis profile, same lateral convention as
    % plot_axial_profile. beam_sign matches the resample mapping
    % (resample_babelbrain_vol: prestus_k = beam_sign*(z - anchor)/res +
    % focus_pos(3)), so depth-from-focus increases toward the far side of the
    % beam regardless of which way the beam runs along grid dim 3. Omitting it
    % (the old bug) mirror-labelled the axis for -z beams, reporting the
    % truncated side of the domain backwards.
    dims = size(p_amp_resampled);
    ix_p = max(1, min(round(focus_pos(1)), dims(1)));
    iy_p = max(1, min(round(focus_pos(2)), dims(2)));
    Nz_p = dims(3);
    z_rel_focus_prestus = beam_sign * ((1:Nz_p)' - focus_pos(3)) * res_mm;
    profile_p = squeeze(p_amp_resampled(ix_p, iy_p, :));

    % BabelBrain domain extent, relative to focal spot [mm] (z_mm may be decreasing).
    bb_z_lo = min(z_rel_focus_bb);
    bb_z_hi = max(z_rel_focus_bb);
    prestus_z_lo = min(z_rel_focus_prestus);
    prestus_z_hi = max(z_rel_focus_prestus);
    trunc_proximal_mm = max(0, bb_z_lo - prestus_z_lo);   % PRESTUS extends closer to transducer than BB does
    trunc_distal_mm   = max(0, prestus_z_hi - bb_z_hi);   % PRESTUS extends deeper than BB does

    fprintf(['[read_babelbrain_output] FOV coverage check: BabelBrain domain covers depth-from-focus ' ...
        '[%.1f, %.1f] mm; PRESTUS grid requires [%.1f, %.1f] mm.\n'], bb_z_lo, bb_z_hi, prestus_z_lo, prestus_z_hi);
    if trunc_proximal_mm > 0 || trunc_distal_mm > 0
        fprintf(['[read_babelbrain_output] -> zero-filled (no BabelBrain data): %.1f mm proximal ' ...
            '(transducer side) / %.1f mm distal (far side) of the PRESTUS grid.\n'], ...
            trunc_proximal_mm, trunc_distal_mm);
    else
        fprintf('[read_babelbrain_output] -> BabelBrain domain fully covers the PRESTUS grid; no zero-fill.\n');
    end

    out_dir = parameters.io.dir_img;
    if isfield(parameters, 'io') && isfield(parameters.io, 'dir_debug_babelbrain') ...
            && isfolder(parameters.io.dir_debug_babelbrain)
        out_dir = parameters.io.dir_debug_babelbrain;
    end

    n_panels = 2 + has_material;
    h = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 800 220*n_panels]);

    % Transducer element axial positions relative to the focal spot [mm], along
    % the central axis. The bowl apex sits focal_length_mm proximal of the focus
    % (z = -focal_length_mm); each annular ring is set back toward the focus by
    % its sagitta on the R = curv_radius_mm cap, so ring i is at
    %   z_i = -focal_length_mm + (R - sqrt(R^2 - r_i^2)),   r_i = mean ring radius.
    % (Independent of any electronic steering — the aim only moves the focus,
    % not the physical bowl geometry.)
    apex_z = NaN; elem_z = [];
    if ~isnan(focal_length_mm) && isfield(parameters, 'transducer') && ~isempty(parameters.transducer)
        tr_p = parameters.transducer(1);
        if isfield(tr_p, 'annular') && ~isempty(tr_p.annular) ...
                && all(isfield(tr_p.annular, {'elem_id_mm', 'elem_od_mm', 'curv_radius_mm'}))
            R_mm = tr_p.annular.curv_radius_mm;
            r_i  = (tr_p.annular.elem_id_mm(:) + tr_p.annular.elem_od_mm(:)) / 4;  % mean radius
            r_i  = r_i(r_i <= R_mm);   % guard against numerical/definition overshoot
            sagitta = R_mm - sqrt(R_mm.^2 - r_i.^2);
            apex_z = -focal_length_mm;
            elem_z = -focal_length_mm + sagitta;
        end
    end

    panel = 1;
    ax1 = subplot(n_panels,1,panel); panel = panel + 1;
    plot(ax1, z_rel_focus_bb, profile_bb / 1e3, 'k-', 'LineWidth', 1.5);
    xline(ax1, 0, '--b', 'LineWidth', 1.0, 'Label', 'anchor (/TargetLocation tissue marker)');
    if has_target && abs(anchor_mismatch_mm) > 1e-6
        xline(ax1, -anchor_mismatch_mm, '--r', 'LineWidth', 1.0, 'Label', 'fallback anchor (unused)');
    end
    ylabel(ax1, 'Pressure (kPa)');
    ttl = sprintf('BabelBrain FOV debug — sub-%03d%s (native grid)', ...
        parameters.subject_id, strrep(parameters.io.output_affix, '_', ' '));
    title(ax1, ttl, 'Interpreter', 'none');
    xlim(ax1, sort([z_rel_focus_bb(1) z_rel_focus_bb(end)]));
    overlay_tx_elements(ax1, apex_z, elem_z);
    grid(ax1, 'on'); box(ax1, 'on');

    if has_material
        ax_mat = subplot(n_panels,1,panel); panel = panel + 1;
        plot(ax_mat, z_rel_focus_bb, material_line_bb, 'k-', 'LineWidth', 1.5);
        xline(ax_mat, 0, '--b', 'LineWidth', 1.0);
        if has_target && abs(anchor_mismatch_mm) > 1e-6
            xline(ax_mat, -anchor_mismatch_mm, '--r', 'LineWidth', 1.0);
        end
        ylabel(ax_mat, 'BB material idx');
        ylim(ax_mat, [-0.5, max(4.5, max(material_line_bb)+0.5)]);
        title(ax_mat, 'BabelBrain MaterialMap (0=water 1=skin 2=cortical 3=trabecular 4=brain), native grid');
        xlim(ax_mat, sort([z_rel_focus_bb(1) z_rel_focus_bb(end)]));
        overlay_tx_elements(ax_mat, apex_z, elem_z);
        grid(ax_mat, 'on'); box(ax_mat, 'on');
    end

    ax2 = subplot(n_panels,1,panel);
    hold(ax2, 'on');
    yl = [0, max(1, max(profile_p) / 1e3 * 1.1)];
    if trunc_proximal_mm > 0
        patch(ax2, [prestus_z_lo, bb_z_lo, bb_z_lo, prestus_z_lo], [yl(1) yl(1) yl(2) yl(2)], ...
            [1 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.6, 'DisplayName', 'zero-filled');
    end
    if trunc_distal_mm > 0
        patch(ax2, [bb_z_hi, prestus_z_hi, prestus_z_hi, bb_z_hi], [yl(1) yl(1) yl(2) yl(2)], ...
            [1 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.6, 'HandleVisibility', 'off');
    end
    plot(ax2, z_rel_focus_prestus, profile_p / 1e3, 'k-', 'LineWidth', 1.5);
    xline(ax2, 0, '--b', 'LineWidth', 1.0);
    hold(ax2, 'off');
    xlabel(ax2, 'Depth relative to focal spot (mm)');
    ylabel(ax2, 'Pressure (kPa)');
    title(ax2, 'Resampled onto PRESTUS grid (shaded = outside BabelBrain domain, zero-filled)');
    xlim(ax2, sort([z_rel_focus_prestus(1) z_rel_focus_prestus(end)]));
    ylim(ax2, yl);
    overlay_tx_elements(ax2, apex_z, elem_z);
    if trunc_proximal_mm > 0 || trunc_distal_mm > 0
        legend(ax2, 'Location', 'northeast', 'FontSize', 7);
    end
    grid(ax2, 'on'); box(ax2, 'on');

    png_path = fullfile(out_dir, sprintf('sub-%03d_babelbrain_fov_debug%s.png', ...
        parameters.subject_id, parameters.io.output_affix));
    set(h, 'InvertHardcopy', 'off');
    if debug_on(parameters, 'babelbrain', 'figures')
        saveas(h, png_path, 'png');
        fprintf('[read_babelbrain_output] Saved: %s\n', png_path);
        debug_register(parameters, 'babelbrain', png_path, ...
            'BabelBrain domain against the PRESTUS grid along the beam (zero-filled depths shaded)', 'figure');
    end
    close(h);

    T = table(bb_z_lo, bb_z_hi, prestus_z_lo, prestus_z_hi, trunc_proximal_mm, trunc_distal_mm, anchor_mismatch_mm, ...
        'VariableNames', {'bb_depth_rel_focus_min_mm', 'bb_depth_rel_focus_max_mm', ...
                          'prestus_depth_rel_focus_min_mm', 'prestus_depth_rel_focus_max_mm', ...
                          'zero_filled_proximal_mm', 'zero_filled_distal_mm', 'fallback_anchor_mismatch_mm'});
    csv_path = fullfile(out_dir, sprintf('sub-%03d_babelbrain_fov_debug%s.csv', ...
        parameters.subject_id, parameters.io.output_affix));
    writetable(T, csv_path);
    fprintf('[read_babelbrain_output] Saved: %s\n', csv_path);
    debug_register(parameters, 'babelbrain', csv_path, ...
        'BabelBrain and PRESTUS depth ranges relative to the focus, and zero-filled depths', 'other');

    % Save the resampled MaterialMap (the BabelBrain INPUT tissue map) as a NIfTI
    % in T1w space so it overlays the T1/medium_masks in an external viewer. See
    % write_bb_t1w_debug for the alignment rationale (back-transform into T1 voxel
    % space + reuse planimg.t1_header, the same path nifti_to_t1w uses).
    if has_material && debug_on(parameters, 'babelbrain', 'volumes')
        try
            mat_path = fullfile(out_dir, sprintf('sub-%03d_babelbrain_material_map%s', ...
                parameters.subject_id, parameters.io.output_affix));
            write_bb_t1w_debug(uint8(round(material_resampled)), mat_path, planimg, fov_offset, ...
                'uint8', 8, 'nearest');
            debug_register(parameters, 'babelbrain', [mat_path '.nii.gz'], ...
                'BabelBrain input material map in T1w space', 'volume');
        catch ME
            warning('read_babelbrain_output:materialMapSaveFailed', ...
                'Could not save resampled MaterialMap NIfTI: %s', ME.message);
        end
    end
end

% -------------------------------------------------------------------------
function write_bb_t1w_debug(vol, out_path_noext, planimg, fov_offset, datatype, bits, resampler)
% Back-transform a (possibly FOV-cropped) acoustic-grid volume into T1 voxel
% space and write it as a compressed NIfTI reusing planimg.t1_header, so it
% overlays the T1/medium_masks exactly. Used for debug-mode overlays of the
% BabelBrain input (material map) and output (pressure) volumes.
%
% ALIGNMENT: rather than store a freshly-composed grid->RAS affine (fragile — a
% single voxel-centre / 0-vs-1-based / flip convention slip and it no longer
% lands on the T1), we back-transform into T1 VOXEL space and reuse
% planimg.t1_header verbatim — byte-for-byte the path nifti_to_t1w.m uses for
% medium_masks, so overlay is guaranteed by construction. The volume lives on
% the (possibly FOV-cropped) acoustic grid, so we prepend the crop offset to
% inv_transf:  crop-grid voxel -> full-grid voxel (offset) -> T1 voxel.
    if nargin < 7 || isempty(resampler), resampler = 'linear'; end
    have_t1 = ~isempty(planimg) && isfield(planimg, 't1_header') && ~isempty(planimg.t1_header) ...
        && isfield(planimg, 'inv_transf') && ~isempty(planimg.inv_transf) ...
        && isfield(planimg, 't1_image_orig') && ~isempty(planimg.t1_image_orig);
    if ~have_t1
        warning('read_babelbrain_output:noT1forDebug', ...
            ['planimg lacks t1_header/inv_transf/t1_image_orig; skipping T1w debug NIfTI ' ...
             '%s (would not align with the T1).'], out_path_noext);
        return;
    end
    M_crop_to_full = eye(4);
    M_crop_to_full(4, 1:3) = fov_offset(:)' - 1;         % crop voxel -> full-grid voxel
    M_crop_to_t1 = M_crop_to_full * planimg.inv_transf;  % -> T1 voxel (row-vector convention)
    vol_t1 = affine_resample_3d(single(vol), M_crop_to_t1, size(planimg.t1_image_orig), resampler, 0);
    if strcmp(datatype, 'uint8')
        vol_t1 = uint8(round(vol_t1));
    else
        vol_t1 = cast(vol_t1, datatype);
    end
    hdr = planimg.t1_header;
    hdr.Datatype = datatype;
    if ~isempty(bits), hdr.BitsPerPixel = bits; end
    niftiwrite(vol_t1, out_path_noext, hdr, 'Compressed', true);
    fprintf('[read_babelbrain_output] Saved T1w debug: %s.nii.gz\n', out_path_noext);
end

% -------------------------------------------------------------------------
function overlay_tx_elements(ax, apex_z, elem_z)
% Overlay the transducer geometry on a depth-profile axis (x = depth relative
% to the focal spot, mm): a magenta dashed line at the bowl apex and magenta
% markers at each annular ring's central-axis position. No-op when geometry is
% unavailable (apex_z NaN / elem_z empty). The apex typically sits well proximal
% of a partial-FOV domain, so the x-axis is extended leftward as needed so the
% transducer (and the water standoff between it and the skull) stays visible.
    if isnan(apex_z) || isempty(elem_z)
        return;
    end
    washold = ishold(ax);
    hold(ax, 'on');

    % Extend the left x-limit to include the transducer, with a small margin.
    xl = xlim(ax);
    left_needed = min([apex_z; elem_z(:)]) - 3;
    if left_needed < xl(1)
        xlim(ax, [left_needed, xl(2)]);
    end

    yl = ylim(ax);
    y_mark = yl(1) + 0.06 * (yl(2) - yl(1));
    % Ring markers along the near-bottom of the axis.
    plot(ax, elem_z, y_mark * ones(size(elem_z)), 'v', ...
        'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'm', 'MarkerSize', 5, ...
        'DisplayName', 'Tx elements');
    % Apex reference line.
    xline(ax, apex_z, ':m', 'LineWidth', 1.0, 'Label', 'Tx apex', ...
        'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right', ...
        'HandleVisibility', 'off');
    if ~washold
        hold(ax, 'off');
    end
end
