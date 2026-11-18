function combine_async_intensity(files_in, file_out, targets_wcm2, mode, groups, pd_s)
% COMBINE_ASYNC_INTENSITY  Scale and combine N independent acoustic results for async (non-simultaneous) transducer firing
%
% For asynchronously firing transducers the pressure fields are independent.
% Two combined fields are produced from the per-transducer scaled pressures:
%
%   p_i_scaled = p_i · sqrt(target_i / baseline_i)
%
%   Thermal field  — incoherent intensity sum, drives the bioheat equation:
%     I_total   = Σ  p_i_scaled² / (2·ρ·c)
%     p_thermal = sqrt(2·I_total·ρ·c) = sqrt(Σ p_i_scaled²)   [ρ·c cancels]
%                                                → sensor_data.p_max_all
%
%   Acoustic safety field  — per-voxel pressure maximum, reflects the
%     instantaneous peak when only one transducer is active at a time:
%     p_safety = max_i( p_i_scaled )            → sensor_data.p_max_async
%
%   Per-transducer fields  — p_i_scaled for each transducer on the full
%     reference grid                           → sensor_data.p_async_tx
%     [Nx Ny Nz N]. THERMAL_SIMULATION heats with these, each switched on its
%     own schedule (THERMAL_SCHEDULE, timing.async_schedule).
%
% Rationale: in async mode the instantaneous peak intensity at any voxel
% equals the maximum over individual transducer fields, not their sum.
% p_max_all, the intensity sum, is the heat source of all transducers on at
% once; the thermal stage uses p_async_tx and the schedule instead. Acoustic
% safety metrics (ISPPA, MI, peak pressure) must be derived from p_max_async,
% not from p_max_all.
%
% Firing groups (groups argument, FIRING_GROUPS)
% ----------------------------------------------
% Transducers with the same group index fire together. Acoustics is linear,
% so a group's pressure is p(t) = Re(Σ_f P_f·exp(1i·ω_f·t)), with P_f the
% coherent sum of the scaled complex fields of its transducers at frequency f
% (they interfere; each is scaled to its own target first). The cross terms
% between frequencies average out over their common period. Per group:
%   time-averaged intensity  Σ_f |P_f|² / (2·ρ·c)
%   heat sources             one field |P_f| per frequency, all in the group's
%                            windows (p_async_tx, with p_async_freqs_hz and
%                            p_async_source_tx: the group's first transducer)
%   peak pressure            |P_f| at one frequency; at several, max / min over
%                            t of p(t) (SIMULTANEOUS_PEAKS), t over the pulse or
%                            over the common period 1/gcd(f), whichever is shorter
% Groups then combine like async transducers: p_max_async / p_min_async are
% the per-voxel extremes over groups, p_max_all the intensity sum of every
% source, and sensor_data.p_intensity_equiv the largest time-averaged
% intensity of any group (ACOUSTIC_ANALYSIS takes Ipa/Ita from it). Group
% members need their complex fields (io.save_p_complex of the acoustic runs,
% phases referred to t = 0) and one shared grid (async_layout 'shared').
%
% Simultaneous mode (mode = 'simultaneous') is one group of all transducers:
% transducers at more than one frequency firing at the same time
% (IS_SIMULTANEOUS_MULTIFREQ). When the common period needs more than 20000
% time samples, the peaks are the upper bounds ±Σ_f|P_f|.
%
% Partial-FOV support
% -------------------
% A cache produced with ACOUSTIC_GRID_FOV (acoustic_provenance.fov_offset_ac
% non-empty) is zero-padded (back-projected) onto its own full grid
% (full_ac_dims) first. In shared grid mode that full grid is the reference
% grid; in per_transducer mode it is then remapped like any other grid (below).
% FOV and full-grid caches can be mixed. The combined output is always on the
% full grid; fov_offset_ac is cleared in the combined provenance, and kgrid and
% acoustic_grid_dims describe the full grid, so the downstream thermal
% pipeline needs no changes (kwave_medium stays file 1's crop).
%
% Grid remapping
% --------------
% When a cache (after any FOV back-projection) is on a different grid than the reference
% (file 1) — e.g. grid.async_layout = 'per_transducer', where each
% transducer is simulated on its own aligned/cropped grid — the field is
% resampled onto the reference physical space (transducer 1's grid). Two
% remap strategies are used, in order of preference:
%   1. World-space affine (rotation + translation aware): each run's planimg
%      sidecar (sim_voxel -> world RAS mm) is composed into a direct
%      sim_voxel_k -> sim_voxel_ref transform and applied via
%      AFFINE_RESAMPLE_3D. Correct for grids that differ in orientation,
%      crop offset, and/or extent (the per_transducer case). Before use, the
%      composed transform is checked against a physical invariant (see
%      CHECK_RIGID_DISTANCE): it must preserve each transducer's own
%      trans_pos-to-focus_pos distance, since both grids share the same
%      voxel resolution and trans_pos/focus_pos describe a fixed physical
%      transducer. A failed check falls back to strategy 2 with a warning
%      rather than silently combining a misaligned field.
%   2. Legacy origin-centred spacing-only remap (used only if a planimg
%      sidecar cannot be recovered, e.g. a legacy cache or a deleted
%      sidecar, or if strategy 1 fails its rigid-distance check): assumes
%      both grids are centred at the origin and differ only in
%      spacing/extent, not rotation.
% Regions of the source field that fall outside the reference extent are
% zero-padded.  This avoids forcing a re-run when the grid change is
% recoverable by interpolation. Complex fields of group members are not
% remapped (a group needs one shared grid).
%
% The target for each transducer is resolved in this order:
%   1. targets_wcm2(i) argument (if supplied and finite)
%   2. acoustic_info.parameters.transducer(1).target_isppa_wcm2 in the cache
%   3. No scaling (use the simulated drive level as-is)
%
% The baseline (free-water ISPPA at the simulated drive level) is read from
% each cache's acoustic_provenance (ISPPA_SCALE_FACTOR: freefield_isppa_wcm2,
% else the library calibration value). Without a baseline, a target that
% needs scaling is an error.
%
% Grid/medium metadata for the combined cache are taken from files_in{1}.
%
% Use as:
%   combine_async_intensity(files_in, file_out)
%   combine_async_intensity(files_in, file_out, targets_wcm2)
%   combine_async_intensity(files_in, file_out, targets_wcm2, mode)
%   combine_async_intensity(files_in, file_out, targets_wcm2, 'async', groups)
%   combine_async_intensity(files_in, file_out, targets_wcm2, 'async', groups, pd_s)
%
% Input:
%   files_in     - cell array of paths to acoustic cache .mat files {N×1}
%   file_out     - path where the combined cache will be written
%   targets_wcm2 - (optional) 1×N vector of target ISPPAs [W/cm²], one per
%                  transducer.  NaN or [] entries fall back to the cached
%                  transducer config.
%   mode         - (optional) 'async' (default): separate firing windows;
%                  'simultaneous': all transducers fire together (above)
%   groups       - (optional) 1×N firing group of each transducer
%                  (FIRING_GROUPS.group); default 1:N ('async') or all ones
%                  ('simultaneous')
%   pd_s         - (optional) 1×N pulse duration of each transducer [s], over
%                  which a group's peak pressure is sampled; NaN or [] = read
%                  from the cache's parameters. files_in may name one file
%                  several times (transducers with the same field).
%
% Notes:
%   - Axisymmetric (2-D) results are not supported.
%   - All inputs must share the same full acoustic grid (full_ac_dims for
%     FOV runs, size(p_max_all) for non-FOV runs).
%   - acoustic_provenance in the combined file records combined_from,
%     combination_mode, per-transducer targets, baselines, groups, and FOV offsets.
%
% See also: ASYNC_TRANSDUCER_PIPELINE, FIRING_GROUPS, ASSEMBLE_LIMITED_FOV_FIELDS,
%           APPLY_ISPPA_SCALING, ACOUSTIC_WRAPPER

arguments
    files_in
    file_out     (1,:) char
    targets_wcm2 (1,:) double = []
    mode         (1,:) char {mustBeMember(mode, {'async', 'simultaneous'})} = 'async'
    groups       (1,:) double = []
    pd_s         (1,:) double = []
end

if ischar(files_in) || isstring(files_in)
    files_in = cellstr(files_in);
end
if ~iscell(files_in)
    error('combine_async_intensity:badInput', 'files_in must be a cell array of file paths.');
end
files_in = files_in(:);
N = numel(files_in);
if N < 2
    error('combine_async_intensity:tooFewFiles', 'At least two acoustic cache files are required.');
end
for k = 1:N
    if ~isfile(files_in{k})
        error('combine_async_intensity:missingFile', 'Acoustic cache not found: %s', files_in{k});
    end
end

% Firing groups: renumber in order of first appearance (group of file 1 = 1)
if isempty(groups)
    if strcmp(mode, 'simultaneous')
        groups = ones(1, N);
    else
        groups = 1:N;
    end
elseif numel(groups) ~= N
    error('combine_async_intensity:badGroups', ...
        'groups must hold one group index per acoustic cache (%d), got %d.', N, numel(groups));
end
[~, ~, groups] = unique(groups, 'stable');
groups   = groups(:)';
G        = max(groups);
n_in_g   = accumarray(groups(:), 1)';
coherent = n_in_g > 1;                       % group of several transducers: complex sum
any_coherent = any(coherent);

% Pad or trim the external targets vector to length N.
if numel(targets_wcm2) < N
    targets_wcm2(end+1:N) = NaN;
end

% =========================================================================
%% Load first file — provides grid, medium, and all metadata
% =========================================================================

fprintf('Loading acoustic result 1/%d:\n  %s\n', N, files_in{1});
data_ref = load(files_in{1}, 'sensor_data', 'kgrid', 'kwave_medium', ...
                'acoustic_provenance', 'acoustic_info', 'acoustic_grid_dims');
p_size_ref = size(data_ref.sensor_data.p_max_all);
if numel(p_size_ref) < 3
    error('combine_async_intensity:axisymmetricNotSupported', ...
        ['Axisymmetric (2-D) results cannot be combined in async mode. ' ...
         'Re-run acoustic simulations with a 3-D grid.']);
end
% The full grid of file 1 is the reference grid. For non-FOV runs the full
% grid IS the sensor_data grid and the offset is [1,1,1].
[full_ac_dims, ~, has_fov_ref] = resolve_grid(data_ref.acoustic_provenance, p_size_ref, 1);

p_sq_sum     = zeros(full_ac_dims);      % Σ over sources of field²: intensity sum
p_max_safety = zeros(full_ac_dims);      % per-voxel maximum over groups
p_min_safety = zeros(full_ac_dims);      % per-voxel minimum over groups
p_ieq_sq     = zeros(full_ac_dims);      % largest time-averaged intensity of any group (× 2ρc)
src_field = cell(1, 0);  src_freq = zeros(1, 0);  src_tx = zeros(1, 0);  src_group = zeros(1, 0);
P_acc = cell(1, G);  F_acc = cell(1, G);  pd_grp = NaN(1, G);  first_tx = zeros(1, G);
for g = 1:G
    P_acc{g} = cell(1, 0);
    F_acc{g} = zeros(1, 0);
    first_tx(g) = find(groups == g, 1);
end

targets_log     = zeros(1, N);
baselines_log   = zeros(1, N);
fov_offsets_log = cell(1, N);
isppa_per_t     = zeros(1, N);   % per-transducer ISPPA diagnostics (W/cm²)

% =========================================================================
%% Scale each transducer and put it on the reference grid
% =========================================================================

for k = 1:N
    if k == 1
        data_k = data_ref;
    else
        fprintf('Loading acoustic result %d/%d:\n  %s\n', k, N, files_in{k});
        data_k = load(files_in{k}, 'sensor_data', 'kgrid', 'acoustic_provenance', 'acoustic_info', ...
                      'acoustic_grid_dims');
    end
    g   = groups(k);
    p_k = double(abs(data_k.sensor_data.p_max_all));
    if ndims(p_k) < 3 %#ok<ISMAT>
        error('combine_async_intensity:axisymmetricNotSupported', ...
            ['Axisymmetric (2-D) results cannot be combined in async mode. ' ...
             'Re-run acoustic simulations with a 3-D grid.']);
    end

    [full_ac_dims_k, fov_offset_k, has_fov_k] = ...
        resolve_grid(data_k.acoustic_provenance, size(p_k), k);

    [p_k_scaled, target_used, baseline_used] = scale_pressure( ...
        p_k, data_k.acoustic_provenance, data_k.acoustic_info, targets_wcm2(k), k);

    scale_k = p_k_scaled ./ max(p_k, eps('double'));
    s_k     = max(p_k_scaled(:)) / max(max(p_k(:)), eps('double'));   % uniform scale
    clear p_k;
    if coherent(g)
        if ~isequal(full_ac_dims_k, full_ac_dims)
            error('combine_async_intensity:simultaneousGridMismatch', ...
                ['File %d is on another grid; transducers firing together need one shared grid ' ...
                 '(grid.async_layout: shared).'], k);
        end
        [P_k, f_k] = load_complex_field(files_in{k});
        P_k = backproject_complex(P_k * s_k, fov_offset_k, full_ac_dims);
        i_f = find(abs(F_acc{g} - f_k) <= 1e-6 * f_k, 1);
        if isempty(i_f)
            P_acc{g}{end+1} = P_k;
            F_acc{g}(end+1) = f_k;
        else
            P_acc{g}{i_f} = P_acc{g}{i_f} + P_k;
        end
        clear P_k;
        if isnan(pd_grp(g))
            % The caller's value: a field shared by several entries
            % (ACOUSTIC_STORE_KEY) carries the timing of the entry that simulated it
            if k <= numel(pd_s) && isfinite(pd_s(k))
                pd_grp(g) = pd_s(k);
            else
                pd_grp(g) = pulse_duration(data_k.acoustic_info);
            end
        end
    end
    p_min_k = double(data_k.sensor_data.p_min_all) .* scale_k;
    clear scale_k;

    % A partial-FOV field is first put back onto its own full grid (the
    % grid its planimg affine describes), so it can be combined like a
    % full-grid run below: directly when that grid is the reference grid
    % (shared mode), else after remapping (per_transducer mode).
    if has_fov_k
        p_k_scaled = backproject(p_k_scaled, fov_offset_k, full_ac_dims_k);
        p_min_k    = backproject(p_min_k,    fov_offset_k, full_ac_dims_k);
        fov_offset_k = [1, 1, 1];
    end
    % Grid remapping: if this file was simulated on a different grid than the
    % reference (e.g. grid.async_layout = 'per_transducer', where
    % each transducer has its own aligned/cropped grid), resample it onto the
    % reference physical space. Prefer the true world-space affine (handles
    % rotation + translation between grids, e.g. per_transducer mode);
    % fall back to the legacy origin-centred spacing-only remap when the
    % affine cannot be recovered (e.g. deleted planimg sidecar, legacy cache).
    if ~isequal(full_ac_dims_k, full_ac_dims)
        fprintf(['[combine_async_intensity] File %d grid [%s] differs from ' ...
                 'reference [%s] — remapping via interpolation.\n'], ...
            k, num2str(full_ac_dims_k), num2str(full_ac_dims));
        [M_k2ref, have_affine] = resolve_affine_remap(data_k.acoustic_info, data_ref.acoustic_info, ...
                                                       data_k.kgrid);
        if have_affine
            p_k_scaled_new = affine_resample_3d(p_k_scaled, M_k2ref, full_ac_dims, 'linear', 0);
            if any(p_k_scaled_new(:) ~= 0)
                p_k_scaled = p_k_scaled_new;
                p_min_k    = affine_resample_3d(p_min_k, M_k2ref, full_ac_dims, 'linear', 0);
                fprintf(['[combine_async_intensity] File %d remapped via world-space affine ' ...
                         '(rotation/crop aware).\n'], k);
            else
                % Resampled field is all-zero: the affine is almost certainly wrong
                % (e.g. no spatial overlap or a convention mismatch) — fall back to
                % the legacy remap rather than silently combining a blank field.
                have_affine = false;
                warning('combine_async_intensity:emptyAffineRemap', ...
                    ['File %d: world-space affine remap produced an all-zero field — ' ...
                     'falling back to legacy origin-centred remap. Verify grid alignment.'], k);
            end
        end
        if ~have_affine
            p_k_scaled = remap_field_to_ref(p_k_scaled, data_k.kgrid, full_ac_dims_k, ...
                                            data_ref.kgrid, full_ac_dims);
            p_min_k    = remap_field_to_ref(p_min_k, data_k.kgrid, full_ac_dims_k, ...
                                            data_ref.kgrid, full_ac_dims);
        end
        fov_offset_k = [1, 1, 1];
    end
    p_k_full = backproject(p_k_scaled, fov_offset_k, full_ac_dims);

    if ~coherent(g)
        % A transducer firing on its own is one heat source and one safety field
        p_sq_sum     = p_sq_sum + p_k_full .^ 2;
        p_ieq_sq     = max(p_ieq_sq, p_k_full .^ 2);
        p_max_safety = max(p_max_safety, p_k_full);
        p_min_safety = min(p_min_safety, backproject(p_min_k, fov_offset_k, full_ac_dims));
        src_field{end+1} = single(p_k_full);                    %#ok<AGROW>
        src_freq(end+1)  = transducer_freq(data_k.acoustic_info); %#ok<AGROW>
        src_tx(end+1)    = k;                                   %#ok<AGROW>
        src_group(end+1) = g;                                   %#ok<AGROW>
    end

    isppa_per_t(k)     = max(p_k_full(:).^2, [], 'all') / (2 * 1000 * 1500) * 1e-4;
    targets_log(k)     = target_used;
    baselines_log(k)   = baseline_used;
    fov_offsets_log{k} = fov_offset_k;
    clear data_k p_k_scaled p_min_k p_k_full;
end

% =========================================================================
%% Firing groups: coherent sum per frequency
% =========================================================================
peak_info = struct('group', {}, 'info', {});
for g = find(coherent)
    [F_acc{g}, order] = sort(F_acc{g});
    P_g = P_acc{g}(order);
    P_acc{g} = [];
    I_g = zeros(full_ac_dims);
    for i_f = 1:numel(P_g)
        I_g = I_g + abs(P_g{i_f}) .^ 2;
        src_field{end+1} = single(abs(P_g{i_f}));   %#ok<AGROW>
        src_freq(end+1)  = F_acc{g}(i_f);          %#ok<AGROW>
        src_tx(end+1)    = first_tx(g);            %#ok<AGROW>
        src_group(end+1) = g;                      %#ok<AGROW>
    end
    if numel(P_g) == 1
        pk_pos = abs(P_g{1});                      % one tone: the peak is its amplitude
        pk_neg = -pk_pos;
    else
        [pk_pos, pk_neg, info_g] = simultaneous_peaks(P_g, F_acc{g}, pd_grp(g));
        peak_info(end+1) = struct('group', g, 'info', info_g); %#ok<AGROW>
    end
    clear P_g;
    p_sq_sum     = p_sq_sum + I_g;
    p_ieq_sq     = max(p_ieq_sq, I_g);
    p_max_safety = max(p_max_safety, pk_pos);
    p_min_safety = min(p_min_safety, pk_neg);
    fprintf('  Group %d: transducers %s fire together (%s kHz), summed coherently\n', ...
        g, mat2str(find(groups == g)), num2str(F_acc{g} / 1e3));
    clear I_g pk_pos pk_neg;
end

% Heat sources ordered by group, then frequency (FIRING_GROUPS)
[~, order]  = sortrows([src_group(:), src_freq(:)]);
src_field   = src_field(order);
src_freq    = src_freq(order);
src_tx      = src_tx(order);
p_async_tx  = cat(4, src_field{:});
clear src_field;

% =========================================================================
%% Per-transducer ISPPA summary
% =========================================================================
fprintf('Per-transducer ISPPAs at the combined (target) scale:\n');
for ti = 1:N
    fprintf('  Transducer %d: %.4f W/cm²\n', ti, isppa_per_t(ti));
end
if any_coherent
    fprintf('  Joint ISPPA (largest time-averaged intensity of any group): %.4f W/cm²\n', ...
        max(p_ieq_sq(:)) / (2 * 1000 * 1500) * 1e-4);
else
    fprintf('  Expected joint ISPPA (= max of above): %.4f W/cm²\n', max(isppa_per_t));
end

% =========================================================================
%% Thermal field
% =========================================================================

% p_thermal = sqrt(Σ p_i²)  —  ρ·c cancels in I = p²/(2ρc) → p = sqrt(2Iρc)
p_thermal = sqrt(p_sq_sum);
clear p_sq_sum;

% =========================================================================
%% Assemble combined sensor_data and provenance
% =========================================================================

sensor_data              = data_ref.sensor_data;
sensor_data.p_max_all    = p_thermal;    % thermal heat source (intensity sum)
sensor_data.p_max_async  = p_max_safety; % acoustic safety: peak positive
sensor_data.p_min_async  = p_min_safety; % acoustic safety: peak negative
sensor_data.p_async_tx   = p_async_tx;   % scaled field per heat source [Nx Ny Nz S]: scheduled heat source
if any_coherent
    % time-averaged intensity of the transducers firing together (ACOUSTIC_ANALYSIS: Ipa, Ita)
    sensor_data.p_intensity_equiv = single(sqrt(p_ieq_sq));
    sensor_data.p_async_freqs_hz  = src_freq;  % frequency of each heat source
    sensor_data.p_async_source_tx = src_tx;    % its windows: those of this transducer (THERMAL_SCHEDULE)
end
clear p_thermal p_max_safety p_min_safety p_async_tx p_ieq_sq;

acoustic_provenance                                        = data_ref.acoustic_provenance;
acoustic_provenance.combined_from                          = files_in;
acoustic_provenance.combination_mode                       = 'async_dual_field';
if any_coherent
    if G == 1
        acoustic_provenance.combination_mode               = 'async_dual_field_simultaneous';
    else
        acoustic_provenance.combination_mode               = 'async_dual_field_grouped';
    end
    acoustic_provenance.combination_groups                 = groups;   % per transducer
    acoustic_provenance.combination_frequency_groups_hz    = src_freq; % fields of p_async_tx
    acoustic_provenance.combination_peak                   = peak_info;
end
acoustic_provenance.combination_targets                    = targets_log;
acoustic_provenance.combination_baselines                  = baselines_log;
acoustic_provenance.combination_isppa_per_transducer_wcm2  = isppa_per_t;
acoustic_provenance.combination_per_transducer             = true;   % sensor_data.p_async_tx present
acoustic_provenance.combination_time                       = datestr(now); %#ok<TNOW1,DATST>
% The reference file's free-water axial profile describes that single run at
% its own drive level; it does not match the combined, per-source scaled field.
if isfield(acoustic_provenance, 'freefield_axial')
    acoustic_provenance = rmfield(acoustic_provenance, 'freefield_axial');
end
if has_fov_ref
    % Combined output is on the full grid — clear per-file FOV crop metadata
    % so the downstream thermal pipeline does not attempt a second back-projection.
    acoustic_provenance.combination_mode       = [acoustic_provenance.combination_mode, '_fov_backprojected'];
    acoustic_provenance.combination_fov_offsets = fov_offsets_log;
    acoustic_provenance.combination_full_ac_dims = full_ac_dims;
    acoustic_provenance.fov_offset_ac          = [];
    acoustic_provenance.full_ac_dims           = full_ac_dims;
end

% =========================================================================
%% Save
% =========================================================================

fprintf('Saving combined acoustic cache:\n  %s\n', file_out);

kgrid         = data_ref.kgrid;
if ~isequal([kgrid.Nx, kgrid.Ny, kgrid.Nz], full_ac_dims)
    % File 1 was a partial-FOV run: its kgrid is the crop. The thermal stage
    % uses this kgrid with the full-grid field, so rebuild it at full size.
    kgrid_fov = kgrid;
    kgrid = kWaveGrid(full_ac_dims(1), kgrid_fov.dx, full_ac_dims(2), kgrid_fov.dy, ...
                      full_ac_dims(3), kgrid_fov.dz);
    kgrid.setTime(kgrid_fov.Nt, kgrid_fov.dt);
    clear kgrid_fov;
end
kwave_medium  = data_ref.kwave_medium;  %#ok<NASGU>
acoustic_info = data_ref.acoustic_info; %#ok<NASGU>
% The grid the combined field is on (the full grid): CACHE_ACOUSTIC compares it
% with the loading run's grid. File 1's own value is its partial-FOV crop.
acoustic_grid_dims = full_ac_dims; %#ok<NASGU>
save(file_out, ...
    'sensor_data', 'kgrid', 'kwave_medium', 'acoustic_info', ...
    'acoustic_provenance', 'acoustic_grid_dims', '-v7.3');

fprintf('Done — combined cache written from %d transducers.\n', N);
end

% =========================================================================
%% Local helpers
% =========================================================================

function f = transducer_freq(acoustic_info)
% Frequency of a single-transducer run [Hz], NaN when the cache does not record it
f = NaN;
if isstruct(acoustic_info) && isfield(acoustic_info, 'parameters') && ...
        isfield(acoustic_info.parameters, 'transducer') && ~isempty(acoustic_info.parameters.transducer) && ...
        isfield(acoustic_info.parameters.transducer(1), 'freq_hz') && ...
        ~isempty(acoustic_info.parameters.transducer(1).freq_hz)
    f = double(acoustic_info.parameters.transducer(1).freq_hz);
end
end

function pd = pulse_duration(acoustic_info)
% Pulse duration of a run [s]: transducer(1).timing.pd, else timing.pd; NaN if unknown
pd = NaN;
if ~isstruct(acoustic_info) || ~isfield(acoustic_info, 'parameters')
    return
end
p = acoustic_info.parameters;
if isfield(p, 'transducer') && ~isempty(p.transducer) && isfield(p.transducer(1), 'timing') && ...
        isstruct(p.transducer(1).timing) && isfield(p.transducer(1).timing, 'pd') && ...
        ~isempty(p.transducer(1).timing.pd)
    pd = double(p.transducer(1).timing.pd);
elseif isfield(p, 'timing') && isfield(p.timing, 'pd') && ~isempty(p.timing.pd)
    pd = double(p.timing.pd);
end
end

function [full_dims, fov_offset, has_fov] = resolve_grid(prov, sensor_size, file_idx)
% Return the full acoustic grid dimensions and the FOV start offset for one
% per-transducer cache.  For runs without FOV cropping, full_dims equals
% sensor_size and offset is [1,1,1] (i.e. the field already covers the full
% grid).  For FOV runs, full_dims and fov_offset come from provenance.
has_fov = isfield(prov, 'fov_offset_ac') && ~isempty(prov.fov_offset_ac);
if has_fov
    if ~isfield(prov, 'full_ac_dims') || isempty(prov.full_ac_dims)
        error('combine_async_intensity:missingFullDims', ...
            'File %d has fov_offset_ac but full_ac_dims is missing from provenance.', ...
            file_idx);
    end
    fov_offset = prov.fov_offset_ac(:)';   % [1×3]
    full_dims  = prov.full_ac_dims(:)';    % [1×3]

    fov_end = fov_offset + sensor_size - 1;
    if any(fov_end > full_dims)
        error('combine_async_intensity:fovExceedsFullGrid', ...
            ['File %d: FOV field [size %s, offset %s] extends beyond full_ac_dims [%s]. ' ...
             'Provenance may be corrupt — re-run the acoustic stage.'], ...
            file_idx, num2str(sensor_size), num2str(fov_offset), num2str(full_dims));
    end
else
    fov_offset = [1, 1, 1];
    full_dims  = sensor_size;
end
end

function vol = backproject(p_fov, fov_offset, full_dims)
% Zero-pad p_fov into a full_dims volume starting at fov_offset.
% When fov_offset=[1,1,1] and size(p_fov)==full_dims this is a no-op copy.
if isequal(fov_offset, [1,1,1]) && isequal(size(p_fov), full_dims)
    vol = p_fov;
    return;
end
vol = zeros(full_dims, 'like', p_fov);
fov_end = fov_offset + size(p_fov) - 1;
vol(fov_offset(1):fov_end(1), fov_offset(2):fov_end(2), fov_offset(3):fov_end(3)) = p_fov;
end

function [M, ok] = resolve_affine_remap(info_k, info_ref, kgrid_k)
% Compose the forward affine (row-vector convention: [out 1] = [in 1] * M)
% mapping file k's sim voxel grid onto the reference (file 1) sim voxel grid,
% via each run's planimg sidecar. Returns ok=false if either sidecar cannot
% be recovered, OR if the composed transform fails CHECK_RIGID_DISTANCE (a
% real physical sanity check, not a heuristic) -- in either case the caller
% falls back to the legacy spacing-only remap rather than trust a possibly
% -wrong affine.
    M  = eye(4);
    ok = false;
    A_k   = grid_world_affine(info_k);
    A_ref = grid_world_affine(info_ref);
    if isempty(A_k) || isempty(A_ref)
        return;
    end
    M = A_k / A_ref;   % sim_voxel_k -> world_RAS -> sim_voxel_ref

    dx_mm = kgrid_k.dx * 1e3;   % k-Wave kgrid.dx is in metres
    [rigid_ok, residual_mm, d_native_mm] = check_rigid_distance(info_k, M, dx_mm);
    if ~rigid_ok
        if isnan(residual_mm)
            % No trans_pos/focus_pos available to check against -- proceed
            % without the extra guard (affine itself was still recoverable).
            ok = true;
            return;
        end
        warning('combine_async_intensity:affineNotRigid', ...
            ['Composed grid-to-grid affine failed a rigid-transform sanity check: ' ...
             'transducer trans_pos-focus_pos distance changed by %.2f mm after mapping ' ...
             '(native %.2f mm). Both grids share the same voxel resolution, so a true ' ...
             'rigid transform must preserve this distance -- falling back to the legacy ' ...
             'origin-centred remap instead of trusting a likely-wrong affine.'], ...
            residual_mm, d_native_mm);
        return;
    end
    ok = true;
end

function [ok, residual_mm, d_native_mm] = check_rigid_distance(info_k, M, dx_mm)
% Verify M preserves the trans_pos-to-focus_pos distance of transducer k --
% a geometric invariant that MUST hold for any true rigid (rotation +
% translation, no scale/shear) transform between two grids of identical
% voxel resolution, since trans_pos/focus_pos describe one fixed physical
% transducer, not something that legitimately differs between per-transducer
% grids. A large residual after mapping through M indicates the composed
% affine itself is wrong (e.g. an indexing/convention bug), not a genuine
% geometric difference.
%
% Returns ok=false with residual_mm=NaN if trans_pos/focus_pos are
% unavailable (nothing to check against); callers should treat that as
% "cannot validate" rather than "failed validation".
    ok = false;
    residual_mm = NaN;
    d_native_mm = NaN;
    if ~isstruct(info_k) || ~isfield(info_k, 'parameters') || ...
            ~isfield(info_k.parameters, 'transducer') || isempty(info_k.parameters.transducer)
        return;
    end
    tr = info_k.parameters.transducer(1);
    if ~isfield(tr, 'trans_pos') || ~isfield(tr, 'focus_pos') || ...
            isempty(tr.trans_pos) || isempty(tr.focus_pos)
        return;
    end
    pts_native  = [tr.trans_pos(:)'; tr.focus_pos(:)'];
    d_native_mm = norm(pts_native(2,:) - pts_native(1,:)) * dx_mm;
    if d_native_mm < eps
        d_native_mm = NaN;
        return;
    end

    pts_mapped = affine_apply_pts(pts_native, M);
    d_mapped_mm = norm(pts_mapped(2,:) - pts_mapped(1,:)) * dx_mm;

    residual_mm = abs(d_mapped_mm - d_native_mm);
    ok = residual_mm < max(2 * dx_mm, 0.05 * d_native_mm);  % ~2 voxels or 5% tolerance
end

function p_out = remap_field_to_ref(p_src, kgrid_src, dims_src, kgrid_ref, dims_ref)
% Resample p_src (size dims_src, physical spacing from kgrid_src) onto the
% reference physical space (size dims_ref, spacing from kgrid_ref) via
% trilinear interpolation.  Both grids are assumed to be centred at the
% origin.  Regions of p_src that fall outside the reference extent are zero.
%
% Physical coordinate vectors for each axis: centred at 0 with the grid's
% voxel spacing.  kgrid.dx/.dy/.dz give the spacing; dims_src/ref give the
% number of voxels (the kgrid stored in the cache is the PML-expanded grid,
% so its N* properties are not used — only the spacings).
x_src = ((0:dims_src(1)-1) - (dims_src(1)-1)/2) * kgrid_src.dx;
y_src = ((0:dims_src(2)-1) - (dims_src(2)-1)/2) * kgrid_src.dy;
z_src = ((0:dims_src(3)-1) - (dims_src(3)-1)/2) * kgrid_src.dz;

x_ref = ((0:dims_ref(1)-1) - (dims_ref(1)-1)/2) * kgrid_ref.dx;
y_ref = ((0:dims_ref(2)-1) - (dims_ref(2)-1)/2) * kgrid_ref.dy;
z_ref = ((0:dims_ref(3)-1) - (dims_ref(3)-1)/2) * kgrid_ref.dz;

[Xr, Yr, Zr] = ndgrid(x_ref, y_ref, z_ref);
p_out = interpn(x_src, y_src, z_src, double(p_src), Xr, Yr, Zr, 'linear', 0);
end

function [p_scaled, target_used, baseline_used] = scale_pressure( ...
        p, acoustic_provenance, acoustic_info, target_arg, transducer_idx)
% Scale p to target_arg [W/cm²], reading baseline from provenance and
% falling back to the cached transducer config for the target when needed.

    % Resolve target
    target = NaN;
    if isfinite(target_arg)
        target = target_arg;
    elseif isstruct(acoustic_info) && isfield(acoustic_info, 'parameters') && ...
            isfield(acoustic_info.parameters, 'transducer') && ...
            ~isempty(acoustic_info.parameters.transducer) && ...
            isfield(acoustic_info.parameters.transducer(1), 'target_isppa_wcm2')
        t = acoustic_info.parameters.transducer(1).target_isppa_wcm2;
        if isscalar(t) && isfinite(t)
            target = t;
        end
    end

    % Baseline lookup, the no-baseline error and the large-scale warning are
    % shared with APPLY_ISPPA_SCALING (ISPPA_SCALE_FACTOR). A library calibration
    % value is accepted as the baseline here, as these runs record one per transducer.
    [scale_p, baseline, baseline_src] = isppa_scale_factor(target, acoustic_provenance, ...
        AllowCalibrationBaseline=true, Context=sprintf('transducer %d', transducer_idx));
    target_used   = target;
    baseline_used = baseline;

    if ~isfinite(target)
        fprintf('  Intensity scaling (transducer %d): no target_isppa_wcm2 specified — using simulated drive level.\n', ...
            transducer_idx);
        p_scaled = p;
        return;
    end
    if scale_p == 1
        p_scaled = p;
        return;
    end
    fprintf('  Transducer %d: %.2f → %.2f W/cm² (scale_p = %.4f, baseline source: %s)\n', ...
        transducer_idx, baseline, target, scale_p, baseline_src);
    p_scaled = p * scale_p;
end

function [P, f] = load_complex_field(results_file)
% Complex field written next to an acoustic results file (io.save_p_complex).
    [d, n, e] = fileparts(results_file);
    pc_file = fullfile(d, [regexprep(n, '_results', '_p_complex', 'once'), e]);
    if ~isfile(pc_file)
        error('combine_async_intensity:noComplexField', ...
            ['Transducers firing together at different frequencies need each run''s complex field: ' ...
             '%s is missing (acoustic runs with io.save_p_complex).'], pc_file);
    end
    s = load(pc_file, 'p_complex', 'p_complex_info');
    if ~isfield(s, 'p_complex_info') || ~strcmp(s.p_complex_info.phase_reference, 'simulation_t0')
        error('combine_async_intensity:phaseReference', ...
            '%s: phases must refer to t = 0 of the simulation (k-Wave or jwave runs).', pc_file);
    end
    P = double(s.p_complex);
    f = double(s.p_complex_info.freq_hz);
end

function vol = backproject_complex(P_fov, fov_offset, full_dims)
    vol = complex(backproject(real(P_fov), fov_offset, full_dims), ...
                  backproject(imag(P_fov), fov_offset, full_dims));
end
