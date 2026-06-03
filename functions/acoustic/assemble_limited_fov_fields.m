function assemble_limited_fov_fields(files_in, file_out, targets_wcm2)
% ASSEMBLE_LIMITED_FOV_FIELDS  Back-project and combine N limited-FOV acoustic results.
%
% Each input file was produced by a PRESTUS acoustic simulation that ran on
% a cropped sub-grid (via ACOUSTIC_GRID_FOV).  The crop offset is stored in
% acoustic_provenance.fov_offset_ac and the full-grid dimensions in
% acoustic_provenance.full_ac_dims.  This function:
%
%   1. Optionally scales each pressure field to its target ISPPA.
%   2. Zero-pads (back-projects) each FOV pressure field to the shared full
%      acoustic grid.
%   3. Combines the back-projected fields using an incoherent intensity sum
%      and a per-voxel safety maximum, exactly as COMBINE_ASYNC_INTENSITY
%      does for async transducer runs:
%
%        I_total   = Σ  p_i² / (2·ρ·c)
%        p_thermal = sqrt(2·I_total·ρ·c)     → sensor_data.p_max_all
%        p_safety  = max_i(p_i)              → sensor_data.p_max_async
%
%   4. Writes the combined result as a standard PRESTUS acoustic cache file
%      so the existing thermal pipeline needs no modification.
%
% For simulations where all placements fire simultaneously, the incoherent
% intensity sum is the correct thermal field (each transducer contributes
% independently to heat deposition at any voxel).  See COMBINE_ASYNC_INTENSITY
% for the rationale and the async (non-FOV) equivalent.
%
% Use as:
%   assemble_limited_fov_fields(files_in, file_out)
%   assemble_limited_fov_fields(files_in, file_out, targets_wcm2)
%
% Input:
%   files_in     - cell array of paths to limited-FOV acoustic cache .mat
%                  files {N×1}.  Each must contain acoustic_provenance with
%                  fov_offset_ac ([1×3]) and full_ac_dims ([1×3]).
%   file_out     - path where the combined cache file will be written.
%   targets_wcm2 - (optional) [1×N] target ISPPAs [W/cm²] for pressure
%                  scaling.  NaN entries fall back to the cached transducer
%                  config (same logic as COMBINE_ASYNC_INTENSITY).
%
% Notes:
%   - All input files must share the same full_ac_dims (full acoustic grid).
%   - Grid, medium, source, sensor, and segmentation metadata are taken from
%     files_in{1} and written to the combined cache verbatim; they must
%     therefore correspond to the full acoustic grid, not the FOV sub-grid.
%     The caller should supply files whose metadata already reflects the full
%     grid (e.g. by running acoustic_grid_fov on a copy of the full pipeline
%     state and re-attaching the full metadata when saving, or by loading the
%     full-grid cache separately and passing its metadata as reference).
%     In practice, use the 'ref_metadata_file' workflow described below.
%   - acoustic_provenance in the combined file records combined_from,
%     combination_mode, fov_offset_ac per file, per-transducer targets, and
%     baselines.
%
% Typical workflow
% ----------------
%   % 1. Run limited-FOV acoustic sims for each placement:
%   %      prestus_pipeline(params_t1)  % saves sub-001_..._results_t1.mat
%   %      prestus_pipeline(params_t2)  % saves sub-001_..._results_t2.mat
%   %
%   % 2. Assemble:
%   assemble_limited_fov_fields( ...
%       {'sub-001_..._results_t1.mat', 'sub-001_..._results_t2.mat'}, ...
%       'sub-001_..._results_combined.mat', ...
%       [target_t1_wcm2, target_t2_wcm2]);
%   %
%   % 3. Point the thermal pipeline at the combined file:
%   %      parameters.io.acoustic_cache_affix = '_combined';
%   %      parameters.modules.run_acoustic_sims = 0;
%   %      prestus_pipeline(parameters);   % runs thermal only
%
% See also: ACOUSTIC_GRID_FOV, COMBINE_ASYNC_INTENSITY, PRESTUS_PIPELINE

arguments
    files_in
    file_out     (1,:) char
    targets_wcm2 (1,:) double = []
end

if ischar(files_in) || isstring(files_in)
    files_in = cellstr(files_in);
end
if ~iscell(files_in)
    error('assemble_limited_fov_fields:badInput', ...
        'files_in must be a cell array of file paths.');
end
files_in = files_in(:);
N = numel(files_in);
if N < 2
    error('assemble_limited_fov_fields:tooFewFiles', ...
        'At least two limited-FOV acoustic cache files are required.');
end
for k = 1:N
    if ~isfile(files_in{k})
        error('assemble_limited_fov_fields:missingFile', ...
            'Acoustic cache not found: %s', files_in{k});
    end
end

if numel(targets_wcm2) < N
    targets_wcm2(end+1:N) = NaN;
end

% =========================================================================
%% Load first file — reference grid and metadata
% =========================================================================

fprintf('Loading limited-FOV acoustic result 1/%d:\n  %s\n', N, files_in{1});
data_ref = load(files_in{1}, 'sensor_data', 'kgrid', 'kwave_medium', 'source', ...
    'sensor', 'segmentation', 'source_labels', 'medium_masks', ...
    'acoustic_provenance', 'acoustic_info');

prov_ref = data_ref.acoustic_provenance;
[full_ac_dims, fov_offset_ref] = parse_fov_provenance(prov_ref, 1);

p_fov_ref = double(abs(data_ref.sensor_data.p_max_all));

% For the incoherent combination we need ρ and c on the FULL grid.
% These come from whichever reference captures the complete medium.  If the
% first file was run on a FOV sub-grid its kwave_medium is cropped — use the
% density and sound speed from the full combined accumulator (see below).
% We defer the ρ·c read to after all files are assembled.

[p_fov_ref_scaled, target_used, baseline_used] = ...
    scale_pressure(p_fov_ref, prov_ref, data_ref.acoustic_info, targets_wcm2(1), 1);
clear p_fov_ref;

% Back-project first field into full grid
p_full_sq_sum  = backproject_squared(p_fov_ref_scaled, fov_offset_ref, full_ac_dims);
p_full_max     = backproject(p_fov_ref_scaled, fov_offset_ref, full_ac_dims);

fov_offsets_log = {fov_offset_ref};
targets_log     = target_used;
baselines_log   = baseline_used;
clear p_fov_ref_scaled;

% =========================================================================
%% Accumulate remaining files
% =========================================================================

for k = 2:N
    fprintf('Loading limited-FOV acoustic result %d/%d:\n  %s\n', k, files_in{k});
    data_k  = load(files_in{k}, 'sensor_data', 'acoustic_provenance', 'acoustic_info');
    prov_k  = data_k.acoustic_provenance;

    [full_ac_dims_k, fov_offset_k] = parse_fov_provenance(prov_k, k);

    if ~isequal(full_ac_dims_k, full_ac_dims)
        error('assemble_limited_fov_fields:fullDimsMismatch', ...
            ['full_ac_dims mismatch: file 1 has [%s], file %d has [%s]. ' ...
             'All runs must share the same full acoustic grid.'], ...
            num2str(full_ac_dims), k, num2str(full_ac_dims_k));
    end

    p_fov_k = double(abs(data_k.sensor_data.p_max_all));

    [p_fov_k_scaled, target_used, baseline_used] = ...
        scale_pressure(p_fov_k, prov_k, data_k.acoustic_info, targets_wcm2(k), k);
    clear data_k p_fov_k;

    p_full_sq_sum = p_full_sq_sum + backproject_squared(p_fov_k_scaled, fov_offset_k, full_ac_dims);
    p_full_max    = max(p_full_max,  backproject(p_fov_k_scaled, fov_offset_k, full_ac_dims));

    fov_offsets_log{k} = fov_offset_k;           %#ok<AGROW>
    targets_log         = [targets_log,  target_used];   %#ok<AGROW>
    baselines_log       = [baselines_log, baseline_used]; %#ok<AGROW>
    clear p_fov_k_scaled;
end

% =========================================================================
%% Back-convert to effective pressure (intensity-sum → thermal field)
%
% Incoherent intensity:  I = p² / (2·ρ·c)
% Sum of intensities:    I_total = Σ p_i² / (2·ρ·c)
%                                = p_sq_sum / (2·ρ·c)
% Effective pressure:    p_thermal = sqrt(2·I_total·ρ·c)
%                                  = sqrt(p_sq_sum)       [ρ,c cancel]
%
% =========================================================================
p_thermal = sqrt(p_full_sq_sum);
clear p_full_sq_sum;

% =========================================================================
%% Assemble combined sensor_data and provenance
% =========================================================================

sensor_data             = data_ref.sensor_data;
sensor_data.p_max_all   = single(p_thermal);   % thermal heat source
sensor_data.p_max_async = single(p_full_max);  % instantaneous safety maximum
clear p_thermal p_full_max;

acoustic_provenance                          = data_ref.acoustic_provenance;
acoustic_provenance.combined_from            = files_in;
acoustic_provenance.combination_mode         = 'limited_fov_incoherent_sum';
acoustic_provenance.combination_fov_offsets  = fov_offsets_log;
acoustic_provenance.combination_full_ac_dims = full_ac_dims;
acoustic_provenance.combination_targets      = targets_log;
acoustic_provenance.combination_baselines    = baselines_log;
acoustic_provenance.combination_time         = datestr(now); %#ok<TNOW1,DATST>
% Clear per-FOV crop metadata — the combined file represents the full grid
acoustic_provenance.fov_offset_ac  = [];
acoustic_provenance.full_ac_dims   = full_ac_dims;

% =========================================================================
%% Save
% =========================================================================

fprintf('Saving combined acoustic cache:\n  %s\n', file_out);

kgrid         = data_ref.kgrid;         %#ok<NASGU>
kwave_medium  = data_ref.kwave_medium;  %#ok<NASGU>
source        = data_ref.source;        %#ok<NASGU>
sensor        = data_ref.sensor;        %#ok<NASGU>
segmentation  = data_ref.segmentation;  %#ok<NASGU>
source_labels = data_ref.source_labels; %#ok<NASGU>
medium_masks  = data_ref.medium_masks;  %#ok<NASGU>
acoustic_info = data_ref.acoustic_info; %#ok<NASGU>

save(file_out, ...
    'sensor_data', 'kgrid', 'kwave_medium', 'source', 'sensor', ...
    'segmentation', 'source_labels', 'medium_masks', 'acoustic_info', ...
    'acoustic_provenance', '-v7.3');

fprintf('Done — combined cache written from %d limited-FOV fields.\n', N);
end

% =========================================================================
%% Local helpers
% =========================================================================

function [full_dims, fov_offset] = parse_fov_provenance(prov, file_idx)
    if ~isfield(prov, 'fov_offset_ac') || isempty(prov.fov_offset_ac)
        error('assemble_limited_fov_fields:missingFovOffset', ...
            ['File %d acoustic_provenance.fov_offset_ac is missing or empty. ' ...
             'Only files produced by an acoustic_grid_fov run can be assembled. ' ...
             'Run the acoustic simulation with grid.acoustic_fov_diameter_mm set.'], ...
            file_idx);
    end
    if ~isfield(prov, 'full_ac_dims') || isempty(prov.full_ac_dims)
        error('assemble_limited_fov_fields:missingFullDims', ...
            'File %d acoustic_provenance.full_ac_dims is missing.', file_idx);
    end
    fov_offset = prov.fov_offset_ac(:)';   % [1×3]
    full_dims  = prov.full_ac_dims(:)';    % [1×3]
end

function vol_full = backproject(p_fov, fov_offset, full_dims)
% Insert p_fov into a zero-padded full-grid volume.
    vol_full = zeros(full_dims, 'like', p_fov);
    fov_dims = size(p_fov);
    fov_end  = fov_offset + fov_dims - 1;
    vol_full(fov_offset(1):fov_end(1), ...
             fov_offset(2):fov_end(2), ...
             fov_offset(3):fov_end(3)) = p_fov;
end

function vol_sq_full = backproject_squared(p_fov, fov_offset, full_dims)
% Insert p_fov.^2 into a zero-padded full-grid volume.
    vol_sq_full = zeros(full_dims, 'like', p_fov);
    fov_dims = size(p_fov);
    fov_end  = fov_offset + fov_dims - 1;
    vol_sq_full(fov_offset(1):fov_end(1), ...
                fov_offset(2):fov_end(2), ...
                fov_offset(3):fov_end(3)) = p_fov .^ 2;
end

function [p_scaled, target_used, baseline_used] = scale_pressure( ...
        p, acoustic_provenance, acoustic_info, target_arg, transducer_idx)
% Mirrors the scale_pressure helper in COMBINE_ASYNC_INTENSITY.
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

    baseline = NaN;
    if isfield(acoustic_provenance, 'freefield_isppa_wcm2') && ...
            isscalar(acoustic_provenance.freefield_isppa_wcm2) && ...
            acoustic_provenance.freefield_isppa_wcm2 > 0
        baseline = acoustic_provenance.freefield_isppa_wcm2;
    end

    target_used   = target;
    baseline_used = baseline;

    if ~isfinite(target) || ~isfinite(baseline)
        if ~isfinite(target) && ~isfinite(baseline)
            fprintf('  Placement %d: no target or baseline — using simulated drive level.\n', ...
                transducer_idx);
        elseif ~isfinite(baseline)
            warning('assemble_limited_fov_fields:noBaseline', ...
                ['Placement %d: freefield_isppa_wcm2 missing from provenance — ' ...
                 'pressure scaling skipped. Re-run with modules.run_water_baseline = 1.'], ...
                transducer_idx);
        else
            fprintf('  Placement %d: no target specified — using simulated drive level.\n', ...
                transducer_idx);
        end
        p_scaled = p;
        return;
    end

    scale_p = sqrt(target / baseline);
    if abs(scale_p - 1) <= 0.01
        p_scaled = p;
        return;
    end
    if scale_p > 4 || scale_p < 0.25
        warning('assemble_limited_fov_fields:largeScale', ...
            ['Placement %d: large pressure scale factor (%.2fx) from %.1f to %.1f W/cm². ' ...
             'Verify linear acoustic regime.'], transducer_idx, scale_p, baseline, target);
    end
    fprintf('  Placement %d: %.2f → %.2f W/cm² (scale_p = %.4f)\n', ...
        transducer_idx, baseline, target, scale_p);
    p_scaled = p * scale_p;
end
