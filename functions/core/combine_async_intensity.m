function combine_async_intensity(files_in, file_out, targets_wcm2)
% COMBINE_ASYNC_INTENSITY  Scale and combine N independent acoustic results
%                          for async (non-simultaneous) transducer firing.
%
% For asynchronously firing transducers the pressure fields are independent.
% Two combined fields are produced from the per-transducer scaled pressures:
%
%   p_i_scaled = p_i · sqrt(target_i / baseline_i)
%
%   Thermal field  — incoherent intensity sum, drives the bioheat equation:
%     I_total  = Σ  p_i_scaled² / (2·ρ·c)
%     p_thermal = sqrt(2·I_total·ρ·c)          → sensor_data.p_max_all
%
%   Acoustic safety field  — per-voxel pressure maximum, reflects the
%     instantaneous peak when only one transducer is active at a time:
%     p_safety = max_i( p_i_scaled )            → sensor_data.p_max_async
%
% Rationale: in async mode the instantaneous peak intensity at any voxel
% equals the maximum over individual transducer fields, not their sum.
% Summing is correct only for the time-averaged heat deposition used in
% thermal modelling (assuming equal, non-overlapping duty cycles that
% together fill the total protocol duty cycle).  Acoustic safety metrics
% (ISPPA, MI, peak pressure) must therefore be derived from p_max_async,
% not from p_max_all.
%
% The target for each transducer is resolved in this order:
%   1. targets_wcm2(i) argument (if supplied and finite)
%   2. acoustic_info.parameters.transducer(1).target_isppa_wcm2 in the cache
%   3. No scaling (use the simulated drive level as-is)
%
% The baseline (free-water ISPPA at the simulated drive level) is read from
% acoustic_provenance.freefield_isppa_wcm2 in each cache file.  If absent,
% scaling is skipped for that transducer with a warning.
%
% Grid/medium metadata for the combined cache are taken from files_in{1};
% all files must share the same simulation grid.
%
% Use as:
%   combine_async_intensity(files_in, file_out)
%   combine_async_intensity(files_in, file_out, targets_wcm2)
%
% Input:
%   files_in     - cell array of paths to acoustic cache .mat files {N×1}
%   file_out     - path where the combined cache will be written
%   targets_wcm2 - (optional) 1×N vector of target ISPPAs [W/cm²], one per
%                  transducer.  NaN or [] entries fall back to the cached
%                  transducer config.
%
% Notes:
%   - Axisymmetric (2-D) results are not supported.
%   - acoustic_provenance in the combined file records combined_from,
%     combination_mode, per-transducer targets, and baselines.
%
% See also: ASYNC_TRANSDUCER_PIPELINE, APPLY_ISPPA_SCALING, ACOUSTIC_WRAPPER

arguments
    files_in
    file_out     (1,:) char
    targets_wcm2 (1,:) double = []
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

% Pad or trim the external targets vector to length N.
if numel(targets_wcm2) < N
    targets_wcm2(end+1:N) = NaN;
end

% =========================================================================
%% Load first file — provides grid, medium, and all metadata
% =========================================================================

fprintf('Loading acoustic result 1/%d:\n  %s\n', N, files_in{1});
data_ref = load(files_in{1}, 'sensor_data', 'kgrid', 'kwave_medium', 'source', ...
                'sensor', 'segmentation', 'source_labels', 'medium_masks', ...
                'acoustic_provenance', 'acoustic_info');

p_ref = double(abs(data_ref.sensor_data.p_max_all));

if ndims(p_ref) < 3 %#ok<ISMAT>
    error('combine_async_intensity:axisymmetricNotSupported', ...
        ['Axisymmetric (2-D) results cannot be combined in async mode. ' ...
         'Re-run acoustic simulations with a 3-D grid.']);
end

rho = data_ref.kwave_medium.density;
c   = data_ref.kwave_medium.sound_speed;

[p_ref_scaled, target_used, baseline_used] = scale_pressure( ...
    p_ref, data_ref.acoustic_provenance, data_ref.acoustic_info, targets_wcm2(1), 1);
clear p_ref;

I_total        = p_ref_scaled.^2 ./ (2 .* rho .* c);
p_max_safety   = p_ref_scaled;   % per-voxel pressure maximum across transducers
targets_log    = target_used;
baselines_log  = baseline_used;
clear p_ref_scaled;

% =========================================================================
%% Accumulate remaining transducers
% =========================================================================

for k = 2:N
    fprintf('Loading acoustic result %d/%d:\n  %s\n', k, N, files_in{k});
    data_k = load(files_in{k}, 'sensor_data', 'acoustic_provenance', 'acoustic_info');
    p_k    = double(abs(data_k.sensor_data.p_max_all));

    if ~isequal(size(p_k), size(I_total))
        error('combine_async_intensity:gridMismatch', ...
            ['p_max_all size mismatch between file 1 [%s] and file %d [%s]. ' ...
             'All runs must use the same simulation grid.'], ...
            num2str(size(I_total)), k, num2str(size(p_k)));
    end

    [p_k_scaled, target_used, baseline_used] = scale_pressure( ...
        p_k, data_k.acoustic_provenance, data_k.acoustic_info, targets_wcm2(k), k);
    clear data_k p_k;

    I_total      = I_total + p_k_scaled.^2 ./ (2 .* rho .* c);
    p_max_safety = max(p_max_safety, p_k_scaled);  % instantaneous peak
    targets_log   = [targets_log,  target_used];   %#ok<AGROW>
    baselines_log = [baselines_log, baseline_used]; %#ok<AGROW>
    clear p_k_scaled;
end

% =========================================================================
%% Back-convert to effective pressures
% =========================================================================

% Thermal field: intensity sum → drives bioheat equation
p_thermal = sqrt(2 .* I_total .* rho .* c);
clear I_total rho c;

% =========================================================================
%% Assemble combined sensor_data and provenance
% =========================================================================

sensor_data              = data_ref.sensor_data;
sensor_data.p_max_all    = p_thermal;   % thermal heat source (intensity sum)
sensor_data.p_max_async  = p_max_safety; % acoustic safety (instantaneous max)
clear p_thermal p_max_safety;

acoustic_provenance                      = data_ref.acoustic_provenance;
acoustic_provenance.combined_from        = files_in;
acoustic_provenance.combination_mode     = 'async_dual_field';  % p_max_all=thermal sum; p_max_async=safety max
acoustic_provenance.combination_targets  = targets_log;
acoustic_provenance.combination_baselines = baselines_log;
acoustic_provenance.combination_time     = datestr(now); %#ok<TNOW1,DATST>

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

fprintf('Done — combined cache written from %d transducers.\n', N);
end

% =========================================================================
%% Local helper
% =========================================================================

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

    % Resolve baseline
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
            fprintf('  Transducer %d: no target or baseline — using simulated drive level.\n', ...
                transducer_idx);
        elseif ~isfinite(baseline)
            warn('combine_async_intensity:noBaseline', ...
                ['Transducer %d: freefield_isppa_wcm2 missing from provenance — ' ...
                 'pressure scaling skipped. Re-run the acoustic stage with ' ...
                 'modules.run_water_baseline = 1.'], transducer_idx);
        else
            fprintf('  Transducer %d: no target specified — using simulated drive level.\n', ...
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
        warn('combine_async_intensity:largeScale', ...
            ['Transducer %d: large pressure scale factor (%.2fx) from %.1f to %.1f W/cm². ' ...
             'Verify linear acoustic regime.'], transducer_idx, scale_p, baseline, target);
    end
    fprintf('  Transducer %d: %.2f → %.2f W/cm² (scale_p = %.4f)\n', ...
        transducer_idx, baseline, target, scale_p);
    p_scaled = p * scale_p;
end
