function [s, dist] = skull_property_stats(parameters, medium_masks, sound_speed, density, alpha_coeff, alpha_power)
% SKULL_PROPERTY_STATS  Modelled skull acoustic properties, in all skull voxels and on the beam path
%
% Summarises the sound speed, density and attenuation at the source frequency
% of the skull as the simulation models it: every medium label starting with
% 'skull' that is a layer of the model (a layered skull, cortical and
% trabecular layers, or a pseudo-CT/CT-mapped skull). Two summaries:
%   - all skull voxels: 10th, 50th and 90th percentile of each property;
%   - the skull on each transducer's beam axis, from the transducer to its
%     focus (SAMPLE_BEAM_AXIS, half-voxel steps): the path length through
%     skull per transducer, and the median of each property over all path
%     samples in skull.
% RUN_STATS records s; telemetry receives its medians only. dist holds the full
% distributions (DISTRIBUTION_STATS: percentiles and histogram) that MEDIUM_SETUP
% writes to the medium_stats JSON for the Nexus.
%
% Use as:
%   [s, dist] = skull_property_stats(parameters, medium_masks, sound_speed, density, alpha_coeff, alpha_power)
%
% Input:
%   parameters   - PRESTUS config with medium_properties, layers,
%                  grid.resolution_mm [mm], transducer(1).freq_hz [Hz], and
%                  transducer(ti).trans_pos / focus_pos [grid voxels] (the beam
%                  path is skipped without them)
%   medium_masks - medium label map (values = index into medium_properties)
%   sound_speed  - sound speed map [m/s], same size as medium_masks
%   density      - density map [kg/m^3]
%   alpha_coeff  - attenuation coefficient map [dB/(cm MHz^y)], before any
%                  power-law refit
%   alpha_power  - power-law exponent map y (or a scalar)
%                  Attenuation at the source frequency [dB/cm] = alpha_coeff *
%                  f_MHz^alpha_power (first transducer's frequency)
%
% Output:
%   s - struct; empty without skull voxels. Fields:
%       skull_c_p10, skull_c_median, skull_c_p90         [m/s]
%       skull_rho_p10, skull_rho_median, skull_rho_p90   [kg/m^3]
%       skull_atten_p10, skull_atten_median, skull_atten_p90 [dB/cm]
%       skull_path_mm               [mm], one value per transducer with positions
%       skull_path_c_median, skull_path_rho_median, skull_path_atten_median
%                                   (only when the beam axes cross skull)
%   dist - struct; empty without skull voxels. Fields: layers (skull layer
%          names), n_voxels, properties.(density|sound_speed|alpha_coeff|
%          attenuation_db_cm) over all skull voxels, and beam_path.path_mm (per
%          transducer), beam_path.n_samples, beam_path.properties.(same names)
%          over the skull samples on the beam paths
%
% See also: MEDIUM_SETUP, SAMPLE_BEAM_AXIS, RUN_STATS

s = struct();
dist = struct();
medium_labels = fieldnames(parameters.medium_properties);
layer_labels = fieldnames(parameters.layers);
ids = find(startsWith(medium_labels, 'skull') & ismember(medium_labels, layer_labels));
if isempty(ids)
    return
end
mask = ismember(medium_masks, ids);
if ~any(mask(:))
    return
end

f_mhz = parameters.transducer(1).freq_hz / 1e6;
atten_db_cm = alpha_coeff .* f_mhz .^ alpha_power;
maps = {'c', sound_speed; 'rho', density; 'atten', atten_db_cm; 'alpha', alpha_coeff};
dist_names = {'sound_speed', 'density', 'attenuation_db_cm', 'alpha_coeff'};
dist.layers = medium_labels(ids)';
dist.n_voxels = nnz(mask);
for k = 1:size(maps, 1)
    dist.properties.(dist_names{k}) = distribution_stats(maps{k, 2}(mask));
end
for k = 1:3                                 % c, rho, atten (alpha only in dist)
    v = double(maps{k, 2}(mask));
    v = v(isfinite(v));
    if isempty(v)
        continue
    end
    q = quantiles(v, [0.1 0.5 0.9]);
    s.(['skull_' maps{k, 1} '_p10'])    = q(1);
    s.(['skull_' maps{k, 1} '_median']) = q(2);
    s.(['skull_' maps{k, 1} '_p90'])    = q(3);
end

% Beam path: transducer to focus, per transducer
res = parameters.grid.resolution_mm;
path_mm = [];
on_path = cell(size(maps, 1), 1);
for ti = 1:numel(parameters.transducer)
    tr = parameters.transducer(ti);
    if ~isfield(tr, 'trans_pos') || ~isfield(tr, 'focus_pos') || ...
            numel(tr.trans_pos) ~= ndims(mask) || numel(tr.focus_pos) ~= ndims(mask)
        continue
    end
    t = double(tr.trans_pos(:)');
    f = double(tr.focus_pos(:)');
    [dist_mm, in_skull] = sample_beam_axis(double(mask), t, f, res);
    if numel(dist_mm) < 2
        continue
    end
    keep = dist_mm <= norm(f - t) * res & in_skull > 0.5;
    path_mm(end + 1) = nnz(keep) * (dist_mm(2) - dist_mm(1)); %#ok<AGROW>
    for k = 1:size(maps, 1)
        [~, v] = sample_beam_axis(double(maps{k, 2}), t, f, res);
        v = v(keep);
        on_path{k} = [on_path{k}, v(isfinite(v))];
    end
end
if ~isempty(path_mm)
    s.skull_path_mm = path_mm;
    dist.beam_path.path_mm = path_mm;
    dist.beam_path.n_samples = numel(on_path{1});
    for k = 1:size(maps, 1)
        if k <= 3 && ~isempty(on_path{k})
            s.(['skull_path_' maps{k, 1} '_median']) = median(on_path{k});
        end
        dist.beam_path.properties.(dist_names{k}) = distribution_stats(on_path{k});
    end
end
end

function q = quantiles(v, p)
% Linear-interpolation quantiles (numpy's default definition; no Statistics Toolbox needed).
v = sort(v(:));
n = numel(v);
if n == 1
    q = repmat(v, size(p));
    return
end
pos = 1 + p * (n - 1);
lo = floor(pos);
hi = min(lo + 1, n);
q = v(lo)' + (pos - lo) .* (v(hi)' - v(lo)');
end
