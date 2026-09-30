function suggestions = neuronav_suggest_target_names(results, candidates, tolerance_mm, depth_mm, target_labels)
% NEURONAV_SUGGEST_TARGET_NAMES  Heuristic target names for ingested series.
%
% Names each series (as returned by NEURONAV_INGEST_MARKERS or
% NEURONAV_COMPUTE_SERIES_STATISTICS) after the planned target it was aimed
% at, instead of a generic 'series_N':
%
%   * Pose series (TriggerMarkers, InstrumentMarker, and InstrumentMarker
%     elements inside GUMMarkers files): the beam axis runs from the mean
%     position (matrix column 4) along the mean beam direction (column 1,
%     coil centre → head, as in LOCALITE_MATRIX_TO_POSITIONS). The candidate
%     target with the smallest perpendicular distance to that axis, with its
%     depth inside depth_mm (default: anywhere in front of the transducer),
%     is chosen if the distance is within tolerance_mm.
%   * GUMMarkers EntryTargetPair series: the series *is* a named target
%     point, so the pair's own description is used.
%
% A GUMMarkers file that holds both a planned pair and a recorded pose for
% the same label gives the pose the bare name and the pair '<label>_target'.
% Several poses matching the same target (repeated blocks) get '_rep2',
% '_rep3', … appended to all but the first so names stay unique.
%
% Use as:
%   s = neuronav_suggest_target_names(results, candidates)
%   s = neuronav_suggest_target_names(results, candidates, tolerance_mm)
%   s = neuronav_suggest_target_names(results, candidates, tolerance_mm, depth_mm)
%   s = neuronav_suggest_target_names(results, candidates, tolerance_mm, depth_mm, target_labels)
%
% Input:
%   results      - 1×N cell of series stats (empty cells allowed)
%   candidates   - struct array from NEURONAV_READ_LABEL_CANDIDATES
%                  (.description, .xyz [1×3 mm], .source)
%   tolerance_mm - max. perpendicular distance target ↔ beam axis [mm]
%                  (default 10)
%   depth_mm     - [min max] accepted target depth along the beam axis [mm]
%                  (default [0 Inf]); e.g. a window around the focal distance
%                  rejects targets the axis only passes on the far side
%   target_labels - (optional) label map applied to the matched label
%                  (NEURONAV_MAP_LABEL); [] = use labels as found
%
% Output:
%   suggestions  - 1×N struct array with fields
%                    .name     suggested name, '' if no match
%                    .label    label as recorded (before target_labels)
%                    .source   where it came from ('' if no match)
%                    .dist_mm  target distance from the beam axis [mm] (NaN if n/a)
%                    .depth_mm target depth along the beam axis [mm] (NaN if n/a)
%
% Note: both the series positions and the candidates are in the Localite
% planning-image space of the same session, so no transform is needed.
%
% See also: NEURONAV_READ_LABEL_CANDIDATES, NEURONAV_INGEST_MARKERS

    if nargin < 3 || isempty(tolerance_mm), tolerance_mm = 10; end
    if nargin < 4 || isempty(depth_mm), depth_mm = [0 Inf]; end
    if nargin < 5, target_labels = []; end
    is_pair = false(1, numel(results));

    n = numel(results);
    suggestions = repmat(struct('name', '', 'label', '', 'source', '', ...
        'dist_mm', NaN, 'depth_mm', NaN), 1, n);

    for s = 1:n
        r = results{s};
        if isempty(r), continue; end

        % GUMMarkers EntryTargetPair: the point carries its own label.
        mk = r.markers(1);
        if isfield(mk, 'element_type') && strcmp(mk.element_type, 'EntryTargetPair')
            % An unset placeholder pair sits at the origin; its label means nothing.
            if isfield(mk, 'description') && strlength(string(mk.description)) > 0 ...
                    && any(abs(r.position_mean) > 1e-9)
                suggestions(s).label  = strtrim(char(string(mk.description)));
                suggestions(s).name   = neuronav_map_label(suggestions(s).label, target_labels);
                suggestions(s).source = 'GUMMarkers pair label';
                is_pair(s) = true;
            end
            continue;
        end

        if isempty(candidates), continue; end
        M = reshape(r.matrix4d_mean, 4, 4);
        p = M(1:3, 4)';
        d = M(1:3, 1)';
        if all(abs(p) < 1e-9) || norm(d) < 1e-9, continue; end   % untracked
        d = d / norm(d);

        xyz   = vertcat(candidates.xyz);
        v     = xyz - p;
        depth = v * d';
        perp  = sqrt(max(sum(v.^2, 2) - depth.^2, 0));
        perp(depth <= depth_mm(1) | depth > depth_mm(2)) = Inf;

        [dmin, k] = min(perp);
        if dmin <= tolerance_mm
            suggestions(s).label    = candidates(k).description;
            suggestions(s).name     = neuronav_map_label(candidates(k).description, target_labels);
            suggestions(s).source   = candidates(k).source;
            suggestions(s).dist_mm  = dmin;
            suggestions(s).depth_mm = depth(k);
        elseif isfinite(dmin)
            % Report the nearest miss for QC, without naming.
            suggestions(s).dist_mm  = dmin;
            suggestions(s).depth_mm = depth(k);
        end
    end

    % Planned pair vs. recorded pose with the same label: the pose keeps it.
    names = {suggestions.name};
    for s = find(is_pair)
        if any(strcmp(names(~is_pair), names{s}))
            suggestions(s).name = [names{s} '_target'];
        end
    end

    % Unique names for repeated blocks on the same target.
    names = {suggestions.name};
    % Index loop: a for-loop over an empty 0×1 cell would still run once.
    u = unique(names(~cellfun(@isempty, names)));
    for i = 1:numel(u)
        idx = find(strcmp(names, u{i}));
        for j = 2:numel(idx)
            suggestions(idx(j)).name = sprintf('%s_rep%d', u{i}, j);
        end
    end
end
