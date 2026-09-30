function N = neuronav_name_series(D, parameters)
% NEURONAV_NAME_SERIES  Name the detected series of a session (shared GUI/script rule).
%
% Applies the same naming rule in neuronav_gui and NEURONAV_INGEST_BATCH, in
% this order of priority:
%
%   1. coil_map name   - a target given in neuronav.coil_map for that series
%   2. heuristic       - the planned target nearest the series' beam axis
%                        (NEURONAV_SUGGEST_TARGET_NAMES), its label mapped
%                        through neuronav.target_labels (NEURONAV_MAP_LABEL);
%                        only if neuronav.suggest_names is true (default)
%   3. default         - 'series_<k>'
%
% Names are then made unique: repeats get '_rep2', '_rep3', ….
%
% Settings (parameters.neuronav):
%   suggest_names      true   use the heuristic (2)
%   axis_tolerance_mm  10     max. target distance from the beam axis [mm]
%   target_depth_mm    [0 150] accepted target depth along the axis [mm]
%   target_labels      []     label map, see NEURONAV_MAP_LABEL
%
% Use as:
%   N = neuronav_name_series(D, parameters)      % D from NEURONAV_DETECT_SERIES
%
% Output:
%   N - 1×numel(D.results) struct array with fields
%         .name      final target name
%         .rule      'coil_map' / 'heuristic' / 'default'
%         .label     planned-target label as recorded ('' if none)
%         .source    file the label came from ('' if none)
%         .dist_mm   target distance from the beam axis (NaN if n/a)
%         .depth_mm  target depth along the beam axis (NaN if n/a)
%         .tracked   false for empty or untracked (all positions 0) series
%
% See also: NEURONAV_DETECT_SERIES, NEURONAV_SUGGEST_TARGET_NAMES, NEURONAV_MAP_LABEL

    opt = naming_options(parameters);
    n = numel(D.results);
    N = repmat(struct('name', '', 'rule', 'default', 'label', '', 'source', '', ...
        'dist_mm', NaN, 'depth_mm', NaN, 'tracked', false), 1, n);

    if opt.suggest_names && n > 0
        cand = D.candidates;
        if isempty(cand), cand = struct('description', {}, 'xyz', {}, 'source', {}); end
        sg = neuronav_suggest_target_names(D.results, cand, opt.axis_tolerance_mm, ...
            opt.target_depth_mm, opt.target_labels);
    else
        sg = repmat(struct('name', '', 'label', '', 'source', '', 'dist_mm', NaN, 'depth_mm', NaN), 1, n);
    end

    for k = 1:n
        r = D.results{k};
        N(k).tracked = ~isempty(r) && ~all(abs(r.matrix4d_all(:, 1:3, 4)) < 1e-9, 'all');
        N(k).label = sg(k).label; N(k).source = sg(k).source;
        N(k).dist_mm = sg(k).dist_mm; N(k).depth_mm = sg(k).depth_mm;
        if k <= numel(D.typed_names) && ~isempty(D.typed_names{k})
            N(k).name = D.typed_names{k};  N(k).rule = 'coil_map';
        elseif ~isempty(sg(k).name)
            N(k).name = sg(k).name;        N(k).rule = 'heuristic';
        else
            N(k).name = sprintf('series_%d', k);
        end
    end

    % Unique names (suggestions are unique already; coil_map/defaults may clash)
    names = {N.name};
    for u = unique(names)
        idx = find(strcmp(names, u{1}));
        for j = 2:numel(idx)
            N(idx(j)).name = sprintf('%s_rep%d', regexprep(u{1}, '_rep\d+$', ''), j);
        end
    end
end


function opt = naming_options(parameters)
    nv = struct();
    if isfield(parameters, 'neuronav') && isstruct(parameters.neuronav)
        nv = parameters.neuronav;
    end
    opt.suggest_names     = logical(get_field(nv, 'suggest_names', true));
    opt.axis_tolerance_mm = get_field(nv, 'axis_tolerance_mm', 10);
    opt.target_depth_mm   = get_field(nv, 'target_depth_mm', [0 150]);
    opt.target_labels     = get_field(nv, 'target_labels', []);
end

function v = get_field(s, f, default)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = default; end
end
