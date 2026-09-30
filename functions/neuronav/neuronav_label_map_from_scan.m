function target_labels = neuronav_label_map_from_scan(labels, mode)
% NEURONAV_LABEL_MAP_FROM_SCAN  Build a neuronav.target_labels map from a label scan.
%
% Turns the label table of NEURONAV_SCAN_SESSIONS into a label map that can
% be edited and stored as parameters.neuronav.target_labels (e.g. saved to
% the study YAML from neuronav_gui).
%
% Use as:
%   target_labels = neuronav_label_map_from_scan(labels)            % by .target
%   target_labels = neuronav_label_map_from_scan(labels, 'stem')    % group by .stem
%
% Input:
%   labels - label table from NEURONAV_SCAN_SESSIONS (columns label, stem,
%            target); edit .target first to set names by hand
%   mode   - 'target' (default): one entry per distinct .target, listing its
%            labels; 'stem': one entry per .stem
%
% Output:
%   target_labels - struct array with .target and .labels (cell of labels)
%
% See also: NEURONAV_SCAN_SESSIONS, NEURONAV_MAP_LABEL

    if nargin < 2 || isempty(mode), mode = 'target'; end
    key = labels.(mode);
    target_labels = struct('target', {}, 'labels', {});
    [u, ~, idx] = unique(key, 'stable');
    for k = 1:numel(u)
        target_labels(k).target = u{k};
        target_labels(k).labels = labels.label(idx == k)';
    end
end
