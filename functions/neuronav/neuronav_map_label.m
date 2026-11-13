function [target, matched] = neuronav_map_label(label, target_labels)
% NEURONAV_MAP_LABEL  Map a Localite target label to the study's target name.
%
% Localite target labels (EntryTarget / GUMMarkers descriptions) follow
% whatever convention the operator used; target_labels maps them to the
% names a study uses (e.g. 'LEFT V1', 'L_V' -> 'L_V1'). Configured as
% parameters.neuronav.target_labels, shared by neuronav_gui and scripts:
%
%   neuronav:
%     target_labels:
%       - target: 'L_V1'
%         labels: ['LEFT V1', 'L_V*']     # wildcards * and ?, case-insensitive
%       - target: 'R_V1'
%         labels: ['RIGHT V1', 'R_V*']
%
% Entries are tried in order; the first matching pattern wins. A label that
% matches no entry is returned unchanged (so an empty map uses the labels as
% found in the data; see NEURONAV_SCAN_SESSIONS to list them).
%
% Use as:
%   [target, matched] = neuronav_map_label(label, target_labels)
%
% Input:
%   label         - char/string label as recorded
%   target_labels - struct array (.target, .labels) or cell of such structs
%                   (as loaded from YAML); [] = no mapping
%
% Output:
%   target        - mapped name (char), or label itself if unmatched
%   matched       - true if an entry matched
%
% See also: NEURONAV_NAME_SERIES, NEURONAV_SCAN_SESSIONS

    target  = strtrim(char(string(label)));
    matched = false;
    if isempty(target_labels) || isempty(target), return; end
    if isstruct(target_labels), target_labels = num2cell(target_labels); end

    for i = 1:numel(target_labels)
        e = target_labels{i};
        pats = e.labels;
        if ischar(pats) || isstring(pats), pats = cellstr(pats); end
        for j = 1:numel(pats)
            rx = ['^' regexptranslate('wildcard', strtrim(char(pats{j}))) '$'];
            if ~isempty(regexpi(target, rx, 'once'))
                target  = char(e.target);
                matched = true;
                return;
            end
        end
    end
end
