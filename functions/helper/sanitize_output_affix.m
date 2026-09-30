function affix = sanitize_output_affix(label)
% SANITIZE_OUTPUT_AFFIX  Normalise io.output_affix to the BIDS desc entity
%
% Replaces characters other than letters, digits, '_' and '-' with '_' and
% wraps a non-empty label as '_desc-<label>' (leading underscores stripped).
% An affix that already starts with '_desc-' is returned unchanged. Used by
% LOAD_PARAMETERS and by the study tools, so both derive the same output
% file names from the same label.
%
% Use as:
%   affix = sanitize_output_affix(label)
%
% Input:
%   label - char/string label; [] or '' yield ''
%
% Output:
%   affix - '' or '_desc-<label>'
%
% See also: LOAD_PARAMETERS, STUDY_RUN_TABLE

    if isempty(label) || ~(ischar(label) || isstring(label))
        affix = '';
        return
    end
    affix = regexprep(char(label), '[^a-zA-Z0-9_-]', '_');
    if ~isempty(affix) && ~startsWith(affix, '_desc-')
        affix = ['_desc-' regexprep(affix, '^_+', '')];
    end
end
