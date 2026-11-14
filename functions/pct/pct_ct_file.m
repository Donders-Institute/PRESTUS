function ct_file = pct_ct_file(parameters)
% PCT_CT_FILE  This subject's measured CT, or '' when the pseudo-CT is used
%
% The CT-based skull mapping reads one CT image per subject. A measured CT
% is given by path.ct_pattern (relative to path.anat, like path.t1_pattern);
% when that is empty, the pipeline uses the pseudo-CT generated from the UTE
% image (PCT_CREATE_PCT_RUN) instead.
%
% Use as:
%   ct_file = pct_ct_file(parameters)
%
% Input:
%   parameters - PRESTUS parameter struct (subject_id, path.anat, path.ct_pattern)
%
% Output:
%   ct_file - full path of the measured CT; '' when path.ct_pattern is empty.
%             Errors when the pattern is set but matches no file.
%
% See also: PCT_LOAD_CT, FIND_PSEUDOCT_DIR

arguments
    parameters (1,1) struct
end

    ct_file = '';
    if ~(isfield(parameters, 'path') && isfield(parameters.path, 'ct_pattern') && ...
            (ischar(parameters.path.ct_pattern) || isstring(parameters.path.ct_pattern)) && ...
            strlength(string(parameters.path.ct_pattern)) > 0)
        return
    end
    % sprintf specifiers (e.g. %03d) and shell-style globs, as for path.t1_pattern
    ct_glob = fullfile(char(parameters.path.anat), ...
        sprintf(char(parameters.path.ct_pattern), parameters.subject_id));
    d = dir(ct_glob);
    if isempty(d)
        error('PRESTUS:pct:ctNotFound', 'No CT file matched path.ct_pattern: %s', ct_glob);
    end
    if numel(d) > 1
        warning('PRESTUS:pct:ctAmbiguous', '%d files match path.ct_pattern (%s); using %s.', ...
            numel(d), ct_glob, d(1).name);
    end
    ct_file = fullfile(d(1).folder, d(1).name);
end
