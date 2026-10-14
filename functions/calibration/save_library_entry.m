function lib = save_library_entry(lib, out_file, mode, report, extra_meta)
% SAVE_LIBRARY_ENTRY  Write a user-calibrated library entry with provenance
%
% Stamps the provenance of a CALIBRATE_LIBRARY_ENTRY result into lib.meta
% and writes the entry to out_file (YAML), creating the folder if needed.
%
% Use as:
%   lib = save_library_entry(lib, out_file, mode, report)
%   lib = save_library_entry(lib, out_file, mode, report, extra_meta)
%
% Input:
%   lib        - library struct
%   out_file   - output YAML path ({library_path}/{combo_name}.yaml)
%   mode       - calibration mode ('per_depth' | 'geometric')
%   report     - report table from CALIBRATE_LIBRARY_ENTRY
%   extra_meta - optional struct of further meta fields (e.g. preset, source)
%
% Output:
%   lib - the entry as written
%
% See also: CALIBRATE_LIBRARY_ENTRY, REFINE_TRANSDUCER_LIBRARY, CALIBRATION_TRANSDUCER

    if nargin < 5; extra_meta = struct(); end
    lib.meta.calibration_type = 'user_empirical';
    lib.meta.last_updated     = datestr(now, 'yyyy-mm-dd'); %#ok<TNOW1,DATST>
    try; [~, lib.meta.prestus_version] = prestus_version(); catch; lib.meta.prestus_version = 'unknown'; end
    lib.meta.refinement_mode  = mode;
    lib.meta.refinement_report = struct('depths_ep_mm', report.depth_ep_mm', ...
        'r2_library', report.r2_library', 'r2_refined', report.r2_refined');
    fn = fieldnames(extra_meta);
    for k = 1:numel(fn); lib.meta.(fn{k}) = extra_meta.(fn{k}); end
    out_dir = fileparts(out_file);
    if ~isempty(out_dir) && ~isfolder(out_dir); mkdir(out_dir); end
    yaml.dumpFile(out_file, lib);
end
