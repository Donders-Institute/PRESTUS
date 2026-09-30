function parameters = study_row_parameters(study_config, run, extra_overrides)
% STUDY_ROW_PARAMETERS  Load the full config of one study run (one manifest row)
%
% Merge order, later wins:
%   config_default.yaml  <  study config  <  row config (manifest 'config'
%   column)  <  manifest key columns  <  extra_overrides  <  subject_id and
%   io.output_affix of the run
%
% All of it goes through LOAD_PARAMETERS in one call, so derived values
% (transducer geometry, thermal timing, ...) are computed from the row's
% own settings. The run's condition and manifest row are recorded in
% parameters.study.condition / parameters.study.row.
%
% Use as:
%   parameters = study_row_parameters(study_config, run)
%   parameters = study_row_parameters(study_config, run, extra_overrides)
%
% Input:
%   study_config    - path to the study YAML (the file that sets study.manifest)
%   run             - one row of the table from STUDY_RUN_TABLE
%   extra_overrides - (optional) N×2 cell {key, value} applied after the
%                     manifest columns (e.g. set by PRESTUS_STUDY_START)
%
% Output:
%   parameters - loaded PRESTUS config for this run
%
% See also: STUDY_RUN_TABLE, PRESTUS_STUDY_START, LOAD_PARAMETERS

arguments
    study_config {mustBeTextScalar}
    run table
    extra_overrides cell = cell(0, 2)
end

    assert(height(run) == 1, 'study_row_parameters:oneRow', 'Pass exactly one run (table row).');

    [sdir, sname, sext] = fileparts(char(study_config));
    if isempty(sdir); sdir = pwd; end
    args = {[sname sext], sdir};
    if strlength(run.config) > 0
        [cdir, cname, cext] = fileparts(char(run.config));
        args = [args, {[cname cext], cdir}];
    end

    ov = [run.overrides{1}; extra_overrides; ...
          {'subject_id', run.subject_id; 'io.output_affix', char(run.affix)}];
    parameters = load_parameters(args{:}, ov);

    parameters.study.condition = char(run.condition);
    parameters.study.row       = run.row;
end
