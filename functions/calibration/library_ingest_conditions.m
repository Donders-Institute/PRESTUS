function [lib, report] = library_ingest_conditions(lib_file, csv_path, equip_param, options)
% LIBRARY_INGEST_CONDITIONS  Add dephased conditions from a profile CSV + sidecar [EXPERIMENTAL]
%
% EXPERIMENTAL — NOT VALIDATED. Reads a dephased axial-profile CSV
% (tran-{TRAN}_ds-{DS}_acq-{LABEL}_axial.csv) and its YAML sidecar, whose
% 'columns' list describes each CSV column in order:
%
%   columns:
%     - label: geom_180              # condition name used in transducer.condition
%       phase_shift_deg: [0, 180, 0, 180]   # or elem_phase_deg: absolute phases
%
% The CSV depth row holds each column's nominal focal setpoint [mm from EP].
% Columns are passed to LIBRARY_ADD_CONDITIONS, which predicts each profile
% from the entry's correction at that setpoint and reports the agreement;
% the conditions are written to lib_file. Run without a measured CSV by
% calling LIBRARY_ADD_CONDITIONS directly.
%
% Use as:
%   [lib, report] = library_ingest_conditions(lib_file, csv_path)
%   [lib, report] = library_ingest_conditions(lib_file, csv_path, equip_param, options)
%
% Input:
%   lib_file    - library YAML to update ({tran}_{ds}.yaml)
%   csv_path    - dephased profile CSV (sidecar: same stem, .yaml)
%   equip_param - (optional) LOAD_EQUIPMENT_CONFIG output
%   options     - (optional) LIBRARY_ADD_CONDITIONS options; .DryRun = true
%                 returns the result without writing lib_file
%
% Output:
%   lib    - updated library entry
%   report - per-condition agreement table (LIBRARY_ADD_CONDITIONS)
%
% See also: LIBRARY_ADD_CONDITIONS, LIBRARY_CONDITION, READ_AXIAL_PROFILES

arguments
    lib_file {mustBeTextScalar}
    csv_path {mustBeTextScalar}
    equip_param = []
    options (1,1) struct = struct()
end

    if nargin < 3 || isempty(equip_param); equip_param = load_equipment_config(); end
    if nargin < 4; options = struct(); end
    dry_run = isfield(options, 'DryRun') && options.DryRun;
    if isfield(options, 'DryRun'); options = rmfield(options, 'DryRun'); end

    lib = yaml.loadFile(lib_file, 'ConvertToArray', true);
    tran_serial = char(lib.meta.tran_serial);
    tran = equip_param.trans.(tran_serial);

    [d, stem] = fileparts(csv_path);
    sidecar = fullfile(d, [stem, '.yaml']);
    if ~isfile(sidecar)
        error('library_ingest_conditions:noSidecar', 'Sidecar not found: %s', sidecar);
    end
    sc = yaml.loadFile(sidecar, 'ConvertToArray', true);
    if ~isfield(sc, 'columns') || isempty(sc.columns)
        error('library_ingest_conditions:noColumns', ...
            'Sidecar %s has no ''columns'' list describing the dephased conditions.', sidecar);
    end
    cols = sc.columns;
    if iscell(cols); cols = [cols{:}]; end

    prof = read_axial_profiles(csv_path);
    if numel(cols) ~= numel(prof.depths_ep_mm)
        error('library_ingest_conditions:columnCount', ...
            '%s lists %d columns but %s has %d profile columns.', sidecar, numel(cols), csv_path, ...
            numel(prof.depths_ep_mm));
    end

    conds = struct('label', {}, 'elem_phase_deg', {}, 'phase_shift_deg', {}, ...
                   'setpoint_ep_mm', {}, 'profile_column', {});
    for k = 1:numel(cols)
        c = cols(k);
        conds(k).label           = char(c.label);
        conds(k).elem_phase_deg  = field_or_empty(c, 'elem_phase_deg');
        conds(k).phase_shift_deg = field_or_empty(c, 'phase_shift_deg');
        conds(k).setpoint_ep_mm  = prof.depths_ep_mm(k);
        conds(k).profile_column  = k;
    end

    options.Profile   = prof;
    options.SourceCsv = csv_path;
    [lib, report] = library_add_conditions(lib, tran, conds, options);
    disp(report);
    warning('library_ingest_conditions:experimental', ...
        'Dephased conditions are EXPERIMENTAL; check the predicted R² before using them.');

    if ~dry_run
        lib.meta.last_updated = datestr(now, 'yyyy-mm-dd'); %#ok<TNOW1,DATST>
        yaml.dumpFile(lib_file, lib);
        fprintf('Wrote %d condition(s) to %s\n', numel(conds), lib_file);
    end
end

function v = field_or_empty(s, f)
    if isfield(s, f) && ~isempty(s.(f)) && ~isa(s.(f), 'yaml.Null'); v = double(s.(f)(:)'); else; v = []; end
end
