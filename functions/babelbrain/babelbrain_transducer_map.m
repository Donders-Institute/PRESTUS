function [preset_name, elem_csv_path] = babelbrain_transducer_map(parameters, equip_param, out_dir)
% BABELBRAIN_TRANSDUCER_MAP  Map a PRESTUS transducer to a BabelBrain preset
%
% Attempts to identify a BabelBrain named preset matching the PRESTUS
% transducer by (1) serial-number lookup against a known-device table and
% (2) geometry matching by {elem_n, curv_radius_mm, max_OD_mm, freq_Hz}.
%
% If a known preset is found, preset_name is returned and elem_csv_path is
% empty. If no preset matches, a custom element CSV is generated from the
% PRESTUS geometry (annular arrays only) and written to out_dir; the caller
% must place this CSV in BabelBrain's TranscranialModeling/ directory and
% implement a matching RUN_SIM subclass (see NOTES below).
%
% BabelBrain preset names (as of BabelBrain v2, 2026):
%   Single, BSonix, CTX_500, CTX_250, CTX_250_2ch, DPX_500, DPXPC_300,
%   R15287, R15473, R15646, H317, H246, IGT64_500, H301, REMOPD,
%   I12378, ATAC, R15148, DomeTx
%
% Use as:
%   [preset_name, elem_csv_path] = babelbrain_transducer_map(parameters, equip_param, out_dir)
%
% Input:
%   parameters  - PRESTUS config struct; uses parameters.transducer(1)
%   equip_param - struct from load_equipment_config (may be empty)
%   out_dir     - directory for writing custom element CSV (if needed)
%
% Output:
%   preset_name   - BabelBrain preset string, e.g. 'CTX_500'; empty if no
%                   match and custom CSV was generated
%   elem_csv_path - path to generated custom element CSV, or '' if a named
%                   preset was matched
%
% NOTES ON CUSTOM TRANSDUCERS:
%   When no preset matches, a CSV is written in the BabelBrain element format
%   for annular arrays (columns: RingIndex, InnerRadius_mm, OuterRadius_mm,
%   CurvRadius_mm). To use it:
%     1. Copy the CSV to BabelBrain/TranscranialModeling/
%     2. Create a subclass in BabelIntegrationANNULAR_ARRAY.py (annular) or
%        BabelIntegrationCONCAVE_PHASEDARRAY.py (matrix/phased) with a
%        _GetTransducerName() method returning the preset name and
%        _LoadTransducerInfo() reading the CSV.
%     3. Register the class in CalculateFieldProcess.py's TxSystem dispatch.
%   See BabelBrain/TranscranialModeling/BabelIntegrationANNULAR_ARRAY.py
%   and OfflineBatchExamples/ for the expected subclass structure.
%
% See also: BABELBRAIN_WRAPPER, WRITE_BABELBRAIN_INPUT

arguments
    parameters  (1,1) struct
    equip_param       struct  = struct()
    out_dir     (1,:) char    = tempdir()
end

preset_name   = '';
elem_csv_path = '';

tr = parameters.transducer(1);
freq_hz = tr.freq_hz;

%% Known-preset lookup table (serial → preset)
% Add entries here as new devices are validated.
% Keys are tran.serial strings from the equipment YAML.
SERIAL_MAP = struct( ...
    'IS_PCD15287_01001', 'R15287', ...   % Sonic Concepts H-115 / R15287
    'IS_PCD15473_01001', 'R15473', ...   % Sonic Concepts 500 kHz annular 4-ring
    'IS_PCD15646_01001', 'R15646', ...   % Sonic Concepts R15646
    'CTX500_default',    'CTX_500', ...
    'CTX250_default',    'CTX_250', ...
    'H317_default',      'H317', ...
    'H246_default',      'H246', ...
    'H301_default',      'H301', ...
    'IGT64_500',         'IGT64_500' ...
);

% Try serial lookup: prefer the active transducer's own serial field so that
% multi-transducer equipment configs pick the correct entry rather than
% always falling back to the first fieldname.
tran_serial = '';
if isfield(tr, 'serial') && ~isempty(tr.serial)
    tran_serial = tr.serial;
elseif ~isempty(equip_param) && isfield(equip_param, 'trans')
    fnames = fieldnames(equip_param.trans);
    if ~isempty(fnames)
        tran_serial = fnames{1};
    end
end

if ~isempty(tran_serial) && isfield(SERIAL_MAP, tran_serial)
    preset_name = SERIAL_MAP.(tran_serial);
    fprintf('[babelbrain_transducer_map] Serial ''%s'' → preset ''%s''\n', ...
        tran_serial, preset_name);
    return
end

%% Geometry-based matching for annular arrays
% Known preset specs: {preset_name, n_elem, curv_radius_mm, max_od_mm, freq_Hz}
% Geometry values from BabelBrain/BabelBrain/Babel_<preset>/default.yaml
% Tolerances: ±2 mm radius, ±3 mm OD, ±5% frequency
GEOM_TABLE = {
    % preset        n_elem  curv_r_mm  max_od_mm  freq_hz
    'CTX_500',    7,   62.94,  64.0,  500e3;   % CTX500 default.yaml
    'CTX_250',    7,   62.94,  64.0,  250e3;
    'R15287',     4,   62.94,  64.0,  500e3;   % R15287 default.yaml (4-ring, 500 kHz)
    'R15473',    10,  100.00,  65.0,  300e3;   % R15473 default.yaml (10-ring, 300 kHz)
    'R15646',    10,   72.00,  72.0,  500e3;   % R15646 (10-ring, 500 kHz; verify from yaml)
    'H317',       1,   62.94,  64.0,  500e3;
    'H246',       1,   62.94,  32.0,  500e3;
    'H301',       1,   40.00,  32.0, 1000e3;
    'Single',     1,     [],    [],    [];      % fallback single-element
};

if isfield(tr, 'annular') && isfield(tr.annular, 'elem_n')
    n_elem      = tr.annular.elem_n;
    curv_r      = tr.annular.curv_radius_mm;
    max_od      = max(tr.annular.elem_od_mm);

    for k = 1:size(GEOM_TABLE, 1)
        ref_n    = GEOM_TABLE{k,2};
        ref_r    = GEOM_TABLE{k,3};
        ref_od   = GEOM_TABLE{k,4};
        ref_freq = GEOM_TABLE{k,5};

        n_match    = (n_elem == ref_n);
        r_match    = isempty(ref_r)    || abs(curv_r - ref_r) < 2.0;
        od_match   = isempty(ref_od)   || abs(max_od - ref_od) < 3.0;
        freq_match = isempty(ref_freq) || abs(freq_hz - ref_freq) / ref_freq < 0.05;

        if n_match && r_match && od_match && freq_match
            preset_name = GEOM_TABLE{k,1};
            fprintf(['[babelbrain_transducer_map] Geometry match: n=%d  R=%.1fmm  ' ...
                'OD=%.1fmm  f=%.0fkHz → preset ''%s''\n'], ...
                n_elem, curv_r, max_od, freq_hz/1e3, preset_name);
            return
        end
    end
end

%% No match — generate custom element CSV
fprintf(['[babelbrain_transducer_map] No known preset matched. ' ...
    'Generating custom element CSV.\n']);

if ~isfield(tr, 'annular') || ~isfield(tr.annular, 'elem_n')
    error(['babelbrain_transducer_map: Custom CSV generation is only implemented ' ...
        'for annular arrays (tr.type = ''annular''). Matrix/phased arrays require ' ...
        'a manual BabelBrain RUN_SIM subclass.']);
end

freq_kHz = round(freq_hz / 1e3);
csv_name = sprintf('CustomAnnular_%delem_R%.0fmm_%dkHz.csv', ...
    tr.annular.elem_n, tr.annular.curv_radius_mm, freq_kHz);
elem_csv_path = fullfile(out_dir, csv_name);

% BabelBrain annular array element CSV format (BabelIntegrationANNULAR_ARRAY):
%   Column 1: ring index (0-based)
%   Column 2: inner radius [mm]
%   Column 3: outer radius [mm]
%   Column 4: radius of curvature [mm]  (same for all rings)
fid = fopen(elem_csv_path, 'w');
fprintf(fid, 'RingIndex,InnerRadius_mm,OuterRadius_mm,CurvRadius_mm\n');
n  = tr.annular.elem_n;
id = tr.annular.elem_id_mm(:)';
od = tr.annular.elem_od_mm(:)';
rc = tr.annular.curv_radius_mm;
for k = 1:n
    fprintf(fid, '%d,%.4f,%.4f,%.4f\n', k-1, id(k)/2, od(k)/2, rc);
end
fclose(fid);

fprintf('[babelbrain_transducer_map] Element CSV written: %s\n', elem_csv_path);
fprintf(['  Next steps:\n' ...
    '  1. Copy CSV to BabelBrain/TranscranialModeling/\n' ...
    '  2. Create a RUN_SIM subclass in BabelIntegrationANNULAR_ARRAY.py\n' ...
    '     with _GetTransducerName() returning a unique preset string.\n' ...
    '  3. Register the preset in BabelBrain/CalculateFieldProcess.py.\n' ...
    '  4. Set parameters.babelbrain.tx_system to the new preset name.\n']);

end
