function build_geometric_library_entry(tran_serial, equip_path, library_path, varargin)
% BUILD_GEOMETRIC_LIBRARY_ENTRY  Create a library entry from transducer geometry alone
%
% Computes geometric phases at a set of standard depths using ring-geometry
% path-length delays (elem_id_mm, elem_od_mm, curv_radius_mm from the equipment
% YAML) and stores them as a global_model entry in the library.
%
% No external phase table file (.ini or steer CSV) is required.  Results are
% suitable for users who do not have access to manufacturer calibration files.
%
% This produces a DS-agnostic baseline entry ({tran_serial}.yaml) that is
% immediately usable in study configs, and serves as an automatic fallback
% when a DS-specific library file is absent.
%
% Use as:
%   build_geometric_library_entry(tran_serial, equip_path, library_path)
%   build_geometric_library_entry(tran_serial, equip_path, library_path, ...
%       'Depths', [40 50 60 70 80], 'SoundSpeed', 1500, 'Combo', 'MY_DS')
%
% Input:
%   tran_serial  - equipment serial matching a YAML in equip_path
%   equip_path   - path to config/equipment/
%   library_path - path to config/calibration/  (output folder)
%   'Depths'     - (optional) focal depths from exit plane [mm]; default auto
%   'SoundSpeed' - (optional) water sound speed [m/s]; default 1500
%   'Combo'      - (optional) DS serial; if set writes {tran_serial}_{ds_serial}.yaml
%
% Output:
%   Writes/updates {library_path}/{tran_serial}.yaml (or combo key) with a
%   global_model section containing geometric phases. Amplitude scaling is
%   not populated — add it via calibrate_transducer.
%
% See also: SET_REAL_PHASES, GENERATE_TRAN_INI_FROM_GEOMETRY,
%           UPDATE_TRANSDUCER_LIBRARY, VALIDATE_TRANSDUCER_LIBRARY

    p = inputParser();
    p.addParameter('Depths',      []);
    p.addParameter('SoundSpeed',  1500);
    p.addParameter('Combo',       '');
    p.parse(varargin{:});
    opts = p.Results;

    if nargin < 3 || isempty(library_path)
        library_path = fullfile(get_prestus_path(), 'config', 'calibration');
    end
    if nargin < 2 || isempty(equip_path)
        equip_path = fullfile(get_prestus_path(), 'config', 'equipment');
    end

    equip = load_equipment_config(equip_path);
    if ~isfield(equip.trans, tran_serial)
        error('build_geometric_library_entry: serial ''%s'' not found in equipment config at %s.', ...
            tran_serial, equip_path);
    end
    tran = equip.trans.(tran_serial);

    %% Resolve output key
    ds_serial = opts.Combo;
    if ~isempty(ds_serial)
        combo_key = [tran_serial '_' ds_serial];
        if ~isfield(tran, 'combos') || ~isfield(tran.combos, ds_serial)
            warning('build_geometric_library_entry:missingCombo', ...
                ['DS ''%s'' not found under %s.combos — writing DS-agnostic entry ''%s.yaml'' instead.\n' ...
                'For a DS-specific entry, add a combos.%s entry to config/equipment/%s.yaml.'], ...
                ds_serial, tran_serial, tran_serial, ds_serial, tran_serial);
            combo_key = tran_serial;
        end
    else
        combo_key = tran_serial;
    end

    %% Determine depth grid
    depths_mm = opts.Depths;
    if isempty(depths_mm)
        % Auto: 10 steps spanning min_foc to max_foc (or curv_radius as fallback)
        r_c = tran.transducer.annular.curv_radius_mm;
        lo  = r_c * 0.15;
        hi  = r_c * 0.95;
        if isfield(tran, 'min_foc') && ~isempty(tran.min_foc); lo = tran.min_foc; end
        if isfield(tran, 'max_foc') && ~isempty(tran.max_foc); hi = tran.max_foc; end
        depths_mm = round(linspace(lo, hi, 10), 1);
        fprintf('Auto depth grid: %s mm\n', mat2str(depths_mm));
    end

    %% Build minimal parameters struct for set_real_phases
    params = struct();
    params.medium_properties.water.sound_speed = opts.SoundSpeed;
    params.medium_properties.water.density     = 994;
    params.transducer = tran.transducer;
    params.transducer.n_elem = tran.transducer.annular.elem_n;
    params = load_parameters(params);

    %% Compute geometric phases at each depth
    n_elem   = tran.transducer.annular.elem_n;
    n_depths = numel(depths_mm);
    phase_mat = zeros(n_depths, n_elem);

    for j = 1:n_depths
        d = depths_mm(j);
        params.transducer.focal_distance_ep   = d;
        params.transducer.focal_distance_bowl = d + ...
            (tran.transducer.annular.curv_radius_mm - tran.transducer.annular.dist_geom_ep_mm);
        phase_mat(j, :) = set_real_phases([], tran, d, params);
    end

    %% Write library YAML
    out_file = fullfile(library_path, [combo_key '.yaml']);

    % Load existing file to preserve any per-depth refinements already present
    if isfile(out_file)
        lib = yaml.loadFile(out_file, 'ConvertToArray', true);
        fprintf('Updating existing library file: %s\n', out_file);
    else
        lib = struct();
    end

    lib.meta.tran_serial       = tran_serial;
    lib.meta.created_at        = datestr(now, 'yyyy-mm-dd'); %#ok<TNOW1,DATST>
    lib.meta.calibration_type  = 'geometric';
    lib.meta.note              = ['Geometric phases from ring-geometry path-length delays. ' ...
        'No external phase table required. ' ...
        'Run build_transducer_library for an empirically calibrated entry.'];
    try
        lib.meta.prestus_version = prestus_version();
    catch
        lib.meta.prestus_version = 'unknown';
    end

    lib.global_model.depths_ep_mm  = depths_mm;
    lib.global_model.elem_phase_deg = num2cell(phase_mat, 2)';
    % Amplitude not set — library lookup will warn and require inline or calibrated value

    yaml.dumpFile(out_file, lib);
    fprintf('Written: %s\n  Depths: %s mm\n  Type: geometric (no DS correction)\n', ...
        out_file, mat2str(depths_mm));
end
