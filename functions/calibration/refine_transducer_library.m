function [out_file, fits, report] = refine_transducer_library(combo_name, profile_csv, out_library_path, options)
% REFINE_TRANSDUCER_LIBRARY  Optional refinement of a library entry to your own hydrophone profiles
%
% Fits the transducer model (FIT_TRANSDUCER_MODEL) to measured axial profiles
% and writes the result to a library YAML in out_library_path — a user
% folder, so the shipped library in config/calibration/ is not modified.
% Use the result by setting transducer.library_path in the study config.
%
% Two modes:
%   'per_depth' — (default) fit phase offsets and gains at each given depth,
%                 regularised towards the entry's global model (geometric
%                 steering, or a legacy entry's static offset), and store them
%                 with the drive that gives Intensity at that depth under
%                 'refinement'. load_transducer_from_library uses them at and
%                 between those depths and fades back to the global model
%                 outside (refinement.blend_mm).
%   'geometric' — no fit: geometric steering with the measured drive per depth.
%
% Optionally (KWave = true) each refined depth is then checked in a
% free-water k-Wave simulation of the runtime library path, and its drive
% corrected to the simulated ISPPA (KWAVE_CALIBRATE_DRIVE).
%
% Use as:
%   out_file = refine_transducer_library(combo_name, profile_csv, out_library_path)
%   [out_file, fits, report] = refine_transducer_library(..., options)
%
% Typical use — measure axial profiles at the depths you will stimulate, then:
%   refine_transducer_library('CTX_250_001_SC_105_010', 'my_profiles.csv', ...
%       '/project/my_study/transducer_library', struct('Plot', true));
% and set transducer.library_path to that folder in the study config.
% See examples/demo_transducer_refinement.m.
%
% Input:
%   combo_name       - '{tran_serial}_{ds_serial}', e.g. 'CTX_250_001_SC_105_010'
%   profile_csv      - axial profile CSV (doc_calibration.md, Appendix)
%   out_library_path - folder to write {combo_name}.yaml into (created if needed).
%                      An existing entry there is updated; otherwise the shipped
%                      entry (config/calibration/) is copied first.
%   options          - struct, all optional:
%     .Mode           - 'per_depth' (default) | 'geometric'
%     .Preset         - per-depth fit settings, CALIBRATION_FIT_OPTIONS:
%                       'rayleigh_phase_gain' (default) | 'oneil_phase_only'
%                       (former PERFORM_GLOBAL_SEARCH behaviour)
%     .MaxGapMm       - warn when a depth is interpolated across a larger gap
%                       between measured depths [mm] (default 10)
%     .Depths         - depths [mm from EP] to calibrate (default: all in the CSV).
%                       A depth within 0.5 mm of a measured one uses that
%                       measurement; any other is interpolated from the
%                       neighbouring measurements (INTERPOLATE_AXIAL_PROFILES);
%                       depths that cannot be derived are skipped with a warning
%     .InterpolateOptions - passed to INTERPOLATE_AXIAL_PROFILES (Align,
%                       MaxExtrapolationMm; default no extrapolation)
%     .ExportProfiles - true | false | CSV path: write the profiles fitted to
%                       (measured and derived, with a provenance sidecar);
%                       default true when any depth is derived, to
%                       {combo_name}_calibration_profiles.csv next to the YAML
%     .ExportSettings - true | CSV path: write the per-depth settings the entry
%                       resolves to (LIBRARY_SETTINGS_TABLE) to
%                       {combo_name}_settings.csv (default false)
%     .SettingsDepths - depths for ExportSettings (default: the calibrated depths)
%     .Intensity      - reference ISPPA for the amplitude [W/cm²] (default 30)
%     .FitOptions     - FIT_TRANSDUCER_MODEL options (FreeParams, ForwardModel,
%                       NearFieldSkipMm, HuberDelta, Lambda, NStarts, ...) on top
%                       of the Preset; Prior defaults to the entry's global
%                       model, with the z shift fixed at 0. Adding 'z_shift' fits an axial
%                       offset as a diagnostic (refinement.fit_z_shift_mm); it is
%                       not applied at run time, so the default keeps it at 0.
%     .BlendMm        - refinement.blend_mm for per_depth (default 5)
%     .KWaveDriveCalibration - true: calibrate the drive of each refined depth
%                       in free-water k-Wave (former name KWave, still read)
%                       (default false)
%     .KWaveDriveOptions - options for KWAVE_CALIBRATE_DRIVE (SimConfig, WorkDir, ...;
%                       former name KWaveOptions, still read)
%     .Plot           - true: save {combo_name}_refinement.png next to the YAML
%                       (measured vs library vs refined, per depth)
%     .RingProfile    - within-ring velocity profile (CALIBRATE_LIBRARY_ENTRY):
%                       'auto' (default: keep the entry's) | 'fit' | 'none' | struct
%     .RingProfileOptions - FIT_RING_PROFILE options (NSub, Smooth, NearFieldSkipMm, HeldOut)
%     .EquipParam     - struct from load_equipment_config (default: loaded)
%     .ShippedLibrary - source of the base entry (default config/calibration/)
%
% Output:
%   out_file - path of the written YAML
%   fits     - FIT_TRANSDUCER_MODEL results, one per depth
%   report   - table per calibrated depth: profile type, R², |focal error| and
%              near-field excess of the entry before (library) and after refinement
%              (on derived profiles these compare with the derived profile)
%
% See also: CALIBRATE_LIBRARY_ENTRY, FIT_TRANSDUCER_MODEL, LOAD_TRANSDUCER_FROM_LIBRARY, KWAVE_CALIBRATE_DRIVE,
%           LIBRARY_ENTRY_FROM_MODEL, READ_AXIAL_PROFILES, CALIBRATION_GUI

    if nargin < 4; options = struct(); end
    mode      = opt(options, 'Mode', 'per_depth');
    intensity = opt(options, 'Intensity', 30);
    fit_opts  = opt(options, 'FitOptions', struct());
    blend_mm  = opt(options, 'BlendMm', 5);
    do_plot   = opt(options, 'Plot', false);
    equip     = opt(options, 'EquipParam', []);
    shipped   = opt(options, 'ShippedLibrary', fullfile(get_prestus_path(), 'config', 'calibration'));
    assert(ismember(mode, {'per_depth', 'geometric'}), 'refine_transducer_library:mode', ...
        'refine_transducer_library: Mode must be ''per_depth'' or ''geometric'' (''static'' has been removed).');
    if isempty(equip); equip = load_equipment_config(); end

    %% Base entry ----------------------------------------------------------
    if ~isfolder(out_library_path); mkdir(out_library_path); end
    out_file  = fullfile(out_library_path, [combo_name '.yaml']);
    base_file = fullfile(shipped, [combo_name '.yaml']);
    if isfile(out_file)
        lib = yaml.loadFile(out_file, 'ConvertToArray', true);
    elseif isfile(base_file)
        lib = yaml.loadFile(base_file, 'ConvertToArray', true);
    else
        error(['refine_transducer_library: no entry for ''%s'' in %s or %s.\n' ...
            'Create one with build_transducer_library first.'], combo_name, out_library_path, shipped);
    end
    tran = equip.trans.(char(lib.meta.tran_serial));

    %% Calibrate (CALIBRATE_LIBRARY_ENTRY) -----------------------------------
    meas = read_axial_profiles(profile_csv);
    co = struct('Mode', mode, 'Preset', opt(options, 'Preset', 'rayleigh_phase_gain'), 'FitOptions', fit_opts, ...
        'Intensity', intensity, 'BlendMm', blend_mm, 'Source', profile_csv, ...
        'InterpolateOptions', opt(options, 'InterpolateOptions', struct()), ...
        'KWaveDriveCalibration', opt(options, 'KWaveDriveCalibration', opt(options, 'KWave', false)), ...
        'KWaveDriveOptions', opt(options, 'KWaveDriveOptions', opt(options, 'KWaveOptions', struct())), ...
        'EquipParam', equip, 'RingProfile', opt(options, 'RingProfile', 'auto'), ...
        'RingProfileOptions', opt(options, 'RingProfileOptions', struct()));
    if isfield(options, 'Depths'); co.Depths = options.Depths; end
    if isfield(options, 'MaxGapMm'); co.MaxGapMm = options.MaxGapMm; end
    [lib, fits, report, prof, diag] = calibrate_library_entry(lib, tran, meas, co);
    derived = ~strcmp(prof.profile_type, 'measured');
    disp(report);

    %% Exports ---------------------------------------------------------------
    % The profiles the entry was fitted to (always when any is derived), so
    % an interpolated calibration can be traced and reproduced
    prof_csv = opt(options, 'ExportProfiles', any(derived));
    if islogical(prof_csv) || isnumeric(prof_csv)
        if prof_csv; prof_csv = fullfile(out_library_path, [combo_name '_calibration_profiles.csv']); else; prof_csv = ''; end
    end
    if ~isempty(prof_csv)
        write_axial_profiles(prof, prof_csv, struct('combo', combo_name));
        fprintf('Calibration profiles: %s\n', prof_csv);
        lib.refinement.calibration_profiles_csv = prof_csv;
    end

    %% Provenance and write ------------------------------------------------
    save_library_entry(lib, out_file, mode, report);
    fprintf('refine_transducer_library (%s): %s\n', mode, out_file);

    % Per-depth settings the written entry resolves to (optional)
    set_csv = opt(options, 'ExportSettings', false);
    if islogical(set_csv) || isnumeric(set_csv)
        if set_csv; set_csv = fullfile(out_library_path, [combo_name '_settings.csv']); else; set_csv = ''; end
    end
    if ~isempty(set_csv)
        sd = opt(options, 'SettingsDepths', prof.depths_ep_mm);
        so = struct('IsppaWcm2', intensity, 'CsvPath', set_csv, ...
                    'NearFieldSkipMm', opt(fit_opts, 'NearFieldSkipMm', 10));
        if isfield(diag.fit_options, 'ForwardModel'); so.ForwardModel = diag.fit_options.ForwardModel; end
        library_settings_table(tran, lib, sd, so);
        fprintf('Settings per depth: %s\n', set_csv);
    end

    if do_plot && ~isempty(diag.ring_profile)
        plot_ring_profile(tran, lib.global_model.ring_profile, diag.ring_profile, ...
            fullfile(out_library_path, [combo_name '_ring_profile.png']));
    end
    if do_plot
        plot_calibration_report(prof, report, diag, ...
            sprintf('%s — %s refinement', strrep(combo_name, '_', '\_'), mode), ...
            fullfile(out_library_path, [combo_name '_refinement.png']));
    end
end

function v = opt(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end
