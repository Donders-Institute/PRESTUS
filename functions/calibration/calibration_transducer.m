function [lib, report, out] = calibration_transducer(...
    profile_empirical, ...
    equipment_name, ...
    desired_intensity, ...
    desired_focal_distance_ep, ...
    parameters, ...
    sim_id)

% CALIBRATION_TRANSDUCER  Calibrate a transducer to hydrophone profiles at chosen depths and intensities
%
% Config-driven front end to CALIBRATE_LIBRARY_ENTRY, the same calibration
% core that builds and refines the transducer library. For each requested
% focal depth the measured axial profile is taken (a measured depth) or
% interpolated from the neighbouring measured depths, and the element phases,
% per-ring gains and drive are fitted to it analytically (FIT_TRANSDUCER_MODEL).
% Optionally, a free-water k-Wave simulation per depth checks the result and
% rescales the drive to the simulated ISPPA (KWAVE_CALIBRATE_DRIVE).
%
% Two endpoints are written:
%   1. Library entry  {library_path}/{combo}.yaml — the calibrated depths as
%      refinements on top of the existing (user or shipped) entry. Studies use
%      it by setting transducer.library_path to library_path.
%   2. Calibration output folder (path_output_profiles) — per depth ×
%      intensity: element phases and amplitude (CSV table and one PRESTUS
%      transducer YAML each, SAVE_OPTIMIZED_VALUES), per-depth settings per
%      intensity (LIBRARY_SETTINGS_TABLE), the profiles fitted to, the fit
%      report and a figure. All values are resolved from the library entry,
%      so both endpoints agree exactly.
%
% Configuration (parameters.calibration, see config_calibration.yaml):
%   phase_mode        - 'single_depth' (default): fit each depth |
%                       'geometric': no fit, geometric steering with the
%                       measured drive
%   fit_preset        - 'rayleigh_phase_gain' (default) | 'oneil_phase_only'
%                       (CALIBRATION_FIT_OPTIONS)
%   fit_options       - struct of FIT_TRANSDUCER_MODEL options on top of the preset
%   forward_model     - 'rayleigh' | 'oneil' (overrides the preset)
%   skip_front_peak_mm, opt_phase_precession, opt_limits [mm from the bowl],
%   opt_weights       - mapped to NearFieldSkipMm, PhaseConstraint,
%                       FitWindowMm, FocalWeight
%   opt_method, opt_seed - optimiser and seed (oneil_phase_only preset only)
%   elem_phase_correction_deg - prior phase offsets (deg) to regularise towards
%   profile_align, max_extrapolation_mm, max_gap_mm - profile interpolation
%   add_FDO           - 1 (default): profile distances are from the exit
%                       plane | 0: from the bowl apex
%   kwave_drive_calibration - free-water k-Wave drive calibration per depth
%                       (default true; former name run_free_water_sim);
%                       force_kwavearray, axisymmetric2D, kwave_grid_dims set up its grid
%   ring_profile      - within-ring velocity profile (FIT_RING_PROFILE): 'auto'
%                       (default: keep the entry's) | 'fit' | 'none';
%                       ring_profile_n_sub (4), ring_profile_smooth (0.3),
%                       ring_profile_near_field_skip_mm (4) set the fit;
%                       ring_profile_offaxis (plane-scan file, READ_OFFAXIS_SCANS;
%                       '' = axial only) and ring_profile_offaxis_kind ('all')
%                       add off-axis data
%   library_path      - library endpoint (default {path_output}/library)
%   path_output_profiles - calibration output folder (default {path_output}/profiles)
%
% Use as:
%   [lib, report, out] = calibration_transducer(profile_empirical, equipment_name, ...
%       desired_intensity, desired_focal_distance_ep, parameters)
%   ... = calibration_transducer(..., sim_id)
%
% Input:
%   profile_empirical         - axial profile CSV path, measured profiles
%                               (READ_AXIAL_PROFILES struct), or a single
%                               profile struct with axial_intensity and
%                               axial_distance_bowl [mm] (former interface)
%   equipment_name            - combo '{tran_serial}_{ds_serial}' or a transducer serial
%   desired_intensity         - target ISPPA(s) [W/cm²]; the first is the
%                               reference drive stored in the library
%   desired_focal_distance_ep - focal depth(s) from the exit plane [mm];
%                               [] = all measured depths
%   parameters                - PRESTUS config with a calibration section
%   sim_id                    - run number for the k-Wave work folder (default 1)
%
% Output:
%   lib    - the calibrated library entry as written
%   report - per-depth table (CALIBRATE_LIBRARY_ENTRY)
%   out    - .library_file, .output_folder, .settings_files, .report_file
%
% See also: CALIBRATE_LIBRARY_ENTRY, CALIBRATION_SETUP, CALIBRATION_PIPELINE_START,
%           REFINE_TRANSDUCER_LIBRARY, CALIBRATION_FIT_OPTIONS

    if nargin < 6 || isempty(sim_id); sim_id = 1; end
    % calibration settings may come from an older calibration_config.yaml loaded directly
    tmp   = upgrade_config_keys(struct('calibration', parameters.calibration), 'parameters.calibration');
    parameters.calibration = tmp.calibration;
    cal   = parameters.calibration;
    own_lib = '';
    if isfield(cal, 'path_output'); own_lib = fullfile(char(cal.path_output), 'library'); end
    own_libs = [{char(get_or(cal, 'library_path', own_lib))}, transducer_library_folders(parameters)];
    equip = load_equipment_config([], own_libs(~cellfun(@isempty, own_libs)));   % + your own equipment entries

    %% Equipment -------------------------------------------------------------
    combo = char(string(equipment_name));
    if isfield(equip.combos, combo)
        tran_serial = char(equip.combos.(combo).tran_serial);
        ds_serial   = combo(numel(tran_serial) + 2:end);
    elseif isfield(equip.trans, combo)
        tran_serial = combo; ds_serial = '';
    else
        error('calibration_transducer:unknownEquipment', ...
            ['calibration_transducer: ''%s'' is neither an equipment combo nor a transducer serial ' ...
             'in config/equipment/. Since v0.7 the transducer geometry comes from the equipment ' ...
             'library, not from parameters.transducer: add an equipment entry for it, or, for a ' ...
             'virtual transducer, scale its drive from a free-water simulation instead.'], combo);
    end
    tran = equip.trans.(tran_serial);
    bowl_offset = annular_bowl_depth_mm(tran.transducer.annular);

    %% Measured profiles -------------------------------------------------------
    depths = reshape(desired_focal_distance_ep, 1, []);
    [meas, source] = measured_profiles(profile_empirical, depths, bowl_offset, get_or(cal, 'add_FDO', 1));
    if isempty(depths); depths = meas.depths_ep_mm; end
    intensities = reshape(desired_intensity, 1, []);

    %% Options -----------------------------------------------------------------
    mode   = phase_mode_to_mode(get_or(cal, 'phase_mode', 'single_depth'));
    preset = char(get_or(cal, 'fit_preset', 'rayleigh_phase_gain'));
    co = struct('Depths', depths, 'Mode', mode, 'Preset', preset, ...
        'FitOptions', fit_overrides(cal, preset, bowl_offset), 'Intensity', intensities(1), ...
        'Source', source, 'EquipParam', equip, 'MaxGapMm', get_or(cal, 'max_gap_mm', 10), ...
        'InterpolateOptions', struct('Align', char(get_or(cal, 'profile_align', 'label')), ...
            'MaxExtrapolationMm', get_or(cal, 'max_extrapolation_mm', 0)));
    co.KWaveDriveCalibration = logical(get_or(cal, 'kwave_drive_calibration', get_or(cal, 'run_free_water_sim', true)));
    if co.KWaveDriveCalibration
        sim_cfg = parameters;
        if get_or(cal, 'force_kwavearray', 0) == 1; sim_cfg.grid.use_kWaveArray = 1; end
        co.KWaveDriveOptions = struct('SimConfig', sim_cfg, 'EquipParam', equip, ...
            'Axisymmetric', logical(get_or(cal, 'axisymmetric2D', 0)), ...
            'GridDims', get_or(cal, 'kwave_grid_dims', []), ...
            'WorkDir', fullfile(cal.path_output, sprintf('sub-%03d', sim_id), 'kwave_drive_calibration'));
    end
    co.RingProfile = char(get_or(cal, 'ring_profile', 'auto'));
    co.RingProfileOptions = struct('NSub', get_or(cal, 'ring_profile_n_sub', 4), ...
        'Smooth', get_or(cal, 'ring_profile_smooth', 0.3), ...
        'NearFieldSkipMm', get_or(cal, 'ring_profile_near_field_skip_mm', 4));
    if ~isempty(get_or(cal, 'ring_profile_offaxis', ''))
        co.RingProfileOptions.OffAxis = char(cal.ring_profile_offaxis);
        co.RingProfileOptions.OffAxisOptions = struct('Kind', char(get_or(cal, 'ring_profile_offaxis_kind', 'all')));
    end
    warn_obsolete(cal);

    %% Base entry: user library, else shipped, else new ----------------------
    lib_path = char(get_or(cal, 'library_path', fullfile(cal.path_output, 'library')));
    out_file = library_entry_file(lib_path, combo);
    shipped  = library_entry_file(fullfile(get_prestus_path(), 'config', 'calibration'), combo);
    if isfile(out_file)
        lib = yaml.loadFile(out_file, 'ConvertToArray', true);
    elseif isfile(shipped)
        lib = yaml.loadFile(shipped, 'ConvertToArray', true);
    else
        lib = struct('meta', struct('tran_serial', tran_serial, ...
            'created_at', datestr(now, 'yyyy-mm-dd'))); %#ok<TNOW1,DATST>
        if ~isempty(ds_serial); lib.meta.driving_system = ds_serial; end
    end

    %% Calibrate -----------------------------------------------------------------
    [lib, ~, report, prof, diag] = calibrate_library_entry(lib, tran, meas, co);
    disp(report);

    %% Endpoint 1: library entry ---------------------------------------------
    lib = save_library_entry(lib, out_file, mode, report, ...
        struct('calibration_preset', preset, 'calibration_source', source));
    fprintf('Library entry: %s\n', out_file);

    %% Endpoint 2: calibration output folder ---------------------------------
    out_dir = char(get_or(cal, 'path_output_profiles', fullfile(cal.path_output, 'profiles')));
    if ~isfolder(out_dir); mkdir(out_dir); end
    out = struct('library_file', out_file, 'output_folder', out_dir, 'settings_files', {{}});
    fm  = diag.fit_options.ForwardModel;

    csv_name = char(get_or(cal, 'filename_calibrated_CSV', ''));
    if isempty(csv_name)
        prefix = 'PRESTUS_virtual_';
        if isfield(equip.gen, 'prestus_virt_name'); prefix = char(equip.gen.prestus_virt_name); end
        csv_name = [prefix combo '.csv'];
    end
    for d = prof.depths_ep_mm
        for I = intensities
            m = library_model_at_depth(lib, tran, d, I, out_file);
            p = parameters;
            tl = library_transducer(tran, lib);   % sub-rings if the entry has a ring profile
            p.transducer = tl.transducer;
            p.transducer.annular = ring_profile_drive(p.transducer.annular, m.elem_phase_deg, m.elem_p_amp);
            p.calibration.path_output_profiles  = out_dir;
            p.calibration.filename_calibrated_CSV = csv_name;
            p.calibration.equipment_name        = combo;
            p.calibration.desired_focal_distance_ep = round(d, 2);
            p.calibration.desired_intensity     = I;
            p.calibration.phase_mode            = mode;
            p.calibration.forward_model         = fm;
            p.calibration.calibration_type      = 'user_empirical';
            save_optimized_values(p);
        end
    end
    for I = intensities
        f = fullfile(out_dir, sprintf('%s_settings_I%g.csv', combo, I));
        library_settings_table(tran, lib, prof.depths_ep_mm, struct('IsppaWcm2', I, 'CsvPath', f, ...
            'ProfileType', {prof.profile_type}, 'ForwardModel', fm));
        out.settings_files{end+1} = f;
    end
    write_axial_profiles(prof, fullfile(out_dir, [combo '_calibration_profiles.csv']), struct('combo', combo));
    out.report_file = fullfile(out_dir, [combo '_calibration_report.csv']);
    writetable(report, out.report_file);
    if ~isempty(diag.ring_profile)
        plot_ring_profile(tran, lib.global_model.ring_profile, diag.ring_profile, ...
            fullfile(out_dir, [combo '_ring_profile.png']));
    end
    if ~isempty(diag.kwave)
        writetable(diag.kwave, fullfile(out_dir, [combo '_kwave_drive_calibration.csv']));
    end
    plot_calibration_report(prof, report, diag, ...
        sprintf('%s — calibration (%s, %s preset)', strrep(combo, '_', '\_'), mode, preset), ...
        fullfile(out_dir, [combo '_calibration.png']));
    fprintf('Calibration outputs: %s\n', out_dir);
end

% -------------------------------------------------------------------------
function [meas, source] = measured_profiles(src, depths, bowl_offset, add_fdo)
% Measured profiles in READ_AXIAL_PROFILES form, distances from the exit plane
    if ischar(src) || isstring(src)
        source = char(src);
        meas = read_axial_profiles(source);
        if ~add_fdo   % distances in the CSV are from the bowl apex
            meas.z_ep_mm = meas.z_ep_mm - bowl_offset;
        end
    elseif isfield(src, 'depths_ep_mm')
        meas = src; source = get_or(src, 'source', 'profiles struct');
    elseif isfield(src, 'axial_intensity')
        % Former interface: one profile at the (single) requested depth
        assert(numel(depths) == 1, 'calibration_transducer:singleProfile', ...
            'A single axial profile can only calibrate one focal depth.');
        z  = src.axial_distance_bowl(:) - bowl_offset;
        I  = src.axial_intensity(:);
        ok = z >= 0;                       % drop the padding before the exit plane
        meas = struct('depths_ep_mm', depths, 'z_ep_mm', z(ok), 'I', I(ok), 'source', 'profile struct');
        source = meas.source;
    else
        error('calibration_transducer:profiles', ...
            'profile_empirical must be a CSV path, a READ_AXIAL_PROFILES struct or an axial profile struct.');
    end
end

function mode = phase_mode_to_mode(phase_mode)
    switch char(phase_mode)
        case 'single_depth'; mode = 'per_depth';
        case 'geometric';    mode = 'geometric';
        case 'multi_depth'
            error('calibration_transducer:phaseMode', ...
                'calibration.phase_mode ''multi_depth'' (joint depth-independent model) has been removed; use ''single_depth''.');
        otherwise
            error('calibration_transducer:phaseMode', ...
                'Unknown calibration.phase_mode ''%s'' (single_depth | geometric).', char(phase_mode));
    end
end

function fo = fit_overrides(cal, preset, bowl_offset)
% calibration.* keys → FIT_TRANSDUCER_MODEL options on top of the preset
    fo = get_or(cal, 'fit_options', struct());
    if isempty(fo); fo = struct(); end
    if isfield(cal, 'forward_model') && ~isempty(cal.forward_model)
        fo.ForwardModel = lower(char(cal.forward_model));
    end
    if isfield(cal, 'skip_front_peak_mm'); fo.NearFieldSkipMm = cal.skip_front_peak_mm; end
    pp = get_or(cal, 'opt_phase_precession', false);
    if ischar(pp) || isstring(pp); fo.PhaseConstraint = char(pp); end
    lim = get_or(cal, 'opt_limits', []);
    if numel(lim) == 2; fo.FitWindowMm = double(lim(:)') - bowl_offset; end
    w = get_or(cal, 'opt_weights', 0);
    if w > 0; fo.FocalWeight = w; end
    if strcmp(preset, 'oneil_phase_only')
        if isfield(cal, 'opt_method'); fo.Optimizer = lower(char(cal.opt_method)); end
        if isfield(cal, 'opt_seed') && ~isempty(cal.opt_seed); fo.Seed = cal.opt_seed; end
    end
    corr = get_or(cal, 'elem_phase_correction_deg', []);
    if ~isempty(corr)
        fo.Prior = struct('elem_phase_correction_deg', double(corr(:)'), 'z_shift_mm', 0);
    end
end

function warn_obsolete(cal)
    if isfield(cal, 'save_elem_correction') && cal.save_elem_correction
        warning('calibration_transducer:obsolete', ['calibration.save_elem_correction is obsolete: ' ...
            'the calibrated offsets are written to the library entry (library_path).']);
    end
end

function v = get_or(s, f, d)
% Field value, or d when absent or empty (YAML '' loads as the string "")
    if isstruct(s) && isfield(s, f) && ~isempty(s.(f)) && ...
            ~((isstring(s.(f)) || ischar(s.(f))) && strlength(string(s.(f))) == 0)
        v = s.(f);
    else
        v = d;
    end
end
