function prestus_gui(init_params)
%PRESTUS_GUI  Interactive parameter setup and simulation launcher.
%
%   Opens a single-window tabbed GUI that:
%     - Loads config_default.yaml and highlights mandatory/optional fields
%     - Allows loading and saving study-specific YAML configs
%     - Submits simulations (local via backgroundPool, or HPC via scheduler)
%     - Displays simulation images and the HTML report inline
%
%   Usage:
%     prestus_gui()
%     prestus_gui(parameters)   % pre-load a parameter struct
%
%   Requirements: MATLAB R2023b, PRESTUS on the MATLAB path.

    if nargin < 1
        init_params = [];
    end

%% ── Figure ────────────────────────────────────────────────────────────────

st = sty();

fig = uifigure( ...
    'Name',            'PRESTUS', ...
    'Position',        [80 60 820 680], ...
    'Color',           st.bg_fig, ...
    'Resize',          'on', ...
    'AutoResizeChildren', 'off', ...
    'Visible',         'off');

% App state stored in UserData
app.dir_output  = '';
app.config_path = '';
fig.UserData    = app;

%% ── Header bar ────────────────────────────────────────────────────────────

hdr = uipanel(fig, ...
    'Position',       [0 fig.Position(4)-42 fig.Position(3) 42], ...
    'BackgroundColor', st.bg_dark, ...
    'BorderType',     'none');

uilabel(hdr, ...
    'Text',       'PRESTUS', ...
    'Position',   [12 4 120 32], ...
    'FontName',   st.font, ...
    'FontSize',   st.fs_xl, ...
    'FontWeight', 'bold', ...
    'FontColor',  [1 1 1]);

uilabel(hdr, ...
    'Text',       'PREprocessing & Simulations for Transcranial Ultrasound Stimulation', ...
    'Position',   [130 10 480 20], ...
    'FontName',   st.font, ...
    'FontSize',   st.fs_sm, ...
    'FontColor',  [0.65 0.75 0.78]);   % teal-tinted grey

lbl_config = uilabel(hdr, ...
    'Text',       'Config: defaults', ...
    'Position',   [130 -4 480 16], ...
    'FontName',   st.font, ...
    'FontSize',   10, ...
    'FontColor',  [0.5 0.6 0.63]);

% Toolbar buttons (top-right)
btn_load = uibutton(hdr, ...
    'Text',            '⬆ Load', ...
    'Position',        [fig.Position(3)-300 6 88 28], ...
    'FontName',        st.font, 'FontSize', st.fs_sm, ...
    'BackgroundColor', [0.043 0.482 0.600], ...
    'FontColor',       [1 1 1], ...
    'ButtonPushedFcn', @(~,~) cb_load_yaml()); %#ok<NASGU>

btn_save = uibutton(hdr, ...
    'Text',            '⬇ Save', ...
    'Position',        [fig.Position(3)-204 6 88 28], ...
    'FontName',        st.font, 'FontSize', st.fs_sm, ...
    'BackgroundColor', [0.043 0.482 0.600], ...
    'FontColor',       [1 1 1], ...
    'ButtonPushedFcn', @(~,~) cb_save_yaml()); %#ok<NASGU>

btn_reset = uibutton(hdr, ...
    'Text',            '↺ Defaults', ...
    'Position',        [fig.Position(3)-108 6 98 28], ...
    'FontName',        st.font, 'FontSize', st.fs_sm, ...
    'BackgroundColor', [0.027 0.322 0.404], ...
    'FontColor',       [0.85 0.92 0.95], ...
    'ButtonPushedFcn', @(~,~) load_defaults()); %#ok<NASGU>

%% ── Tab group ─────────────────────────────────────────────────────────────

tg = uitabgroup(fig, ...
    'Position',      [0 0 fig.Position(3) fig.Position(4)-42], ...
    'TabLocation',   'left');

tabs.io          = uitab(tg, 'Title', '  I/O & Paths   ');
tabs.simulation  = uitab(tg, 'Title', '  Simulation    ');
tabs.transducer  = uitab(tg, 'Title', '  Transducer    ');
tabs.grid        = uitab(tg, 'Title', '  Grid          ');
tabs.medium      = uitab(tg, 'Title', '  Medium        ');
tabs.thermal     = uitab(tg, 'Title', '  Thermal       ');
tabs.hpc         = uitab(tg, 'Title', '  HPC           ');
tabs.advanced         = uitab(tg, 'Title', '  Advanced      ');
tabs.placement        = uitab(tg, 'Title', '  Placement     ');
tabs.multitransducer  = uitab(tg, 'Title', '  Multi-Transducer  ');
tabs.calibration      = uitab(tg, 'Title', '  Calibration   ');
tabs.run         = uitab(tg, 'Title', '  ▶  Run        ');
tabs.results     = uitab(tg, 'Title', '  Results       ');

for f = fieldnames(tabs)'
    tabs.(f{1}).BackgroundColor = st.bg_panel;
end

%% ── Build each tab ────────────────────────────────────────────────────────

build_tab_io(tabs.io);
build_tab_simulation(tabs.simulation);
build_tab_transducer(tabs.transducer);
build_tab_grid(tabs.grid);
build_tab_medium(tabs.medium);
build_tab_thermal(tabs.thermal);
build_tab_hpc(tabs.hpc);
build_tab_advanced(tabs.advanced);
build_tab_placement(tabs.placement);
build_tab_multitransducer(tabs.multitransducer);
build_tab_calibration(tabs.calibration);
build_tab_run(tabs.run);
build_tab_results(tabs.results);

%% ── Build widget tag → handle cache ──────────────────────────────────────

% Collect all tagged descendants once so set_widget/findobj callers can use
% the cache instead of walking the full component tree on each lookup.
widget_cache = build_widget_cache();

%% ── Wire cross-tab callbacks ──────────────────────────────────────────────

% Thermal timing fields are only meaningful when thermal simulation is on.
% Disable them initially (run_heating_sims defaults to false) and toggle
% whenever the checkbox changes.
h_heating = findobj(fig, 'Tag', 'modules.run_heating_sims');
if ~isempty(h_heating)
    h_heating.ValueChangedFcn = @(cb,~) cb_toggle_thermal_timing(cb.Value);
    cb_toggle_thermal_timing(h_heating.Value);  % apply initial state
end

%% ── Show window, then load defaults ───────────────────────────────────────

fig.Visible = 'on';
drawnow;

load_defaults();

if ~isempty(init_params) && isstruct(init_params)
    apply_params_to_gui(init_params);
    if isfield(init_params, 'config_path') && ~isempty(init_params.config_path)
        set_config_label(init_params.config_path);
    end
end

%% ════════════════════════════════════════════════════════════════════════
%%  TAB BUILDERS
%% ════════════════════════════════════════════════════════════════════════

    %% ── Tab 1: I/O & Paths ────────────────────────────────────────────
    function build_tab_io(t)
        gl = tab_grid(t, 39, {200, '1x', 55, 80});

        sec(gl, 1, 'Subject');
        lbl(gl, 2, 1, '* Subject ID');
        nedt(gl, 2, 2, 'subject_id', 1);
        sp(gl, 3);

        sec(gl, 4, 'Data Paths');
        lbl(gl, 5, 1, '* Structural MRI path');
        edt(gl, 5, 2, 'path.anat', '', 1);
        lbl(gl, 5, 3, '');
        brw(gl, 5, 4, 'path.anat', 'dir');

        lbl(gl, 6, 1, '* Simulation output path');
        edt(gl, 6, 2, 'path.sim', '', 1);
        brw(gl, 6, 4, 'path.sim', 'dir');

        lbl(gl, 7, 1, '* Segmentation path');
        edt(gl, 7, 2, 'path.seg', '', 1);
        brw(gl, 7, 4, 'path.seg', 'dir');

        lbl(gl, 8, 1, 'Localite path');
        edt(gl, 8, 2, 'path.localite', '', 0);
        brw(gl, 8, 4, 'path.localite', 'dir');

        lbl(gl, 9, 1, 'T1 filename pattern');
        edt(gl, 10, 2, 'path.t1_pattern', 'sub-%1$03d_T1w.nii*', 0);
        note_lbl(gl, 11, 'Relative to path.anat. Use %1$03d for subject ID substitution.');

        lbl(gl, 12, 1, 'T2 filename pattern');
        edt(gl, 12, 2, 'path.t2_pattern', '', 0);
        note_lbl(gl, 13, 'Optional. Leave empty to run SimNIBS charm with T1 only.');

        sec(gl, 14, 'Environment');
        lbl(gl, 15, 1, '* SimNIBS bin path');
        edt(gl, 15, 2, 'startup.simnibs_bin_path', '', 1);
        brw(gl, 15, 4, 'startup.simnibs_bin_path', 'dir');

        lbl(gl, 16, 1, 'Extra paths (addpath)');
        edt(gl, 16, 2, 'startup.paths_to_add', '', 0);
        note_lbl(gl, 17, 'Semicolon-separated absolute paths.');

        lbl(gl, 18, 1, 'Extra subpaths (genpath)');
        edt(gl, 18, 2, 'startup.subpaths_to_add', '', 0);
        note_lbl(gl, 19, 'Semicolon-separated paths; each added recursively.');
        sp(gl, 20);

        sec(gl, 21, 'Output Settings');
        lbl(gl, 22, 1, 'Output affix');
        edt(gl, 22, 2, 'io.output_affix', '', 0);

        lbl(gl, 23, 1, 'Overwrite files');
        drp(gl, 23, 2, 'io.overwrite_files', {'always','never','ask'}, 'always');

        lbl(gl, 24, 1, 'Overwrite SimNIBS');
        chk(gl, 24, 2, 'io.overwrite_simnibs', false, 'Re-run SimNIBS segmentation');
        sp(gl, 25);

        sec(gl, 26, 'Matrix Saving');
        note_lbl(gl, 27, 'Global flag is the fallback; individual flags take precedence when set.');
        lbl(gl, 28, 1, 'Global save flag');
        chk(gl, 28, 2, 'io.save_matrices', false, 'Global fallback: save all intermediate matrices');

        lbl(gl, 29, 1, 'Save source matrices');
        chk(gl, 29, 2, 'io.save_source_matrices', true, 'kwave_source.mat (shared across variants)');

        lbl(gl, 30, 1, 'Save acoustic matrices');
        chk(gl, 30, 2, 'io.save_acoustic_matrices', false, 'sensor_data + kgrid + medium (large)');

        lbl(gl, 31, 1, 'Save thermal matrices');
        chk(gl, 31, 2, 'io.save_thermal_matrices', true, 'heating_res.mat — temperature timeseries');

        lbl(gl, 32, 1, 'Save heating video');
        chk(gl, 32, 2, 'io.save_heatingvideo', false, 'Save MP4 of incremental heating');
        sp(gl, 33);

        sec(gl, 34, 'Advanced I/O');
        lbl(gl, 35, 1, 'External acoustic NIfTI');
        edt(gl, 35, 2, 'io.external_acoustic_nifti', '', 0);
        brw(gl, 35, 4, 'io.external_acoustic_nifti', 'file');
        note_lbl(gl, 36, 'When set, skips acoustic simulation. NIfTI must contain p-max [Pa] in T1 space (e.g. from BabelBrain or a prior run).');

        lbl(gl, 37, 1, 'Save complex pressure');
        chk(gl, 37, 2, 'io.save_p_complex', false, 'Save magnitude + phase of steady-state pressure field');

        lbl(gl, 38, 1, 'Save property maps');
        chk(gl, 38, 2, 'io.save_property_maps', false, 'Save per-tissue acoustic/thermal property NIfTIs to nii/properties/');

        lbl(gl, 39, 1, 'Save MNI outputs');
        chk(gl, 39, 2, 'io.save_MNI', true, 'Save simulation outputs in MNI space in addition to native T1w');
    end

    %% ── Tab 2: Simulation ─────────────────────────────────────────────
    function build_tab_simulation(t)
        gl = tab_grid(t, 52, {200, '1x', 55, 80});

        sec(gl, 1, 'Simulation Type');
        lbl(gl, 2, 1, '* Medium');
        drp(gl, 2, 2, 'simulation.medium', {'layered','water','phantom'}, 'layered');

        lbl(gl, 3, 1, '* Code backend');
        drp(gl, 3, 2, 'simulation.code_type', ...
            {'matlab_gpu','matlab_cpu','cpp_gpu','cpp_cpu'}, 'matlab_gpu');

        lbl(gl, 4, 1, 'Precision');
        drp(gl, 4, 2, 'simulation.precision', {'single','double'}, 'single');

        lbl(gl, 5, 1, 'Platform');
        drp(gl, 5, 2, 'platform', {'auto','matlab','slurm','qsub'}, 'auto');

        lbl(gl, 6, 1, 'Interactive mode');
        chk(gl, 6, 2, 'simulation.interactive', false, ...
            'Ask for user input and plot evolving figures (desktop only)');

        lbl(gl, 7, 1, 'Debug mode');
        chk(gl, 7, 2, 'simulation.debug', false, ...
            'Enable additional intermediate diagnostic outputs');

        lbl(gl, 8, 1, 'Uncertainty mode');
        chk(gl, 8, 2, 'simulation.uncertainty', false, ...
            'Run default / liberal / conservative variants + combined report');
        sp(gl, 9);

        sec(gl, 10, 'Pipeline Modules');
        note = uilabel(gl, 'Text', ...
            '  Toggle individual pipeline stages. Disable stages only if you know what you are doing.', ...
            'FontName', st.font, 'FontSize', st.fs_sm, 'FontColor', st.text_sub);
        note.Layout.Row    = 11;
        note.Layout.Column = [1 4];

        modules = { ...
            'run_grid_setup',          'Grid setup & head preprocessing',     true;  ...
            'run_medium_setup',        'Medium acoustic property mapping',     true;  ...
            'run_source_setup',        'Acoustic source setup',                true;  ...
            'run_acoustic_sims',       'Acoustic simulation (k-Wave)',         true;  ...
            'run_acoustic_analysis',   'Acoustic analysis (ISPPA, MI, etc.)', true;  ...
            'run_heating_sims',        'Thermal simulation',                   false; ...
            'run_thermal_analysis',    'Thermal analysis (CEM43, maxT)',       true;  ...
            'run_nifti_creation',      'NIfTI export',                         true;  ...
            'run_posthoc_water_sims',  'Post-hoc free-water simulation',       true;  ...
            'generate_report',         'Generate HTML report',                 true   ...
        };

        col_span = [1 2];
        for i = 1:size(modules,1)
            mtag = ['modules.' modules{i,1}];
            chk(gl, 11+i, col_span, mtag, modules{i,3}, modules{i,2});
        end

        chk(gl, 11+size(modules,1)+1, col_span, 'modules.segmentation_only', false, ...
            'Segmentation only — stop after SimNIBS (skip grid and simulations)');
        sp(gl, 23);

        sec(gl, 24, 'Multi-ISPPA Intensity Sweep');
        lbl(gl, 25, 1, 'Target Isppa (W/cm²)');
        edt(gl, 25, 2, 'transducer.target_isppa_wcm2', '', 0);
        note_lbl(gl, 26, 'One value for a single run; comma-separated list triggers intensity sweep (separate thermal sim per value, shared acoustic sim).');
        sp(gl, 27);

        sec(gl, 28, 'Sequential Simulations');
        chk(gl, 29, [1 2], 'options.sequential.enabled', false, ...
            'Chain follow-up sonications (thermal state inherited between runs)');
        note_lbl(gl, 30, 'Enter one YAML config path per line. Each follow-up inherits the temperature and CEM43 maps from the preceding run.');

        lbl(gl, 31, 1, 'Follow-up configs');
        seq_ta = uitextarea(gl, ...
            'Value',       {''}, ...
            'FontName',    st.font, 'FontSize', st.fs_sm, ...
            'Placeholder', 'path/to/followup1.yaml', ...
            'Tag',         'sequential_configs_textarea');
        seq_ta.Layout.Row = [32 38]; seq_ta.Layout.Column = [1 3];

        btn_seq_brw = uibutton(gl, 'Text', '+ Add', ...
            'FontName', st.font, 'FontSize', st.fs_sm, ...
            'BackgroundColor', [0.88 0.91 0.96], ...
            'ButtonPushedFcn', @(~,~) cb_seq_browse_append());
        btn_seq_brw.Layout.Row = 32; btn_seq_brw.Layout.Column = 4;
    end

    %% ── Tab 3: Transducer ─────────────────────────────────────────────
    function build_tab_transducer(t)
        gl = tab_grid(t, 35, {200, '1x', 55, 80});

        % ── Library-based selection (rows 1-4) ────────────────────────
        sec(gl, 1, 'Transducer Library');
        lbl(gl, 2, 1, 'Serial');

        % Populate serial list from equipment config; fall back to empty.
        serial_items = {'(manual)'};
        try
            eq = load_equipment_config();
            serial_items = [{'(manual)'}; fieldnames(eq.trans)];
        catch
        end
        h_serial = uidropdown(gl, ...
            'Items',    serial_items, ...
            'Value',    serial_items{1}, ...
            'Tag',      'transducer.serial', ...
            'FontName', st.font, 'FontSize', st.fs);
        h_serial.Layout.Row = 2; h_serial.Layout.Column = [2 3];
        h_serial.ValueChangedFcn = @(dd,~) cb_serial_changed(dd);

        lbl(gl, 3, 1, 'Driving system serial');
        edt(gl, 3, 2, 'transducer.combo.ds_serial', '', 0);
        note_lbl(gl, 4, 'Optional. Omit for a generic (DS-agnostic) calibration. Specify to use a driving-system-specific calibration.');

        h_lib_cov = uilabel(gl, ...
            'Text',      'Library: —', ...
            'Tag',       'lbl_library_coverage', ...
            'FontName',  st.font, ...
            'FontSize',  st.fs_sm, ...
            'FontColor', [0.5 0.5 0.5], ...
            'WordWrap',  'on');
        h_lib_cov.Layout.Row = 5; h_lib_cov.Layout.Column = [1 3];

        sp(gl, 6);

        % ── General (rows 7-12) ───────────────────────────────────────
        sec(gl, 7, 'General');
        lbl(gl, 8, 1, '* Type');
        h_type = drp(gl, 8, 2, 'transducer.type', {'annular','matrix'}, 'annular');
        h_type.ValueChangedFcn = @(dd,~) cb_transducer_type(dd);

        lbl(gl, 9, 1, '* Frequency');
        nedt(gl, 9, 2, 'transducer.freq_hz', 500000);
        lbl(gl, 9, 3, 'Hz');

        lbl(gl, 10, 1, 'Focal distance (exit plane)');
        nedt(gl, 10, 2, 'transducer.focal_distance_ep', NaN);
        lbl(gl, 10, 3, 'mm');

        lbl(gl, 11, 1, 'Focal distance (bowl)');
        nedt(gl, 11, 2, 'transducer.focal_distance_bowl', NaN);
        lbl(gl, 11, 3, 'mm');
        sp(gl, 12);

        sec(gl, 13, 'Transducer Position (T1 voxels)');
        lbl(gl, 14, 1, 'Transducer position');
        xyz_panel(gl, 14, 'transducer.trans_pos', [NaN NaN NaN]);

        lbl(gl, 15, 1, 'Focus position');
        xyz_panel(gl, 15, 'transducer.focus_pos', [NaN NaN NaN]);
        sp(gl, 16);

        % ── Annular panel (rows 17-30) ────────────────────────────────
        pnl_ann = uipanel(gl, ...
            'Title',           '', ...
            'BackgroundColor', st.bg_panel, ...
            'BorderType',      'none', ...
            'Tag',             'panel_annular');
        pnl_ann.Layout.Row    = [17 30];
        pnl_ann.Layout.Column = [1 4];

        gl_ann = uigridlayout(pnl_ann, ...
            'RowHeight',    repmat({32},1,14), ...
            'ColumnWidth',  {200,'1x',55,80}, ...
            'Padding',      [8 4 8 4], ...
            'RowSpacing',   4, ...
            'BackgroundColor', st.bg_panel);

        sec(gl_ann, 1, 'Annular Array');
        lbl(gl_ann, 2, 1, '* Elements');
        nedt(gl_ann, 2, 2, 'transducer.annular.elem_n', 4);

        lbl(gl_ann, 3, 1, '* Inner diameters (comma-sep)');
        edt(gl_ann, 3, 2, 'transducer.annular.elem_id_mm', '0,34,53,70', 0);
        lbl(gl_ann, 3, 3, 'mm');

        lbl(gl_ann, 4, 1, '* Outer diameters (comma-sep)');
        edt(gl_ann, 4, 2, 'transducer.annular.elem_od_mm', '32,51,67,80', 0);
        lbl(gl_ann, 4, 3, 'mm');

        lbl(gl_ann, 5, 1, '* Radius of curvature');
        nedt(gl_ann, 5, 2, 'transducer.annular.curv_radius_mm', 63.2);
        lbl(gl_ann, 5, 3, 'mm');

        lbl(gl_ann, 6, 1, '* Pressure amplitude');
        nedt(gl_ann, 6, 2, 'transducer.annular.elem_amp', 1);
        lbl(gl_ann, 6, 3, 'Pa');

        lbl(gl_ann, 7, 1, 'Phase per element (°)');
        edt(gl_ann, 7, 2, 'transducer.annular.elem_phase_deg', '0', 0);
        lbl(gl_ann, 7, 3, 'deg');

        lbl(gl_ann, 8, 1, 'Geometric focus-to-EP dist.');
        nedt(gl_ann, 8, 2, 'transducer.annular.dist_geom_ep_mm', NaN);
        lbl(gl_ann, 8, 3, 'mm');

        lbl(gl_ann, 9, 1, 'Visualization depth');
        nedt(gl_ann, 9, 2, 'transducer.annular.depth_mm', 16);
        lbl(gl_ann, 9, 3, 'mm');

        % ── Matrix panel (rows 17-30, hidden by default) ──────────────
        pnl_mat = uipanel(gl, ...
            'Title',           '', ...
            'BackgroundColor', st.bg_panel, ...
            'BorderType',      'none', ...
            'Tag',             'panel_matrix', ...
            'Visible',         'off');
        pnl_mat.Layout.Row    = [17 30];
        pnl_mat.Layout.Column = [1 4];

        gl_mat = uigridlayout(pnl_mat, ...
            'RowHeight',   repmat({32},1,14), ...
            'ColumnWidth', {200,'1x',55,80}, ...
            'Padding',     [8 4 8 4], ...
            'RowSpacing',  4, ...
            'BackgroundColor', st.bg_panel);

        sec(gl_mat, 1, 'Matrix Array');
        lbl(gl_mat, 2, 1, 'Element shape');
        drp(gl_mat, 2, 2, 'transducer.matrix.elem_shape', {'rect','disc','bowl'}, 'rect');

        lbl(gl_mat, 3, 1, 'Element height');
        nedt(gl_mat, 3, 2, 'transducer.matrix.elem_height_mm', 1);
        lbl(gl_mat, 3, 3, 'mm');

        lbl(gl_mat, 4, 1, 'Element width');
        nedt(gl_mat, 4, 2, 'transducer.matrix.elem_width_mm', 1);
        lbl(gl_mat, 4, 3, 'mm');

        lbl(gl_mat, 5, 1, 'Outer diameter');
        nedt(gl_mat, 5, 2, 'transducer.matrix.outer_diameter_mm', 70);
        lbl(gl_mat, 5, 3, 'mm');

        lbl(gl_mat, 6, 1, 'Pressure amplitude');
        nedt(gl_mat, 6, 2, 'transducer.matrix.elem_amp', 1);
        lbl(gl_mat, 6, 3, 'Pa');

        lbl(gl_mat, 7, 1, 'Curved surface');
        chk(gl_mat, 7, 2, 'transducer.matrix.is_curved', false, 'Curved aperture');

        lbl(gl_mat, 8, 1, 'Radius of curvature');
        nedt(gl_mat, 8, 2, 'transducer.matrix.curv_radius_mm', NaN);
        lbl(gl_mat, 8, 3, 'mm');

        lbl(gl_mat, 9, 1, 'Dist. geom. focus-to-EP');
        nedt(gl_mat, 9, 2, 'transducer.matrix.dist_geom_ep_mm', NaN);
        lbl(gl_mat, 9, 3, 'mm');

        lbl(gl_mat, 10, 1, 'Steering');
        drp(gl_mat, 10, 2, 'transducer.matrix.steering', {'1D','3D'}, '3D');

        lbl(gl_mat, 11, 1, 'Visualization depth');
        nedt(gl_mat, 11, 2, 'transducer.matrix.depth_mm', 16);
        lbl(gl_mat, 11, 3, 'mm');

        lbl(gl_mat, 12, 1, 'Clover multi-aperture');
        chk(gl_mat, 12, 2, 'transducer.matrix.is_clover_setup', false, 'Enable Clover replication');
    end

    %% ── Tab 4: Grid ───────────────────────────────────────────────────
    function build_tab_grid(t)
        gl = tab_grid(t, 29, {200, '1x', 55, 80});

        sec(gl, 1, 'Spatial Resolution');
        lbl(gl, 2, 1, '* Grid resolution');
        nedt(gl, 2, 2, 'grid.resolution_mm', 0.5);
        lbl(gl, 2, 3, 'mm');

        lbl(gl, 3, 1, 'Minimum PPW');
        nedt(gl, 3, 2, 'grid.min_ppw', 6);
        lbl(gl, 3, 3, '');
        note_lbl(gl, 4, 'Minimum points per wavelength at transducer frequency. Warning raised if violated.');
        sp(gl, 5);

        sec(gl, 6, 'Temporal Resolution');
        lbl(gl, 7, 1, 'CFL number');
        nedt(gl, 7, 2, 'grid.source_cfl', 0.15);
        note_lbl(gl, 8, 'Courant-Friedrichs-Lewy fraction. Default 0.15 (half of k-Wave default) for additional stability in heterogeneous skull.');

        lbl(gl, 9, 1, 'PPW override');
        nedt(gl, 9, 2, 'grid.source_ppw', NaN);
        note_lbl(gl, 10, 'Leave NaN to compute PPW from resolution and max sound speed.');

        lbl(gl, 11, 1, 'Stability limit fraction');
        nedt(gl, 11, 2, 'grid.source_limit_fraction', 0.9);
        note_lbl(gl, 12, 'If dt exceeds stability limit, rescale to this fraction of the limit. Set 0 to disable.');
        sp(gl, 13);

        sec(gl, 14, 'Grid Dimensions');
        lbl(gl, 15, 1, 'Default dims (water/phantom)');
        xyz_panel(gl, 15, 'grid.default_dims', [144 144 400]);

        lbl(gl, 16, 1, 'PML size');
        nedt(gl, 16, 2, 'grid.pml_size', 10);
        lbl(gl, 16, 3, 'voxels');

        lbl(gl, 17, 1, 'Max grid expand');
        nedt(gl, 17, 2, 'grid.max_expand', 40);
        lbl(gl, 17, 3, 'voxels');

        lbl(gl, 18, 1, 'Axisymmetric mode');
        chk(gl, 18, 2, 'grid.axisymmetric', false, '2D axisymmetric (kspaceFirstOrderAS)');

        lbl(gl, 19, 1, 'Use kWaveArray');
        chk(gl, 19, 2, 'grid.use_kWaveArray', true, 'Recommended for accurate transducer modelling');
        sp(gl, 20);

        sec(gl, 21, 'Grid Orientation & Thermal Grid');
        lbl(gl, 22, 1, 'Grid orientation mode');
        drp(gl, 22, 2, 'grid.mode', {'transducer_axis','ras_plus'}, 'transducer_axis');
        note_lbl(gl, 23, 'transducer_axis: rotate volume so focal axis aligns with z (compact grid). ras_plus: keep scanner RAS+ orientation without rotation.');

        lbl(gl, 24, 1, 'Thermal resolution');
        nedt(gl, 24, 2, 'grid.thermal_resolution_mm', NaN);
        lbl(gl, 24, 3, 'mm');
        note_lbl(gl, 25, 'Leave NaN to use the acoustic grid resolution for thermal simulation.');

        lbl(gl, 26, 1, 'Thermal FOV');
        nedt(gl, 26, 2, 'grid.thermal_fov_mm', NaN);
        lbl(gl, 26, 3, 'mm');
        note_lbl(gl, 27, 'Leave NaN to match the acoustic simulation field-of-view.');
    end

    %% ── Tab 5: Medium Properties ──────────────────────────────────────
    function build_tab_medium(t)
        gl = uigridlayout(t, ...
            'RowHeight',   {'fit','1x','fit'}, ...
            'ColumnWidth', {'1x'}, ...
            'Padding',     [16 16 16 16], ...
            'RowSpacing',  10, ...
            'BackgroundColor', st.bg_panel);

        uilabel(gl, ...
            'Text', 'Acoustic and thermal tissue properties. Red cells override PRESTUS defaults. Leave cells blank to use defaults.', ...
            'FontName', st.font, 'FontSize', st.fs_sm, 'FontColor', st.text_sub, ...
            'WordWrap', 'on');

        tissues = {'water','brain','skin','skull','skull_cortical','skull_trabecular'};
        props   = {'sound_speed','density','alpha_coeff','alpha_power', ...
                   'thermal_conductivity','specific_heat_capacity','perfusion','absorption_fraction'};
        col_names = {'Sound Speed [m/s]','Density [kg/m³]','α coeff [dB/cm/MHz]', ...
                     'α power','Therm. cond.','Heat cap.','Perfusion','Absorb. frac.'};

        % Default values (from PRESTUS defaults)
        defaults = { ...
            1500,   1000, 0.002,  2.0, 0.6,  4000, 0,  0.02;  ... % water
            1546,   1046, 0.6,    1.3, 0.51, 3630, 559, 0.02; ... % brain
            1537,   1116, 0.8,    1.3, 0.37, 3391, 106, 0.02; ... % skin
            2800,   1912, 11.6,   1.3, 0.32, 1313, 0,   0.02; ... % skull
            2800,   1912, 11.6,   1.3, 0.32, 1313, 0,   0.02; ... % skull_cortical
            2300,   1178, 6.9,    1.3, 0.32, 1313, 0,   0.02  ... % skull_trabecular
        };

        tbl = uitable(gl, ...
            'Data',                defaults, ...
            'RowName',             tissues, ...
            'ColumnName',          col_names, ...
            'ColumnEditable',      true(1,8), ...
            'FontName',            st.font, ...
            'FontSize',            st.fs_sm, ...
            'Tag',                 'medium_table', ...
            'ColumnWidth',         {115, 100, 140, 80, 100, 90, 100, 110});

        uilabel(gl, ...
            'Text', sprintf('Tip: layers can be excluded in the Advanced tab by removing them from the active layers list.'), ...
            'FontName', st.font, 'FontSize', st.fs_sm, 'FontColor', st.text_sub);
    end

    %% ── Tab 6: Thermal ────────────────────────────────────────────────
    function build_tab_thermal(t)
        gl = tab_grid(t, 30, {220, '1x', 55, 80});

        sec(gl, 1, 'Sonication Protocol Timing');
        lbl(gl, 2, 1, 'Pulse duration (PD)');
        nedt(gl, 2, 2, 'timing.pd', NaN); lbl(gl, 2, 3, 's');

        lbl(gl, 3, 1, 'Pulse repetition interval (PRI)');
        nedt(gl, 3, 2, 'timing.pri', NaN); lbl(gl, 3, 3, 's');

        lbl(gl, 4, 1, 'Pulse train duration (PTD)');
        nedt(gl, 4, 2, 'timing.ptd', NaN); lbl(gl, 4, 3, 's');

        lbl(gl, 5, 1, 'Pulse train rep. interval (PTRI)');
        nedt(gl, 5, 2, 'timing.ptri', NaN); lbl(gl, 5, 3, 's');

        lbl(gl, 6, 1, 'Pulse train rep. duration (PTRD)');
        nedt(gl, 6, 2, 'timing.ptrd', NaN); lbl(gl, 6, 3, 's');

        lbl(gl, 7, 1, 'Post-PTRI steady-state');
        nedt(gl, 7, 2, 'timing.post_ptri_dur', NaN); lbl(gl, 7, 3, 's');

        lbl(gl, 8, 1, 'In-train timestep');
        nedt(gl, 8, 2, 'timing.pt_timestep', 0.02); lbl(gl, 8, 3, 's');

        lbl(gl, 9, 1, 'Post-train timestep');
        nedt(gl, 9, 2, 'timing.post_pt_timestep', 1); lbl(gl, 9, 3, 's');

        lbl(gl, 10, 1, 'Equal on/off step durations');
        chk(gl, 10, 2, 'timing.equal_step_duration', false, ...
            'Force equal durations for on and off cycles');
        sp(gl, 11);

        sec(gl, 12, 'CEM43 Thermal Dose');
        lbl(gl, 13, 1, 'ISO CEM43');
        chk(gl, 13, 2, 'thermal.cem43_iso', false, ...
            'Use ISO definition (R=0.5 above 43°C; zero below 39°C; ∞ above 57°C)');
        note_lbl(gl, 14, 'Default uses k-Wave kWaveDiffusion built-in formula (R=0.25 for T<43, R=0.5 for T≥43, Sapareto-Dewey).');
        sp(gl, 15);

        sec(gl, 16, 'Sensor & Recording');
        lbl(gl, 17, 1, 'Sensor half-window size');
        nedt(gl, 17, 2, 'thermal.sensor_xy_halfsize', 100); lbl(gl, 17, 3, 'voxels');
        note_lbl(gl, 18, 'Max half-size of temperature recording window. Reduce if memory is limited.');

        lbl(gl, 19, 1, 'Record T at every step');
        chk(gl, 19, 2, 'thermal.record_t_at_every_step', false, ...
            'Record full sensor window at every timestep (memory intensive)');
        sp(gl, 20);

        sec(gl, 21, 'Initial Temperature (°C)');
        tissues_t = {'water','brain','skin','skull','skull_cortical','skull_trabecular'};
        for i = 1:numel(tissues_t)
            lbl(gl, 21+i, 1, tissues_t{i});
            nedt(gl, 21+i, 2, ['thermal.temp_0.' tissues_t{i}], 37);
            lbl(gl, 21+i, 3, '°C');
        end
    end

    %% ── Tab 7: HPC ────────────────────────────────────────────────────
    function build_tab_hpc(t)
        gl = tab_grid(t, 20, {200, '1x', 55, 80});

        sec(gl, 1, 'Scheduler');
        lbl(gl, 2, 1, 'HPC profile');
        drp(gl, 2, 2, 'hpc.name', {'default','snellius'}, 'default');

        lbl(gl, 3, 1, 'Partition');
        edt(gl, 3, 2, 'hpc.partition', '', 0);

        lbl(gl, 4, 1, 'GPU request');
        edt(gl, 4, 2, 'hpc.gpu', '', 0);
        note_lbl(gl, 5, 'e.g.  nvidia_a100-sxm4-40gb:1');

        lbl(gl, 6, 1, 'Reservation');
        edt(gl, 6, 2, 'hpc.reservation', '', 0);

        lbl(gl, 7, 1, 'Job prefix');
        edt(gl, 7, 2, 'hpc.job_prefix', 'PRESTUS', 0);
        sp(gl, 8);

        sec(gl, 9, 'Resources');
        lbl(gl, 10, 1, 'Wall time limit');
        edt(gl, 10, 2, 'hpc.timelimit', '04:00:00', 0);

        lbl(gl, 11, 1, 'Memory limit');
        nedt(gl, 11, 2, 'hpc.memorylimit', 20);
        lbl(gl, 11, 3, 'GB');
        sp(gl, 12);

        sec(gl, 13, 'Behaviour');
        lbl(gl, 14, 1, 'Wait for job');
        chk(gl, 14, 2, 'hpc.wait_for_job', false, 'Block until HPC job completes');

        lbl(gl, 15, 1, 'Max wait checks');
        nedt(gl, 15, 2, 'hpc.max_wait_checks', 540);
        note_lbl(gl, 16, 'At ~20 s per check; 540 ≈ 3 hours.');

        lbl(gl, 17, 1, 'LD_LIBRARY_PATH');
        edt(gl, 17, 2, 'hpc.ld_library_path', '', 0);
        note_lbl(gl, 18, 'Set if SimNIBS shows "undefined symbol" errors on HPC.');
    end

    %% ── Tab 8: Advanced ───────────────────────────────────────────────
    function build_tab_advanced(t)
        gl = tab_grid(t, 55, {220, '1x', 55, 80});

        sec(gl, 1, 'Segmentation');
        lbl(gl, 2, 1, 'Force qform reorientation');
        chk(gl, 2, 2, 'segmentation.use_qform', false, 'Fix qform/sform mismatch before charm');
        lbl(gl, 3, 1, 'Debug segmentation');
        chk(gl, 3, 2, 'segmentation.debug', false, 'Pass --debug flag to charm');
        sp(gl, 4);

        sec(gl, 5, 'Head Model Processing');
        lbl(gl, 6, 1, 'Head pad');
        nedt(gl, 6, 2, 'headmodel.head_pad_mm', 0); lbl(gl, 6, 3, 'mm');

        lbl(gl, 7, 1, 'CSF expansion');
        nedt(gl, 7, 2, 'headmodel.csf_expansion', 40); lbl(gl, 7, 3, 'voxels');

        lbl(gl, 8, 1, 'Smoothing method');
        drp(gl, 8, 2, 'headmodel.smooth_method', {'gaussian','box','off'}, 'gaussian');

        lbl(gl, 9, 1, 'Smoothing FWHM');
        nedt(gl, 9, 2, 'headmodel.smooth_fwhm_mm', 1); lbl(gl, 9, 3, 'mm');

        lbl(gl, 10, 1, 'Skull threshold');
        nedt(gl, 10, 2, 'headmodel.smooth_threshold_skull', 0.5);
        note_lbl(gl, 11, 'Higher = thinner skull mask after smoothing.');

        lbl(gl, 12, 1, 'Other tissue threshold');
        nedt(gl, 12, 2, 'headmodel.smooth_threshold_other', 0.5);

        lbl(gl, 13, 1, 'Smooth acoustic properties');
        chk(gl, 13, 2, 'headmodel.smooth_properties', false, ...
            'Apply smoothing kernel to acoustic property maps');

        lbl(gl, 14, 1, 'Skull fill method');
        drp(gl, 14, 2, 'headmodel.skull_fill_method', {'rubberwrap','imclose'}, 'rubberwrap');

        lbl(gl, 15, 1, 'Rubber-wrap radius');
        nedt(gl, 15, 2, 'headmodel.skull_wrap_radius', 10); lbl(gl, 15, 3, 'voxels');

        lbl(gl, 16, 1, 'Visualize rubber-wrap');
        chk(gl, 16, 2, 'headmodel.skull_wrap_visualize', false, ...
            'Show rubber-wrap result (disable on HPC)');
        sp(gl, 17);

        sec(gl, 18, 'Pseudo-CT Skull Mapping');
        lbl(gl, 19, 1, 'Enable pCT');
        h_pct_enabled = chk(gl, 19, 2, 'pct.enabled', false, 'Use CT/pseudo-CT to inform skull properties');

        lbl(gl, 20, 1, 'UTE→HU skull mapping');
        h_skull_map = drp(gl, 20, 2, 'pct.skull_mapping', ...
            {'kosciessa','miscouridou','carpino','wiesinger','treeby'}, 'kosciessa');

        lbl(gl, 21, 1, 'Debug (keep intermediates)');
        h_debug = chk(gl, 21, 2, 'pct.debug', true, ...
            'Keep intermediate NIfTI files in pseudoCT/ subfolder');

        lbl(gl, 22, 1, 'Density mapping');
        h_dens = drp(gl, 22, 2, 'pct.mapping_density', {'k-plan','k-wave','marsac','aubry','none'}, 'k-plan');

        lbl(gl, 23, 1, 'Sound speed mapping');
        h_ss = drp(gl, 23, 2, 'pct.mapping_soundspeed', {'k-plan','marsac','aubry','none'}, 'k-plan');

        lbl(gl, 24, 1, 'Attenuation mapping');
        h_att = drp(gl, 24, 2, 'pct.mapping_attenuation', {'k-plan','mueller','aubry','none'}, 'k-plan');
        sp(gl, 25);

        % Enable/disable pCT sub-controls based on the checkbox state
        pct_sub_controls = [h_skull_map, h_debug, h_dens, h_ss, h_att];
        cb_pct_toggle(h_pct_enabled, pct_sub_controls);  % set initial state
        h_pct_enabled.ValueChangedFcn = @(src,~) cb_pct_toggle(src, pct_sub_controls);

        sec(gl, 26, 'Analysis');
        lbl(gl, 27, 1, 'Focus area radius');
        nedt(gl, 27, 2, 'analysis.focus_area_radius', 5); lbl(gl, 27, 3, 'mm');
        note_lbl(gl, 28, 'Radius around focus for ISPPA averaging in output metrics.');
        sp(gl, 29);

        sp(gl, 30);
        note_lbl(gl, 31, 'Transducer placement settings (manual / Localite / heuristic / PlanTUS) are in the Placement tab.');
        sp(gl, 32);

        sec(gl, 33, 'Neuronavigation Ingestion');
        lbl(gl, 34, 1, 'Subject ID');
        h_sub = uieditfield(gl, 'text', 'Tag', 'neuronav.sub_id', 'Value', '', ...
            'FontName', st.font, 'FontSize', st.fs, 'Placeholder', 'sub-001');
        h_sub.Layout.Row = 34; h_sub.Layout.Column = 2;

        lbl(gl, 35, 1, 'Session ID');
        h_ses = uieditfield(gl, 'text', 'Tag', 'neuronav.ses_id', 'Value', '', ...
            'FontName', st.font, 'FontSize', st.fs, 'Placeholder', 'ses-01');
        h_ses.Layout.Row = 35; h_ses.Layout.Column = 2;

        note_lbl(gl, 36, 'Marker type (TriggerMarkers / GUMMarkers) is auto-detected from the folder.');

        lbl(gl, 37, 1, 'Raw Localite path');
        h_raw = uieditfield(gl, 'text', 'Tag', 'neuronav.raw_path', 'Value', '', ...
            'FontName', st.font, 'FontSize', st.fs);
        h_raw.Layout.Row = 37; h_raw.Layout.Column = 2;
        brw(gl, 37, 4, 'neuronav.raw_path', 'dir');

        lbl(gl, 38, 1, 'Output path');
        h_out = uieditfield(gl, 'text', 'Tag', 'neuronav.out_path', 'Value', '', ...
            'FontName', st.font, 'FontSize', st.fs);
        h_out.Layout.Row = 38; h_out.Layout.Column = 2;
        brw(gl, 38, 4, 'neuronav.out_path', 'dir');

        h_ingest_btn = uibutton(gl, 'Text', '▶  Run Ingestion', ...
            'Tag', 'btn_ingest', ...
            'FontName', st.font, 'FontSize', st.fs, 'FontWeight', 'bold', ...
            'BackgroundColor', st.accent, 'FontColor', [1 1 1], ...
            'ButtonPushedFcn', @(~,~) cb_ingest());
        h_ingest_btn.Layout.Row = 39; h_ingest_btn.Layout.Column = [1 2];

        h_ingest_status = uilabel(gl, ...
            'Text', '', 'Tag', 'lbl_ingest_status', ...
            'FontName', st.font, 'FontSize', st.fs_sm, ...
            'FontColor', st.text_sub, 'WordWrap', 'on', ...
            'HorizontalAlignment', 'left');
        h_ingest_status.Layout.Row = 40; h_ingest_status.Layout.Column = [1 4];
        sp(gl, 41);

        sec(gl, 42, 'Simulation Layers');
        note_lbl(gl, 43, 'Comma-separated SimNIBS label indices assigned to each tissue compartment.');

        layer_names = {'water','brain','skin','skull','skull_cortical','skull_trabecular'};
        layer_defaults = {'0,3,6,9,10', '1,2', '5', '4', '7', '8'};
        for i = 1:numel(layer_names)
            lbl(gl, 43+i, 1, layer_names{i});
            edt(gl, 43+i, 2, ['layers.' layer_names{i}], layer_defaults{i}, 0);
        end
    end

    %% ── Tab 9: Placement ─────────────────────────────────────────────
    function build_tab_placement(t)
        gl = tab_grid(t, 60, {220, '1x', 55, 80});

        sec(gl, 1, 'Transducer Placement');
        lbl(gl, 2, 1, 'Placement mode');
        h_placement_mode = drp(gl, 2, 2, 'placement.mode', ...
            {'manual','localite','heuristic','plantus'}, 'manual');
        note_lbl(gl, 3, 'manual: use trans_pos/focus_pos as-is.  localite: read from XML.  heuristic: sphere-expansion search.  plantus: multi-objective optimisation.');
        sp(gl, 4);

        % ── Localite sub-panel ────────────────────────────────────────────
        pnl_loc = uipanel(gl, 'Title', '', 'BackgroundColor', st.bg_panel, ...
            'BorderType', 'none', 'Tag', 'panel_placement_localite', 'Visible', 'off');
        pnl_loc.Layout.Row = [5 16]; pnl_loc.Layout.Column = [1 4];
        gl_loc = uigridlayout(pnl_loc, 'RowHeight', repmat({26},1,12), ...
            'ColumnWidth', {220,'1x',55,80}, 'Padding', [0 2 0 2], ...
            'RowSpacing', 4, 'BackgroundColor', st.bg_panel);

        sec(gl_loc, 1, 'Localite Config');
        lbl(gl_loc, 2, 1, 'Localite XML file');
        edt(gl_loc, 2, 2, 'placement.localite.file', '', 0);
        brw(gl_loc, 2, 4, 'placement.localite.file', 'file');
        lbl(gl_loc, 3, 1, 'Session index');
        nedt(gl_loc, 3, 2, 'placement.localite.session', 1);
        lbl(gl_loc, 4, 1, 'Marker type');
        drp(gl_loc, 4, 2, 'placement.localite.markertype', ...
            {'TriggerMarkers','NavigationMarkers'}, 'TriggerMarkers');
        lbl(gl_loc, 5, 1, 'Marker position index');
        nedt(gl_loc, 5, 2, 'placement.localite.position', 1);
        lbl(gl_loc, 6, 1, 'Reference distance');
        nedt(gl_loc, 6, 2, 'placement.localite.reference_distance_mm', 15);
        lbl(gl_loc, 6, 3, 'mm');
        note_lbl(gl_loc, 7, 'Distance between IR trackers and transducer exit plane.');
        lbl(gl_loc, 8, 1, 'Save Localite-aligned T1');
        chk(gl_loc, 8, 2, 'placement.localite.save_localite_t1', false, ...
            'Save T1 aligned to Localite header for QC');

        % ── Heuristic sub-panel ───────────────────────────────────────────
        pnl_heu = uipanel(gl, 'Title', '', 'BackgroundColor', st.bg_panel, ...
            'BorderType', 'none', 'Tag', 'panel_placement_heuristic', 'Visible', 'off');
        pnl_heu.Layout.Row = [5 16]; pnl_heu.Layout.Column = [1 4];
        gl_heu = uigridlayout(pnl_heu, 'RowHeight', repmat({26},1,12), ...
            'ColumnWidth', {220,'1x',55,80}, 'Padding', [0 2 0 2], ...
            'RowSpacing', 4, 'BackgroundColor', st.bg_panel);

        sec(gl_heu, 1, 'Heuristic (Sphere-Expansion) Config');
        lbl(gl_heu, 2, 1, 'MNI target (mm)');
        xyz_panel(gl_heu, 2, 'placement.heuristic.mni_target_mm', [NaN NaN NaN]);
        lbl(gl_heu, 3, 1, 'Target name');
        edt(gl_heu, 3, 2, 'placement.heuristic.target_name', '', 0);
        lbl(gl_heu, 4, 1, 'Dist. close');
        nedt(gl_heu, 4, 2, 'placement.heuristic.dist_close', NaN);
        lbl(gl_heu, 4, 3, 'mm');
        lbl(gl_heu, 5, 1, 'Ear exclusion radius');
        nedt(gl_heu, 5, 2, 'placement.heuristic.ear_radius', 35);
        lbl(gl_heu, 5, 3, 'mm');
        lbl(gl_heu, 6, 1, 'Save Localite-aligned T1');
        chk(gl_heu, 6, 2, 'placement.heuristic.save_localite_t1', false, ...
            'Save T1 aligned to heuristic result for QC');

        % ── PlanTUS sub-panel ─────────────────────────────────────────────
        pnl_ptu = uipanel(gl, 'Title', '', 'BackgroundColor', st.bg_panel, ...
            'BorderType', 'none', 'Tag', 'panel_placement_plantus', 'Visible', 'off');
        pnl_ptu.Layout.Row = [5 57]; pnl_ptu.Layout.Column = [1 4];
        gl_ptu = uigridlayout(pnl_ptu, 'RowHeight', repmat({26},1,26), ...
            'ColumnWidth', {220,'1x',55,80}, 'Padding', [0 2 0 2], ...
            'RowSpacing', 4, 'BackgroundColor', st.bg_panel, 'Scrollable', 'on');

        sec(gl_ptu, 1, 'PlanTUS Config');
        lbl(gl_ptu, 2, 1, '* PlanTUS script path');
        edt(gl_ptu, 2, 2, 'placement.plantus.script_path', '', 1);
        brw(gl_ptu, 2, 4, 'placement.plantus.script_path', 'dir');
        lbl(gl_ptu, 3, 1, '* SimNIBS env path');
        edt(gl_ptu, 3, 2, 'placement.plantus.env_path', '', 1);
        brw(gl_ptu, 3, 4, 'placement.plantus.env_path', 'dir');
        lbl(gl_ptu, 4, 1, '* MNI target (mm)');
        xyz_panel(gl_ptu, 4, 'placement.plantus.mni_target_mm', [NaN NaN NaN]);
        lbl(gl_ptu, 5, 1, '* Target name');
        edt(gl_ptu, 5, 2, 'placement.plantus.target_name', '', 1);
        lbl(gl_ptu, 6, 1, '* Focal distance list');
        edt(gl_ptu, 6, 2, 'placement.plantus.focal_distance_list', '', 1);
        lbl(gl_ptu, 6, 3, 'mm');
        note_lbl(gl_ptu, 7, 'Comma-separated calibration focal distances from the exit plane.');
        lbl(gl_ptu, 8, 1, '* FLHM list');
        edt(gl_ptu, 8, 2, 'placement.plantus.flhm_list', '', 1);
        lbl(gl_ptu, 8, 3, 'mm');
        note_lbl(gl_ptu, 9, 'Comma-separated full-width half-max focal lengths matching the focal distance list.');
        lbl(gl_ptu, 10, 1, 'Max steering angle');
        nedt(gl_ptu, 10, 2, 'placement.plantus.max_angle_deg', 10);
        lbl(gl_ptu, 10, 3, 'deg');
        lbl(gl_ptu, 11, 1, 'Additional offset');
        nedt(gl_ptu, 11, 2, 'placement.plantus.additional_offset_mm', 0);
        lbl(gl_ptu, 11, 3, 'mm');
        lbl(gl_ptu, 12, 1, 'Mask radius');
        nedt(gl_ptu, 12, 2, 'placement.plantus.mask_radius_mm', 2);
        lbl(gl_ptu, 12, 3, 'mm');
        lbl(gl_ptu, 13, 1, 'Connectome WB path');
        edt(gl_ptu, 13, 2, 'placement.plantus.connectome_wb_path', '', 0);
        brw(gl_ptu, 13, 4, 'placement.plantus.connectome_wb_path', 'dir');
        note_lbl(gl_ptu, 14, 'Optional. Required for interactive GUI visualisation.');
        sp(gl_ptu, 15);
        sec(gl_ptu, 16, 'Optimisation Weights (must sum to 1.0)');
        lbl(gl_ptu, 17, 1, 'Skin-target distance');
        nedt(gl_ptu, 17, 2, 'placement.plantus.weights.skin_target_distances', 0.2);
        lbl(gl_ptu, 18, 1, 'Skin-target angle');
        nedt(gl_ptu, 18, 2, 'placement.plantus.weights.skin_target_angles', 0.2);
        lbl(gl_ptu, 19, 1, 'Skin-target intersection');
        nedt(gl_ptu, 19, 2, 'placement.plantus.weights.skin_target_intersections', 0.2);
        lbl(gl_ptu, 20, 1, 'Skin-skull angle');
        nedt(gl_ptu, 20, 2, 'placement.plantus.weights.skin_skull_angles', 0.2);
        lbl(gl_ptu, 21, 1, 'Skull thickness');
        nedt(gl_ptu, 21, 2, 'placement.plantus.weights.skull_thickness', 0.2);

        % Wire placement mode dropdown to show/hide sub-panels
        h_placement_mode.ValueChangedFcn = @(dd,~) cb_placement_mode(dd.Value);
        cb_placement_mode('manual');  % apply initial state
    end

    %% ── Tab 10: Multi-Transducer ──────────────────────────────────────
    function build_tab_multitransducer(t)
        gl = tab_grid(t, 22, {220, '1x', 55, 80});

        sec(gl, 1, 'Multi-Transducer Mode');
        note_lbl(gl, 2, 'Primary transducer is defined in the Transducer tab. Secondary transducers are loaded from YAML config files listed below.');

        lbl(gl, 3, 1, 'Coupling mode');
        drp(gl, 3, 2, 'simulation.transducer_coupling', {'coherent','async'}, 'coherent');
        note_lbl(gl, 4, 'coherent: all transducers simulated in one k-Wave run (interference modelled).  async: each simulated independently, intensities summed incoherently.');
        sp(gl, 5);

        sec(gl, 6, 'Secondary Transducer Configs');
        note_lbl(gl, 7, 'Each YAML fully defines one additional transducer (type, frequency, position, focus, intensity). They are loaded and appended to parameters.transducer(2..N).');

        lbl(gl, 8, 1, 'Config paths');
        mt_ta = uitextarea(gl, ...
            'Value',       {''}, ...
            'FontName',    st.font, 'FontSize', st.fs_sm, ...
            'Placeholder', 'path/to/transducer2.yaml', ...
            'Tag',         'multitrans_configs_textarea');
        mt_ta.Layout.Row = [9 16]; mt_ta.Layout.Column = [1 3];

        btn_mt_brw = uibutton(gl, 'Text', '+ Add', ...
            'FontName', st.font, 'FontSize', st.fs_sm, ...
            'BackgroundColor', [0.88 0.91 0.96], ...
            'ButtonPushedFcn', @(~,~) cb_multitrans_browse_append());
        btn_mt_brw.Layout.Row = 9; btn_mt_brw.Layout.Column = 4;

        note_lbl(gl, 17, 'One YAML path per line. Leave empty for single-transducer mode.');
    end

    %% ── Tab 10: Calibration ───────────────────────────────────────────
    function build_tab_calibration(t)
        gl = tab_grid(t, 52, {200, '1x', 55, 80});

        sec(gl, 1, 'Input Data');
        lbl(gl, 2, 1, 'Equipment name');
        edt(gl, 2, 2, 'calibration.equipment_name', '', 0);
        note_lbl(gl, 3, 'Identifier matched against characterisation CSV files (e.g. CTX500).');

        lbl(gl, 4, 1, 'Axial profile folder');
        edt(gl, 4, 2, 'calibration.path_input_axial', '', 0);
        brw(gl, 4, 3, 'calibration.path_input_axial', 'dir');
        note_lbl(gl, 5, 'Directory containing per-equipment axial hydrophone CSVs.');

        lbl(gl, 6, 1, 'Phase table folder');
        edt(gl, 6, 2, 'calibration.path_input_phase', '', 0);
        brw(gl, 6, 3, 'calibration.path_input_phase', 'dir');
        note_lbl(gl, 7, 'Directory with manufacturer phase tables (Sonic Concepts / Imasonic).');
        sp(gl, 8);

        sec(gl, 9, 'Calibration Targets');
        lbl(gl, 10, 1, 'Focal depths (ep)');
        edt(gl, 10, 2, 'calibration.focal_depths_wrt_exit_plane', '', 0);
        lbl(gl, 10, 3, 'mm');
        note_lbl(gl, 11, 'Comma-separated target focal distances from the exit plane.');

        lbl(gl, 12, 1, 'Desired intensities');
        edt(gl, 12, 2, 'calibration.desired_intensities', '', 0);
        lbl(gl, 12, 3, 'W/cm²');
        note_lbl(gl, 13, 'Comma-separated target Isppa values; one calibration run per value.');

        lbl(gl, 14, 1, 'Add focal-distance offset');
        chk(gl, 14, 2, 'calibration.add_FDO', false, 'Add bowl-offset to focal distance');
        note_lbl(gl, 15, 'Shifts focal target by (curvature radius − exit-plane distance).');
        sp(gl, 16);

        sec(gl, 17, 'Optimisation Settings');
        lbl(gl, 18, 1, 'Calibration mode');
        drp(gl, 18, 2, 'calibration.Mode', {'single_ref','multi_depth'}, 'single_ref');
        note_lbl(gl, 19, '"single_ref" optimises at one reference depth and extrapolates; "multi_depth" fits all depths jointly (BabelBrain-style, slower).');

        lbl(gl, 20, 1, 'Forward model');
        drp(gl, 20, 2, 'calibration.ForwardModel', {'oneil','rayleigh'}, 'oneil');
        note_lbl(gl, 21, '"oneil" is faster; "rayleigh" matches BabelBrain and is more accurate for ring gaps. Default for multi_depth is rayleigh.');

        lbl(gl, 22, 1, 'Reference depth');
        nedt(gl, 22, 2, 'calibration.RefDepth', NaN);
        lbl(gl, 22, 3, 'mm');
        note_lbl(gl, 23, 'Reference depth for single_ref and for packaging output phases. Leave blank for median of available depths.');

        lbl(gl, 24, 1, 'Initial velocity');
        nedt(gl, 24, 2, 'calibration.initial_velocity', 0.05);
        lbl(gl, 24, 3, 'm/s');
        note_lbl(gl, 25, 'Starting guess for the element normal velocity during fitting.');

        lbl(gl, 26, 1, 'Upper velocity');
        nedt(gl, 26, 2, 'calibration.opt_upper_velocity', 0.2);
        lbl(gl, 26, 3, 'm/s');
        note_lbl(gl, 27, 'Upper bound on velocity during optimisation.');

        lbl(gl, 28, 1, 'Phase precession');
        drp(gl, 28, 2, 'calibration.opt_phase_precession', ...
            {'none','linear','monotonic'}, 'none');
        note_lbl(gl, 29, '"none" = unconstrained per-element phases (required for single_ref geo-correction); "linear"/"monotonic" add phase ramp constraints.');

        lbl(gl, 30, 1, 'Regularisation λ');
        nedt(gl, 30, 2, 'calibration.opt_regularization_lambda', 0);
        note_lbl(gl, 31, 'L2 penalty weight on hardware correction. 0 = none. Recommended for both modes to keep optimised phases near geometric steering (typical: 1e-4–1e-2). Larger values reduce the correction magnitude.');

        lbl(gl, 32, 1, 'Opt. limits');
        edt(gl, 32, 2, 'calibration.opt_limits', '', 0);
        lbl(gl, 32, 3, 'mm');
        note_lbl(gl, 33, 'Two-element [min max] range (from bowl) for profile error computation. Leave blank to use full non-NaN range.');

        lbl(gl, 34, 1, 'Profile weights');
        nedt(gl, 34, 2, 'calibration.opt_weights', 0);
        note_lbl(gl, 35, '0 = uniform weighting; ≥1 = Gaussian centred on FLHM (higher = narrower peak emphasis).');

        lbl(gl, 36, 1, 'Opt. method');
        drp(gl, 36, 2, 'calibration.opt_method', {'FEXminimize','GlobalSearch'}, 'FEXminimize');
        note_lbl(gl, 37, '"FEXminimize" (default, no toolbox required); "GlobalSearch" requires the Global Optimisation Toolbox.');

        lbl(gl, 38, 1, 'Random seed');
        nedt(gl, 38, 2, 'calibration.opt_seed', NaN);
        note_lbl(gl, 39, 'Integer seed for reproducible optimisation. Leave blank for a non-deterministic run.');
        sp(gl, 40);

        sec(gl, 41, 'Simulation Control');
        lbl(gl, 42, 1, 'Run free-water sim');
        chk(gl, 42, 2, 'calibration.run_free_water_sim', true, 'Simulate correction in free water');

        lbl(gl, 43, 1, 'Force kWaveArray');
        chk(gl, 43, 2, 'calibration.force_kwavearray', false, 'Use kWaveArray transducer backend');

        lbl(gl, 44, 1, 'Axisymmetric 2D');
        chk(gl, 44, 2, 'calibration.axisymmetric2D', false, 'Use 2-D axisymmetric grid (annular only)');

        lbl(gl, 45, 1, 'Amplitude validation');
        drp(gl, 45, 2, 'calibration.opt_amp_validation', ...
            {'always','initial','final','none'}, 'final');
        note_lbl(gl, 46, 'When to run a validation simulation after amplitude calibration.');
        sp(gl, 47);

        sec(gl, 48, 'Per-Element Phase Corrections');
        lbl(gl, 49, 1, 'Correction offsets');
        edt(gl, 49, 2, 'calibration.elem_phase_correction_deg', '', 0);
        lbl(gl, 49, 3, 'deg');
        note_lbl(gl, 50, 'Comma-separated hardware phase offsets per element (optional; loaded from transducer YAML if present).');

        lbl(gl, 51, 1, 'Save recovered offsets');
        chk(gl, 51, 2, 'calibration.save_elem_correction', false, 'Write recovered corrections back to transducer YAML');
        sp(gl, 52);

        sec(gl, 53, 'Output');
        lbl(gl, 54, 1, 'Output folder');
        edt(gl, 54, 2, 'calibration.path_output', '', 0);
        brw(gl, 54, 3, 'calibration.path_output', 'dir');

        lbl(gl, 55, 1, 'Profiles output folder');
        edt(gl, 55, 2, 'calibration.path_output_profiles', '', 0);
        brw(gl, 55, 3, 'calibration.path_output_profiles', 'dir');

        lbl(gl, 56, 1, 'Calibrated CSV filename');
        edt(gl, 56, 2, 'calibration.filename_calibrated_CSV', '', 0);
        note_lbl(gl, 57, 'Output filename for the calibrated intensity profile CSV.');

        lbl(gl, 58, 1, 'Save in calibration folder');
        chk(gl, 58, 2, 'calibration.save_in_calibration_folder', false, 'Route all outputs to the calibration output folder');

        lbl(gl, 59, 1, 'BabelBrain export path');
        edt(gl, 59, 2, 'calibration.ExportBabelBrain', '', 0);
        brw(gl, 59, 3, 'calibration.ExportBabelBrain', 'file');
        note_lbl(gl, 60, 'Optional .h5 output path for BabelBrain element weights (delta-from-geo phasors).');
        sp(gl, 61);

        sec(gl, 62, 'Run Calibration');
        btn_cal = uibutton(gl, ...
            'Text',       '▶   Run Calibration', ...
            'FontName',   st.font, 'FontSize', st.fs_lg, 'FontWeight', 'bold', ...
            'BackgroundColor', [0.18 0.55 0.34], ...
            'FontColor',  [1 1 1]);
        btn_cal.Layout.Row    = 63;
        btn_cal.Layout.Column = [1 4];
        btn_cal.ButtonPushedFcn = @(~,~) cb_run_calibration();

        lbl_cal_status = uilabel(gl, ...
            'Text',       'Ready.', ...
            'FontName',   st.font, 'FontSize', st.fs, ...
            'FontColor',  st.text_sub, ...
            'HorizontalAlignment', 'left');
        lbl_cal_status.Layout.Row    = 64;
        lbl_cal_status.Layout.Column = [1 4];
        sp(gl, 65);
    end

    function cb_run_calibration()
        % Gather calibration fields and call calibration_pipeline_start.
        % Minimal stub — expand once the full config-to-struct mapping is wired.
        params = struct();
        fields = findobj(ancestor(gcf,'figure'), '-not', 'Type', 'figure', 'Tag', '-regexp', '^calibration\.');
        for k = 1:numel(fields)
            tag = fields(k).Tag;
            parts = strsplit(tag, '.');
            val   = fields(k).Value;
            params.(parts{1}).(parts{2}) = val;
        end
        disp('calibration parameters collected:');
        disp(params);
        % TODO: call calibration_pipeline_start(params) once config path is wired.
    end

    %% ── Tab 10: Run ────────────────────────────────────────────────────
    function build_tab_run(t)
        gl = uigridlayout(t, ...
            'RowHeight',    {44, 28, '1x'}, ...
            'ColumnWidth',  {'1x'}, ...
            'Padding',      [16 16 16 16], ...
            'RowSpacing',   8, ...
            'BackgroundColor', st.bg_panel);

        btn_run = uibutton(gl, ...
            'Text',            '▶   Run Simulation', ...
            'FontName',        st.font, 'FontSize', st.fs_lg, 'FontWeight', 'bold', ...
            'BackgroundColor', st.success, 'FontColor', [1 1 1], ...
            'Tag',             'btn_run', ...
            'ButtonPushedFcn', @(~,~) cb_run());
        btn_run.Layout.Row = 1; btn_run.Layout.Column = 1;

        % Stage tracker / status label
        lbl_status = uilabel(gl, ...
            'Text',       'Ready.', ...
            'FontName',   st.font, 'FontSize', st.fs, ...
            'FontColor',  st.text_sub, ...
            'Tag',        'lbl_status', ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment',   'center');
        lbl_status.Layout.Row = 2; lbl_status.Layout.Column = 1;

        % Log area — terminal-style
        log_area = uitextarea(gl, ...
            'Value',           {'Simulation log will appear here.'}, ...
            'FontName',        'Courier New', ...
            'FontSize',        st.fs_sm + 1, ...
            'FontColor',       [0.80 0.95 0.75], ...
            'BackgroundColor', [0.10 0.12 0.10], ...
            'Editable',        'off', ...
            'Tag',             'log_area');
        log_area.Layout.Row = 3; log_area.Layout.Column = 1;
    end

    %% ── Tab 11: Results ───────────────────────────────────────────────
    function build_tab_results(t)
        gl = uigridlayout(t, ...
            'RowHeight',   {36, '1x'}, ...
            'ColumnWidth', {'1x'}, ...
            'Padding',     [16 16 16 16], ...
            'RowSpacing',  8, ...
            'BackgroundColor', st.bg_panel);

        % Toolbar
        btn_gl = uigridlayout(gl, ...
            'RowHeight',   {'1x'}, ...
            'ColumnWidth', {160, 160, 160, '1x'}, ...
            'Padding',     [0 0 0 0], 'RowSpacing', 0, 'ColumnSpacing', 8, ...
            'BackgroundColor', st.bg_panel);
        btn_gl.Layout.Row = 1; btn_gl.Layout.Column = 1;

        uibutton(btn_gl, 'Text', '🔄 Refresh', ...
            'FontName', st.font, 'FontSize', st.fs, ...
            'ButtonPushedFcn', @(~,~) refresh_results());
        uibutton(btn_gl, 'Text', '🌐 Open HTML Report', ...
            'FontName', st.font, 'FontSize', st.fs, ...
            'ButtonPushedFcn', @(~,~) open_html_report());
        uibutton(btn_gl, 'Text', '📁 Open Output Folder', ...
            'FontName', st.font, 'FontSize', st.fs, ...
            'ButtonPushedFcn', @(~,~) open_output_folder());

        % Sub-tab group for result types
        res_tg = uitabgroup(gl, 'Tag', 'result_tabgroup');
        res_tg.Layout.Row = 2; res_tg.Layout.Column = 1;

        % Acoustic maps
        tab_ac = uitab(res_tg, 'Title', 'Acoustic');
        tab_ac.BackgroundColor = st.bg_panel;
        ac_gl = uigridlayout(tab_ac, ...
            'RowHeight',   {'1x'}, ...
            'ColumnWidth', {'1x','1x','1x'}, ...
            'Padding',     [8 8 8 8], 'RowSpacing', 8, 'ColumnSpacing', 12, ...
            'BackgroundColor', st.bg_panel);
        ax_labels = {'X','Y','Z'};
        for k = 1:3
            ax = uiaxes(ac_gl);
            ax.Layout.Row = 1; ax.Layout.Column = k;
            ax.Tag = sprintf('ax_acoustic_%d', k);
            ax.XTick = []; ax.YTick = [];
            title(ax, ax_labels{k}, 'FontName', st.font);
        end

        % Thermal maps
        tab_th = uitab(res_tg, 'Title', 'Thermal');
        tab_th.BackgroundColor = st.bg_panel;
        th_gl = uigridlayout(tab_th, ...
            'RowHeight',   {'1x'}, ...
            'ColumnWidth', {'1x','1x','1x'}, ...
            'Padding',     [8 8 8 8], 'RowSpacing', 8, 'ColumnSpacing', 12, ...
            'BackgroundColor', st.bg_panel);
        for k = 1:3
            ax = uiaxes(th_gl);
            ax.Layout.Row = 1; ax.Layout.Column = k;
            ax.Tag = sprintf('ax_thermal_%d', k);
            ax.XTick = []; ax.YTick = [];
            title(ax, ax_labels{k}, 'FontName', st.font);
        end

        % HTML report viewer
        tab_rp = uitab(res_tg, 'Title', 'Report');
        tab_rp.BackgroundColor = st.bg_panel;
        rp_gl = uigridlayout(tab_rp, ...
            'RowHeight',   {'1x'}, ...
            'ColumnWidth', {'1x'}, ...
            'Padding',     [0 0 0 0], ...
            'BackgroundColor', st.bg_panel);
        html_viewer = uihtml(rp_gl, 'Tag', 'html_viewer');
        html_viewer.Layout.Row = 1; html_viewer.Layout.Column = 1;

        % CSV table
        tab_csv = uitab(res_tg, 'Title', 'Data Table');
        tab_csv.BackgroundColor = st.bg_panel;
        csv_gl = uigridlayout(tab_csv, ...
            'RowHeight',   {'1x'}, ...
            'ColumnWidth', {'1x'}, ...
            'Padding',     [8 8 8 8], ...
            'BackgroundColor', st.bg_panel);
        uitable(csv_gl, 'Tag', 'csv_table', 'FontName', st.font, 'FontSize', st.fs_sm);
    end

%% ════════════════════════════════════════════════════════════════════════
%%  CALLBACKS
%% ════════════════════════════════════════════════════════════════════════

    function cb_pct_toggle(src, sub_controls)
        % Enable/disable pCT sub-controls based on the pct.enabled checkbox
        if src.Value
            en = 'on';
        else
            en = 'off';
        end
        for k = 1:numel(sub_controls)
            sub_controls(k).Enable = en;
        end
    end

    function cb_load_yaml()
        [f, d] = uigetfile({'*.yaml;*.yml','YAML config'}, 'Load YAML');
        if isequal(f,0), return; end
        load_yaml_to_gui(fullfile(d,f));
    end

    function cb_save_yaml()
        % Guard: prevent dialog from opening twice if button is hit rapidly
        if isappdata(fig, 'save_dialog_open') && getappdata(fig, 'save_dialog_open')
            return;
        end
        setappdata(fig, 'save_dialog_open', true);
        drawnow;
        [f, d] = uiputfile({'*.yaml','YAML config'}, 'Save parameters as YAML', 'config_study.yaml');
        % Restore focus — uiputfile can push the figure behind on macOS
        figure(fig);
        setappdata(fig, 'save_dialog_open', false);
        if isequal(f,0), return; end
        try
            params = collect_params();
            yaml.dumpFile(fullfile(d,f), params);
            uialert(fig, sprintf('Saved to:\n%s', fullfile(d,f)), 'Saved', 'Icon', 'success');
        catch ME
            fprintf('\n=== PRESTUS GUI Save Error ===\n%s\n==============================\n', ME.getReport('extended'));
            uialert(fig, ME.message, 'Save Error', 'Icon', 'error');
        end
    end

    function cb_run()
        try
            params = collect_params();
        catch ME
            uialert(fig, ME.message, 'Parameter Error', 'Icon', 'error'); return;
        end

        % Write temporary YAML for the worker
        tmp = [tempname '.yaml'];
        try
            yaml.dumpFile(tmp, params);
        catch ME
            uialert(fig, ['Could not write config: ' ME.message], 'YAML Error', 'Icon', 'error'); return;
        end

        % Resolve output dir and log path
        app = fig.UserData;
        out_dir = '';
        if isfield(params, 'io') && isfield(params.io, 'dir_output')
            out_dir = params.io.dir_output;
        end
        app.dir_output = out_dir;
        app.log_path   = '';
        if ~isempty(out_dir)
            sub_id = 1;
            if isfield(params, 'subject_id') && isnumeric(params.subject_id) && ~isnan(params.subject_id)
                sub_id = params.subject_id;
            end
            medium = 'sim';
            if isfield(params, 'simulation') && isfield(params.simulation, 'medium')
                medium = params.simulation.medium;
            end
            affix = '';
            if isfield(params, 'io') && isfield(params.io, 'output_affix')
                affix = params.io.output_affix;
            end
            app.log_path = fullfile(out_dir, sprintf('sub-%03d_%s%s.txt', sub_id, medium, affix));
        end

        % Update UI before run
        set_status('⟳  Running...', st.accent);
        h_run = findobj(fig, 'Tag', 'btn_run');
        if ~isempty(h_run), h_run.Enable = 'off'; end
        clear_log();
        append_log('▶ Starting simulation...');
        drawnow;

        % Run simulation synchronously.
        % diary captures all disp/fprintf output; populated into log after run.
        % For HPC platforms (slurm/qsub) the pipeline submits a job and returns
        % immediately, so blocking is not an issue in practice.
        % The pipeline sets up its own diary via path_log_setup
        % (parameters.io.log_file). We run synchronously and read that file
        % afterwards — no need to open a separate diary here.
        try
            parameters = load_parameters(tmp);
            prestus_pipeline_start(parameters);
            % Read the pipeline's own log file
            if isfield(parameters, 'io') && isfield(parameters.io, 'log_file')
                populate_log_from_file(parameters.io.log_file);
            end
            set_status('✓  Done', st.success);
            append_log('✓ Simulation completed.');
            app = fig.UserData;
            app.dir_output = out_dir;
            fig.UserData = app;
            refresh_results();
        catch ME
            if exist('parameters', 'var') && isfield(parameters, 'io') && isfield(parameters.io, 'log_file')
                populate_log_from_file(parameters.io.log_file);
            end
            fprintf('\n=== PRESTUS GUI Error ===\n%s\n========================\n', ME.getReport('extended'));
            set_status(['✗  ' ME.message], st.mandatory);
            append_log(['✗ Error: ' ME.message]);
            uialert(fig, sprintf('%s\n\n(Full stack trace printed to MATLAB terminal)', ME.message), ...
                'Simulation Error', 'Icon', 'error');
        end
        h_run = findobj(fig, 'Tag', 'btn_run');
        if ~isempty(h_run), h_run.Enable = 'on'; end
    end

    function cb_ingest()
        % Collect required fields from the Neuronavigation Ingestion widgets
        sub_id_h   = findobj(fig, 'Tag', 'neuronav.sub_id');
        ses_id_h   = findobj(fig, 'Tag', 'neuronav.ses_id');
        raw_path_h = findobj(fig, 'Tag', 'neuronav.raw_path');
        out_path_h = findobj(fig, 'Tag', 'neuronav.out_path');
        status_h   = findobj(fig, 'Tag', 'lbl_ingest_status');

        sub_id   = strtrim(sub_id_h.Value);
        ses_id   = strtrim(ses_id_h.Value);
        raw_path = strtrim(raw_path_h.Value);
        out_path = strtrim(out_path_h.Value);

        if isempty(sub_id) || isempty(ses_id) || isempty(raw_path) || isempty(out_path)
            status_h.Text      = '✗ Subject ID, Session ID, Raw path and Output path are all required.';
            status_h.FontColor = st.mandatory;
            return;
        end

        % Build a minimal parameters struct sufficient for neuronav_ingest_markers
        params_ingest = struct();
        params_ingest.path.localite_raw  = raw_path;
        params_ingest.path.localite_post = out_path;
        params_ingest.io.overwrite_files = 'always';

        % Build a minimal single-target map (series 1, transducer 1)
        target_map.series_index  = 1;
        target_map.transducer_id = 1;
        target_map.target_name   = 'target';

        btn_h = findobj(fig, 'Tag', 'btn_ingest');
        if ~isempty(btn_h), btn_h.Enable = 'off'; end
        status_h.Text      = '⟳ Running ingestion…';
        status_h.FontColor = st.accent;
        drawnow;

        try
            [xml_path, json_path] = neuronav_ingest_markers(params_ingest, sub_id, ses_id, target_map);
            status_h.Text      = sprintf('✓ Done.  XML: %s', xml_path);
            status_h.FontColor = st.success;
        catch ME
            status_h.Text      = sprintf('✗ %s', ME.message);
            status_h.FontColor = st.mandatory;
            fprintf('\n=== Neuronavigation Ingestion Error ===\n%s\n=======================================\n', ME.getReport('extended'));
        end

        if ~isempty(btn_h), btn_h.Enable = 'on'; end
    end

    function populate_log_from_file(diary_file)
        if ~exist(diary_file, 'file'), return; end
        fid = fopen(diary_file, 'r');
        if fid == -1, return; end
        la = findobj(fig, 'Tag', 'log_area');
        while ~feof(fid)
            line = fgetl(fid);
            if ~ischar(line), continue; end
            line = strtrim(line);
            if isempty(line), continue; end
            if ~isempty(la), la.Value{end+1} = line; end
            update_stage(line);
        end
        fclose(fid);
        try; scroll(la, 'bottom'); catch; end
        drawnow;
    end

    function update_stage(line)
        % Detect pipeline section banners and update status label
        persistent last_was_banner;
        if isempty(last_was_banner), last_was_banner = false; end
        if contains(line, '========')
            last_was_banner = true;
        elseif last_was_banner && ~isempty(line)
            set_status(['⟳  ' strtrim(line)], st.accent);
            last_was_banner = false;
        else
            last_was_banner = false;
        end
    end


    function cb_toggle_thermal_timing(enabled)
    % Enable or disable the sonication timing fields in the Thermal tab.
    % These fields are only meaningful when thermal simulation is requested.
        en = 'off';
        if enabled; en = 'on'; end
        timing_tags = { ...
            'timing.pd', 'timing.pri', 'timing.ptd', 'timing.ptri', ...
            'timing.ptrd', 'timing.post_ptri_dur', ...
            'timing.pt_timestep', 'timing.post_pt_timestep', ...
            'timing.equal_step_duration' };
        for k = 1:numel(timing_tags)
            h = findobj(fig, 'Tag', timing_tags{k});
            if ~isempty(h); set(h, 'Enable', en); end
        end
    end

    function cb_transducer_type(dd)
        pnl_ann = findobj(fig, 'Tag', 'panel_annular');
        pnl_mat = findobj(fig, 'Tag', 'panel_matrix');
        if strcmp(dd.Value, 'annular')
            pnl_ann.Visible = 'on';
            pnl_mat.Visible = 'off';
        else
            pnl_ann.Visible = 'off';
            pnl_mat.Visible = 'on';
        end
    end

    function cb_serial_changed(dd)
        serial = dd.Value;
        if strcmp(serial, '(manual)'), return; end
        try
            eq   = load_equipment_config();
            if ~isfield(eq.trans, serial), return; end
            geom = eq.trans.(serial).transducer;

            % Populate type
            if isfield(geom, 'annular')
                set_widget('transducer.type', 'annular');
                cb_transducer_type(findobj(fig, 'Tag', 'transducer.type'));
            elseif isfield(geom, 'matrix')
                set_widget('transducer.type', 'matrix');
                cb_transducer_type(findobj(fig, 'Tag', 'transducer.type'));
            end

            % Populate frequency
            if isfield(geom, 'freq_hz')
                set_widget('transducer.freq_hz', geom.freq_hz);
            end

            % Populate annular geometry fields
            if isfield(geom, 'annular')
                ann = geom.annular;
                if isfield(ann, 'elem_n')
                    set_widget('transducer.annular.elem_n', ann.elem_n);
                end
                if isfield(ann, 'elem_id_mm') && isnumeric(ann.elem_id_mm)
                    h = findobj(fig, 'Tag', 'transducer.annular.elem_id_mm');
                    if ~isempty(h)
                        h.Value = strjoin(arrayfun(@num2str, ann.elem_id_mm(:)', 'UniformOutput', false), ',');
                    end
                end
                if isfield(ann, 'elem_od_mm') && isnumeric(ann.elem_od_mm)
                    h = findobj(fig, 'Tag', 'transducer.annular.elem_od_mm');
                    if ~isempty(h)
                        h.Value = strjoin(arrayfun(@num2str, ann.elem_od_mm(:)', 'UniformOutput', false), ',');
                    end
                end
                if isfield(ann, 'curv_radius_mm')
                    set_widget('transducer.annular.curv_radius_mm', ann.curv_radius_mm);
                end
                if isfield(ann, 'dist_geom_ep_mm')
                    set_widget('transducer.annular.dist_geom_ep_mm', ann.dist_geom_ep_mm);
                end
            end
        catch ME
            append_log(sprintf('Could not load geometry for %s: %s', serial, ME.message));
        end

        % Update library coverage label
        update_library_coverage_label(serial);
    end

    function update_library_coverage_label(serial_val)
        lbl_cov = findobj(fig, 'Tag', 'lbl_library_coverage');
        if isempty(lbl_cov), return; end
        if strcmp(serial_val, '(manual)')
            lbl_cov.Text = 'Library: —';
            return;
        end
        try
            lib_path  = fullfile(get_prestus_path(), 'config', 'transducer');
            yaml_path = fullfile(lib_path, [serial_val '.yaml']);
            if isfile(yaml_path)
                lib = yaml.loadFile(yaml_path, 'ConvertToArray', true);
                if isfield(lib, 'global_model') && isfield(lib.global_model, 'depths_ep_mm')
                    depths = lib.global_model.depths_ep_mm(:)';
                elseif isfield(lib, 'calibration') && isfield(lib.calibration, 'focal_depths')
                    dk     = fieldnames(lib.calibration.focal_depths);
                    depths = sort(cellfun(@(k) str2double(strrep(strrep(k,'f',''),'p','.')), dk));
                else
                    depths = [];
                end
                if isempty(depths)
                    lbl_cov.Text = 'Library: no calibrated depths';
                else
                    depth_strs   = arrayfun(@(d) sprintf('%.0f', d), depths(:)', 'UniformOutput', false);
                    lbl_cov.Text = ['Library: ' strjoin(depth_strs, ', ') ' mm'];
                end
            else
                lbl_cov.Text = 'Library: not found';
            end
        catch
            lbl_cov.Text = 'Library: (error reading)';
        end
    end

    function cb_placement_mode(mode)
        % Show the sub-panel matching the selected placement mode; hide others.
        tags   = {'panel_placement_localite','panel_placement_heuristic','panel_placement_plantus'};
        modes  = {'localite','heuristic','plantus'};
        for k = 1:numel(tags)
            h = findobj(fig, 'Tag', tags{k});
            if ~isempty(h)
                if strcmp(mode, modes{k})
                    h.Visible = 'on';
                else
                    h.Visible = 'off';
                end
            end
        end
    end

    function cb_seq_browse_append()
        [f, d] = uigetfile({'*.yaml;*.yml','YAML files';'*.*','All files'}, ...
            'Select follow-up config YAML');
        if isequal(f, 0), return; end
        p = fullfile(d, f);
        h = findobj(fig, 'Tag', 'sequential_configs_textarea');
        if ~isempty(h)
            lines = h.Value;
            if numel(lines) == 1 && isempty(strtrim(lines{1}))
                h.Value = {p};
            else
                h.Value = [lines; {p}];
            end
        end
    end

    function cb_multitrans_browse_append()
        [f, d] = uigetfile({'*.yaml;*.yml','YAML files';'*.*','All files'}, ...
            'Select secondary transducer config YAML');
        if isequal(f, 0), return; end
        p = fullfile(d, f);
        h = findobj(fig, 'Tag', 'multitrans_configs_textarea');
        if ~isempty(h)
            lines = h.Value;
            if numel(lines) == 1 && isempty(strtrim(lines{1}))
                h.Value = {p};
            else
                h.Value = [lines; {p}];
            end
        end
    end

    function cb_browse(target_tag, browse_type)
        if strcmp(browse_type, 'dir')
            p = uigetdir();
            if isequal(p, 0), return; end
        elseif strcmp(browse_type, 'file')
            [f, d] = uigetfile({'*.*','All files'});
            if isequal(f, 0), return; end
            p = fullfile(d, f);
        else
            [f, d] = uigetfile({'*.yaml;*.yml','YAML files';'*.*','All files'});
            if isequal(f, 0), return; end
            p = fullfile(d, f);
        end
        h = findobj(fig, 'Tag', target_tag);
        if ~isempty(h), h.Value = p; end
    end

    function refresh_results()
        app = fig.UserData;
        if isempty(app.dir_output) || ~isfolder(app.dir_output)
            try
                params = collect_params();
                if isfield(params, 'io') && isfield(params.io, 'dir_output') && ~isempty(params.io.dir_output)
                    app.dir_output = params.io.dir_output;
                    fig.UserData   = app;
                end
            catch; end
        end
        if ~isempty(app.dir_output) && isfolder(app.dir_output)
            try; show_results(app.dir_output); catch; end
        end
    end

    function open_html_report()
        app = fig.UserData;
        if isempty(app.dir_output), return; end
        reports = dir(fullfile(app.dir_output, '*report*.html'));
        if isempty(reports)
            uialert(fig, 'No HTML report found in output directory.', 'Not found');
            return;
        end
        web(fullfile(app.dir_output, reports(end).name), '-browser');
    end

    function open_output_folder()
        app = fig.UserData;
        if isempty(app.dir_output) || ~isfolder(app.dir_output)
            uialert(fig, 'Output directory not set or does not exist.', 'Not found');
            return;
        end
        if ispc
            system(['explorer "' app.dir_output '"']);
        elseif ismac
            system(['open "' app.dir_output '"']);
        else
            system(['xdg-open "' app.dir_output '"']);
        end
    end

%% ════════════════════════════════════════════════════════════════════════
%%  I/O: LOAD / SAVE / COLLECT
%% ════════════════════════════════════════════════════════════════════════

    function load_defaults()
        try
            prestus_root  = get_prestus_path();
            default_yaml  = fullfile(prestus_root, 'config', 'config_default.yaml');
            params        = yaml.loadFile(default_yaml, 'ConvertToArray', true);
            apply_params_to_gui(params);
            set_config_label(default_yaml);
        catch ME
            append_log(sprintf('Could not load defaults: %s', ME.message));
        end
    end

    function load_yaml_to_gui(filepath)
        try
            params = yaml.loadFile(filepath, 'ConvertToArray', true);
            apply_params_to_gui(params);
            set_config_label(filepath);
            set_status(sprintf('Loaded: %s', filepath), st.success);
        catch ME
            fprintf('\n=== PRESTUS GUI Load Error ===\n%s\n==============================\n', ME.getReport('extended'));
            uialert(fig, ME.message, 'Load Error', 'Icon', 'error');
        end
    end

    function set_config_label(filepath)
        lbl_config.Text = sprintf('Config: %s', filepath);
        app = fig.UserData;
        app.config_path = filepath;
        fig.UserData = app;
    end

    function apply_params_to_gui(params)
        % Flatten struct to tag->value map and set each widget
        flat = flatten_struct(params, '');
        keys = fieldnames(flat);
        for i = 1:numel(keys)
            tag = strrep(keys{i}, '__', '.');
            set_widget(tag, flat.(keys{i}));
        end
        % XYZ array fields stored as separate _1/_2/_3 widgets
        xyz_tags = {'transducer.trans_pos', 'transducer.focus_pos', 'grid.default_dims', ...
                    'placement.heuristic.mni_target_mm', 'placement.plantus.mni_target_mm'};
        for i = 1:numel(xyz_tags)
            tag     = xyz_tags{i};
            safe    = strrep(tag, '.', '__');
            if isfield(flat, safe) && isnumeric(flat.(safe))
                v = double(flat.(safe));
                if numel(v) >= 3
                    % 3-D: load all three components directly
                    for k = 1:3
                        hw = cached_findobj(sprintf('%s_%d', tag, k));
                        if ~isempty(hw), hw.Value = v(k); end
                    end
                elseif numel(v) == 2
                    % 2-D (axisymmetric): map to X(_1) and Z(_3), leave Y(_2) as-is
                    for k = [1 3]
                        hw = cached_findobj(sprintf('%s_%d', tag, k));
                        if ~isempty(hw), hw.Value = v(ceil(k/2)); end
                    end
                end
            end
        end
        % Numeric array fields displayed as comma-separated strings
        csv_array_tags = {'transducer.annular.elem_id_mm', ...
                          'transducer.annular.elem_od_mm', ...
                          'transducer.annular.elem_phase_deg', ...
                          'transducer.target_isppa_wcm2', ...
                          'placement.plantus.focal_distance_list', ...
                          'placement.plantus.flhm_list', ...
                          'calibration.elem_phase_correction_deg', ...
                          'layers.water', 'layers.brain', 'layers.skin', ...
                          'layers.skull', 'layers.skull_cortical', 'layers.skull_trabecular'};
        for i = 1:numel(csv_array_tags)
            tag  = csv_array_tags{i};
            safe = strrep(tag, '.', '__');
            if isfield(flat, safe) && isnumeric(flat.(safe)) && ~isempty(flat.(safe))
                hw = cached_findobj(tag);
                if ~isempty(hw) && isa(hw, 'matlab.ui.control.EditField')
                    hw.Value = strjoin(arrayfun(@num2str, flat.(safe)(:)', 'UniformOutput', false), ',');
                end
            end
        end
        % Medium properties table
        if isfield(params, 'medium_properties')
            set_medium_table(params.medium_properties);
        end

        % Sequential configs textarea
        seq_ta_h = cached_findobj('sequential_configs_textarea');
        if ~isempty(seq_ta_h) && isfield(params, 'options') && isfield(params.options, 'sequential_configs')
            sc = params.options.sequential_configs;
            lines = {};
            for k = 1:100
                fname = sprintf('config_%d', k);
                if isfield(sc, fname), lines{end+1} = sc.(fname); else, break; end %#ok<AGROW>
            end
            if ~isempty(lines), seq_ta_h.Value = lines; end
        end

        % Multi-transducer configs textarea
        mt_ta_h = cached_findobj('multitrans_configs_textarea');
        if ~isempty(mt_ta_h) && isfield(params, 'options') && isfield(params.options, 'additional_transducer_configs')
            ac = params.options.additional_transducer_configs;
            lines = {};
            for k = 1:100
                fname = sprintf('config_%d', k);
                if isfield(ac, fname), lines{end+1} = ac.(fname); else, break; end %#ok<AGROW>
            end
            if ~isempty(lines), mt_ta_h.Value = lines; end
        end

        % Update placement sub-panel visibility to match loaded mode
        h_pm = cached_findobj('placement.mode');
        if ~isempty(h_pm)
            cb_placement_mode(h_pm.Value);
        end

        % Re-sync thermal timing enable state to match loaded run_heating_sims value.
        h_heat = cached_findobj('modules.run_heating_sims');
        if ~isempty(h_heat)
            cb_toggle_thermal_timing(h_heat.Value);
        end

        % Sync transducer serial dropdown: add serial to items list if not
        % present (e.g. config was written without launching GUI), then select it.
        h_serial = cached_findobj('transducer.serial');
        if ~isempty(h_serial) && isfield(params, 'transducer')
            tr1 = params.transducer;
            if isstruct(tr1) && numel(tr1) >= 1; tr1 = tr1(1); end
            if isstruct(tr1) && isfield(tr1, 'serial') && ~isempty(tr1.serial)
                serial_val = char(tr1.serial);
                if ~any(strcmp(h_serial.Items, serial_val))
                    h_serial.Items = [h_serial.Items, {serial_val}];
                end
                h_serial.Value = serial_val;
                update_library_coverage_label(serial_val);
            end
        end
    end

    function params = collect_params()
        params = struct();
        % Collect all tagged input widgets
        all_widgets = [findobj(fig, '-isa', 'matlab.ui.control.EditField'); ...
                       findobj(fig, '-isa', 'matlab.ui.control.NumericEditField'); ...
                       findobj(fig, '-isa', 'matlab.ui.control.DropDown');   ...
                       findobj(fig, '-isa', 'matlab.ui.control.CheckBox')];
        for i = 1:numel(all_widgets)
            tag = all_widgets(i).Tag;
            if isempty(tag) || strcmp(tag, ''), continue; end
            val = get_widget_value(all_widgets(i));
            % Skip empty strings: yaml.dumpFile serialises '' as [] (null),
            % which breaks pipeline code that expects char. Empty fields
            % intentionally fall back to the config_default.yaml value.
            if ischar(val) && isempty(val), continue; end
            % Skip NaN-sentinel numerics (fields created with nedt(…, NaN)
            % that the user left at their "not set" state).  These fall back
            % to config_default.yaml / pipeline defaults.
            if isnumeric(val) && isscalar(val) && isnan(val), continue; end
            % Skip the sentinel value for the serial dropdown — it means
            % "not selected" so no serial should be written to the config.
            if strcmp(tag, 'transducer.serial') && strcmp(val, '(manual)'), continue; end
            params = set_nested(params, tag, val);
        end
        % XYZ triplets stored as separate tagged fields
        xyz_tags = {'transducer.trans_pos', 'transducer.focus_pos', 'grid.default_dims', ...
                    'placement.heuristic.mni_target_mm', 'placement.plantus.mni_target_mm'};
        for i = 1:numel(xyz_tags)
            tag = xyz_tags{i};
            vals = zeros(1,3);
            for ax_i = 1:3
                hw = cached_findobj(sprintf('%s_%d', tag, ax_i));
                if ~isempty(hw), vals(ax_i) = hw.Value; end
            end
            params = set_nested(params, tag, vals);
        end
        % For 2D axisymmetric simulations drop the Y (2nd) component from
        % position and dimension fields so they match the 2-element grid.dims.
        ax_h = cached_findobj('grid.axisymmetric');
        if ~isempty(ax_h) && ax_h.Value
            for pos_tag = {'transducer.trans_pos', 'transducer.focus_pos', 'grid.default_dims'}
                t = pos_tag{1};
                v = get_nested(params, t);
                if isnumeric(v) && numel(v) == 3
                    params = set_nested(params, t, v([1 3]));
                end
            end
        end
        % Parse comma-separated fields into numeric arrays
        csv_array_tags = {'transducer.annular.elem_id_mm', ...
                          'transducer.annular.elem_od_mm', ...
                          'transducer.annular.elem_phase_deg', ...
                          'transducer.target_isppa_wcm2', ...
                          'placement.plantus.focal_distance_list', ...
                          'placement.plantus.flhm_list', ...
                          'calibration.elem_phase_correction_deg', ...
                          'layers.water', 'layers.brain', 'layers.skin', ...
                          'layers.skull', 'layers.skull_cortical', 'layers.skull_trabecular'};
        for i = 1:numel(csv_array_tags)
            tag = csv_array_tags{i};
            hw = cached_findobj(tag);
            if ~isempty(hw) && isa(hw, 'matlab.ui.control.EditField')
                raw = strtrim(hw.Value);
                if ~isempty(raw)
                    nums = str2num(raw); %#ok<ST2NM>
                    if ~isempty(nums)
                        params = set_nested(params, tag, nums);
                    end
                end
            end
        end

        % Medium properties table
        tbl = cached_findobj('medium_table');
        if ~isempty(tbl)
            params = read_medium_table(params, tbl);
        end

        % Sequential configs textarea (one YAML path per line)
        seq_ta_h = cached_findobj('sequential_configs_textarea');
        if ~isempty(seq_ta_h)
            cfg_idx = 0;
            for k = 1:numel(seq_ta_h.Value)
                p = strtrim(seq_ta_h.Value{k});
                if ~isempty(p)
                    cfg_idx = cfg_idx + 1;
                    params = set_nested(params, ...
                        sprintf('options.sequential_configs.config_%d', cfg_idx), p);
                end
            end
        end

        % Multi-transducer configs textarea (one YAML path per line)
        mt_ta_h = cached_findobj('multitrans_configs_textarea');
        if ~isempty(mt_ta_h)
            cfg_idx = 0;
            for k = 1:numel(mt_ta_h.Value)
                p = strtrim(mt_ta_h.Value{k});
                if ~isempty(p)
                    cfg_idx = cfg_idx + 1;
                    params = set_nested(params, ...
                        sprintf('options.additional_transducer_configs.config_%d', cfg_idx), p);
                end
            end
        end

        % Derived: output_affix and output_dir
        if ~isfield(params, 'io') || ~isfield(params.io, 'output_affix')
            params.io.output_affix = '';
        end
        if isfield(params, 'path') && isfield(params.path, 'sim') && ~isempty(params.path.sim)
            sub_id = 1;
            if isfield(params, 'subject_id') && isnumeric(params.subject_id) && ~isnan(params.subject_id)
                sub_id = params.subject_id;
            end
            params.io.dir_output = fullfile(params.path.sim, sprintf('sub-%03d', sub_id));
        end
    end

    function ok = validate_params(params)
        ok = true;
        missing = {};
        required = {'subject_id', 'transducer.freq_hz'};
        labels   = {'Subject ID (I/O tab)', 'Transducer frequency (Transducer tab)'};
        for i = 1:numel(required)
            val = get_nested(params, required{i});
            if isempty(val) || (ischar(val) && isempty(strtrim(val))) || (isnumeric(val) && (isnan(val) || val == 0))
                missing{end+1} = labels{i}; %#ok<AGROW>
            end
        end
        if ~isempty(missing)
            uialert(fig, ['Please fill in: ' strjoin(missing, ', ')], ...
                'Required fields missing', 'Icon', 'warning');
            ok = false;
        end
    end

    function show_results(output_dir)
        % Acoustic PNG maps
        dims = {'x','y','z'};
        img_dir = fullfile(output_dir, 'img');
        for d = 1:3
            pngs = dir(fullfile(img_dir, sprintf('*intensity_%s*.png', dims{d})));
            ax   = findobj(fig, 'Tag', sprintf('ax_acoustic_%d', d));
            if ~isempty(pngs) && ~isempty(ax)
                img = imread(fullfile(img_dir, pngs(end).name));
                imshow(img, 'Parent', ax);
            end
        end
        % Thermal PNG maps
        for d = 1:3
            pngs = dir(fullfile(img_dir, sprintf('*maxT_%s*.png', dims{d})));
            ax   = findobj(fig, 'Tag', sprintf('ax_thermal_%d', d));
            if ~isempty(pngs) && ~isempty(ax)
                img = imread(fullfile(img_dir, pngs(end).name));
                imshow(img, 'Parent', ax);
            end
        end
        % HTML report
        reports = dir(fullfile(output_dir, '*report*.html'));
        viewer  = findobj(fig, 'Tag', 'html_viewer');
        if ~isempty(reports) && ~isempty(viewer)
            viewer.HTMLSource = fullfile(output_dir, reports(end).name);
        end
        % CSV table
        csvs  = dir(fullfile(output_dir, 'sub-*_medium-*.csv'));
        tbl_h = findobj(fig, 'Tag', 'csv_table');
        if ~isempty(csvs) && ~isempty(tbl_h)
            try
                T = readtable(fullfile(output_dir, csvs(end).name));
                tbl_h.Data        = table2cell(T);
                tbl_h.ColumnName  = T.Properties.VariableNames;
            catch
            end
        end
    end

%% ════════════════════════════════════════════════════════════════════════
%%  WIDGET HELPERS
%% ════════════════════════════════════════════════════════════════════════

    function gl = tab_grid(parent, n_rows, col_widths)
        row_h = repmat({26}, 1, n_rows);
        gl = uigridlayout(parent, ...
            'RowHeight',   row_h, ...
            'ColumnWidth', col_widths, ...
            'Padding',     [20 16 20 16], ...
            'RowSpacing',  4, ...
            'ColumnSpacing', 10, ...
            'BackgroundColor', st.bg_panel, ...
            'Scrollable',  'on');
    end

    function lbl(parent, row, col, text)
        mandatory = startsWith(text, '*');
        disp_text = strtrim(strrep(text, '*', ''));
        if mandatory
            disp_text = ['<font color="#D12B1E">*</font> ' disp_text];
        end
        h = uilabel(parent, ...
            'Text',                disp_text, ...
            'Interpreter',         'html', ...
            'FontName',            st.font, ...
            'FontSize',            st.fs, ...
            'FontColor',           st.text, ...
            'HorizontalAlignment', 'right', ...
            'VerticalAlignment',   'center');
        if isnumeric(col) && numel(col) == 2
            h.Layout.Row = row; h.Layout.Column = col;
        else
            h.Layout.Row = row; h.Layout.Column = col;
        end
    end

    function h = edt(parent, row, col, tag, default_val, mandatory)
        h = uieditfield(parent, 'text', ...
            'Value',               char(default_val), ...
            'Tag',                 tag, ...
            'HorizontalAlignment', 'right', ...
            'FontName', st.font, 'FontSize', st.fs);
        if mandatory
            h.BackgroundColor = [1.0 0.95 0.95];
        end
        h.Layout.Row = row; h.Layout.Column = col;
    end

    function h = nedt(parent, row, col, tag, default_val)
        if isnan(default_val)
            h = uieditfield(parent, 'numeric', ...
                'Value',              0, ...
                'Placeholder',        'NaN', ...
                'HorizontalAlignment','right', ...
                'Tag',                tag, ...
                'FontName',           st.font, 'FontSize', st.fs);
        else
            h = uieditfield(parent, 'numeric', ...
                'Value',               default_val, ...
                'HorizontalAlignment', 'right', ...
                'Tag',                 tag, ...
                'FontName',            st.font, 'FontSize', st.fs);
        end
        h.Layout.Row = row; h.Layout.Column = col;
    end

    function h = drp(parent, row, col, tag, items, default_val)
        h = uidropdown(parent, ...
            'Items',    items, ...
            'Tag',      tag, ...
            'FontName', st.font, 'FontSize', st.fs);
        if any(strcmp(items, default_val))
            h.Value = default_val;
        end
        h.Layout.Row = row; h.Layout.Column = col;
    end

    function h = chk(parent, row, col, tag, default_val, label_text)
        h = uicheckbox(parent, ...
            'Text',       label_text, ...
            'Value',      default_val, ...
            'Tag',        tag, ...
            'FontName',   st.font, 'FontSize', st.fs, ...
            'FontColor',  st.text);
        if numel(col) == 2
            h.Layout.Row = row; h.Layout.Column = col;
        else
            h.Layout.Row = row; h.Layout.Column = col;
        end
    end

    function brw(parent, row, col, target_tag, browse_type)
        h = uibutton(parent, 'Text', '···', ...
            'FontName',        st.font, 'FontSize', st.fs, ...
            'BackgroundColor', [0.88 0.91 0.96], ...
            'ButtonPushedFcn', @(~,~) cb_browse(target_tag, browse_type));
        h.Layout.Row = row; h.Layout.Column = col;
    end

    function sec(parent, row, text)
        h = uilabel(parent, ...
            'Text',                ['  ' upper(text)], ...
            'FontName',            st.font, ...
            'FontSize',            st.fs_sm, ...
            'FontWeight',          'bold', ...
            'FontColor',           st.accent, ...
            'BackgroundColor',     st.accent_lt, ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment',   'center');
        h.Layout.Row = row; h.Layout.Column = [1 4];
    end

    function note_lbl(parent, row, text)
        h = uilabel(parent, ...
            'Text',       text, ...
            'FontName',   st.font, ...
            'FontSize',   st.fs_sm, ...
            'FontColor',  st.text_sub, ...
            'WordWrap',   'on', ...
            'HorizontalAlignment', 'left');
        h.Layout.Row = row; h.Layout.Column = [2 4];
    end

    function sp(parent, row)
        h = uilabel(parent, 'Text', '');
        h.Layout.Row = row; h.Layout.Column = 1;
    end

    function xyz_panel(parent, row, base_tag, defaults)
        % 6-column layout: [lbl ef lbl ef lbl ef] for X Y Z
        gl_xyz = uigridlayout(parent, ...
            'RowHeight',   {'1x'}, ...
            'ColumnWidth', {22,'1x', 22,'1x', 22,'1x'}, ...
            'Padding',     [0 0 0 0], 'RowSpacing', 0, 'ColumnSpacing', 2, ...
            'BackgroundColor', st.bg_panel);
        gl_xyz.Layout.Row = row; gl_xyz.Layout.Column = 2;
        labels = {'X','Y','Z'};
        for k = 1:3
            lbl_h = uilabel(gl_xyz, 'Text', labels{k}, 'FontName', st.font, ...
                'FontSize', st.fs_sm, 'FontColor', st.text_sub, ...
                'HorizontalAlignment', 'center');
            lbl_h.Layout.Row = 1; lbl_h.Layout.Column = 2*k - 1;
            if isnan(defaults(k))
                ef = uieditfield(gl_xyz, 'numeric', ...
                    'Value',       0, ...
                    'Placeholder', 'NaN', ...
                    'Tag',         sprintf('%s_%d', base_tag, k), ...
                    'FontName', st.font, 'FontSize', st.fs);
            else
                ef = uieditfield(gl_xyz, 'numeric', ...
                    'Value', defaults(k), ...
                    'Tag',   sprintf('%s_%d', base_tag, k), ...
                    'FontName', st.font, 'FontSize', st.fs);
            end
            ef.Layout.Row = 1; ef.Layout.Column = 2*k;
        end
    end

%% ════════════════════════════════════════════════════════════════════════
%%  STYLE CONSTANTS
%% ════════════════════════════════════════════════════════════════════════

    function s = sty()
        % Palette derived from the PRESTUS logo
        % Deep red: #a91d25  Teal: #0b7b99  Greys from logo background
        s.bg_fig    = [0.945 0.945 0.947];   % near-white warm grey
        s.bg_panel  = [1.000 1.000 1.000];
        s.bg_dark   = [0.122 0.122 0.125];   % near-black (logo text dark)
        s.accent    = [0.043 0.482 0.600];   % teal  #0b7b99
        s.accent_lt = [0.851 0.933 0.945];   % teal tint
        s.accent_dk = [0.027 0.322 0.404];   % teal dark
        s.mandatory = [0.663 0.114 0.145];   % deep red #a91d25
        s.success   = [0.180 0.631 0.239];   % keep green for success
        s.text      = [0.082 0.082 0.090];
        s.text_sub  = [0.400 0.400 0.420];
        s.border    = [0.820 0.820 0.831];
        s.font      = 'Helvetica Neue';
        s.fs        = 13;
        s.fs_sm     = 11;
        s.fs_lg     = 15;
        s.fs_xl     = 20;
    end

%% ════════════════════════════════════════════════════════════════════════
%%  UTILITY: STRUCT ↔ WIDGET
%% ════════════════════════════════════════════════════════════════════════

    function cache = build_widget_cache()
        % Returns a containers.Map: tag -> handle for all tagged descendants.
        all_h = findobj(fig, '-not', 'Tag', '');
        cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
        for ki = 1:numel(all_h)
            t = all_h(ki).Tag;
            if ~isempty(t) && ~isKey(cache, t)
                cache(t) = all_h(ki);
            end
        end
    end

    function h = cached_findobj(tag)
        if isKey(widget_cache, tag)
            h = widget_cache(tag);
        else
            h = findobj(fig, 'Tag', tag);
            if ~isempty(h), widget_cache(tag) = h(1); end
        end
    end

    function set_widget(tag, val)
        h = cached_findobj(tag);
        if isempty(h), return; end
        h = h(1);
        try
            if isa(h, 'matlab.ui.control.DropDown')
                if ischar(val) || isstring(val)
                    if any(strcmp(h.Items, char(val)))
                        h.Value = char(val);
                    end
                end
            elseif isa(h, 'matlab.ui.control.CheckBox')
                h.Value = logical(val);
            elseif isa(h, 'matlab.ui.control.NumericEditField')
                if isnumeric(val) && isscalar(val) && ~isnan(val)
                    h.Value = double(val);
                end
            elseif isa(h, 'matlab.ui.control.EditField')
                if ischar(val) || isstring(val)
                    h.Value = char(val);
                elseif isnumeric(val)
                    h.Value = num2str(val);
                end
            end
        catch
        end
    end

    function val = get_widget_value(h)
        if isa(h, 'matlab.ui.control.CheckBox')
            val = h.Value;
        elseif isa(h, 'matlab.ui.control.NumericEditField')
            % Fields created with nedt(…, NaN) use Value=0 as a sentinel
            % for "not set" (NumericEditField cannot store NaN at creation
            % time).  Return NaN so collect_params can skip them, keeping
            % the downstream code's isempty() / isfield() guards intact.
            if strcmp(h.Placeholder, 'NaN') && h.Value == 0
                val = NaN;
            else
                val = h.Value;
            end
        elseif isa(h, 'matlab.ui.control.DropDown')
            val = h.Value;
        else
            val = h.Value;
        end
    end

    function s = set_nested(s, path, val)
        parts = strsplit(path, '.');
        s = rebuild(s, parts, val);
    end

    function s = rebuild(s, parts, val)
        if numel(parts) == 1
            s.(parts{1}) = val;
        else
            if ~isfield(s, parts{1}), s.(parts{1}) = struct(); end
            s.(parts{1}) = rebuild(s.(parts{1}), parts(2:end), val);
        end
    end

    function val = get_nested(s, path)
        try
            parts = strsplit(path, '.');
            val   = s;
            for i = 1:numel(parts)
                val = val.(parts{i});
            end
        catch
            val = [];
        end
    end

    function flat = flatten_struct(s, prefix)
        flat  = struct();
        flds  = fieldnames(s);
        for i = 1:numel(flds)
            key = flds{i};
            val = s.(key);
            full_key = key;
            if ~isempty(prefix), full_key = [prefix '.' key]; end
            safe_key  = strrep(full_key, '.', '__');
            if isstruct(val) && ~isempty(val) && numel(val) == 1
                sub = flatten_struct(val, full_key);
                sub_flds = fieldnames(sub);
                for j = 1:numel(sub_flds)
                    flat.(sub_flds{j}) = sub.(sub_flds{j});
                end
            elseif isnumeric(val) || islogical(val) || ischar(val) || isstring(val)
                if numel(safe_key) <= 63
                    flat.(safe_key) = val;
                end
            end
        end
    end

    function set_medium_table(mp)
        tbl = findobj(fig, 'Tag', 'medium_table');
        if isempty(tbl), return; end
        tissues = {'water','brain','skin','skull','skull_cortical','skull_trabecular'};
        props   = {'sound_speed','density','alpha_coeff','alpha_power', ...
                   'thermal_conductivity','specific_heat_capacity','perfusion','absorption_fraction'};
        data = tbl.Data;
        for ti = 1:numel(tissues)
            if isfield(mp, tissues{ti})
                for pi = 1:numel(props)
                    if isfield(mp.(tissues{ti}), props{pi})
                        data{ti, pi} = mp.(tissues{ti}).(props{pi});
                    end
                end
            end
        end
        tbl.Data = data;
    end

    function params = read_medium_table(params, tbl)
        tissues = {'water','brain','skin','skull','skull_cortical','skull_trabecular'};
        props   = {'sound_speed','density','alpha_coeff','alpha_power', ...
                   'thermal_conductivity','specific_heat_capacity','perfusion','absorption_fraction'};
        for ti = 1:numel(tissues)
            for pi = 1:numel(props)
                val = tbl.Data{ti, pi};
                if isnumeric(val)
                    params.medium_properties.(tissues{ti}).(props{pi}) = val;
                end
            end
        end
    end

    function set_status(text, color)
        h = findobj(fig, 'Tag', 'lbl_status');
        if ~isempty(h)
            h.Text      = text;
            h.FontColor = color;
        end
    end

    function append_log(text)
        la = findobj(fig, 'Tag', 'log_area');
        if ~isempty(la)
            la.Value{end+1} = text;
            try; scroll(la, 'bottom'); catch; end
        end
    end

    function clear_log()
        la = findobj(fig, 'Tag', 'log_area');
        if ~isempty(la)
            la.Value = {''};
        end
    end


    function p = get_prestus_path()
        % Return the PRESTUS root directory (two levels above this file).
        p = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    end

end % prestus_gui
