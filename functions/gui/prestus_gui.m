function prestus_gui()
%PRESTUS_GUI  Interactive parameter setup and simulation launcher.
%
%   Opens a single-window tabbed GUI that:
%     - Loads default_config.yaml and highlights mandatory/optional fields
%     - Allows loading and saving study-specific YAML configs
%     - Submits simulations (local via backgroundPool, or HPC via scheduler)
%     - Displays simulation images and the HTML report inline
%
%   Usage:
%     prestus_gui()
%
%   Requirements: MATLAB R2023b, PRESTUS on the MATLAB path.

%% ── Figure ────────────────────────────────────────────────────────────────

st = sty();

fig = uifigure( ...
    'Name',            'PRESTUS', ...
    'Position',        [80 60 820 680], ...
    'Color',           st.bg_fig, ...
    'Resize',          'on', ...
    'AutoResizeChildren', 'off');

% App state stored in UserData
app.output_dir  = '';
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
tabs.advanced    = uitab(tg, 'Title', '  Advanced      ');
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
build_tab_run(tabs.run);
build_tab_results(tabs.results);

%% ── Wire cross-tab callbacks ──────────────────────────────────────────────

% Thermal timing fields are only meaningful when thermal simulation is on.
% Disable them initially (run_heating_sims defaults to false) and toggle
% whenever the checkbox changes.
h_heating = findobj(fig, 'Tag', 'modules.run_heating_sims');
if ~isempty(h_heating)
    h_heating.ValueChangedFcn = @(cb,~) cb_toggle_thermal_timing(cb.Value);
    cb_toggle_thermal_timing(h_heating.Value);  % apply initial state
end

%% ── Load defaults ─────────────────────────────────────────────────────────

load_defaults();

%% ════════════════════════════════════════════════════════════════════════
%%  TAB BUILDERS
%% ════════════════════════════════════════════════════════════════════════

    %% ── Tab 1: I/O & Paths ────────────────────────────────────────────
    function build_tab_io(t)
        gl = tab_grid(t, 36, {200, '1x', 55, 80});

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

        lbl(gl, 9, 1, 'Subject subfolder');
        chk(gl, 9, 2, 'path.subject_subfolder', true, 'Create per-subject output subfolder');

        lbl(gl, 10, 1, 'T1 filename pattern');
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
    end

    %% ── Tab 2: Simulation ─────────────────────────────────────────────
    function build_tab_simulation(t)
        gl = tab_grid(t, 28, {200, '1x', 55, 80});

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
    end

    %% ── Tab 3: Transducer ─────────────────────────────────────────────
    function build_tab_transducer(t)
        gl = tab_grid(t, 30, {200, '1x', 55, 80});

        sec(gl, 1, 'General');
        lbl(gl, 2, 1, '* Type');
        h_type = drp(gl, 2, 2, 'transducer.type', {'annular','matrix'}, 'annular');
        h_type.ValueChangedFcn = @(dd,~) cb_transducer_type(dd);

        lbl(gl, 3, 1, '* Frequency');
        nedt(gl, 3, 2, 'transducer.freq_hz', 500000);
        lbl(gl, 3, 3, 'Hz');

        lbl(gl, 4, 1, 'Focal distance (exit plane)');
        nedt(gl, 4, 2, 'transducer.focal_distance_ep', NaN);
        lbl(gl, 4, 3, 'mm');

        lbl(gl, 5, 1, 'Focal distance (bowl)');
        nedt(gl, 5, 2, 'transducer.focal_distance_bowl', NaN);
        lbl(gl, 5, 3, 'mm');
        sp(gl, 6);

        sec(gl, 7, 'Transducer Position (T1 voxels)');
        lbl(gl, 8, 1, 'Transducer position');
        xyz_panel(gl, 8, 'transducer.trans_pos', [NaN NaN NaN]);

        lbl(gl, 9, 1, 'Focus position');
        xyz_panel(gl, 9, 'transducer.focus_pos', [NaN NaN NaN]);
        sp(gl, 10);

        % ── Annular panel ─────────────────────────────────────────────
        pnl_ann = uipanel(gl, ...
            'Title',           '', ...
            'BackgroundColor', st.bg_panel, ...
            'BorderType',      'none', ...
            'Tag',             'panel_annular');
        pnl_ann.Layout.Row    = [11 24];
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

        % ── Matrix panel ───────────────────────────────────────────────
        pnl_mat = uipanel(gl, ...
            'Title',           '', ...
            'BackgroundColor', st.bg_panel, ...
            'BorderType',      'none', ...
            'Tag',             'panel_matrix', ...
            'Visible',         'off');
        pnl_mat.Layout.Row    = [11 24];
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
        gl = tab_grid(t, 20, {200, '1x', 55, 80});

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
        gl = tab_grid(t, 56, {220, '1x', 55, 80});

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
        drp(gl, 8, 2, 'headmodel.smooth_method', {'gaussian','box'}, 'gaussian');

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
        chk(gl, 19, 2, 'pct.enabled', false, 'Use CT/pseudo-CT to inform skull properties');

        lbl(gl, 20, 1, 'Density mapping');
        drp(gl, 20, 2, 'pct.mapping_density', {'k-plan','k-wave','marsac','aubry','none'}, 'k-plan');

        lbl(gl, 21, 1, 'Sound speed mapping');
        drp(gl, 21, 2, 'pct.mapping_soundspeed', {'k-plan','marsac','aubry','none'}, 'k-plan');

        lbl(gl, 22, 1, 'Attenuation mapping');
        drp(gl, 22, 2, 'pct.mapping_attenuation', {'k-plan','mueller','aubry','none'}, 'k-plan');
        sp(gl, 23);

        sec(gl, 24, 'Analysis');
        lbl(gl, 25, 1, 'Focus area radius');
        nedt(gl, 25, 2, 'analysis.focus_area_radius', 5); lbl(gl, 25, 3, 'mm');
        note_lbl(gl, 26, 'Radius around focus for ISPPA averaging in output metrics.');
        sp(gl, 27);

        sec(gl, 28, 'Transducer Placement');
        lbl(gl, 29, 1, 'Localite placement');
        chk(gl, 29, 2, 'placement.localite.enabled', false, ...
            'Use position_transducer_localite for placement');

        lbl(gl, 30, 1, 'Localite reference dist.');
        nedt(gl, 30, 2, 'placement.localite.reference_distance_mm', 15); lbl(gl, 30, 3, 'mm');
        note_lbl(gl, 31, 'Correct for distance between IR trackers and transducer exit plane.');

        lbl(gl, 32, 1, 'Save Localite-aligned T1');
        chk(gl, 32, 2, 'placement.heuristic.save_localite_t1', false, ...
            'Save T1 aligned to Localite header for correction');

        lbl(gl, 33, 1, 'Heuristic dist. close');
        nedt(gl, 33, 2, 'placement.heuristic.dist_close', NaN); lbl(gl, 33, 3, 'mm');

        lbl(gl, 34, 1, 'Ear radius');
        nedt(gl, 34, 2, 'placement.heuristic.ear_radius', 35); lbl(gl, 34, 3, 'mm');
        sp(gl, 35);

        sec(gl, 36, 'Simulation Layers');
        note_lbl(gl, 37, 'Comma-separated SimNIBS label indices assigned to each tissue compartment.');

        layer_names = {'water','brain','skin','skull','skull_cortical','skull_trabecular'};
        layer_defaults = {'0,3,6,9,10', '1,2', '5', '4', '7', '8'};
        for i = 1:numel(layer_names)
            lbl(gl, 37+i, 1, layer_names{i});
            edt(gl, 37+i, 2, ['layers.' layer_names{i}], layer_defaults{i}, 0);
        end
    end

    %% ── Tab 9: Run ────────────────────────────────────────────────────
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

    %% ── Tab 10: Results ───────────────────────────────────────────────
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
        if isfield(params, 'io') && isfield(params.io, 'output_dir')
            out_dir = params.io.output_dir;
        end
        app.output_dir = out_dir;
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
            app.output_dir = out_dir;
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

    function cb_browse(target_tag, browse_type)
        if strcmp(browse_type, 'dir')
            p = uigetdir();
            if isequal(p, 0), return; end
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
        if isempty(app.output_dir) || ~isfolder(app.output_dir)
            try
                params = collect_params();
                if isfield(params, 'io') && isfield(params.io, 'output_dir') && ~isempty(params.io.output_dir)
                    app.output_dir = params.io.output_dir;
                    fig.UserData   = app;
                end
            catch; end
        end
        if ~isempty(app.output_dir) && isfolder(app.output_dir)
            try; show_results(app.output_dir); catch; end
        end
    end

    function open_html_report()
        app = fig.UserData;
        if isempty(app.output_dir), return; end
        reports = dir(fullfile(app.output_dir, '*.html'));
        if isempty(reports)
            uialert(fig, 'No HTML report found in output directory.', 'Not found');
            return;
        end
        web(fullfile(app.output_dir, reports(end).name), '-browser');
    end

    function open_output_folder()
        app = fig.UserData;
        if isempty(app.output_dir) || ~isfolder(app.output_dir)
            uialert(fig, 'Output directory not set or does not exist.', 'Not found');
            return;
        end
        if ispc
            system(['explorer "' app.output_dir '"']);
        elseif ismac
            system(['open "' app.output_dir '"']);
        else
            system(['xdg-open "' app.output_dir '"']);
        end
    end

%% ════════════════════════════════════════════════════════════════════════
%%  I/O: LOAD / SAVE / COLLECT
%% ════════════════════════════════════════════════════════════════════════

    function load_defaults()
        try
            prestus_root  = get_prestus_path();
            default_yaml  = fullfile(prestus_root, 'configs', 'default_config.yaml');
            params        = yaml.loadFile(default_yaml, 'ConvertToArray', true);
            apply_params_to_gui(params);
        catch ME
            append_log(sprintf('Could not load defaults: %s', ME.message));
        end
    end

    function load_yaml_to_gui(filepath)
        try
            params = yaml.loadFile(filepath, 'ConvertToArray', true);
            apply_params_to_gui(params);
            set_status(sprintf('Loaded: %s', filepath), st.success);
        catch ME
            fprintf('\n=== PRESTUS GUI Load Error ===\n%s\n==============================\n', ME.getReport('extended'));
            uialert(fig, ME.message, 'Load Error', 'Icon', 'error');
        end
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
        xyz_tags = {'transducer.trans_pos', 'transducer.focus_pos', 'grid.default_dims'};
        for i = 1:numel(xyz_tags)
            tag     = xyz_tags{i};
            safe    = strrep(tag, '.', '__');
            if isfield(flat, safe) && isnumeric(flat.(safe)) && numel(flat.(safe)) >= 3
                for k = 1:3
                    hw = findobj(fig, 'Tag', sprintf('%s_%d', tag, k));
                    if ~isempty(hw), hw.Value = double(flat.(safe)(k)); end
                end
            end
        end
        % Numeric array fields displayed as comma-separated strings
        csv_array_tags = {'transducer.annular.elem_id_mm', ...
                          'transducer.annular.elem_od_mm', ...
                          'transducer.annular.elem_phase_deg', ...
                          'layers.water', 'layers.brain', 'layers.skin', ...
                          'layers.skull', 'layers.skull_cortical', 'layers.skull_trabecular'};
        for i = 1:numel(csv_array_tags)
            tag  = csv_array_tags{i};
            safe = strrep(tag, '.', '__');
            if isfield(flat, safe) && isnumeric(flat.(safe)) && ~isempty(flat.(safe))
                hw = findobj(fig, 'Tag', tag);
                if ~isempty(hw) && isa(hw, 'matlab.ui.control.EditField')
                    hw.Value = strjoin(arrayfun(@num2str, flat.(safe)(:)', 'UniformOutput', false), ',');
                end
            end
        end
        % Medium properties table
        if isfield(params, 'medium_properties')
            set_medium_table(params.medium_properties);
        end
    end

    function params = collect_params()
        params = struct();
        % Collect all tagged widgets
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
            % intentionally fall back to the default_config.yaml value.
            if ischar(val) && isempty(val), continue; end
            params = set_nested(params, tag, val);
        end
        % XYZ triplets stored as separate tagged fields
        xyz_tags = {'transducer.trans_pos', 'transducer.focus_pos', 'grid.default_dims'};
        for i = 1:numel(xyz_tags)
            tag = xyz_tags{i};
            vals = zeros(1,3);
            for ax_i = 1:3
                hw = findobj(fig, 'Tag', sprintf('%s_%d', tag, ax_i));
                if ~isempty(hw), vals(ax_i) = hw.Value; end
            end
            params = set_nested(params, tag, vals);
        end
        % Parse comma-separated fields into numeric arrays
        csv_array_tags = {'transducer.annular.elem_id_mm', ...
                          'transducer.annular.elem_od_mm', ...
                          'transducer.annular.elem_phase_deg', ...
                          'layers.water', 'layers.brain', 'layers.skin', ...
                          'layers.skull', 'layers.skull_cortical', 'layers.skull_trabecular'};
        for i = 1:numel(csv_array_tags)
            tag = csv_array_tags{i};
            hw = findobj(fig, 'Tag', tag);
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
        tbl = findobj(fig, 'Tag', 'medium_table');
        if ~isempty(tbl)
            params = read_medium_table(params, tbl);
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
            if isfield(params.path, 'subject_subfolder') && params.path.subject_subfolder
                params.io.output_dir = fullfile(params.path.sim, sprintf('sub-%03d', sub_id));
            else
                params.io.output_dir = params.path.sim;
            end
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
        for d = 1:3
            pngs = dir(fullfile(output_dir, sprintf('*intensity_%s*.png', dims{d})));
            ax   = findobj(fig, 'Tag', sprintf('ax_acoustic_%d', d));
            if ~isempty(pngs) && ~isempty(ax)
                img = imread(fullfile(output_dir, pngs(end).name));
                imshow(img, 'Parent', ax);
            end
        end
        % Thermal PNG maps
        for d = 1:3
            pngs = dir(fullfile(output_dir, sprintf('*maxT_%s*.png', dims{d})));
            ax   = findobj(fig, 'Tag', sprintf('ax_thermal_%d', d));
            if ~isempty(pngs) && ~isempty(ax)
                img = imread(fullfile(output_dir, pngs(end).name));
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
        csvs  = dir(fullfile(output_dir, '*output_table*.csv'));
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

    function set_widget(tag, val)
        h = findobj(fig, 'Tag', tag);
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
            val = h.Value;
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
