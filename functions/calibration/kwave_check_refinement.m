function [lib, report, sims] = kwave_check_refinement(lib, prof, options)
% KWAVE_CHECK_REFINEMENT  Check per-depth refinements in free-water k-Wave and correct their drive
%
% The analytical fit (FIT_TRANSDUCER_MODEL, Rayleigh) sets phases and the
% drive that should give the reference ISPPA at each refined depth. This
% function runs, for every refined depth, one free-water k-Wave simulation
% of exactly what a study would get: the entry is written to a scratch
% library and the transducer is resolved from it by serial, driving system,
% depth and target ISPPA (RESOLVE_TRANSDUCER_FROM_SERIAL). In-pipeline ISPPA
% targeting is off, so the simulated drive is the library drive.
%
% Per depth it reports the simulated ISPPA, the agreement of the simulated
% on-axis profile with the measured one (R², FLHM-centre error; same
% definitions as FIT_TRANSDUCER_MODEL) and, with ApplyCorrection, rescales
% the stored drive so the simulation would hit the reference ISPPA:
%   elem_p_amp <- elem_p_amp * sqrt(isppa_ref / isppa_kwave)
% Phases are not changed: k-Wave is a check of the analytical phase fit and
% a correction of its amplitude, not an optimiser.
%
% Use as:
%   [lib, report] = kwave_check_refinement(lib, prof)
%   [lib, report, sims] = kwave_check_refinement(lib, prof, options)
%
% Input:
%   lib     - library struct with refinement.depths_ep_mm and refinement.elem_p_amp
%             (LIBRARY_REFINEMENT_FROM_FITS)
%   prof    - measured profiles (READ_AXIAL_PROFILES) for the comparison
%   options - struct, all optional:
%     .Depths          - depths to check [mm from EP] (default: every refined
%                        depth that is also in prof, ±0.5 mm)
%     .SimConfig       - study/simulation config merged under the check's own
%                        settings: YAML path or struct (grid, code_type,
%                        platform/hpc, ...). Its transducer block is ignored.
%                        Without platform: SLURM/qsub if detected, else local.
%     .WorkDir         - scratch folder (default tempdir/prestus_kwave_check/<key>)
%     .Axisymmetric    - use a 2D axisymmetric grid unless SimConfig sets
%                        grid.axisymmetric (default false: the axisymmetric
%                        free-water grid currently puts the focus far short
%                        of the target; see doc_transducer_library.md)
%     .ApplyCorrection - rescale the stored drives (default true)
%     .NearFieldSkipMm - ignore z < this [mm from EP] (default 10)
%     .EquipParam      - struct from load_equipment_config (default: loaded)
%     .ProgressFcn     - @(k, n, depth_ep_mm) called before each simulation
%
% Output:
%   lib    - entry with, per checked depth, refinement.kwave_isppa_wcm2 (at the
%            uncorrected drive), kwave_r2, kwave_focal_err_mm and (if applied)
%            the corrected elem_p_amp with amp_source 'kwave'
%   report - table per depth: depth_ep_mm, isppa_ref_wcm2, isppa_kwave_wcm2,
%            amp_scale, r2_kwave, focal_err_kwave_mm
%   sims   - struct array per depth: depth_ep_mm, z_ep_mm, I_wcm2 (simulated
%            on-axis intensity at the uncorrected drive)
%
% See also: REFINE_TRANSDUCER_LIBRARY, LIBRARY_REFINEMENT_FROM_FITS,
%           CALIBRATION_GUI, PRESTUS_PIPELINE_START

    if nargin < 3; options = struct(); end
    assert(isfield(lib, 'refinement') && isfield(lib.refinement, 'depths_ep_mm') && ...
        isfield(lib.refinement, 'elem_p_amp'), ...
        'kwave_check_refinement:noRefinement', ...
        'The entry has no per-depth refinement with drives; fit per-depth refinements first.');
    r     = lib.refinement;
    rd    = double(r.depths_ep_mm(:)');
    ramp  = double(r.elem_p_amp(:)');
    I_ref = double(r.isppa_wcm2);

    equip = opt(options, 'EquipParam', []);
    if isempty(equip); equip = load_equipment_config(); end
    tran_serial = char(lib.meta.tran_serial);
    ds_serial   = '';
    if isfield(lib.meta, 'driving_system'); ds_serial = char(lib.meta.driving_system); end
    key = tran_serial;
    if ~isempty(ds_serial); key = [tran_serial '_' ds_serial]; end
    tran = equip.trans.(tran_serial);
    bowl_offset = tran.transducer.annular.curv_radius_mm - tran.transducer.annular.dist_geom_ep_mm;

    in_prof = @(d) any(abs(prof.depths_ep_mm - d) < 0.5);
    depths  = opt(options, 'Depths', rd(arrayfun(in_prof, rd) & isfinite(ramp)));
    work    = opt(options, 'WorkDir', fullfile(tempdir, 'prestus_kwave_check', key));
    skip_mm = opt(options, 'NearFieldSkipMm', 10);
    apply   = opt(options, 'ApplyCorrection', true);
    progress = opt(options, 'ProgressFcn', []);
    base    = sim_config(opt(options, 'SimConfig', []), opt(options, 'Axisymmetric', false));

    lib_dir = fullfile(work, 'library');
    if ~isfolder(lib_dir); mkdir(lib_dir); end
    n = numel(depths);
    rep  = nan(n, 5);
    sims = struct('depth_ep_mm', cell(1, n), 'z_ep_mm', [], 'I_wcm2', []);
    for f = {'kwave_isppa_wcm2', 'kwave_r2', 'kwave_focal_err_mm'}
        if ~isfield(r, f{1}) || numel(r.(f{1})) ~= numel(rd); r.(f{1}) = nan(1, numel(rd)); end
    end
    src = repmat("analytical", 1, numel(rd));
    if isfield(r, 'amp_source') && numel(r.amp_source) == numel(rd); src = string(r.amp_source(:)'); end

    for k = 1:n
        d = depths(k);
        m = find(abs(rd - d) < 0.5, 1);
        assert(~isempty(m) && isfinite(ramp(m)), 'kwave_check_refinement:depth', ...
            'Depth %.2f mm has no refinement drive.', d);
        if ~isempty(progress); progress(k, n, d); end

        % Simulate what a study resolves from the uncorrected entry
        yaml.dumpFile(fullfile(lib_dir, [key '.yaml']), lib);
        s = base;
        s.subject_id = k;
        s.path.sim = fullfile(work, 'sims'); s.path.anat = s.path.sim;
        s.path.seg = s.path.sim;             s.path.data = s.path.sim;
        s.transducer = struct('serial', tran_serial, 'library_path', lib_dir, ...
            'focal_distance_ep', d, 'target_isppa_wcm2', I_ref);
        if ~isempty(ds_serial); s.transducer.combo = struct('ds_serial', ds_serial); end
        sub_dir = fullfile(s.path.sim, sprintf('sub-%03d', k));
        if isfolder(sub_dir); rmdir(sub_dir, 's'); end   % never reuse a cached field
        p = load_parameters(s);
        prestus_pipeline_start(p);

        res = dir(fullfile(sub_dir, 'cache', 'sub-*_water_results*.mat'));
        assert(~isempty(res), 'kwave_check_refinement:noResult', ...
            'No free-water result in %s.', fullfile(sub_dir, 'cache'));
        [z_ep, I_sim] = axial_intensity(load(fullfile(res(1).folder, res(1).name)), bowl_offset);
        sims(k) = struct('depth_ep_mm', d, 'z_ep_mm', z_ep, 'I_wcm2', I_sim);

        % Compare with the measured profile (peak-normalised, near field excluded)
        j   = find(abs(prof.depths_ep_mm - d) < 0.5, 1);
        zm  = prof.z_ep_mm(:);
        sel = zm >= skip_mm & zm >= min(z_ep) & zm <= max(z_ep);
        yt  = prof.I(sel, j);            yt = yt / max(max(yt), eps);
        ym  = interp1(z_ep, I_sim, zm(sel), 'linear');  ym = ym / max(max(ym), eps);
        r2  = 1 - sum((yt - ym).^2) / max(sum((yt - mean(yt)).^2), eps);
        fe  = compute_flhm_center(ym, zm(sel)) - compute_flhm_center(yt, zm(sel));
        I_kw  = max(I_sim(z_ep >= skip_mm));
        scale = sqrt(I_ref / I_kw);
        rep(k, :) = [I_ref, I_kw, scale, r2, fe];

        r.kwave_isppa_wcm2(m) = I_kw; r.kwave_r2(m) = r2; r.kwave_focal_err_mm(m) = fe;
        if apply
            ramp(m) = ramp(m) * scale;
            src(m)  = "kwave";
        end
        r.elem_p_amp = round(ramp); r.amp_source = cellstr(src);
        lib.refinement = r;
        fprintf('k-Wave check %.1f mm: ISPPA %.1f W/cm² (target %.1f) → drive ×%.3f | R² %.3f | focal %+.2f mm\n', ...
            d, I_kw, I_ref, scale, r2, fe);
    end

    lib.refinement.kwave_config = struct('code_type', char(string(get_field(base, {'simulation', 'code_type'}, 'default'))), ...
        'axisymmetric', get_field(base, {'grid', 'axisymmetric'}, 0), ...
        'resolution_mm', get_field(base, {'grid', 'resolution_mm'}, []), ...
        'checked', datestr(now, 'yyyy-mm-dd')); %#ok<TNOW1,DATST>
    report = array2table([depths(:), rep], 'VariableNames', {'depth_ep_mm', 'isppa_ref_wcm2', ...
        'isppa_kwave_wcm2', 'amp_scale', 'r2_kwave', 'focal_err_kwave_mm'});
end

%% ==========================================================================
function s = sim_config(cfg, axisym)
% Simulation settings for a free-water check; the user's config fills the gaps.
    if isempty(cfg)
        s = struct();
    elseif ischar(cfg) || isstring(cfg)
        s = yaml.loadFile(char(cfg), 'ConvertToArray', true);
    else
        s = cfg;
    end
    if isfield(s, 'transducer'); s = rmfield(s, 'transducer'); end
    s.simulation.medium      = 'water';
    s.simulation.interactive = 0;
    if ismac && ~has_field(s, {'simulation', 'code_type'})
        s.simulation.code_type = 'matlab_cpu';   % GPU is not supported by MATLAB on macOS
    end
    s.modules.run_heating_sims   = 0;
    s.modules.run_nifti_creation = 0;
    s.modules.run_water_baseline = 0;            % no post-hoc ISPPA rescaling
    s.io.save_matrices          = 1;
    s.io.save_acoustic_matrices = 1;
    s.hpc.wait_for_job          = true;
    if ~has_field(s, {'platform'}) || strcmp(s.platform, 'auto')
        try; s.platform = hpc_detect_system(); catch; s.platform = 'matlab'; end   % no scheduler: run here
    end
    if ~has_field(s, {'grid', 'axisymmetric'}); s.grid.axisymmetric = double(axisym); end
    if s.grid.axisymmetric && ~has_field(s, {'grid', 'default_dims'})
        dflt = yaml.loadFile(fullfile(get_prestus_path(), 'config', 'config_default.yaml'), 'ConvertToArray', true);
        dims = dflt.grid.default_dims;
        if numel(dims) == 3; dims(2) = []; end
        s.grid.default_dims = dims;
    end
end

function [z_ep, I] = axial_intensity(res, bowl_offset)
% On-axis intensity [W/cm²] from the transducer position along the beam axis.
    p  = res.acoustic_info.parameters;
    pm = gather(res.sensor_data.p_max_all);
    tp = p.transducer(1).trans_pos;
    if ismatrix(pm)
        line = pm(tp(1), tp(2):end);
    else
        line = pm(tp(1), tp(2), tp(3):end);
    end
    w  = p.medium_properties.water;
    I  = double(line(:)).^2 / (2 * w.sound_speed * w.density) * 1e-4;
    z_ep = (0:numel(I) - 1)' * p.grid.resolution_mm - bowl_offset;
end

function tf = has_field(s, path)
    tf = true;
    for k = 1:numel(path)
        if ~isstruct(s) || ~isfield(s, path{k}) || isempty(s.(path{k})); tf = false; return; end
        s = s.(path{k});
    end
end

function v = get_field(s, path, d)
    v = d;
    if has_field(s, path)
        for k = 1:numel(path); s = s.(path{k}); end
        v = s;
    end
end

function v = opt(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end
