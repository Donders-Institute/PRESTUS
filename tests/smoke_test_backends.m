function smoke_test_backends(varargin)
% SMOKE_TEST_BACKENDS  Fast environment check for BabelBrain and jwave backends.
%
% Verifies that both backends can be imported, key dependencies are available,
% and the PRESTUS bridge scripts are present. Does NOT run a full simulation.
% Completes in <30 s and is safe to call after environment changes (new conda env,
% BabelBrain update, JAX upgrade, etc.).
%
% Use as:
%   smoke_test_backends                % test both backends
%   smoke_test_backends('jwave')       % test jwave only
%   smoke_test_backends('babelbrain')  % test BabelBrain only
%   smoke_test_backends('babelbrain', 'preset', 'CTX_500')  % with specific preset
%
% Environment variables (override auto-detection):
%   JWAVE_PYTHON          CPU Python path / env root for jwave
%   JWAVE_PYTHON_GPU      GPU Python path / env root for jwave (optional)
%   BABELBRAIN_PYTHON     Python path / env root for BabelBrain
%   BABELBRAIN_ROOT       BabelBrain source root (override auto-detect)
%
% Each check is printed with PASS / FAIL / WARN status. A non-zero exit code
% is returned (via error) if any FAIL-level check does not pass.

p = inputParser;
addOptional(p,  'backends',  'both', @ischar);
addParameter(p, 'preset',    'CTX_500', @ischar);
parse(p, varargin{:});

run_bb  = ismember(p.Results.backends, {'both', 'babelbrain'});
run_jw  = ismember(p.Results.backends, {'both', 'jwave'});
preset  = p.Results.preset;

n_fail = 0;
n_pass = 0;
n_warn = 0;

fprintf('\n=== PRESTUS Backend Smoke Test  (%s) ===\n\n', datestr(now, 'yyyy-mm-dd HH:MM:SS')); %#ok<DATST,TNOW1>

%% -------------------------------------------------------------------------
%% 1. Common: locate PRESTUS root and bridge scripts
%% -------------------------------------------------------------------------
section('Common');

prestus_root = get_prestus_path();
check('PRESTUS root accessible', ~isempty(prestus_root) && isfolder(prestus_root));

bb_bridge = fullfile(prestus_root, 'functions', 'babelbrain', 'run_babelbrain_case.py');
check('BabelBrain bridge script exists', isfile(bb_bridge));

bb_thermal_bridge = fullfile(prestus_root, 'functions', 'babelbrain', 'run_babelbrain_thermal.py');
check('BabelBrain thermal bridge exists', isfile(bb_thermal_bridge));

jw_bridge = fullfile(prestus_root, 'functions', 'jwave', 'run_jwave_sim.py');
check('jwave bridge script exists', isfile(jw_bridge));

%% -------------------------------------------------------------------------
%% 2. BabelBrain
%% -------------------------------------------------------------------------
if run_bb
    section('BabelBrain');

    % Locate Python
    bb_py_env = getenv('BABELBRAIN_PYTHON');
    if isempty(bb_py_env)
        bb_py = resolve_python_env('', 'python');
    else
        bb_py = resolve_python_env(bb_py_env, 'python');
    end
    check('Python for BabelBrain resolves', ~isempty(bb_py));

    % Python version
    [st, ver] = system(sprintf('"%s" --version 2>&1', bb_py));
    check_warn('Python version readable', st == 0, strtrim(ver));

    % BabelBrain root
    bb_root_env = getenv('BABELBRAIN_ROOT');
    if ~isempty(bb_root_env) && isfolder(fullfile(bb_root_env, 'TranscranialModeling'))
        bb_root = bb_root_env;
    else
        bb_root = fullfile(prestus_root, 'external', 'BabelBrain');
    end
    check('BabelBrain submodule initialised', isfolder(fullfile(bb_root, 'TranscranialModeling')));

    % Key Python imports
    for pkg = {'numpy', 'h5py', 'yaml', 'nibabel'}
        [st, ~] = system(sprintf('"%s" -c "import %s" 2>/dev/null', bb_py, pkg{1}));
        check(sprintf('Python: import %s', pkg{1}), st == 0);
    end

    % BabelBrain version query
    [st, out] = system(sprintf('"%s" "%s" --get_version --babelbrain_root "%s" 2>/dev/null', ...
        bb_py, bb_bridge, bb_root));
    check('BabelBrain --get_version returns 0', st == 0);
    if st == 0
        ver_line = regexp(out, 'BABELBRAIN_VERSION_TXT=(.+)', 'tokens', 'once');
        if ~isempty(ver_line)
            fprintf('  BabelBrain version: %s\n', strtrim(ver_line{1}));
        end
    end

    % Preset default.yaml
    [~, preset_folder] = babelbrain_transducer_name_to_folder(preset);
    if ~isempty(preset_folder)
        yaml_path = fullfile(bb_root, 'BabelBrain', ['Babel_' preset_folder], 'default.yaml');
        check_warn(sprintf('Preset %s default.yaml exists', preset), isfile(yaml_path), ...
            sprintf('Missing: %s', yaml_path));
    end
end

%% -------------------------------------------------------------------------
%% 3. jwave
%% -------------------------------------------------------------------------
if run_jw
    section('jwave');

    % Locate Python
    jw_py_env = getenv('JWAVE_PYTHON');
    if isempty(jw_py_env)
        jw_py = resolve_python_env('', 'python3');
    else
        jw_py = resolve_python_env(jw_py_env, 'python3');
    end
    check('Python for jwave resolves', ~isempty(jw_py));

    % Python version
    [st, ver] = system(sprintf('"%s" --version 2>&1', jw_py));
    check_warn('Python version readable', st == 0, strtrim(ver));

    % Key Python imports
    for pkg = {'numpy', 'h5py', 'jax', 'jwave', 'jaxdf'}
        [st, ~] = system(sprintf('"%s" -c "import %s" 2>/dev/null', jw_py, pkg{1}));
        check(sprintf('Python: import %s', pkg{1}), st == 0);
    end

    % JAX version and backend
    [st, jax_info] = system(sprintf('"%s" -c "import jax; print(jax.__version__, jax.default_backend())" 2>/dev/null', jw_py));
    if st == 0
        fprintf('  JAX version + backend: %s\n', strtrim(jax_info));
    end

    % Run bridge --help (fast import check)
    [st, ~] = system(sprintf('"%s" "%s" --help 2>/dev/null', jw_py, jw_bridge));
    check('jwave bridge --help returns 0', st == 0);

    % GPU Python env (warn-only if absent)
    jw_gpu_env = getenv('JWAVE_PYTHON_GPU');
    if ~isempty(jw_gpu_env)
        jw_gpu_py = resolve_python_env(jw_gpu_env, 'python3');
        check_warn('GPU Python env resolves', ~isempty(jw_gpu_py), ...
            sprintf('JWAVE_PYTHON_GPU=%s could not be resolved', jw_gpu_env));
    else
        warn('JWAVE_PYTHON_GPU not set (GPU simulations will not be available)');
    end
end

%% -------------------------------------------------------------------------
%% Summary
%% -------------------------------------------------------------------------
fprintf('\n');
section('Summary');
fprintf('  PASS: %d   WARN: %d   FAIL: %d\n\n', n_pass, n_warn, n_fail);
if n_fail > 0
    error('smoke_test_backends: %d check(s) FAILED. See output above.', n_fail);
end

%% =========================================================================
%% Local helpers
%% =========================================================================

    function section(name)
        fprintf('--- %s ---\n', name);
    end

    function check(label, ok)
        if ok
            fprintf('  [PASS] %s\n', label);
            n_pass = n_pass + 1;
        else
            fprintf('  [FAIL] %s\n', label);
            n_fail = n_fail + 1;
        end
    end

    function check_warn(label, ok, detail)
        if ok
            if nargin >= 3 && ~isempty(detail)
                fprintf('  [PASS] %s  (%s)\n', label, detail);
            else
                fprintf('  [PASS] %s\n', label);
            end
            n_pass = n_pass + 1;
        else
            if nargin >= 3 && ~isempty(detail)
                fprintf('  [WARN] %s  (%s)\n', label, detail);
            else
                fprintf('  [WARN] %s\n', label);
            end
            n_warn = n_warn + 1;
        end
    end

    function warn(msg)
        fprintf('  [WARN] %s\n', msg);
        n_warn = n_warn + 1;
    end

end


function [module_path, folder] = babelbrain_transducer_name_to_folder(preset)
% Minimal copy of the preset map for smoke test use.
map = struct( ...
    'CTX_500',  'CTX500', ...
    'CTX_250',  'CTX250', ...
    'R15287',   'R15287', ...
    'R15473',   'R15473', ...
    'R15646',   'R15646', ...
    'DPX_500',  'DPX500', ...
    'Single',   'SingleTx', ...
    'H317',     'H317', ...
    'H246',     'H246', ...
    'H301',     'H301' ...
);
module_path = '';
if isfield(map, preset)
    folder = map.(preset);
else
    folder = '';
end
end
