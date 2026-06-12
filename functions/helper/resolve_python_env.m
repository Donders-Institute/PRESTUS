function python_bin = resolve_python_env(env_path, exe_name)
% RESOLVE_PYTHON_ENV  Locate a Python executable inside a venv/conda environment.
%
%   python_bin = resolve_python_env(env_path)
%   python_bin = resolve_python_env(env_path, exe_name)
%
%   env_path  - path to the conda/venv environment root or its bin/ dir,
%               or '' / [] to use the system executable.
%   exe_name  - executable name to look for (default: 'python3').
%               Pass 'python' for BabelBrain, 'python3' for jwave, etc.
%
%   The function tries, in order:
%     1. <env_path>/bin/<exe_name>   (Linux/macOS env root)
%     2. <env_path>/<exe_name>       (Windows, or caller passed bin/ directly)
%     3. <exe_name>                  (system PATH fallback when env_path empty)

if nargin < 2 || isempty(exe_name)
    exe_name = 'python3';
end

if isempty(env_path)
    python_bin = exe_name;
    return
end

candidate_root = fullfile(env_path, 'bin', exe_name);
candidate_bin  = fullfile(env_path, exe_name);

if exist(candidate_root, 'file')
    python_bin = candidate_root;
elseif exist(candidate_bin, 'file')
    python_bin = candidate_bin;
else
    % Neither found — return the root-layout path so the OS error names the
    % expected location rather than silently falling back to system python.
    python_bin = candidate_root;
end
end
