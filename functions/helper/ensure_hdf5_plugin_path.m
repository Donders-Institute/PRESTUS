function ensure_hdf5_plugin_path(python_bin)
% ENSURE_HDF5_PLUGIN_PATH  Set HDF5_PLUGIN_PATH so MATLAB's h5read/h5write don't fail
%
% MATLAB's HDF5 C library searches for filter plugins (e.g. blosc) at a
% compiled-in default path.  When that directory does not exist or is empty,
% h5read raises "required filter 'blosc' is not registered" for files written
% with blosc compression.
%
% This function locates a plugin directory that actually contains the required
% filter libraries and registers it via HDF5_PLUGIN_PATH.  Priority:
%   1. HDF5_PLUGIN_PATH already set, directory exists, and contains *.so/dylib/dll → done
%   2. python_bin supplied: query hdf5plugin Python package for its plugin path
%   3. Default system path (/usr/local/hdf5/lib/plugin) exists with plugins
%   4. Default system path, created if needed (used without a plugin check)
%   5. Fall back to tempdir; warns when it holds no plugins
%
% Called automatically by babelbrain_wrapper (which passes python_bin) and by
% jwave_wrapper / read_babelbrain_output (no python_bin, best-effort).
%
% Use as:
%   ensure_hdf5_plugin_path()
%   ensure_hdf5_plugin_path(python_bin)   % preferred — resolves blosc filter
%
% See also: BABELBRAIN_WRAPPER, JWAVE_WRAPPER, SETENV, GETENV

persistent done
if nargin < 1
    python_bin = '';
end

% If already set and the directory has plugin files, we're done.
current = getenv('HDF5_PLUGIN_PATH');
if ~isempty(current) && isfolder(current) && dir_has_plugins(current)
    done = true;
    return
end

% Re-run if python_bin is now provided (may resolve plugins that weren't found before)
if ~isempty(done) && done && isempty(python_bin)
    return
end

% ---- Try hdf5plugin Python package (most reliable for BabelBrain output) ----
if ~isempty(python_bin) && isfile(python_bin)
    % run_command invokes python directly via argv, no shell at all.
    [status, out] = run_command({python_bin, '-c', 'import hdf5plugin; print(hdf5plugin.FILTERS_PATH)'});
    if status == 0
        plugin_path = strtrim(out);
        if isfolder(plugin_path) && dir_has_plugins(plugin_path)
            setenv('HDF5_PLUGIN_PATH', plugin_path);
            fprintf('[ensure_hdf5_plugin_path] HDF5_PLUGIN_PATH → %s (hdf5plugin)\n', plugin_path);
            done = true;
            return
        end
    end
end

% ---- Try the compiled-in default path ----
default_path = '/usr/local/hdf5/lib/plugin';
if ispc
    default_path = fullfile(getenv('ALLUSERSPROFILE'), 'HDF Group', 'HDF5', 'lib', 'plugin');
end

if isfolder(default_path) && dir_has_plugins(default_path)
    setenv('HDF5_PLUGIN_PATH', default_path);
    done = true;
    return
end

% Try to create it (requires write access to parent)
[ok, ~] = mkdir(default_path);
if ok && isfolder(default_path)
    setenv('HDF5_PLUGIN_PATH', default_path);
    fprintf('[ensure_hdf5_plugin_path] Created HDF5 plugin dir: %s\n', default_path);
    done = true;
    return
end

% ---- Fall back to tempdir ----
fallback = fullfile(tempdir(), 'hdf5_plugins');
if ~isfolder(fallback)
    mkdir(fallback);
end
setenv('HDF5_PLUGIN_PATH', fallback);
fprintf('[ensure_hdf5_plugin_path] HDF5_PLUGIN_PATH set to: %s\n', fallback);
if ~dir_has_plugins(fallback)
    fprintf(['  WARNING: no blosc/filter plugins found in this directory.\n' ...
             '  HDF5 files written with blosc compression (e.g. BabelBrain output)\n' ...
             '  will fail to read.  Install hdf5plugin in the BabelBrain Python\n' ...
             '  environment and set parameters.babelbrain.python_env.\n']);
end
done = true;
end

function tf = dir_has_plugins(dir_path)
% Return true if dir_path contains at least one HDF5 filter plugin file.
d = [dir(fullfile(dir_path, '*.so')); ...
     dir(fullfile(dir_path, '*.dylib')); ...
     dir(fullfile(dir_path, '*.dll'))];
tf = ~isempty(d);
end
