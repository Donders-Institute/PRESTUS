function info = simnibs_install_info(bin_path)
% SIMNIBS_INSTALL_INFO  Version and default charm settings of a SimNIBS install
%
% Reads the SimNIBS install whose bin folder holds charm, without running it:
% the simnibs package that charm calls (from the charm wrapper, or the
% site-packages of the conda environment), its version (simnibs-X.Y.Z.dist-info
% or _version.py) and the charm.ini shipped with it. Mirrors
% nexus/prestus_nexus/simnibs.py.
%
% Use as:
%   info = simnibs_install_info(bin_path)
%
% Input:
%   bin_path - SimNIBS bin folder (parameters.startup.simnibs_bin_path)
%
% Output:
%   info - struct with fields
%          .version          version string, '' if unknown
%          .major            major version, NaN if unknown
%          .package          simnibs package folder, '' if not found
%          .default_settings charm.ini of the install, '' if not found
%          .python           Python interpreter charm runs ('' if not found);
%                            PlanTUS uses it when placement.plantus.env_path is empty
%          .env_root         environment folder holding the simnibs package
%          .ld_library_path  the install's own libraries (<package>/external/lib/linux,
%                            '' if absent): hpc.ld_library_path when that is empty
%
% See also: SIMNIBS_VERSION, SEGMENTATION_RUN, SEGMENTATION_CHECK_EXISTING

arguments
    bin_path (1,:) char
end

    info = struct('version', '', 'major', NaN, 'package', '', 'default_settings', '', ...
        'python', '', 'env_root', '', 'ld_library_path', '');

    charm = '';
    for n = {'charm', 'charm.cmd', 'charm.exe', 'charm.bat'}
        if isfile(fullfile(bin_path, n{1}))
            charm = fullfile(bin_path, n{1});
            break
        end
    end
    if isempty(charm)
        return
    end

    % installer / postinstall links: a small wrapper naming .../site-packages/simnibs/cli/charm.py
    pkg = '';
    txt = '';
    d = dir(charm);
    if d.bytes < 64000
        txt = fileread(charm);
        tok = regexp(txt, '["'']?([^"''\s]*site-packages[/\\]simnibs)[/\\]cli[/\\]charm\.py', 'tokens', 'once');
        if ~isempty(tok) && isfolder(tok{1})
            pkg = tok{1};
        end
    end
    info.python = charm_python(txt, bin_path);
    % conda environment: charm is a console script of the environment itself
    if isempty(pkg)
        env_root = fileparts(bin_path);
        roots = {env_root, fullfile(env_root, 'simnibs_env')};
        for r = roots
            hits = [dir(fullfile(r{1}, 'lib', 'python3*', 'site-packages', 'simnibs')); ...
                    dir(fullfile(r{1}, 'Lib', 'site-packages', 'simnibs'))];
            hits = hits([hits.isdir] & strcmp({hits.name}, '.'));
            if ~isempty(hits)
                pkg = hits(1).folder;
                break
            end
        end
    end
    if isempty(pkg)
        return
    end
    info.package = pkg;

    di = dir(fullfile(fileparts(pkg), 'simnibs-*.dist-info'));
    if ~isempty(di)
        tok = regexp(di(1).name, '^simnibs-([\d.]+[\w.+-]*)\.dist-info$', 'tokens', 'once');
        if ~isempty(tok), info.version = tok{1}; end
    end
    if isempty(info.version) && isfile(fullfile(pkg, '_version.py'))
        tok = regexp(fileread(fullfile(pkg, '_version.py')), '__version__\s*=\s*[''"]([^''"]+)[''"]', 'tokens', 'once');
        if ~isempty(tok), info.version = tok{1}; end
    end
    tok = regexp(info.version, '^(\d+)', 'tokens', 'once');
    if ~isempty(tok), info.major = str2double(tok{1}); end

    if isfile(fullfile(pkg, 'charm.ini'))
        info.default_settings = fullfile(pkg, 'charm.ini');
    end
    if isfolder(fullfile(pkg, 'external', 'lib', 'linux'))
        info.ld_library_path = fullfile(pkg, 'external', 'lib', 'linux');
    end
    % <env>/lib/python3.x/site-packages/simnibs (Windows: <env>/Lib/site-packages/simnibs)
    root = fileparts(fileparts(pkg));
    [up, name] = fileparts(root);
    if startsWith(name, 'python'), root = up; end
    info.env_root = fileparts(root);
    if isempty(info.python)
        for c = {fullfile(info.env_root, 'bin', 'python'), fullfile(info.env_root, 'python.exe')}
            if isfile(c{1}), info.python = c{1}; break; end
        end
    end
end


function py = charm_python(txt, bin_path)
% The interpreter charm runs: the '#!' line of a conda console script, or the
% quoted python of an installer wrapper; else python next to charm.
    py = '';
    lines = splitlines(txt);
    if ~isempty(lines) && startsWith(lines{1}, '#!') && ~contains(lines{1}, '/env ')
        tok = regexp(lines{1}, '^#!\s*(\S+)', 'tokens', 'once');
        if ~isempty(tok) && contains(tok{1}, 'python') && isfile(tok{1}), py = tok{1}; return; end
    end
    tok = regexp(txt, '"([^"]*python[^"]*)"\s+-E', 'tokens', 'once');
    if ~isempty(tok) && isfile(tok{1}), py = tok{1}; return; end
    if isfile(fullfile(bin_path, 'python')), py = fullfile(bin_path, 'python'); end
end
