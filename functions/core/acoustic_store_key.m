function [key, shareable] = acoustic_store_key(parameters, transducer_idx, p_acoustic)
% ACOUSTIC_STORE_KEY  Content key of one transducer's acoustic field in an async run
%
% Two transducer entries with the same key have the same acoustic field: the
% async pipeline simulates it once and both read it (ASYNC_TRANSDUCER_PIPELINE).
% With shareable true, the field is also kept under its key in the output
% folder's cache (sub-XXX_<medium>_results_ac-<key>.mat), so a later run with
% the same key, e.g. a follow-up, a repeated study, or a sequence that returns
% to a stance, reads it instead of simulating it again.
%
% The key is computed when the run is planned, before placement, from:
%   - CACHE_SIGNATURE(p_acoustic, 'acoustic'): subject, medium and its
%     properties (so every uncertainty variant has its own key), grid
%     settings incl. the acoustic FOV, head model, pseudo-CT, segmentation
%     content, and this transducer without its drive level;
%   - the acoustic backend and its settings, code type and precision;
%   - with grid.async_layout 'shared', every transducer of the run, because
%     the shared grid is sized and oriented for all of them;
%   - the PRESTUS version (VERSION file).
% Not in the key: the target ISPPA (applied when the fields are combined),
% timing_offset_s, transducer(i).timing and the name (the field is a
% steady-state solution that pulse timing does not change; the BabelBrain
% backend reads the timing, so there it is in the key).
%
% Positions are known when the run is planned only when they come from the
% transducer entries (placement.mode 'manual', or no placement stage). Other
% placement modes resolve each transducer from inputs that depend on its
% index or on files, so the key then includes the index (no two entries
% share a field) and shareable is false. shareable is also false for a
% layered medium whose segmentation is not there yet (its content is part of
% the key).
%
% Use as:
%   [key, shareable] = acoustic_store_key(parameters, transducer_idx, p_acoustic)
%
% Input:
%   parameters     - the async run's parameters (all transducers)
%   transducer_idx - index of the transducer
%   p_acoustic     - that transducer's acoustic sub-run parameters
%
% Output:
%   key       - 16-character hexadecimal key
%   shareable - true when the field may be kept for later runs
%
% See also: ASYNC_TRANSDUCER_PIPELINE, CACHE_SIGNATURE, STRUCT_MD5

    s = struct();
    sig = cache_signature(p_acoustic, 'acoustic');
    if isfield(sig, 'transducer')
        sig.transducer = strip_stance_free(sig.transducer);
    end
    if isfield(sig, 'acoustic_solver')
        sig = rmfield(sig, 'acoustic_solver');   % the solver enters the key below, so kept fields keep their keys
    end
    s.signature = sig;

    backend = 'kwave';
    if isfield(parameters, 'acoustic_backend') && ~isempty(parameters.acoustic_backend)
        backend = char(parameters.acoustic_backend);
    end
    s.backend = backend;
    for f = {'jwave', 'babelbrain', 'elastic'}
        if isfield(parameters, f{1})
            s.(f{1}) = parameters.(f{1});
        end
    end
    s.code_type = field_or(parameters, {'simulation', 'code_type'});
    s.precision = field_or(parameters, {'simulation', 'precision'});
    if strcmpi(backend, 'babelbrain')
        s.timing = field_or(parameters, {'timing'});
        s.transducer_timing = field_or(parameters.transducer(transducer_idx), {'timing'});
    end

    layout = field_or(parameters, {'grid', 'async_layout'});
    if isempty(layout)
        layout = 'shared';
    end
    s.layout = char(layout);
    if strcmp(s.layout, 'shared')
        s.all_transducers = arrayfun(@strip_stance_free_drive, parameters.transducer, ...
                                     'UniformOutput', false);
    end
    s.version = prestus_version_file();

    explicit = placement_is_explicit(parameters);
    if explicit
        s.placement = field_or(parameters, {'placement', 'manual'});
    else
        s.transducer_idx = transducer_idx;   % resolved per index: never shared
    end
    shareable = explicit;
    if strcmp(field_or(parameters, {'simulation', 'medium'}), 'layered') && ~isfield(sig, 'segmentation')
        shareable = false;                   % segmentation not there yet: its content is unknown
    end

    hex = struct_md5(s);
    key = hex(1:16);
end

% -------------------------------------------------------------------------
function t = strip_stance_free(t)
% Fields of a transducer entry that do not change its acoustic field
    for f = {'timing_offset_s', 'timing', 'name'}
        if isstruct(t) && isfield(t, f{1})
            t = rmfield(t, f{1});
        end
    end
end

function t = strip_stance_free_drive(t)
    for f = {'target_isppa_wcm2', 'tx_calibration'}
        if isfield(t, f{1})
            t = rmfield(t, f{1});
        end
    end
    t = strip_stance_free(t);
end

function tf = placement_is_explicit(parameters)
% Positions come from the transducer entries themselves
    run_placement = true;
    if isfield(parameters, 'modules') && isfield(parameters.modules, 'run_transducer_placement')
        run_placement = logical(parameters.modules.run_transducer_placement);
    end
    mode = field_or(parameters, {'placement', 'mode'});
    tf = ~run_placement || isempty(mode) || strcmp(char(mode), 'manual');
end

function v = field_or(s, path)
    v = [];
    for i = 1:numel(path)
        if ~isstruct(s) || ~isfield(s, path{i})
            return
        end
        s = s.(path{i});
    end
    v = s;
end

function v = prestus_version_file()
    v = 'unknown';
    try
        f = fullfile(get_prestus_path, 'VERSION');
        if isfile(f)
            v = strtrim(fileread(f));
        end
    catch
    end
end
