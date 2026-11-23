function [parameters, options] = legacy_bridge_entry(caller, parameters, options, legacy_options)
% LEGACY_BRIDGE_ENTRY  Accept pipeline calls and parameter structs written for earlier PRESTUS versions
%
% Called first in PRESTUS_PIPELINE_START and PRESTUS_PIPELINE. Study scripts
% written for earlier versions keep running, with a warning that says what to
% change:
%
%   - The call form of the v0.5 development builds, caller(subject_id,
%     parameters[, options]), is accepted: subject_id goes to
%     parameters.subject_id (prestus:legacyCallForm).
%   - Old key names a script sets on the parameters struct after
%     LOAD_PARAMETERS (e.g. parameters.code_type = 'cpp_gpu') are converted
%     like old keys in a YAML file (LEGACY_BRIDGE_STRUCT), also in the
%     follow-up configs of a sequential chain (options.sequential_configs).
%
% A call in the current form with current key names passes through unchanged.
%
% Use as:
%   [parameters, options] = legacy_bridge_entry(caller, parameters, options, legacy_options)
%
% Input:
%   caller         - name of the calling function (for messages)
%   parameters     - first argument the caller received (struct, or subject_id)
%   options        - second argument (options struct, or parameters)
%   legacy_options - third argument ([] when not given; options in the old form)
%
% Output:
%   parameters - parameters struct with current key names and subject_id
%   options    - options struct
%
% See also: LEGACY_BRIDGE_STRUCT, LEGACY_BRIDGE_LAYER, UPGRADE_CONFIG_KEYS,
%           PRESTUS_PIPELINE_START, PRESTUS_PIPELINE

    if nargin < 4
        legacy_options = [];
    end

    if isnumeric(parameters) && isscalar(parameters) && isstruct(options)
        subject_id = double(parameters);
        parameters = options;
        options = legacy_options;
        if isempty(options)
            options = struct();
        end
        if isfield(parameters, 'subject_id') && ~isempty(parameters.subject_id) && ...
                ~isequal(double(parameters.subject_id), subject_id) && parameters.subject_id ~= 0
            warning('prestus:legacyCallForm', ...
                '%s: subject_id argument (%d) differs from parameters.subject_id (%d); using %d as earlier versions did.', ...
                caller, subject_id, parameters.subject_id, subject_id);
        end
        parameters.subject_id = subject_id;
        warning('prestus:legacyCallForm', ...
            ['%s(subject_id, parameters, ...) is the call form of PRESTUS before v0.6.0; it still works. ' ...
             'Update the script to:\n  parameters.subject_id = %d;\n  %s(parameters%s)'], ...
            caller, subject_id, caller, iff(isempty(fieldnames(options)), '', ', options'));
    elseif ~isempty(legacy_options)
        error('prestus:tooManyInputs', ...
            'Too many input arguments: use %s(parameters) or %s(parameters, options).', caller, caller);
    end

    if ~isstruct(parameters) || ~isscalar(parameters)
        error('prestus:invalidParameters', ...
            '%s: parameters must be the struct from load_parameters (got %s).', caller, class(parameters));
    end
    if ~isstruct(options)
        error('prestus:invalidOptions', '%s: options must be a struct (got %s).', caller, class(options));
    end

    parameters = legacy_bridge_struct(parameters, caller);
    if isfield(options, 'sequential_configs') && isstruct(options.sequential_configs)
        fn = fieldnames(options.sequential_configs);
        for i = 1:numel(fn)
            options.sequential_configs.(fn{i}) = legacy_bridge_struct( ...
                options.sequential_configs.(fn{i}), sprintf('%s (sequential_configs.%s)', caller, fn{i}));
        end
    end
end

function out = iff(c, a, b)
    if c; out = a; else; out = b; end
end
