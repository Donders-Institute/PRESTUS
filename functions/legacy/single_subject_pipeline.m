function parameters = single_subject_pipeline(subject_id, parameters, options)
% SINGLE_SUBJECT_PIPELINE  Entry point of PRESTUS v0.5.0 and earlier; runs PRESTUS_PIPELINE_START in this MATLAB
%
% Kept so that scripts written for v0.5.0 keep running. Warns
% (prestus:legacyFunction) with the current call. Old parameter names are
% converted (LEGACY_BRIDGE_STRUCT), then the run starts on platform 'matlab',
% as v0.5.0 ran it in the calling MATLAB session.
%
% Use as:
%   parameters = single_subject_pipeline(subject_id, parameters)
%   parameters = single_subject_pipeline(subject_id, parameters, sequential_configs=s)
%
% Input:
%   subject_id - subject number
%   parameters - struct from LOAD_PARAMETERS (old key names allowed)
%   options.sequential_configs - follow-up configs of a sequential chain
%   options.adopted_heatmap    - not supported: a v0.5.0 in-memory temperature map;
%                                set io.adopted_heatmap to a T1w-space NIfTI instead
%
% Output:
%   parameters - the parameters the run started with (current key names)
%
% See also: PRESTUS_PIPELINE_START, LEGACY_BRIDGE_ENTRY

    arguments
        subject_id (1,1) double
        parameters struct
        options.adopted_heatmap = []
        options.sequential_configs struct = struct()
    end
    if ~isempty(options.adopted_heatmap)
        error('prestus:legacyFunction', ...
            ['single_subject_pipeline: adopted_heatmap (an in-memory temperature map) is no longer ' ...
             'supported. Set parameters.io.adopted_heatmap to the T1w-space NIfTI of the previous ' ...
             'run''s end temperature, or chain the runs with options.sequential_configs.']);
    end
    warning('prestus:legacyFunction', ...
        ['single_subject_pipeline is the entry point of PRESTUS v0.5.0; it still works. Update the script to:\n' ...
         '  parameters.subject_id = %d;\n  parameters.platform = ''matlab'';\n  prestus_pipeline_start(parameters%s)'], ...
        subject_id, seq_text(options.sequential_configs));

    parameters.subject_id = subject_id;
    parameters = legacy_bridge_struct(parameters, 'single_subject_pipeline');
    parameters.platform = 'matlab';
    prestus_pipeline_start(parameters, seq_options(options.sequential_configs));
end

function o = seq_options(s)
    o = struct();
    if ~isempty(fieldnames(s)); o.sequential_configs = s; end
end

function t = seq_text(s)
    t = '';
    if ~isempty(fieldnames(s)); t = ', struct(''sequential_configs'', sequential_configs)'; end
end
