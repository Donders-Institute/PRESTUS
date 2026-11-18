%% EXAMPLE_MIXED_UNCERTAINTY_ERNIE  Mixed sonication protocol with uncertainty mode on Ernie
% T1 alone, then T2 + T3 firing together on a deep target, then T1's stance again, simulated
% as one async run per tissue-property variant (example_mixed_uncertainty_ernie.yaml).
% See doc/doc_multitransducer_async.md, "Example Protocols".
%
% Before running: set the paths in tutorial_ernie.yaml and run its segmentation step
% (section 1 of tutorial_ernie.m). For your own image, change path.* there; nothing here.
% On a cluster (platform: slurm in the YAML) this submits every job with its dependencies
% and returns. Not yet run end to end.

prestus_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(prestus_root, 'functions')));
addpath(genpath(fullfile(prestus_root, 'config')));
addpath(genpath(fullfile(prestus_root, 'external')));

cfg_dir = fullfile(prestus_root, 'examples');
parameters = load_parameters('tutorial_ernie.yaml', cfg_dir, ...
                             'example_mixed_uncertainty_ernie.yaml', cfg_dir, {'subject_id', 1});

%% 1. What fires when, and which transducers fire together
fg = firing_groups(parameters);
for g = 1:fg.n_groups
    fprintf('Group %d: transducer(s) %s\n', g, mat2str(fg.members{g}));
end
sched = thermal_schedule(parameters, numel(fg.source_tx), ...
    SourceTx=fg.source_tx, SourceFreq=fg.source_freq);
for i = 1:sched.n_sources
    fprintf('Heat source %d (windows of transducer %d): starts at %g s, sonicates %g s\n', ...
        i, fg.source_tx(i), sched.source(i).offset_s, sched.source(i).ptrd);
end
fprintf('Protocol: %g s of sonication and pauses, %g s in total\n', sched.t_sonication, sched.t_end);

%% 2. Run: five variants, each with its own acoustic, combine and joint jobs, then the report
prestus_pipeline_start(parameters);
