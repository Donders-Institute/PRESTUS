%% Hypothetical sonication sequence on SimNIBS' Ernie (example_sequence_ernie.yaml)
% Three stances sonicated one after another: left 20 s, right 20 s, left 40 s, with 30 s pauses
% and a 60 s steady-state at the end. Each stance gets its own acoustic simulation; one thermal
% simulation covers the whole sequence, so every stance starts from the heat the earlier ones
% left. See doc/doc_sequential.md ("Planning a Sequence in the Nexus").
%
% Before running: set the four paths in example_sequence_ernie.yaml and run the segmentation
% step of tutorial_ernie.m (the sequence reuses its m2m_sub-001 folder). The stances are given
% in MNI mm and are not snapped to the scalp: open the config in the Nexus (prestus_nexus)
% and check them in the Viewer first.
%
% Not yet run end to end; outputs are named with the affix '_desc-sequence'.

prestus_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(prestus_root, 'functions')));

parameters = load_parameters('example_sequence_ernie.yaml', fullfile(prestus_root, 'examples'));

% The firing schedule (thermal_schedule): starts, durations, and the end of the protocol.
sched = thermal_schedule(parameters, numel(parameters.transducer));
for i = 1:sched.n_sources
    s = sched.source(i);
    fprintf('Stance %d: starts at %g s, sonicates %g s\n', i, s.offset_s, s.ptrd);
end
fprintf('Protocol: %g s of sonication and pauses, %g s in total\n', sched.t_sonication, sched.t_end);

% Run (on a cluster, set platform: slurm in the config: one job per stance, then the joint
% thermal job).
prestus_pipeline_start(parameters);
