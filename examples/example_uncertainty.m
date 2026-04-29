%% EXAMPLE_UNCERTAINTY  Uncertainty quantification pipeline example
%
% This script demonstrates the PRESTUS uncertainty quantification workflow.
% It runs three parallel acoustic + thermal simulations (default, liberal,
% conservative medium properties) and combines results into a single
% uncertainty HTML report.
%
% Usage:
%   Edit the configuration variables in Section 1, then run the script.
%   On a local MATLAB session all stages run sequentially.
%   On HPC (SLURM/qsub) jobs are submitted with automatic dependencies.
%
% The uncertainty pipeline is activated by setting
%   parameters.simulation.uncertainty = true
% and calling prestus_pipeline as normal. The pipeline detects this flag
% and automatically redirects to the uncertainty workflow.
%
% See also: uncertainty_pipeline, generate_uncertainty_report
% Documentation: documentation/doc_uncertainty.md

clear; close all; clc;

% =========================================================================
%% 1. CONFIGURATION
% =========================================================================

% --- Paths ---
prestus_path  = fileparts(fileparts(mfilename('fullpath')));  % PRESTUS root
config_path   = fullfile(prestus_path, 'config');            % adjust if your config lives elsewhere
config_file   = 'config_study.yaml';                          % your study config

% --- Subject ---
subject_id = 1;

% =========================================================================
%% 2. SETUP
% =========================================================================

addpath(fullfile(prestus_path, 'functions', 'helper'));
safe_addpath(fullfile(prestus_path, 'functions'));
safe_addpath(fullfile(prestus_path, 'external'));

% Load base parameters from your study config.
% Do NOT set io.output_affix here — the pipeline manages it per variant.
parameters = load_parameters(fullfile(config_path, config_file));
parameters.subject_id              = subject_id;
parameters.simulation.medium       = 'layered';
parameters.simulation.uncertainty  = true;   % <-- activates uncertainty mode
parameters.platform                = 'auto'; % 'matlab' | 'slurm' | 'qsub' | 'auto'

% path.sim MUST be an absolute path — the pipeline errors if it is empty.
parameters.path.sim = '/absolute/path/to/sim_outputs';  % <-- set this

% =========================================================================
%% 3. OPTIONAL: OVERRIDE UNCERTAINTY OPTIONS
% =========================================================================
% All options below show their defaults; comment out or remove any line
% to accept the default value.

options = struct();

% Output affixes that identify each simulation variant in file names:
options.affixes.default       = '';
options.affixes.liberal       = '_liberal';
options.affixes.conservative  = '_conservative';

% Medium property overrides (YAML files merged on top of base parameters).
% The built-in files are a reasonable starting point; copy and edit them
% to reflect the uncertainty range relevant to your study population.
options.liberal_config = fullfile(prestus_path, 'config', 'uncertainty', ...
    'config_medium_liberal.yaml');
options.conservative_config = fullfile(prestus_path, 'config', 'uncertainty', ...
    'config_medium_conservative.yaml');

% HPC wall-time limits (ignored when running locally in MATLAB):
options.stage1_timelimit  = '00:30:00';   % preprocessing & source setup
options.sim_timelimit     = '03:00:00';   % acoustic + thermal (per variant)
options.report_timelimit  = '00:30:00';   % uncertainty report generation

% =========================================================================
%% 4. RUN
% =========================================================================

% prestus_pipeline_start detects simulation.uncertainty = true and
% automatically runs the full uncertainty workflow (stages 1–5).
prestus_pipeline_start(parameters, options);
