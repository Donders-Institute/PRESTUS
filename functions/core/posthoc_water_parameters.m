function water_parameters = posthoc_water_parameters(parameters, it_water)
% POSTHOC_WATER_PARAMETERS  Parameters of one post-hoc free-water re-entry run
%
% Builds the config PRESTUS_PIPELINE passes to PRESTUS_PIPELINE_START for the
% post-hoc water simulation (modules.run_posthoc_water_sims): the same run in
% simulation.medium 'water' on the run's grid, without heating and without a
% further post-hoc water run. Everything else, including acoustic_backend and
% thermal_backend, is inherited, so a BabelBrain run re-enters as a BabelBrain
% water-only run (babelbrain_wrapper --water_only).
%
% Use as:
%   water_parameters = posthoc_water_parameters(parameters, it_water)
%
% Input:
%   parameters - PRESTUS config of the finished run (grid.dims set)
%   it_water   - transducer index; with more than one transducer each gets
%                its own water run with a trans_suffix output affix
%
% Output:
%   water_parameters - config of the water run
%
% See also: PRESTUS_PIPELINE, PRESTUS_PIPELINE_START, TRANS_SUFFIX

    n_tx = numel(parameters.transducer);
    water_parameters = parameters;
    water_parameters.simulation.medium              = 'water';
    water_parameters.modules.run_heating_sims       = 0;
    water_parameters.modules.run_posthoc_water_sims = 0;
    water_parameters.simulation.debug               = 0;
    water_parameters.grid.default_dims              = water_parameters.grid.dims;
    water_parameters.io.log_file                    = '';
    % Water runs have their own HPC limits (hpc.posthoc_water_*)
    water_parameters.hpc.timelimit    = posthoc_water_limit(parameters, 'timelimit',   '05:00:00');
    water_parameters.hpc.memorylimit  = posthoc_water_limit(parameters, 'memorylimit', 40);
    water_parameters.hpc.wait_for_job = false;
    if n_tx > 1
        water_parameters.transducer      = parameters.transducer(it_water:it_water);
        water_parameters.io.output_affix = [parameters.io.output_affix, ...
            trans_suffix(it_water, n_tx)];
    end
end

function v = posthoc_water_limit(parameters, field, fallback)
% hpc.posthoc_water_<field> when set, else the built-in water-run default.
    key = ['posthoc_water_' field];
    v = fallback;
    if isfield(parameters, 'hpc') && isfield(parameters.hpc, key) && ~isempty(parameters.hpc.(key))
        v = parameters.hpc.(key);
    end
end
