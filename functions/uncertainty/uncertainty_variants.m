function variants = uncertainty_variants()
% UNCERTAINTY_VARIANTS  The simulation variants of the uncertainty pipeline
%
% Two pairs of medium-property bounds follow the ITRUSST setups: Temp-H /
% Temp-L take, per property, the bound that maximises / minimises brain
% heating; Trans-H / Trans-L take, per acoustic property, the bound that
% maximises / minimises brain intensity and run without heating. The two
% pairs differ because brain sound speed, density and attenuation and skull
% sound speed raise heating and intensity in opposite directions.
%
% Use as:
%   variants = uncertainty_variants()
%
% Output:
%   variants - 1x5 struct array in run order (default first), with fields
%              .name     'default' | 'temp_high' | 'temp_low' | 'trans_high' | 'trans_low'
%              .label    report label ('Default', 'Temp-H', ...)
%              .desc     output affix appended to the base affix ('' for
%                        default, else '_desc-<name without underscore>')
%              .pair     '' | 'temp' | 'trans'
%              .bound    '' | 'high' | 'low'
%              .heating  true when the variant runs the thermal stages
%                        (the default variant follows modules.run_heating_sims)
%
% See also: UNCERTAINTY_PIPELINE, UNCERTAINTY_MEDIUM_LAYER, GENERATE_UNCERTAINTY_REPORT

names   = {'default', 'temp_high', 'temp_low', 'trans_high', 'trans_low'};
labels  = {'Default', 'Temp-H',    'Temp-L',   'Trans-H',    'Trans-L'};
pairs   = {'',        'temp',      'temp',     'trans',      'trans'};
bounds  = {'',        'high',      'low',      'high',       'low'};
heating = {true,      true,        true,       false,        false};

descs = cellfun(@(n) ['_desc-' strrep(n, '_', '')], names, 'UniformOutput', false);
descs{1} = '';

variants = struct('name', names, 'label', labels, 'desc', descs, ...
                  'pair', pairs, 'bound', bounds, 'heating', heating);
end
