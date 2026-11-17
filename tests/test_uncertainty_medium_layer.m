classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_uncertainty_medium_layer < matlab.unittest.TestCase
% TEST_UNCERTAINTY_MEDIUM_LAYER  Medium layer of each uncertainty variant: shipped bounds + edits.
%
%   Run with:   results = runtests('tests/test_uncertainty_medium_layer.m');
%
% See also: UNCERTAINTY_MEDIUM_LAYER, UNCERTAINTY_PIPELINE

    methods (Test)

        function noEditsGiveShippedBounds(tc)
            folder = fullfile(get_prestus_path(), 'config', 'uncertainty');
            p = defaults();
            for name = {'temp_high', 'temp_low', 'trans_high', 'trans_low'}
                shipped = yaml.loadFile(fullfile(folder, ['config_medium_' name{1} '.yaml']), ...
                    'ConvertToArray', true);
                [layer, edits] = uncertainty_medium_layer(p, name{1});
                tc.verifyEqual(layer, shipped, name{1});
                tc.verifyEmpty(fieldnames(edits), name{1});
            end
        end

        function editIsMergedOnTop(tc)
            p = defaults();
            p.simulation.uncertainty_overrides = struct('temp_high', ...
                struct('brain', struct('sound_speed', 1600)));
            [layer, edits] = uncertainty_medium_layer(p, 'temp_high');
            tc.verifyEqual(layer.medium_properties.brain.sound_speed, 1600);
            tc.verifyEqual(edits.brain.sound_speed, 1600);
            other = uncertainty_medium_layer(p, 'temp_low');            % other variants unchanged
            tc.verifyNotEqual(other.medium_properties.brain.sound_speed, 1600);
        end

        function editedAbsorptionFieldDropsTheOther(tc)
            p = defaults();
            p.simulation.uncertainty_overrides = struct('temp_high', ...
                struct('brain', struct('absorption_coeff', 0.5)));
            layer = uncertainty_medium_layer(p, 'temp_high');
            tc.verifyEqual(layer.medium_properties.brain.absorption_coeff, 0.5);
            tc.verifyEmpty(layer.medium_properties.brain.absorption_fraction);
        end

        function badEditsError(tc)
            p = defaults();
            bad = {struct('liberal', struct('brain', struct('density', 1))), ...            % not a variant
                   struct('temp_high', struct('brainz', struct('density', 1))), ...         % not a tissue
                   struct('temp_high', struct('brain', struct('densty', 1))), ...           % not a property
                   struct('temp_high', struct('brain', struct('density', 'x'))), ...        % not a number
                   struct('trans_high', struct('brain', struct('perfusion', 1))), ...       % thermal, acoustic only
                   struct('temp_high', struct('brain', struct('absorption_coeff', 1, ...
                                                              'absorption_fraction', 1)))}; % both absorption fields
            for i = 1:numel(bad)
                p.simulation.uncertainty_overrides = bad{i};
                variant = fieldnames(bad{i});           % edits are checked for the variant built
                variant = variant{1};
                if ~any(strcmp(variant, {'temp_high', 'temp_low', 'trans_high', 'trans_low'}))
                    variant = 'temp_high';
                end
                tc.verifyError(@() uncertainty_medium_layer(p, variant), ...
                    'uncertainty_medium_layer:badOverride', sprintf('case %d', i));
            end
        end

        function transFilesListNoThermalProperties(tc)
            % Trans variants run without heating; their bounds are acoustic only
            thermal = {'thermal_conductivity', 'specific_heat_capacity', 'perfusion', ...
                       'absorption_fraction', 'absorption_coeff'};
            for name = {'trans_high', 'trans_low'}
                m = uncertainty_medium_layer(defaults(), name{1}).medium_properties;
                for t = fieldnames(m)'
                    tc.verifyFalse(any(isfield(m.(t{1}), thermal)), sprintf('%s.%s', name{1}, t{1}));
                end
            end
        end

        function shippedFilesSetOneAbsorptionField(tc)
            for name = {'temp_high', 'temp_low'}
                m = uncertainty_medium_layer(defaults(), name{1}).medium_properties;
                for t = fieldnames(m)'
                    tc.verifyEqual(sum(isfield(m.(t{1}), {'absorption_fraction', 'absorption_coeff'})), 1, ...
                        sprintf('%s.%s', name{1}, t{1}));
                end
            end
        end

        function unknownVariantErrors(tc)
            tc.verifyError(@() uncertainty_medium_layer(defaults(), 'default'), ...
                'uncertainty_medium_layer:badVariant');
            tc.verifyError(@() uncertainty_medium_layer(defaults(), 'liberal'), ...
                'uncertainty_medium_layer:badVariant');
        end
    end
end

function p = defaults()
    p = yaml.loadFile(fullfile(get_prestus_path(), 'config', 'config_default.yaml'), 'ConvertToArray', true);
end
