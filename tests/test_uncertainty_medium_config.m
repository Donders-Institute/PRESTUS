classdef test_uncertainty_medium_config < matlab.unittest.TestCase
% TEST_UNCERTAINTY_MEDIUM_CONFIG  Which medium YAML each uncertainty variant uses.
%
%   Run with:   results = runtests('tests/test_uncertainty_medium_config.m');
%
% See also: UNCERTAINTY_MEDIUM_CONFIG, UNCERTAINTY_PIPELINE

    methods (Test)

        function emptyOrAbsentUsesShippedFile(tc)
            folder = fullfile(get_prestus_path(), 'config', 'uncertainty');
            p = struct('simulation', struct('uncertainty_temp_high', ''));
            tc.verifyEqual(uncertainty_medium_config(p, 'temp_high'), ...
                fullfile(folder, 'config_medium_temp_high.yaml'));
            for name = {'temp_high', 'temp_low', 'trans_high', 'trans_low'}
                tc.verifyEqual(uncertainty_medium_config(struct(), name{1}), ...
                    fullfile(folder, ['config_medium_' name{1} '.yaml']));
                tc.verifyTrue(isfile(uncertainty_medium_config(struct(), name{1})), name{1});
            end
        end

        function fileNameIsLookedUpInConfigFolder(tc)
            p = struct('simulation', struct('uncertainty_trans_low', 'my_bounds.yaml'));
            tc.verifyEqual(uncertainty_medium_config(p, 'trans_low'), ...
                fullfile(get_prestus_path(), 'config', 'uncertainty', 'my_bounds.yaml'));
        end

        function absolutePathIsKept(tc)
            f = fullfile(tempdir, 'study_temp_high.yaml');
            p = struct('simulation', struct('uncertainty_temp_high', f));
            tc.verifyEqual(uncertainty_medium_config(p, 'temp_high'), f);
        end

        function shippedDefaultMatchesConfigDefault(tc)
            % config_default.yaml leaves all keys empty, so a plain run uses the shipped files
            p = yaml.loadFile(fullfile(get_prestus_path(), 'config', 'config_default.yaml'));
            for name = {'temp_high', 'temp_low', 'trans_high', 'trans_low'}
                tc.verifyEmpty(char(p.simulation.(['uncertainty_' name{1}])), name{1});
            end
        end

        function transFilesListNoThermalProperties(tc)
            % Trans variants run without heating; their files bound acoustic properties only
            thermal = {'thermal_conductivity', 'specific_heat_capacity', 'perfusion', ...
                       'absorption_fraction', 'absorption_coeff'};
            for name = {'trans_high', 'trans_low'}
                m = yaml.loadFile(uncertainty_medium_config(struct(), name{1})).medium_properties;
                for t = fieldnames(m)'
                    tc.verifyFalse(any(isfield(m.(t{1}), thermal)), sprintf('%s.%s', name{1}, t{1}));
                end
            end
        end

        function shippedFilesSetOneAbsorptionField(tc)
            for name = {'temp_high', 'temp_low'}
                m = yaml.loadFile(uncertainty_medium_config(struct(), name{1})).medium_properties;
                for t = fieldnames(m)'
                    tc.verifyEqual(sum(isfield(m.(t{1}), {'absorption_fraction', 'absorption_coeff'})), 1, ...
                        sprintf('%s.%s', name{1}, t{1}));
                end
            end
        end

        function unknownVariantErrors(tc)
            tc.verifyError(@() uncertainty_medium_config(struct(), 'default'), ...
                'uncertainty_medium_config:badVariant');
            tc.verifyError(@() uncertainty_medium_config(struct(), 'liberal'), ...
                'uncertainty_medium_config:badVariant');
        end
    end
end
