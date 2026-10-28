classdef test_uncertainty_medium_config < matlab.unittest.TestCase
% TEST_UNCERTAINTY_MEDIUM_CONFIG  Which medium YAML each uncertainty variant uses.
%
%   Run with:   results = runtests('tests/test_uncertainty_medium_config.m');
%
% See also: UNCERTAINTY_MEDIUM_CONFIG, UNCERTAINTY_PIPELINE

    methods (Test)

        function emptyOrAbsentUsesShippedFile(tc)
            folder = fullfile(get_prestus_path(), 'config', 'uncertainty');
            p = struct('simulation', struct('uncertainty_liberal', ''));
            tc.verifyEqual(uncertainty_medium_config(p, 'liberal'), ...
                fullfile(folder, 'config_medium_liberal.yaml'));
            tc.verifyEqual(uncertainty_medium_config(struct(), 'conservative'), ...
                fullfile(folder, 'config_medium_conservative.yaml'));
            tc.verifyTrue(isfile(uncertainty_medium_config(struct(), 'liberal')));
            tc.verifyTrue(isfile(uncertainty_medium_config(struct(), 'conservative')));
        end

        function fileNameIsLookedUpInConfigFolder(tc)
            p = struct('simulation', struct('uncertainty_conservative', 'my_bounds.yaml'));
            tc.verifyEqual(uncertainty_medium_config(p, 'conservative'), ...
                fullfile(get_prestus_path(), 'config', 'uncertainty', 'my_bounds.yaml'));
        end

        function absolutePathIsKept(tc)
            f = fullfile(tempdir, 'study_liberal.yaml');
            p = struct('simulation', struct('uncertainty_liberal', f));
            tc.verifyEqual(uncertainty_medium_config(p, 'liberal'), f);
        end

        function shippedDefaultMatchesConfigDefault(tc)
            % config_default.yaml leaves both keys empty, so a plain run uses the shipped files
            p = yaml.loadFile(fullfile(get_prestus_path(), 'config', 'config_default.yaml'));
            tc.verifyEmpty(char(p.simulation.uncertainty_liberal));
            tc.verifyEmpty(char(p.simulation.uncertainty_conservative));
        end

        function unknownVariantErrors(tc)
            tc.verifyError(@() uncertainty_medium_config(struct(), 'default'), ...
                'uncertainty_medium_config:badVariant');
        end
    end
end
