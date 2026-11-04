classdef test_pct_uncertainty < matlab.unittest.TestCase
% TEST_PCT_UNCERTAINTY  Pseudo-CT in uncertainty mode: only the default variant
% maps the skull, and the skull properties a pseudo-CT run will model can be
% predicted without simulating.
%
%   Run with:   results = runtests('tests/test_pct_uncertainty.m');
%
% See also: PCT_SKULL_MAPPED, PCT_PREDICT_SKULL_PROPERTIES, UNCERTAINTY_PIPELINE

    properties
        TmpDir
    end

    methods (TestMethodSetup)
        function setup(tc)
            tc.TmpDir = tempname;
            mkdir(tc.TmpDir);
            tc.addTeardown(@() rmdir(tc.TmpDir, 's'));
        end
    end

    methods (Test)

        function homogeneousFlagDisablesMapping(tc)
            p = struct('pct', struct('enabled', 1));
            tc.verifyTrue(pct_skull_mapped(p));
            p.uncertainty.homogeneous_skull = true;
            tc.verifyFalse(pct_skull_mapped(p));
            tc.verifyFalse(pct_skull_mapped(struct('pct', struct('enabled', 0))));
            tc.verifyFalse(pct_skull_mapped(struct()));
        end

        function signatureOnlyMarksHomogeneousPctVariants(tc)
            p = make_minimal_parameters(fullfile(tc.TmpDir, 'sims'));
            p.pct.enabled = 1;
            tc.verifyFalse(isfield(cache_signature(p, 'medium'), 'homogeneous_skull'));
            p.uncertainty.homogeneous_skull = true;
            tc.verifyTrue(cache_signature(p, 'medium').homogeneous_skull);
            p.pct.enabled = 0;
            tc.verifyFalse(isfield(cache_signature(p, 'medium'), 'homogeneous_skull'), ...
                'no pseudo-CT: nothing to mark');
        end

        function predictionMatchesMappingFunctions(tc)
            p = native_pct(tc.TmpDir);
            s = pct_predict_skull_properties(p);
            tc.verifyEqual(s.source, 'native');
            tc.verifyEqual(s.n_voxels, 8);
            tc.verifyTrue(isfile(s.file));
            tc.verifyTrue(s.predicted);

            % same calls as medium_setup on the same skull voxels
            hu = 500 + 100 * (0:7)';
            rho = medium_pct_density(p, zeros(8, 1), hu, (1:8)', char(p.pct.mapping_density));
            c   = medium_pct_soundspeed(p, zeros(8, 1), rho, hu, (1:8)', char(p.pct.mapping_soundspeed));
            tc.verifyEqual(s.properties.density.median, median(rho), 'RelTol', 1e-12);
            tc.verifyEqual(s.properties.sound_speed.max, max(c), 'RelTol', 1e-12);
            tc.verifyEqual(s.properties.pseudoCT_HU.min, 500);
            tc.verifyEqual(s.reference.skull.sound_speed, p.medium_properties.skull.sound_speed);
            tc.verifyEqual(sort(fieldnames(s.bounds)), sort({'temp_high'; 'temp_low'; 'trans_high'; 'trans_low'}));
            tc.verifyEqual(s.bounds.trans_high.sound_speed, 2190);
        end

        function predictionErrorsWithoutPseudoCT(tc)
            p = make_minimal_parameters(fullfile(tc.TmpDir, 'sims'));
            p.subject_id = 1;
            p.path.pct = fullfile(tc.TmpDir, 'nopct');
            tc.verifyError(@() pct_predict_skull_properties(p, write=false), ...
                'pct_predict_skull_properties:notFound');
        end
    end
end

function p = native_pct(tmp)
% pCT folder with a 4x4x4 head: 8 skull voxels (label 4) of 500..1200 HU.
    p = make_minimal_parameters(fullfile(tmp, 'sims'));
    p.subject_id = 1;
    p.simulation.medium = 'layered';
    p.path.pct = fullfile(tmp, 'pct');
    d = fullfile(p.path.pct, 'sub-001');
    mkdir(d);
    hu = zeros(4, 4, 4, 'single');
    mask = zeros(4, 4, 4, 'uint8');
    mask(1:2, 1:2, 1:2) = 4;
    hu(mask == 4) = 500 + 100 * (0:7);
    niftiwrite(hu, fullfile(d, 'pseudoCT'), 'Compressed', true);
    niftiwrite(mask, fullfile(d, 'tissues_mask'), 'Compressed', true);
end
