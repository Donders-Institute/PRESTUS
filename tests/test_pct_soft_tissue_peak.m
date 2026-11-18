classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_pct_soft_tissue_peak < matlab.unittest.TestCase
% TEST_PCT_SOFT_TISSUE_PEAK  The soft-tissue peak is the centre of the modal log-intensity bin.
%
% Synthetic UTE and tissue labels: only WM/GM voxels (labels 1 and 2)
% count, and pCT_soft_tissue_value.txt holds exp() of the centre of the
% histogram bin with most counts (not its lower edge). Pure logic.
%
%   results = runtests('tests/test_pct_soft_tissue_peak.m');

    properties
        Dir
    end

    methods (TestMethodSetup)
        function make_dir(tc)
            tc.Dir = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
        end
    end

    methods (Test)
        function test_peak_is_bin_centre(tc)
            rng(1);
            n = 20 * 20 * 20;
            log_ute = 5 + 0.05 * randn(n, 1);           % soft-tissue mode near log(UTE) = 5
            labels = repmat([1; 2], n / 2, 1);          % all WM/GM
            peak = run_peak(tc, log_ute, labels);

            [counts, edges] = histcounts(log_ute);
            [~, k] = max(counts);
            centre = (edges(k) + edges(k + 1)) / 2;
            % txt stores %.2f of exp(centre): compare in log space with that rounding
            tc.verifyEqual(log(peak), centre, 'AbsTol', 0.006 / exp(centre) + 1e-9);
            tc.verifyGreaterThan(abs(log(peak) - edges(k)), (edges(k + 1) - edges(k)) / 4, ...
                'Peak must not be the lower bin edge');
        end

        function test_only_brain_labels_count(tc)
            % Non-brain voxels (label 5) at a far intensity must not move the peak
            n = 20 * 20 * 20;
            log_ute = 5 * ones(n, 1);
            labels = repmat([1; 5], n / 2, 1);
            log_ute(labels == 5) = 8;
            peak = run_peak(tc, log_ute, labels);
            tc.verifyLessThan(abs(log(peak) - 5), 0.5);
        end
    end
end

function peak = run_peak(tc, log_ute, labels)
    sz = [20 20 20];
    m2m = fullfile(tc.Dir, 'm2m_sub-001');
    pct = fullfile(tc.Dir, 'pct');
    mkdir(m2m); mkdir(pct);
    niftiwrite(single(reshape(exp(log_ute), sz)), fullfile(pct, 'UTE_reg_thr0_corr.nii'), 'Compressed', true);
    niftiwrite(int16(reshape(labels, sz)), fullfile(m2m, 'final_tissues.nii'), 'Compressed', true);
    evalc('pct_soft_tissue_peak(m2m, pct)');
    close all
    peak = str2double(strtrim(fileread(fullfile(pct, 'pCT_soft_tissue_value.txt'))));
end
