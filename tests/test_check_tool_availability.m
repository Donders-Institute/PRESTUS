classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_check_tool_availability < matlab.unittest.TestCase
% TEST_CHECK_TOOL_AVAILABILITY  The SimNIBS part of the startup tool check.
%
% Every setting that needs the m2m/toMNI warps (MNI NIfTI outputs, placement.mode 'mni', manual
% placement in MNI coordinates), for each state of the segmentation: none yet (charm runs first and
% writes the warps), a re-run requested (io.overwrite_simnibs), complete, and an existing
% segmentation without toMNI (the negative control, which must fail). A fake SimNIBS bin folder
% holds the tools; the acoustic stages are off, so k-Wave and the GPU are not checked.
%
% Run with:   results = runtests('tests/test_check_tool_availability.m');

    properties
        root
    end

    methods (TestMethodSetup)
        function make_dirs(tc)
            tc.root = tempname;
            mkdir(fullfile(tc.root, 'bin'));
            for t = {'charm', 'subject2mni', 'mni2subject_coords'}
                fclose(fopen(fullfile(tc.root, 'bin', t{1}), 'w'));
            end
            mkdir(fullfile(tc.root, 'seg'));
            tc.addTeardown(@() rmdir(tc.root, 's'));
        end
    end

    methods (Test)
        function tomni_needed_only_when_charm_does_not_run(tc)
            triggers = {'mni_outputs', 'placement_mni', 'manual_mni'};
            states = {'none', 'rerun', 'complete', 'no_tomni'};
            for i = 1:numel(triggers)
                for j = 1:numel(states)
                    p = tc.params(triggers{i}, states{j});
                    what = sprintf('%s, segmentation %s', triggers{i}, states{j});
                    if strcmp(states{j}, 'no_tomni')
                        tc.verifyError(@() run_quiet(p), 'prestus:missingTools', what);
                        try
                            run_quiet(p);
                        catch err
                            tc.verifySubstring(err.message, 'toMNI', what);
                        end
                    else
                        checks = run_quiet(p);
                        tc.verifyTrue(all([checks.ok]), what);
                        has_tomni = any(strcmp({checks.name}, 'm2m toMNI'));
                        tc.verifyEqual(has_tomni, strcmp(states{j}, 'complete'), what);
                    end
                end
            end
        end

        function no_mni_need_no_tomni_check(tc)
            p = tc.params('', 'no_tomni');
            checks = run_quiet(p);
            tc.verifyFalse(any(strcmp({checks.name}, 'm2m toMNI')));
        end
    end

    methods
        function p = params(tc, trigger, state)
            p = struct();
            p.platform = 'matlab';
            p.subject_id = 1;
            p.simulation.medium = 'layered';
            p.simulation.code_type = 'matlab_cpu';
            p.startup.simnibs_bin_path = fullfile(tc.root, 'bin');
            p.path.seg = fullfile(tc.root, 'seg');
            p.io.overwrite_simnibs = strcmp(state, 'rerun');
            p.io.save_MNI = strcmp(trigger, 'mni_outputs');
            p.analysis.mni_warp_method = 'simnibs';
            p.modules = struct('run_nifti_creation', 1, 'run_transducer_placement', 1, ...
                'run_acoustic_sims', 0, 'run_water_baseline', 0, 'run_posthoc_water_sims', 0);
            p.placement.mode = 'manual';
            p.placement.manual.coordinate_space = 'voxel';
            if strcmp(trigger, 'placement_mni')
                p.placement.mode = 'mni';
            elseif strcmp(trigger, 'manual_mni')
                p.placement.manual.coordinate_space = 'mni';
            end
            m2m = fullfile(p.path.seg, 'm2m_sub-001');
            if exist(m2m, 'dir'), rmdir(m2m, 's'); end
            if ~strcmp(state, 'none')
                mkdir(m2m);
                fclose(fopen(fullfile(m2m, 'final_tissues.nii.gz'), 'w'));
                if strcmp(state, 'complete')
                    mkdir(fullfile(m2m, 'toMNI'));
                end
            end
        end
    end
end

function checks = run_quiet(p)
    [~, checks] = evalc('check_tool_availability(p)');
end
