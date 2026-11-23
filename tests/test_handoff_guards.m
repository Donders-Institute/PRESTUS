classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_handoff_guards < matlab.unittest.TestCase
% TEST_HANDOFF_GUARDS  Input checks that stop a run before it writes wrong
% outputs: ISPPA sweep targets that share an output affix, and thermal
% hand-off maps (io.adopted_*) that name a missing file. Each guard has a
% negative control: valid input must get past it.

    properties (TestParameter)
        adopted = {'adopted_heatmap', 'adopted_cem43', 'adopted_cem43_iso'};
    end

    methods (Test)

        function sweepTargetsRoundingTogetherError(tc)
            p = struct('transducer', struct('target_isppa_wcm2', [10.2 10.4]));
            tc.verifyError(@() multi_isppa_pipeline(p), 'multi_isppa_pipeline:duplicateAffix');
        end

        function distinctSweepTargetsPassGuard(tc)
            % 10.4 and 10.6 round to 10 and 11: the run proceeds past the
            % guard and stops later on the missing fields of this stub
            p = struct('transducer', struct('target_isppa_wcm2', [10.4 10.6]));
            id = caught_id(@() multi_isppa_pipeline(p));
            tc.verifyNotEqual(id, 'multi_isppa_pipeline:duplicateAffix');
        end

        function tooFewSweepTargetsError(tc)
            p = struct('transducer', struct('target_isppa_wcm2', 10));
            tc.verifyError(@() multi_isppa_pipeline(p), 'multi_isppa_pipeline:tooFewTargets');
        end

        function missingAdoptedMapErrors(tc, adopted)
            p.io.(adopted) = fullfile(tempdir, 'prestus_no_such_map_T1w.nii.gz');
            tc.verifyError(@() call_thermal(p), 'thermal_simulation:adoptedMissing');
        end

        function existingOrEmptyAdoptedMapPassesGuard(tc, adopted)
            f = [tempname '.nii.gz'];
            fid = fopen(f, 'w'); fclose(fid);
            tc.addTeardown(@() delete(f));
            for value = {f, ''}
                p = struct('io', struct(adopted, value{1}));
                tc.verifyNotEqual(caught_id(@() call_thermal(p)), 'thermal_simulation:adoptedMissing');
            end
        end

    end
end

function call_thermal(p)
% Stub inputs: the guard runs before any of them is used
thermal_simulation(p, struct(), 0, struct(), struct(), struct(), [], []);
end

function id = caught_id(fn)
id = '';
try
    fn();
catch err
    id = err.identifier;
end
end
