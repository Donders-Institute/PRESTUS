classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_placement_positions < matlab.unittest.TestCase
% TEST_PLACEMENT_POSITIONS  JSON sidecar of resolved positions (WRITE_PLACEMENT_POSITIONS).
%
% The placement stage writes sub-XXX_placement<affix>.json for every mode;
% the Nexus placement viewer reads it to show non-manual placements.
%
%   results = runtests('tests/test_placement_positions.m');

    properties
        Dir
    end


    methods (TestMethodSetup)
        function make_dir(tc)
            tc.Dir = tempname;
            mkdir(tc.Dir);
            tc.addTeardown(@() rmdir(tc.Dir, 's'));
        end
    end

    methods (Static, Access = private)
        function p = params(dir, n)
            p = struct('subject_id', 7);
            p.io = struct('dir_output', dir, 'output_affix', '_desc-A');
            p.transducer = struct('trans_pos', {[10 20 30], [40 50 60]}, ...
                                  'focus_pos', {[11 21 90], [41 51 95]});
            p.transducer = p.transducer(1:n);
        end
    end

    methods (Test, TestTags = {'placement'})

        function test_two_transducers(tc)
            f = write_placement_positions(tc.params(tc.Dir, 2), 'heuristic');
            tc.verifyEqual(f, fullfile(tc.Dir, 'sub-007_placement_desc-A.json'));
            j = jsondecode(fileread(f));
            tc.verifyEqual(j.mode, 'heuristic');
            tc.verifyEqual(j.space, 't1_voxel_1based');
            tc.verifyEqual(numel(j.transducers), 2);
            tc.verifyEqual(j.transducers(2).focus_pos(:)', [41 51 95]);
        end

        function test_single_transducer_is_still_a_list(tc)
            f = write_placement_positions(tc.params(tc.Dir, 1), 'manual');
            txt = fileread(f);
            tc.verifyTrue(~isempty(regexp(txt, '"transducers":\s*\[', 'once')));
        end

        function test_missing_output_dir_writes_nothing(tc)
            p = tc.params(tc.Dir, 1);
            p.io = rmfield(p.io, 'dir_output');
            tc.verifyEqual(write_placement_positions(p, 'manual'), '');
        end

    end

end
