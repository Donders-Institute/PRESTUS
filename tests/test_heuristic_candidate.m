classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_heuristic_candidate < matlab.unittest.TestCase
% TEST_HEURISTIC_CANDIDATE  placement.heuristic.candidate in TP_SELECT_HEURISTIC_POSITION.
%
% With a candidate set (the Nexus Heuristic page), that row of the candidate table is
% placed as it is and the criteria are not applied; [] applies them; an idx that is not
% in the table, or no candidate passing the criteria, stops the run.
%
% Run with:   results = runtests('tests/test_heuristic_candidate.m');

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
        function test_criteria_choose_without_candidate(tc)
            [locs, p, hdr] = setup_case(tc.TmpDir);
            [~, ~, best] = evalc_select(locs, p, hdr);
            tc.verifyEqual(best.idx, 9);            % closest of those under criterion_intersection
        end

        function test_chosen_candidate_is_placed(tc)
            [locs, p, hdr] = setup_case(tc.TmpDir);
            p.placement.heuristic.candidate = 7;    % fails the intersection criterion: used anyway
            [trans, targ, best] = evalc_select(locs, p, hdr);
            tc.verifyEqual(best.idx, 7);
            tc.verifyEqual(trans, [160 137 245]);
            tc.verifyEqual(targ, [117 139 161]);
        end

        function test_unknown_candidate_errors(tc)
            [locs, p, hdr] = setup_case(tc.TmpDir);
            p.placement.heuristic.candidate = 99;
            tc.verifyError(@() evalc_select(locs, p, hdr), 'PRESTUS:heuristic:badCandidate');
        end

        function test_no_candidate_left_errors(tc)
            % Every candidate fails the intersection criterion even after
            % expansion to 0.3, or the table is empty: a clear error, not an
            % index error on an empty table.
            [locs, p, hdr] = setup_case(tc.TmpDir);
            all_fail = locs;
            all_fail.prop_intersect(:) = 0.5;
            tc.verifyError(@() evalc_select(all_fail, p, hdr), 'PRESTUS:heuristic:noCandidate');
            tc.verifyError(@() evalc_select(locs([], :), p, hdr), 'PRESTUS:heuristic:noCandidate');
            % Negative control: one candidate within the expanded criterion is chosen
            one_ok = all_fail;
            one_ok.prop_intersect(2) = 0.2;
            [~, ~, best] = evalc_select(one_ok, p, hdr);
            tc.verifyEqual(best.idx, 7);
        end

        function test_exported_columns_agree(tc)
            % Non-integer candidate coordinates: the returned positions, the
            % table's voxel columns, and the Localite columns use the same
            % nearest voxel (no floor/round mismatch).
            [locs, p, hdr] = setup_case(tc.TmpDir);
            cols = ["trans_x","trans_y","trans_z","targ_x","targ_y","targ_z"];
            locs{:, cols} = locs{:, cols} + [0.6 -0.4 0.5 0.7 0.2 -0.5];
            [trans, targ, best] = evalc_select(locs, p, hdr);
            row = locs{locs.idx == best.idx, cols};
            tc.verifyEqual(trans, round(row(1:3)));
            tc.verifyEqual(targ, round(row(4:6)));
            tc.verifyNotEqual(trans, floor(row(1:3)), 'Negative control: floor would differ');
            tc.verifyEqual([best.trans_x best.trans_y best.trans_z], trans);
            tc.verifyEqual([best.targ_x best.targ_y best.targ_z], targ);
            tc.verifyEqual([best.loc_trans_x_mm best.loc_trans_y_mm best.loc_trans_z_mm], ...
                grid_to_localite_mm(trans, hdr), 'AbsTol', 1e-9);
            tc.verifyEqual([best.loc_targ_x_mm best.loc_targ_y_mm best.loc_targ_z_mm], ...
                grid_to_localite_mm(targ, hdr), 'AbsTol', 1e-9);
        end
    end
end

function [trans, targ, best] = evalc_select(locs, p, hdr)
    [~, trans, targ, best] = evalc('tp_select_heuristic_position(locs, 9, ''M1'', p, hdr)');
end

function [locs, p, hdr] = setup_case(root)
    locs = array2table([ ...
        4 161 138 246 117 139 161 85.2 0.02 3.1 0.4 5.2 0.6; ...
        7 160 137 245 117 139 161 84.9 0.31 2.9 0.5 5.0 0.7; ...
        9 150 140 240 117 139 161 80.1 0.04 3.3 0.2 5.5 0.3], ...
        'VariableNames', ["idx","trans_x","trans_y","trans_z","targ_x","targ_y","targ_z", ...
        "dist_to_target","prop_intersect","mean_dist_skin","var_dist_skin","mean_dist_skull","var_dist_skull"]);
    p = struct();
    p.subject_id = 9;
    p.path.sim = fullfile(root, 'sims');
    p.io.output_affix = '';
    p.placement.mode = 'heuristic';
    p.placement.heuristic = struct('target_name', 'M1', 'criterion_intersection', 0.05);
    nii = fullfile(root, 't1.nii');
    niftiwrite(zeros(200, 200, 260, 'single'), nii);
    hdr = niftiinfo(nii);
end
