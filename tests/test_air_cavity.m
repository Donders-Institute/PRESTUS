classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_air_cavity < matlab.unittest.TestCase
% TEST_AIR_CAVITY  Air-cavity beam-path check (headmodel.air_cavities).
%
% AIR_CAVITY_CHARM on a synthetic charm labelling: the label selection and
% the mapping onto the grid. AIR_CAVITY_MASK on synthetic label and HU volumes: every source, the
% exterior-air and grid-edge rules, the majority filter, 2D grids, and the
% cases without a mask. AIR_CAVITY_METRICS on a synthetic Gaussian beam: a
% cavity on the beam axis is flagged, one beside the beam is not (negative
% control), the result is the same for every beam direction along each grid
% axis and on a diagonal, any of several transducers can trigger the flag,
% and each threshold flags on its own. Config validation in LOAD_PARAMETERS
% and the cache signature rule.
%
% Run with:   results = runtests('tests/test_air_cavity.m');

    methods (Test)

        %% ---- AIR_CAVITY_MASK ----------------------------------------------

        function source_none_gives_no_mask(tc)
            p = base_params('none');
            [m, info] = air_cavity_mask(p, ones(10, 10, 10), []);
            tc.verifyEmpty(m);
            tc.verifyEqual(info.source, 'none');
        end

        function segmentation_source_takes_listed_labels(tc)
            seg = randi([0 12], 20, 20, 20);
            p = base_params('segmentation');
            for labels = {11, [11 12]}
                p.headmodel.air_cavities.labels = labels{1};
                [m, info] = air_cavity_mask(p, seg, []);
                tc.verifyEqual(m, ismember(seg, labels{1}));
                tc.verifyEqual(info.n_vox, nnz(ismember(seg, labels{1})));
            end
        end

        function non_layered_medium_gives_no_mask(tc)
            for medium = {'water', 'phantom'}
                p = base_params('segmentation');
                p.headmodel.air_cavities.labels = 11;
                p.simulation.medium = medium{1};
                [m, info] = air_cavity_mask(p, 11 * ones(5, 5, 5), []);
                tc.verifyEmpty(m, medium{1});
                tc.verifyNotEmpty(info.note);
            end
        end

        function simnibs_source_uses_charm_mask(tc)
            seg = zeros(8, 8, 8);
            charm = false(8, 8, 8); charm(3:4, 3:4, 3:4) = true;
            p = base_params('simnibs');
            [m, info] = air_cavity_mask(p, seg, [], charm);
            tc.verifyEqual(m, charm);
            tc.verifyEqual(info.n_vox, 8);
            tc.verifyEqual(info.n_components, 1);
            for bad = {[], false(7, 8, 8)}             % absent or on another grid
                [m, info] = air_cavity_mask(p, seg, [], bad{1});
                tc.verifyEmpty(m);
                tc.verifyNotEmpty(info.note);
            end
            % the charm mask is ignored by every other source (negative control)
            [m, info] = air_cavity_mask(base_params('none'), seg, [], charm);
            tc.verifyEmpty(m); tc.verifyEqual(info.n_vox, 0);
        end

        function charm_labelling_is_mapped_to_the_grid(tc)
            root = tempname; tc.addTeardown(@() rmdir(root, 's'));
            seg_dir = fullfile(root, 'm2m_sub-007', 'segmentation');
            mkdir(seg_dir);
            lab = int16(zeros(12, 12, 12));
            lab(4:5, 4:5, 4:5) = 501;                  % Air-Internal
            lab(9, 9, 9) = 509;                        % another label (mucosa): not air
            niftiwrite(lab, fullfile(seg_dir, 'labeling.nii'), 'Compressed', true);
            p = base_params('simnibs');
            p.headmodel.air_cavities.simnibs_labels = 501;
            p.path.seg = root;
            p.subject_id = 7;
            M = eye(4); M(4, 1:3) = [2 0 -1];           % grid voxel = T1 voxel + [2 0 -1]
            planimg = struct('transf', M);
            m = air_cavity_charm(p, planimg, [12 12 12]);
            expected = false(12, 12, 12); expected(6:7, 4:5, 3:4) = true;
            tc.verifyEqual(m, expected);
            % several labels
            p.headmodel.air_cavities.simnibs_labels = [501 509];
            m = air_cavity_charm(p, planimg, [12 12 12]);
            tc.verifyEqual(nnz(m), 9);
            tc.verifyTrue(m(11, 9, 8));
            % no mask: other source, missing file, non-layered, 2D grid
            [m, note] = air_cavity_charm(base_params('hu'), planimg, [12 12 12]);
            tc.verifyEmpty(m); tc.verifyEmpty(note);
            q = p; q.subject_id = 8;                    % labelling removed: warns
            [m, note] = tc.verifyWarning(@() air_cavity_charm(q, planimg, [12 12 12]), ...
                'prestus:airCavitiesUnavailable');
            tc.verifyEmpty(m); tc.verifySubstring(note, 'not found');
            for medium = {'water', 'phantom'}               % nothing to check: silent
                q = p; q.simulation.medium = medium{1}; q.subject_id = 8;
                [m, note] = tc.verifyWarningFree(@() air_cavity_charm(q, planimg, [12 12 12]));
                tc.verifyEmpty(m); tc.verifyNotEmpty(note);
            end
            [m, note] = air_cavity_charm(p, planimg, [12 12]);
            tc.verifyEmpty(m); tc.verifyNotEmpty(note);
        end

        function hu_source_keeps_internal_air_only(tc)
            [ct, inner, edge_pocket, c_inner] = synthetic_ct();
            p = base_params('hu');
            m = air_cavity_mask(p, zeros(size(ct)), ct);
            % The majority filter trims pocket corners; the core stays
            tc.verifyTrue(m(c_inner(1), c_inner(2), c_inner(3)), 'internal pocket kept');
            tc.verifyFalse(any(m(~inner)), 'nothing outside the internal pocket');
            tc.verifyFalse(any(m(edge_pocket)), 'pocket at the grid edge removed');
            tc.verifyFalse(m(1, 1, 1), 'air around the head removed');
        end

        function hu_keep_boundary_keeps_edge_pockets(tc)
            [ct, ~, edge_pocket, c_inner, c_edge] = synthetic_ct();
            p = base_params('hu');
            p.headmodel.air_cavities.keep_boundary = true;
            m = air_cavity_mask(p, zeros(size(ct)), ct);
            tc.verifyTrue(m(c_inner(1), c_inner(2), c_inner(3)));
            tc.verifyTrue(m(c_edge(1), c_edge(2), c_edge(3)), 'edge pocket kept');
            tc.verifyFalse(any(m(:) & ~edge_pocket(:) & ct(:) > 0), 'tissue never air');
            tc.verifyFalse(m(1, 1, 1), 'largest component (air around the head) still removed');
        end

        function hu_majority_filter_drops_single_voxels(tc)
            ct = 40 * ones(31, 31, 31);                 % all tissue
            ct(16, 16, 16) = -1000;                     % one isolated air voxel
            ct(5:10, 5:10, 5:10) = -1000;               % a real pocket (largest: removed)
            ct(20:24, 20:24, 20:24) = -1000;            % a smaller real pocket: kept
            p = base_params('hu');
            m = air_cavity_mask(p, zeros(size(ct)), ct);
            tc.verifyFalse(m(16, 16, 16));
            tc.verifyTrue(m(22, 22, 22));
        end

        function hu_window_and_nan_padding(tc)
            [ct, ~, ~, c_inner] = synthetic_ct();
            ct(1:3, :, :) = NaN;                        % padding is never air
            p = base_params('hu');
            p.headmodel.air_cavities.hu_window = [-2000, -1500];   % excludes -1000
            m = air_cavity_mask(p, zeros(size(ct)), ct);
            tc.verifyEqual(nnz(m), 0);
            p.headmodel.air_cavities.hu_window = [-1200, -400];
            m = air_cavity_mask(p, zeros(size(ct)), ct);
            tc.verifyTrue(m(c_inner(1), c_inner(2), c_inner(3)));
        end

        function hu_without_ct_gives_no_mask(tc)
            p = base_params('hu');
            for ct = {[], zeros(4, 4, 4)}               % missing, or on another grid
                [m, info] = air_cavity_mask(p, zeros(5, 5, 5), ct{1});
                tc.verifyEmpty(m);
                tc.verifyNotEmpty(info.note);
            end
        end

        function works_on_2d_grids(tc)
            ct = 40 * ones(41, 41);
            ct(1:3, :) = -1000; ct(end-2:end, :) = -1000;   % exterior air at two edges
            ct(:, 1:3) = -1000; ct(:, end-2:end) = -1000;
            ct(18:24, 18:24) = -1000;                       % internal pocket
            p = base_params('hu');
            m = air_cavity_mask(p, zeros(size(ct)), ct);
            tc.verifyEqual(size(m), size(ct));
            tc.verifyTrue(m(21, 21));
            tc.verifyFalse(m(1, 1));
            seg = zeros(41, 41); seg(18:24, 18:24) = 11;
            p = base_params('segmentation'); p.headmodel.air_cavities.labels = 11;
            tc.verifyEqual(air_cavity_mask(p, seg, []), seg == 11);
        end

        %% ---- AIR_CAVITY_METRICS -------------------------------------------

        function cavity_on_axis_is_flagged(tc)
            [pr, beam, cav_on, ~] = beam_scene(3);
            m = air_cavity_metrics(pr, cav_on, 1, beam, 1, -6, 0.1);
            tc.verifyGreaterThan(m.beam_frac_max, 0.5);
            tc.verifyGreaterThan(m.rel_db, -6);
            tc.verifyEqual(m.flag, 1);
            tc.verifyGreaterThan(m.beam_frac_vol, 0);
            % Depth from the exit plane: cavity slab z = 20:22, exit plane at z = 5 + 3
            tc.verifyGreaterThanOrEqual(m.beam_frac_depth_mm, 20 - 8 - 1);
            tc.verifyLessThanOrEqual(m.beam_frac_depth_mm, 22 - 8 + 1);
        end

        function cavity_beside_beam_is_not_flagged(tc)
            [pr, beam, ~, cav_off] = beam_scene(3);
            m = air_cavity_metrics(pr, cav_off, 1, beam, 1, -6, 0.1);
            tc.verifyEqual(m.beam_frac_max, 0);
            tc.verifyLessThan(m.rel_db, -6);
            tc.verifyEqual(m.flag, 0);
        end

        function every_beam_direction_gives_the_same_result(tc)
            [pr, beam, cav_on] = beam_scene(3);
            ref = air_cavity_metrics(pr, cav_on, 1, beam, 1, -6, 0.1);
            % Beam along each grid axis, in both directions
            for perm = {[1 2 3], [3 1 2], [2 3 1]}
                for flip = [false true]
                    q  = permute(pr, perm{1});
                    cv = permute(cav_on, perm{1});
                    b = beam;
                    b.trans_pos = beam.trans_pos(perm{1});
                    b.focus_pos = beam.focus_pos(perm{1});
                    if flip
                        ax = find(perm{1} == 3);
                        q = flip_dim(q, ax); cv = flip_dim(cv, ax);
                        n = size(q, ax);
                        b.trans_pos(ax) = n + 1 - b.trans_pos(ax);
                        b.focus_pos(ax) = n + 1 - b.focus_pos(ax);
                    end
                    m = air_cavity_metrics(q, cv, 1, b, 1, -6, 0.1);
                    label = sprintf('perm [%s], flip %d', num2str(perm{1}), flip);
                    tc.verifyEqual(m.beam_frac_max, ref.beam_frac_max, 'AbsTol', 1e-12, label);
                    tc.verifyEqual(m.rel_db, ref.rel_db, 'AbsTol', 1e-12, label);
                    tc.verifyEqual(m.flag, ref.flag, label);
                end
            end
        end

        function diagonal_beam_flags_cavity_on_its_axis(tc)
            n = 41;
            [x, y, z] = ndgrid(1:n, 1:n, 1:n);
            tp = [6 6 6]; fp = [34 34 34];
            u = (fp - tp) / norm(fp - tp);
            t = (x - tp(1)) * u(1) + (y - tp(2)) * u(2) + (z - tp(3)) * u(3);
            r2 = (x - tp(1)).^2 + (y - tp(2)).^2 + (z - tp(3)).^2 - t.^2;
            pr = exp(-r2 / (2 * 3^2));
            on  = r2 <= 4 & t >= 18 & t <= 21;
            off = r2 >= 15^2 & t >= 18 & t <= 21;
            beam = struct('trans_pos', tp, 'focus_pos', fp, 'exit_offset_vox', 2, 'radius_vox', 10);
            m_on  = air_cavity_metrics(pr, on, 1, beam, 1, -6, 0.1);
            m_off = air_cavity_metrics(pr, off, 1, beam, 1, -6, 0.1);
            tc.verifyEqual(m_on.flag, 1);
            tc.verifyEqual(m_off.flag, 0);
        end

        function any_transducer_can_trigger_the_flag(tc)
            [pr, beam, cav_on] = beam_scene(3);
            elsewhere = beam;                           % a second beam along x, far from the cavity
            elsewhere.trans_pos = [1 8 36]; elsewhere.focus_pos = [31 8 36];
            elsewhere.radius_vox = 5;
            for order = {[beam, elsewhere], [elsewhere, beam]}
                m = air_cavity_metrics(pr, cav_on, 1, order{1}, 1, -6, 0.1);
                tc.verifyEqual(m.flag, 1);
                tc.verifyGreaterThan(m.beam_frac_max, 0.5);
            end
            m = air_cavity_metrics(pr, cav_on, 1, elsewhere, 1, -6, 0.1);
            tc.verifyEqual(m.beam_frac_max, 0, 'second beam alone misses the cavity');
        end

        function each_threshold_flags_on_its_own(tc)
            [pr, beam, cav_on] = beam_scene(3);
            m = air_cavity_metrics(pr, cav_on, 1, beam, 1, Inf, 0.1);     % pressure test off
            tc.verifyEqual(m.flag, 1, 'beam share alone');
            m = air_cavity_metrics(pr, cav_on, 1, beam, 1, -6, 1);        % share test off
            tc.verifyEqual(m.flag, 1, 'pressure alone');
            m = air_cavity_metrics(pr, cav_on, 1, beam, 1, Inf, 1);
            tc.verifyEqual(m.flag, 0, 'neither');
        end

        function empty_mask_reports_nothing(tc)
            [pr, beam] = beam_scene(3);
            m = air_cavity_metrics(pr, false(size(pr)), 1, beam, 1, -6, 0.1);
            tc.verifyEqual(m.n_vox, 0);
            tc.verifyEqual(m.flag, 0);
            tc.verifyTrue(isnan(m.rel_db));
        end

        function size_mismatch_errors(tc)
            [pr, beam] = beam_scene(3);
            tc.verifyError(@() air_cavity_metrics(pr, false(5, 5, 5), 1, beam, 1, -6, 0.1), ...
                'air_cavity_metrics:size');
        end

        %% ---- Config validation and cache signature ------------------------

        function invalid_settings_error_at_load(tc)
            cases = {
                {'headmodel.air_cavities.source', 'sinus'},                                       'prestus:airCavitiesSource'
                {'headmodel.air_cavities.source', 'segmentation'},                                'prestus:airCavitiesLabels'
                {'headmodel.air_cavities.source', 'segmentation'; 'headmodel.air_cavities.labels', 0}, 'prestus:airCavitiesLabels'
                {'headmodel.air_cavities.source', 'segmentation'; 'headmodel.air_cavities.labels', 1.5}, 'prestus:airCavitiesLabels'
                {'headmodel.air_cavities.source', 'simnibs'; 'headmodel.air_cavities.simnibs_labels', []}, 'prestus:airCavitiesSimnibsLabels'
                {'headmodel.air_cavities.source', 'simnibs'; 'headmodel.air_cavities.simnibs_labels', -1}, 'prestus:airCavitiesSimnibsLabels'
                {'headmodel.air_cavities.source', 'hu'; 'pct.enabled', 0},                        'prestus:airCavitiesNeedsCT'
                {'headmodel.air_cavities.warn_beam_fraction', 2},                                 'prestus:airCavitiesFraction'
                };
            for i = 1:size(cases, 1)
                tc.verifyError(@() load_parameters(cases{i, 1}), cases{i, 2}, cases{i, 2});
            end
            p = load_parameters({'headmodel.air_cavities.source', 'segmentation'; ...
                                 'headmodel.air_cavities.labels', [11 12]});
            tc.verifyEqual(char(p.headmodel.air_cavities.source), 'segmentation');
            p = load_parameters();                        % the default
            tc.verifyEqual(char(p.headmodel.air_cavities.source), 'simnibs');
            tc.verifyEqual(double(p.headmodel.air_cavities.simnibs_labels), 501);
        end

        function settings_leave_cache_signatures_unchanged(tc)
            p0 = load_parameters();
            p1 = p0;
            p1.headmodel.air_cavities.source = 'segmentation';
            p1.headmodel.air_cavities.labels = 11;
            p1.headmodel.air_cavities.warn_db = -12;
            for kind = {'grid', 'medium', 'acoustic', 'headmodel'}
                tc.verifyEqual(cache_signature(p1, kind{1}), cache_signature(p0, kind{1}), kind{1});
            end
            p2 = p0;
            p2.headmodel.smooth_fwhm_mm = 2;            % negative control: a real setting still counts
            tc.verifyNotEqual(cache_signature(p2, 'grid'), cache_signature(p0, 'grid'));
        end
    end
end

%% ---- Helpers -------------------------------------------------------------

function p = base_params(source)
    p = struct();
    p.simulation.medium = 'layered';
    p.grid.resolution_mm = 1;
    p.headmodel.air_cavities = struct('source', source, 'simnibs_labels', 501, 'labels', [], ...
        'hu_window', [-1200, -400], 'keep_boundary', false, 'warn_db', -6, 'warn_beam_fraction', 0.1);
end

function [ct, inner, edge_pocket, c_inner, c_edge] = synthetic_ct()
% Head (40 HU) in air (-1000 HU), with an internal air pocket and a pocket
% cut off from the exterior air that touches the grid edge
    n = 41; c = (n + 1) / 2;
    [x, y, z] = ndgrid(1:n, 1:n, 1:n);
    r = sqrt((x - c).^2 + (y - c).^2 + (z - c).^2);
    ct = -1000 * ones(n, n, n);
    ct(r < 16) = 40;
    inner = sqrt((x - c).^2 + (y - c).^2 + (z - 14).^2) < 3.5;
    c_inner = [c c 14];
    ct(inner) = -1000;
    % Tissue slab at the x = n face separating a pocket from the exterior
    ct(n-4:n, 15:27, 15:27) = 40;
    edge_pocket = false(n, n, n);
    edge_pocket(n-1:n, 19:23, 19:23) = true;
    ct(edge_pocket) = -1000;
    c_edge = [n 21 21];
end

function [pr, beam, cav_on, cav_off] = beam_scene(sigma)
% Gaussian beam along +z from (16,16,5) to (16,16,35); cavity slab on the
% axis (z = 20:22) or beside the beam (12 voxels off axis, same depth)
    n = [31 31 41];
    [x, y, ~] = ndgrid(1:n(1), 1:n(2), 1:n(3));
    r2 = (x - 16).^2 + (y - 16).^2;
    pr = exp(-r2 / (2 * sigma^2));
    beam = struct('trans_pos', [16 16 5], 'focus_pos', [16 16 35], 'exit_offset_vox', 3, 'radius_vox', 14);
    cav_on = false(n);  cav_on(13:19, 13:19, 20:22) = true;
    cav_off = false(n); cav_off(27:30, 14:18, 20:22) = true;
end

function a = flip_dim(a, d)
    a = flip(a, d);
end
