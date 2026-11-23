classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_skull_from_hu < matlab.unittest.TestCase
% TEST_SKULL_FROM_HU  HU-defined skull (pct.bone_source 'hu', SKULL_FROM_HU) on a synthetic head.
%
% A spherical shell of bone (HU 1000) with a small hole and a thin low-HU
% rim stands in for the skull; the segmentation's skull is a thicker shell.
% Checks the threshold, the hole closing, the largest-component rule, the HU
% floor and cap, the relabelling of segmentation skull (nearest medium or a
% fixed layer), that every bone source (HU threshold, mask) is reconciled
% alike, the bone mask file, the cache signature rules, the HU reference in
% the max-normalised mappings and on a resampled grid, the common limits
% scheme, and the CT interpolation.
%
% Run with:   results = runtests('tests/test_skull_from_hu.m');

    properties
        p
        masks
        hu
        idx
    end

    methods (TestMethodSetup)
        function make_head(tc)
            p = load_parameters();
            p.grid.resolution_mm = 0.5;
            p.pct.enabled = 1;
            p.pct.bone_source = 'hu';
            p.pct.bone_hu_threshold = 300;
            p.pct.bone_median_mm = 0;
            p.pct.bone_close_mm = 0;
            p.pct.bone_soft_layer = 'brain';
            labels = fieldnames(p.medium_properties);
            tc.idx = struct('water', find(strcmp(labels, 'water')), 'brain', find(strcmp(labels, 'brain')), ...
                'skin', find(strcmp(labels, 'skin')), 'skull', find(strcmp(labels, 'skull')));
            n = 61; c = (n + 1) / 2;
            [x, y, z] = ndgrid(1:n, 1:n, 1:n);
            r = sqrt((x - c).^2 + (y - c).^2 + (z - c).^2);
            m = repmat(tc.idx.water, n, n, n);
            m(r < 26) = tc.idx.skin;
            m(r < 24) = tc.idx.skull;                    % segmentation skull: r 18-24
            m(r < 18) = tc.idx.brain;
            h = zeros(n, n, n, 'single');
            h(r < 26) = 40;
            h(r >= 19 & r < 23) = 1000;                  % CT bone: r 19-23
            h(r >= 18 & r < 19) = 200;                   % low-HU rim inside the segmentation skull
            h(abs(x - c) <= 1 & abs(y - c) <= 1 & z > c) = 40;   % 3-voxel hole through the top of the shell
            h(5, 5, 5) = 1000;                           % isolated bright voxel in water
            tc.p = p; tc.masks = m; tc.hu = h;
        end
    end

    methods (Test)
        function test_threshold_and_relabelling(tc)
            [m, bone, h, info] = skull_from_hu(tc.p, tc.masks, tc.hu);
            seg_skull = tc.masks == tc.idx.skull;
            tc.verifyTrue(all(tc.hu(bone) > 300), 'without closing, bone is a subset of HU > threshold');
            tc.verifyTrue(all(m(bone) == tc.idx.skull), 'HU bone is skull');
            tc.verifyTrue(all(m(seg_skull & ~bone) == tc.idx.brain), 'segmentation skull outside the bone joins bone_soft_layer');
            tc.verifyFalse(bone(5, 5, 5), 'only the largest component is kept');
            tc.verifyEqual(info.n_bone, nnz(bone));
            tc.verifyGreaterThanOrEqual(min(h(bone)), single(300));
        end

        function test_closing_fills_hole_and_floors_hu(tc)
            p = tc.p; p.pct.bone_close_mm = 3;           % 6-voxel cube closes the 3-voxel hole
            [~, bone, h] = skull_from_hu(p, tc.masks, tc.hu);
            c = 31;
            tc.verifyTrue(bone(c, c, c + 21), 'the hole through the shell is closed');
            tc.verifyEqual(h(c, c, c + 21), single(300), 'closed voxels get the threshold HU');
            [~, bone0] = skull_from_hu(tc.p, tc.masks, tc.hu);
            tc.verifyFalse(bone0(c, c, c + 21), 'without closing the hole stays open');
        end

        function test_median_filter(tc)
            p = tc.p; p.pct.bone_median_mm = 1.5;
            [~, bone] = skull_from_hu(p, tc.masks, tc.hu);
            tc.verifyGreaterThan(nnz(bone), 0);
        end

        function test_nearest_soft_layer(tc)
            % low-HU rim on the inner side joins the brain; an outer rim joins the skin;
            % no brain voxel is left between skin and bone
            p = tc.p; p.pct.bone_soft_layer = 'nearest';
            [m, bone, ~, info] = skull_from_hu(p, tc.masks, tc.hu);
            n = 61; c = (n + 1) / 2;
            [x, y, z] = ndgrid(1:n, 1:n, 1:n);
            r = sqrt((x - c).^2 + (y - c).^2 + (z - c).^2);
            seg_skull = tc.masks == tc.idx.skull;
            outer = seg_skull & ~bone & r >= 23;         % segmentation skull r 23-24, outside the CT bone
            inner = seg_skull & ~bone & r < 19;          % low-HU rim r 18-19
            tc.verifyGreaterThan(nnz(outer), 0); tc.verifyGreaterThan(nnz(inner), 0);
            tc.verifyTrue(all(m(outer) == tc.idx.skin), 'outer rim joins the skin');
            tc.verifyTrue(all(m(inner) == tc.idx.brain), 'inner rim joins the brain');
            tc.verifyFalse(any(m(r >= 23 & r < 26) == tc.idx.brain), 'no brain between bone and skin surface');
            tc.verifyFalse(any(m(seg_skull & ~bone) == tc.idx.skull), 'no segmentation skull left outside the bone');
            tc.verifyEqual(info.to_soft.skin + info.to_soft.brain, info.n_to_soft);
            % negative control: a fixed layer puts brain outside the bone
            p.pct.bone_soft_layer = 'brain';
            m2 = skull_from_hu(p, tc.masks, tc.hu);
            tc.verifyTrue(all(m2(outer) == tc.idx.brain));
        end

        function test_bad_soft_layer_errors(tc)
            p = tc.p; p.pct.bone_soft_layer = 'nope';
            tc.verifyError(@() skull_from_hu(p, tc.masks, tc.hu), 'skull_reconcile:badSoftLayer');
        end

        function test_signature_ignores_keys_for_segmentation_skull(tc)
            % the new keys count only with bone_source 'hu', so existing caches keep their signature
            p = load_parameters(); p.pct.enabled = 1;
            p.simulation.medium = 'water';                % no segmentation fingerprint needed
            p.pct.bone_source = 'segmentation';
            s1 = cache_signature(p, 'grid');
            p.pct.bone_close_mm = 7;
            s2 = cache_signature(p, 'grid');
            tc.verifyEqual(s1, s2);
            tc.verifyFalse(isfield(s1.pct, 'bone_source'));
            p.pct.bone_source = 'hu';
            s3 = cache_signature(p, 'grid');
            tc.verifyEqual(s3.pct.bone_close_mm, 7);
        end

        function test_opening_removes_spur(tc)
            % a 1-voxel spur on the shell survives without opening (it joins
            % the largest component) and goes with a 3-voxel opening; the
            % 4-voxel shell stays
            c = 31; h = tc.hu;
            h(c, c, (c - 30):(c - 23)) = 1000;
            [~, bone0] = skull_from_hu(tc.p, tc.masks, h);
            tc.verifyTrue(bone0(c, c, c - 27), 'negative control: without opening the spur is bone');
            p = tc.p; p.pct.bone_open_mm = 1.5;         % 3-voxel cube at 0.5 mm
            [~, bone] = skull_from_hu(p, tc.masks, h);
            tc.verifyFalse(bone(c, c, c - 27), 'the opening removes the spur');
            tc.verifyTrue(bone(c, c + 21, c), 'the shell stays bone');
        end

        function test_soft_tissue_from_hu(tc)
            % pct.soft_hu_threshold: tissue below it becomes water, water at or
            % above it takes the nearest soft tissue, skull is untouched
            c = 31; n = size(tc.masks, 1);
            [x, y, z] = ndgrid(1:n, 1:n, 1:n);
            r = sqrt((x - c).^2 + (y - c).^2 + (z - c).^2);
            h = tc.hu;
            h(r < 3) = -1000;                            % air pocket inside the brain
            h(r >= 26 & r < 27) = 40;                    % tissue just beyond the segmentation's head
            [m0, ~, h0] = skull_from_hu(tc.p, tc.masks, h);
            p = tc.p; p.pct.soft_hu_threshold = -100;
            [m, info] = tissue_from_hu(p, m0, h0);
            tc.verifyTrue(all(m(r < 3) == tc.idx.water), 'the air pocket becomes water');
            tc.verifyTrue(all(m(r >= 26 & r < 27) == tc.idx.skin), 'tissue beyond the head takes the nearest soft tissue (skin)');
            tc.verifyEqual(m(m0 == tc.idx.skull), m0(m0 == tc.idx.skull), 'skull is untouched');
            tc.verifyEqual(info.n_to_water, nnz(r < 3));
            q = tc.p; q.pct.soft_hu_threshold = [];
            tc.verifyEqual(tissue_from_hu(q, m0, h0), m0, 'negative control: off by default');
            tc.verifyError(@() tissue_from_hu(p, m0, []), 'tissue_from_hu:noCT');
        end

        function test_signature_ignores_new_keys_at_defaults(tc)
            % bone_open_mm 0 and soft_hu_threshold [] keep existing signatures;
            % other values count
            p = load_parameters(); p.pct.enabled = 1; p.simulation.medium = 'water';
            p.pct.bone_source = 'hu';
            p.pct = rmfield(p.pct, intersect(fieldnames(p.pct), {'bone_open_mm', 'soft_hu_threshold'}));
            s0 = cache_signature(p, 'grid');
            p.pct.bone_open_mm = 0; p.pct.soft_hu_threshold = [];
            tc.verifyEqual(cache_signature(p, 'grid'), s0);
            p.pct.bone_open_mm = 1; p.pct.soft_hu_threshold = -100;
            s2 = cache_signature(p, 'grid');
            tc.verifyEqual(s2.pct.bone_open_mm, 1);
            tc.verifyEqual(s2.pct.soft_hu_threshold, -100);
        end

        function test_every_bone_source_reconciles_alike(tc)
            % a mask equal to the HU bone gives the HU path's result, for every
            % soft-layer rule and with or without closing
            for soft = {'nearest', 'brain'}
                for clo = [0 3]
                    p = tc.p; p.pct.bone_soft_layer = soft{1}; p.pct.bone_close_mm = clo;
                    [m1, b1, h1, i1] = skull_from_hu(p, tc.masks, tc.hu);
                    p.pct.bone_source = 'mask';
                    [m2, b2, h2, i2] = skull_reconcile(p, tc.masks, tc.hu > 300, tc.hu, 'mask');
                    msg = sprintf('%s, closing %g mm', soft{1}, clo);
                    tc.verifyEqual(m2, m1, msg); tc.verifyEqual(b2, b1, msg); tc.verifyEqual(h2, h1, msg);
                    tc.verifyEqual(i2.n_to_soft, i1.n_to_soft, msg);
                end
            end
        end

        function test_mask_bone_gets_floor_and_relabelling(tc)
            % a mask that includes the low-HU rim (drawn on another grid): the rim
            % becomes skull at the floor HU; outer segmentation skull joins the skin
            n = 61; c = (n + 1) / 2;
            [x, y, z] = ndgrid(1:n, 1:n, 1:n);
            r = sqrt((x - c).^2 + (y - c).^2 + (z - c).^2);
            mask = r >= 18 & r < 23;
            rim = r >= 18 & r < 19;
            p = tc.p; p.pct.bone_source = 'mask'; p.pct.bone_soft_layer = 'nearest';
            [m, bone, h, info] = skull_reconcile(p, tc.masks, mask, tc.hu, 'mask');
            tc.verifyTrue(all(m(rim) == tc.idx.skull));
            tc.verifyEqual(h(rim), repmat(single(300), nnz(rim), 1), 'rim HU raised to the threshold');
            tc.verifyEqual(info.n_floor, nnz(bone & tc.hu < 300));
            outer = tc.masks == tc.idx.skull & r >= 23;
            tc.verifyTrue(all(m(outer) == tc.idx.skin), 'outer segmentation skull joins the skin');
            % negative control: without a CT there is no floor, the labels are the same
            [m0, ~, h0, info0] = skull_reconcile(p, tc.masks, mask, [], 'mask');
            tc.verifyEmpty(h0); tc.verifyEqual(info0.n_floor, 0); tc.verifyEqual(m0, m);
        end

        function test_bone_mask_file_on_other_grid(tc)
            % a 0.5 mm mask with its own origin lands on the 1 mm segmentation grid
            % through the world affines; 'segmentation' ignores the file
            tmp = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            seg = zeros(20, 20, 20, 'uint8'); seg(5:12, 5:12, 5:12) = 4;
            seg_file = write_nii(fullfile(tmp, 'seg'), seg, eye(4));
            bm = zeros(30, 30, 30, 'uint8'); bm(9:16, 9:16, 9:16) = 1;   % 0-based 8-15: 6 to 9.5 mm
            T = diag([0.5 0.5 0.5 1]); T(4, 1:3) = 2;
            write_nii(fullfile(tmp, 'bone_sub-001'), bm, T);
            p = tc.p; p.subject_id = 1; p.pct.bone_source = 'mask';
            p.path.bone_mask_pattern = fullfile(tmp, 'bone_sub-%03d.nii.gz');
            hdr = niftiinfo(seg_file);
            [bone, f] = pct_load_bone_mask(p, seg, hdr);
            tc.verifyEqual(f, fullfile(tmp, 'bone_sub-001.nii.gz'));
            tc.verifyTrue(all(bone(8:9, 8:9, 8:9), 'all'), 'world 7-8 mm is bone');
            tc.verifyFalse(any(bone([6 13], 8, 8)), 'world 5 and 12 mm are not');
            p.pct.bone_source = 'segmentation';
            tc.verifyEqual(pct_load_bone_mask(p, seg, hdr), ismember(seg, charm_seg_labels().bonemask));
            p.pct.bone_source = 'mask'; p.path.bone_mask_pattern = fullfile(tmp, 'none_%03d.nii.gz');
            tc.verifyError(@() pct_load_bone_mask(p, seg, hdr), 'PRESTUS:pct:boneMaskNotFound');
        end

        function test_bone_source_rules(tc)
            p = tc.p; p.pct.enabled = 0;
            tc.verifyEqual(skull_bone_source(p), 'segmentation', '''hu'' needs the CT');
            p.pct.bone_source = 'mask';
            tc.verifyEqual(skull_bone_source(p), 'mask', 'a mask works without a CT');
            p.pct.bone_source = 'bogus';
            tc.verifyError(@() skull_bone_source(p), 'skull_bone_source:badSource');
            % cache signature: bone keys and the file count for 'mask' only
            tmp = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            write_nii(fullfile(tmp, 'bone_sub-001'), uint8(ones(4, 4, 4)), eye(4));
            q = load_parameters(); q.subject_id = 1; q.simulation.medium = 'water';
            q.path.bone_mask_pattern = fullfile(tmp, 'bone_sub-%03d.nii.gz');
            s0 = cache_signature(q, 'grid');
            tc.verifyFalse(isfield(s0, 'bone_mask')); tc.verifyFalse(isfield(s0.pct, 'bone_close_mm'));
            q.pct.bone_source = 'mask';
            s1 = cache_signature(q, 'grid');
            tc.verifyEqual(s1.bone_mask.name, 'bone_sub-001.nii.gz');
            tc.verifyTrue(isfield(s1.pct, 'bone_close_mm'));
        end

        function test_reference_resolution(tc)
            % the reference on a resampled grid follows the interpolation: spline
            % overshoots the sharp bone edge, linear does not; the segmentation grid
            % gives the data maximum; the HU limits cap the overshoot
            tmp = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            d = fullfile(tmp, 'sub-001'); mkdir(d);
            ct = zeros(24, 24, 24, 'single'); ct(:, :, 11:14) = 1000;
            seg = zeros(24, 24, 24, 'uint8'); seg(:, :, 10:15) = 4;
            write_nii(fullfile(d, 'pseudoCT'), ct, eye(4));
            write_nii(fullfile(d, 'tissues_mask'), seg, eye(4));
            p = load_parameters(); p.subject_id = 1; p.pct.enabled = 1;
            p.path.pct = tmp; p.path.ct_pattern = '';
            p.pct.hu_reference_layers = {'skull'}; p.pct.hu_reference_percentile = 100;
            tc.verifyEqual(pct_hu_reference(p), 1000, 'AbsTol', 1e-3);
            p.pct.hu_reference_resolution_mm = 0.4;
            p.pct.interpolation = 'linear';
            tc.verifyLessThanOrEqual(pct_hu_reference(p), 1000 + 1e-3);
            p.pct.interpolation = 'spline';
            [r, info] = pct_hu_reference(p);
            tc.verifyGreaterThan(r, 1000, 'spline overshoot raises the reference');
            tc.verifyEqual(info.resolution_mm, 0.4);
            p.pct.limits.hu_max = 1000;
            tc.verifyEqual(pct_hu_reference(p), 1000, 'AbsTol', 1e-3);
            p.pct.hu_reference_resolution_mm = -1;
            tc.verifyError(@() pct_hu_reference(p), 'pct_hu_reference:badResolution');
        end

        function test_reference_scales_every_normalised_mapping(tc)
            % max-normalised mappings take the given reference, not the data maximum,
            % so a 'cropped' set of voxels maps the same way
            p = load_parameters(); p.pct.enabled = 1;
            hu = single([400; 800; 1000]);
            for m = {{'mapping_density', 'marsac-babelbrain'}, {'mapping_density', 'aubry'}, ...
                     {'mapping_soundspeed', 'aubry'}, {'mapping_attenuation', 'mueller'}, {'mapping_attenuation', 'aubry'}}
                q = p; q.pct.(m{1}{1}) = m{1}{2};
                full = pct_map_properties(q, hu, 1:3, 2000);
                crop = pct_map_properties(q, hu(1:2), 1:2, 2000);
                f = {'density', 'sound_speed', 'alpha_coeff'};
                for k = 1:3
                    tc.verifyEqual(full.(f{k})(1:2), crop.(f{k}), 'AbsTol', 1e-6, ...
                        sprintf('%s %s: independent of the voxels present', m{1}{2}, f{k}));
                end
                % negative control: without a reference, the data maximum sets the scale
                if strcmp(m{1}{1}, 'mapping_density')
                    a = pct_map_properties(q, hu, 1:3); c = pct_map_properties(q, hu(1:2), 1:2);
                    tc.verifyNotEqual(a.density(1), c.density(1));
                end
            end
            q = p; q.pct.mapping_density = 'marsac-babelbrain';
            d = pct_map_properties(q, hu, 1:3, 2000);
            tc.verifyEqual(double(d.density(3)), 1000 + 1700 * 1000 / 2000, 'AbsTol', 1e-3);
            % HU above the reference map as the densest bone (MEDIUM_PCT_POROSITY
            % limits porosity to [0, 1]): 1000 HU maps like the 900 HU reference;
            % attenuation stays real
            q.pct.mapping_attenuation = 'aubry';
            a = pct_map_properties(q, hu, 1:3, 900);
            r = pct_map_properties(q, single([400; 800; 900]), 1:3, 900);
            tc.verifyTrue(isreal(a.alpha_coeff) && all(isfinite(a.alpha_coeff)));
            for f = {'density', 'alpha_coeff'}
                tc.verifyEqual(double(a.(f{1})(3)), double(r.(f{1})(3)), 'AbsTol', 1e-6, ...
                    sprintf('%s above the reference = densest bone', f{1}));
            end
            tc.verifyNotEqual(double(a.density(2)), double(a.density(3)), 'negative control: below the reference HU still matter');
        end

        function test_limits_common_scheme(tc)
            p = load_parameters(); p.pct.enabled = 1;
            p.pct.mapping_density = 'k-plan'; p.pct.mapping_soundspeed = 'k-plan'; p.pct.mapping_attenuation = 'k-plan';
            hu = single([-500; 300; 1500; 3000]);
            % defaults: water floors on density and sound speed, nothing else
            [r, ~, info] = pct_map_properties(p, hu, 1:4);
            tc.verifyGreaterThanOrEqual(min(r.density), p.medium_properties.water.density - 1e-6);
            tc.verifyGreaterThanOrEqual(min(r.sound_speed), p.medium_properties.water.sound_speed - 1e-6);
            tc.verifyEqual(info.hu.n_high, 0);
            tc.verifyGreaterThan(info.density.n_low, 0, '-500 HU lies below water density');
            % HU cap applies before mapping; numbers and 'water' work for every quantity
            p.pct.limits.hu_max = 2100; p.pct.limits.density_max = 1900;
            p.pct.limits.sound_speed_max = 2800; p.pct.limits.alpha_coeff_min = 'water';
            [r2, ~, info2] = pct_map_properties(p, hu, 1:4);
            tc.verifyEqual(max(r2.hu), single(2100));
            tc.verifyEqual(info2.hu.n_high, 1);
            tc.verifyLessThanOrEqual(max(r2.density), 1900 + 1e-6);
            tc.verifyLessThanOrEqual(max(r2.sound_speed), 2800 + 1e-6);
            tc.verifyGreaterThanOrEqual(min(r2.alpha_coeff), p.medium_properties.water.alpha_coeff);
            % no floor when switched off: the mapping's own value comes through
            p.pct.limits.density_min = [];
            r3 = pct_map_properties(p, hu, 1:4);
            tc.verifyLessThan(r3.density(1), p.medium_properties.water.density);
            % invalid limits
            p.pct.limits.hu_min = 'water';
            tc.verifyError(@() pct_map_properties(p, hu, 1:4), 'pct_apply_limits:badLimit');
            p.pct.limits.hu_min = 3000;
            tc.verifyError(@() pct_map_properties(p, hu, 1:4), 'pct_apply_limits:badLimits');
        end

        function test_signature_keeps_default_caches(tc)
            % new keys at their defaults leave the signature unchanged; changes count
            p = load_parameters(); p.pct.enabled = 1; p.simulation.medium = 'water';
            p.pct.mapping_density = 'k-plan'; p.pct.mapping_soundspeed = 'k-plan'; p.pct.mapping_attenuation = 'k-plan';
            s = cache_signature(p, 'grid');
            tc.verifyFalse(any(isfield(s.pct, {'interpolation', 'limits', 'hu_reference_layers', 'hu_reference_percentile'})));
            q = p; q.pct.interpolation = 'spline';
            tc.verifyEqual(cache_signature(q, 'grid').pct.interpolation, 'spline');
            q = p; q.pct.limits.hu_max = 2100;
            tc.verifyTrue(isfield(cache_signature(q, 'grid').pct, 'limits'));
            q = p; q.pct.mapping_density = 'aubry';
            tc.verifyTrue(isfield(cache_signature(q, 'grid').pct, 'hu_reference_percentile'), ...
                'the reference counts for a normalised mapping');
        end

        function test_interpolation_override(tc)
            % the CT honours pct.interpolation; spline overshoots at a sharp edge, linear does not
            v = zeros(12, 12, 12, 'single'); v(:, :, 7:end) = 1000;
            M = diag([2 2 2 1]); M(4, 1:3) = -1;                       % 2x upsampling: [out 1] = [in 1] * M
            lin = affine_resample_3d(v, M, [22 22 22], 'linear', 0);
            spl = affine_resample_3d(v, M, [22 22 22], 'spline', 0);
            tc.verifyLessThanOrEqual(max(lin(:)), single(1000) + 1e-3);
            tc.verifyGreaterThan(max(spl(:)), single(1000));
            p = load_parameters(); p.pct.limits.hu_max = 1000;
            tc.verifyLessThanOrEqual(max(pct_apply_limits(p, 'hu', spl), [], 'all'), single(1000));
        end
    end
end

function file = write_nii(base, img, T)
% NIfTI with affine T (row-vector convention, 0-based voxels -> mm); returns the .nii.gz path
    niftiwrite(img, base, 'Compressed', true);
    file = [base '.nii.gz'];
    info = niftiinfo(file);
    info.Transform = affine3d(T);
    info.TransformName = 'Sform';          % without it niftiwrite keeps the identity affine
    info.PixelDimensions = abs(diag(T(1:3, 1:3)))';
    niftiwrite(img, base, info, 'Compressed', true);
    assert(max(abs(niftiinfo(file).Transform.T - T), [], 'all') < 1e-6, 'write_nii: affine not stored');
end
