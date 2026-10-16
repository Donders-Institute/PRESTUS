classdef test_grid_cache_transducers < matlab.unittest.TestCase
% TEST_GRID_CACHE_TRANSDUCERS  LOAD_GRID_CACHE_FILE keeps a run's own transducer
% settings when it reuses another run's grid: only the grid-space positions
% come from the cache.
%
%   Regression test: a reused cache with the same number of transducers
%   replaced the whole transducer array, so e.g. the ISPPA targets and drives
%   of a coherent run on an async run's grid were those of the async run.
%
%   Run with:   results = runtests('tests/test_grid_cache_transducers.m');

    methods (Test, TestTags = {'cache'})

        function test_caller_settings_cached_positions(tc)
            addpath(fullfile(fileparts(mfilename('fullpath')), 'fixtures'));
            tmp = tempname; mkdir(tmp);
            tc.addTeardown(@() rmdir(tmp, 's'));
            p = make_minimal_parameters(tmp);
            p.transducer(2) = p.transducer(1);

            % the cache: built by another run (no targets, other drive), grid-space positions
            c = p;
            for k = 1:2
                c.transducer(k).target_isppa_wcm2 = NaN;
                c.transducer(k).annular.elem_p_amp = 1000 * k;
                c.transducer(k).trans_pos = [10 20 30] + k;
                c.transducer(k).focus_pos = [40 50 60] + k;
            end
            grid_cache_info = struct('parameters', c); %#ok<NASGU>
            planimg = struct(); medium_masks = []; segmentation = []; bone_mask = []; pseudoCT = []; %#ok<NASGU>
            f = fullfile(tmp, 'grid_cache.mat');
            save(f, 'grid_cache_info', 'planimg', 'medium_masks', 'segmentation', 'bone_mask', 'pseudoCT');

            % this run: its own targets and drives
            p.transducer(1).target_isppa_wcm2 = 20;
            p.transducer(2).target_isppa_wcm2 = 40;
            p.transducer(2).annular.elem_p_amp = 55555;
            [q, ~, ~, ~, ~, ~, tp, fp] = load_grid_cache_file(f, p);

            tc.verifyEqual([q.transducer.target_isppa_wcm2], [20 40]);
            tc.verifyEqual(q.transducer(2).annular.elem_p_amp(1), 55555);
            for k = 1:2
                tc.verifyEqual(q.transducer(k).trans_pos, [10 20 30] + k);
                tc.verifyEqual(q.transducer(k).focus_pos, [40 50 60] + k);
            end
            tc.verifyEqual(tp, [11 21 31]); tc.verifyEqual(fp, [41 51 61]);
        end

        function test_cached_bilateral_positions(tc)
            % An axisymmetric grid cache also restores the bilateral snapshot
            % that the axisymmetric expansion and the thermal stage read.
            addpath(fullfile(fileparts(mfilename('fullpath')), 'fixtures'));
            tmp = tempname; mkdir(tmp);
            tc.addTeardown(@() rmdir(tmp, 's'));
            p = make_minimal_parameters(tmp);
            p.transducer = p.transducer(1);
            c = p;
            c.grid.dims = [240 71];
            c.grid.axisym_bilateral_dims = [240 140];
            c.grid.axisym_transposed = true;
            c.transducer.trans_pos = [1 1];
            c.transducer.focus_pos = [50 1];
            c.transducer.trans_pos_bilateral = [1 70];
            c.transducer.focus_pos_bilateral = [50 70];
            grid_cache_info = struct('parameters', c); %#ok<NASGU>
            planimg = struct(); medium_masks = []; segmentation = []; bone_mask = []; pseudoCT = []; %#ok<NASGU>
            f = fullfile(tmp, 'grid_cache.mat');
            save(f, 'grid_cache_info', 'planimg', 'medium_masks', 'segmentation', 'bone_mask', 'pseudoCT');

            q = load_grid_cache_file(f, p);
            tc.verifyEqual(q.transducer.trans_pos_bilateral, [1 70]);
            tc.verifyEqual(q.transducer.focus_pos_bilateral, [50 70]);
            tc.verifyEqual(q.grid.axisym_bilateral_dims, [240 140]);
            tc.verifyTrue(q.grid.axisym_transposed);
        end

        function test_caller_grid_settings_cached_dims(tc)
            % The grid stage's outputs come from the cache; the run's other
            % grid settings stay (a partial FOV on a full-grid run's cache, or
            % none for a job reading a back-projected field).
            addpath(fullfile(fileparts(mfilename('fullpath')), 'fixtures'));
            tmp = tempname; mkdir(tmp);
            tc.addTeardown(@() rmdir(tmp, 's'));
            p = make_minimal_parameters(tmp);

            c = p;
            c.grid.dims = [64 72 80];
            c.grid.acoustic_fov_diameter_mm = 'auto';
            c.grid.acoustic_fov_length_mm   = 'auto';
            grid_cache_info = struct('parameters', c); %#ok<NASGU>
            planimg = struct(); medium_masks = []; segmentation = []; bone_mask = []; pseudoCT = []; %#ok<NASGU>
            f = fullfile(tmp, 'grid_cache.mat');
            save(f, 'grid_cache_info', 'planimg', 'medium_masks', 'segmentation', 'bone_mask', 'pseudoCT');

            p.grid.dims = [1 1 1];
            p.grid.acoustic_fov_diameter_mm = [];
            p.grid.acoustic_fov_length_mm   = [];
            q = load_grid_cache_file(f, p);
            tc.verifyEqual(q.grid.dims, [64 72 80]);
            tc.verifyEmpty(q.grid.acoustic_fov_diameter_mm);
            tc.verifyEmpty(q.grid.acoustic_fov_length_mm);
        end

    end
end
