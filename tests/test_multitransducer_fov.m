classdef test_multitransducer_fov < matlab.unittest.TestCase
% TEST_MULTITRANSDUCER_FOV  Partial acoustic FOV and field combination with several transducers
%
%   acoustic_grid_fov       - coherent coupling: the FOV box holds every beam and
%                             every transducer moves into FOV voxels
%   fov_restore_transducers - all transducers return to full-grid voxels
%   combine_async_intensity - FOV caches on a shared grid, mixed FOV/full-grid
%                             caches, and FOV caches on different grids
%   superpose_complex_fields - coherent / incoherent / max combination of
%                             complex fields, FOV back-projection, checks
%   p_complex phase reference - the t0 correction used by acoustic_wrapper
%
%   Run with:   results = runtests('tests/test_multitransducer_fov.m');

    methods (Test, TestTags = {'fov'})

        function test_single_transducer_box_unchanged(tc)
            % One transducer: lateral box centred on the focus, axial from the transducer
            p = fov_params(1);
            [pf, ~, ~, t_fov, f_fov, off] = run_fov(p);
            tc.verifyEqual(off, [20 20 5]);              % focus 29 +- 10 mm, z from the transducer
            tc.verifyEqual(pf.grid.dims, [21 21 41]);
            tc.verifyEqual(t_fov, [11 11 1]);
            tc.verifyEqual(f_fov, [11 11 36]);
        end

        function test_oblique_single_beam_keeps_transducer(tc)
            % One transducer whose beam is far off the grid z axis (an async sub-run
            % on a grid aligned to another transducer): the box must follow its beam
            p = fov_params(1);
            p.transducer(1).trans_pos = [3 30 20];
            p.transducer(1).focus_pos = [30 30 40];          % ~53 deg off z
            [pf, ~, ~, t_fov, f_fov, off] = run_fov(p);
            tc.verifyLessThanOrEqual(off(1), 3);
            tc.verifyEqual(t_fov, [3 30 20] - off + 1);
            tc.verifyEqual(f_fov, [30 30 40] - off + 1);
            tc.verifyTrue(all(t_fov >= 1 & t_fov <= pf.grid.dims));
        end

        function test_opposed_transducers_share_box(tc)
            % Two opposed transducers: the box spans both, both move into it
            p = fov_params(2);
            p.transducer(2).trans_pos = [30 30 75];
            p.transducer(2).focus_pos = [30 30 42];
            [pf, ~, ~, ~, ~, off] = run_fov(p);
            tc.verifyLessThanOrEqual(off(3), 5);
            tc.verifyGreaterThanOrEqual(off(3) + pf.grid.dims(3) - 1, 75);
            for ti = 1:2
                for key = {'trans_pos', 'focus_pos'}
                    tc.verifyEqual(pf.transducer(ti).(key{1}), ...
                        p.transducer(ti).(key{1}) - off + 1, ...
                        sprintf('transducer %d %s', ti, key{1}));
                end
            end
        end

        function test_angled_transducer_widens_box(tc)
            % Second beam along x: the box must extend in x to its transducer
            p = fov_params(2);
            p.transducer(2).trans_pos = [3 30 40];
            p.transducer(2).focus_pos = [30 30 40];
            [pf, ~, ~, ~, ~, off] = run_fov(p);
            tc.verifyLessThanOrEqual(off(1), 3);
            t2 = pf.transducer(2).trans_pos;
            tc.verifyTrue(all(t2 >= 1 & t2 <= pf.grid.dims));
            tc.verifyEqual(t2, [3 30 40] - off + 1);
        end

        function test_per_transducer_auto_length(tc)
            % 'auto' length resolved per transducer: the deeper focus of T2 is covered
            p = fov_params(2);
            p.grid.acoustic_fov_length_mm = 'auto';
            p.grid.z_beyond_focus_mm = 5;
            p.transducer(2).trans_pos = [30 30 78];
            p.transducer(2).focus_pos = [30 30 20];   % 58 mm deep vs 35 mm for T1
            [p, has_fov] = resolve_acoustic_fov_params(p, p.transducer(1).trans_pos, ...
                p.transducer(1).focus_pos);
            tc.verifyTrue(has_fov);
            tc.verifyEqual(p.state.acoustic_fov_per_tx.length_mm, [40 63], 'AbsTol', 1e-9);
            [pf, ~, ~, ~, ~, off] = run_fov(p);
            tc.verifyLessThanOrEqual(off(3), 78 - 63);
            tc.verifyEqual(pf.transducer(2).focus_pos, [30 30 20] - off + 1);
        end

        function test_restore_all_transducers(tc)
            p = fov_params(2);
            p.transducer(2).trans_pos = [30 30 75];
            p.transducer(2).focus_pos = [30 30 42];
            [pf, ~, ~, ~, ~, off] = run_fov(p);
            st.trans_pos_full = p.transducer(1).trans_pos;
            st.focus_pos_full = p.transducer(1).focus_pos;
            st.transducer_pos_full = struct('trans_pos', {p.transducer.trans_pos}, ...
                                            'focus_pos', {p.transducer.focus_pos});
            tc.verifyNotEqual(off, [1 1 1]);
            pr = fov_restore_transducers(pf, st);
            tc.verifyEqual(pr.transducer(2).trans_pos, [30 30 75]);
            tc.verifyEqual(pr.transducer(2).focus_pos, [30 30 42]);
            tc.verifyEqual(pr.transducer(1).focus_pos, p.transducer(1).focus_pos);
        end
    end

    methods (Test, TestTags = {'async_combine'})

        function test_combine_fov_shared_grid(tc)
            % Two FOV caches on one full grid: max/sum of back-projected fields
            d = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            full = [12 12 12];
            f1 = write_cache(fullfile(d, 'a.mat'), 2 * ones(4, 4, 4), [2 2 2], full);
            f2 = write_cache(fullfile(d, 'b.mat'), 3 * ones(4, 4, 4), [4 4 4], full);
            out = fullfile(d, 'c.mat');
            evalc('combine_async_intensity({f1, f2}, out);');
            s = load(out, 'sensor_data', 'acoustic_provenance', 'kgrid', 'acoustic_grid_dims');
            tc.verifySize(s.sensor_data.p_max_async, full);
            % grid metadata describe the full grid the combined field is on
            tc.verifyEqual(s.acoustic_grid_dims, full);
            tc.verifyEqual([s.kgrid.Nx, s.kgrid.Ny, s.kgrid.Nz], full);
            tc.verifyEqual([s.kgrid.Nt, s.kgrid.dt], [10, 1e-7]);
            tc.verifyEqual(s.sensor_data.p_max_async(5, 5, 5), 3, 'AbsTol', 1e-12);   % overlap
            tc.verifyEqual(s.sensor_data.p_max_async(2, 2, 2), 2, 'AbsTol', 1e-12);
            tc.verifyEqual(s.sensor_data.p_max_all(5, 5, 5), sqrt(13), 'AbsTol', 1e-12);
            tc.verifyEqual(s.sensor_data.p_max_all(10, 10, 10), 0);
            tc.verifyEmpty(s.acoustic_provenance.fov_offset_ac);
        end

        function test_combine_mixed_fov_and_full(tc)
            % An FOV cache and a full-grid cache on the same grid (used to error)
            d = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            full = [8 8 8];
            f1 = write_cache(fullfile(d, 'a.mat'), ones(full), [], []);
            f2 = write_cache(fullfile(d, 'b.mat'), 2 * ones(3, 3, 3), [3 3 3], full);
            out = fullfile(d, 'c.mat');
            evalc('combine_async_intensity({f1, f2}, out);');
            s = load(out, 'sensor_data');
            tc.verifyEqual(s.sensor_data.p_max_async(4, 4, 4), 2, 'AbsTol', 1e-12);
            tc.verifyEqual(s.sensor_data.p_max_async(1, 1, 1), 1, 'AbsTol', 1e-12);
        end

        function test_combine_fov_different_grids(tc)
            % FOV caches on different full grids (per_transducer mode) are remapped
            % onto the reference grid instead of erroring
            d = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            f1 = write_cache(fullfile(d, 'a.mat'), ones(4, 4, 4), [3 3 3], [10 10 10]);
            f2 = write_cache(fullfile(d, 'b.mat'), ones(4, 4, 4), [3 3 3], [12 12 12]);
            out = fullfile(d, 'c.mat');
            evalc('combine_async_intensity({f1, f2}, out);');
            s = load(out, 'sensor_data');
            tc.verifySize(s.sensor_data.p_max_async, [10 10 10]);
            tc.verifyGreaterThan(nnz(s.sensor_data.p_max_async), 0);
        end
    end

    methods (Test, TestTags = {'p_complex'})

        function test_superpose_coherent_incoherent_max(tc)
            d = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            full = [10 10 10];
            f1 = write_pc(fullfile(d, 'p1.mat'), ones(4, 4, 4), [2 2 2], full, 3e5);
            f2 = write_pc(fullfile(d, 'p2.mat'), -ones(4, 4, 4), [4 4 4], full, 3e5);
            out = superpose_complex_fields({f1, f2});
            tc.verifySize(out.p_coherent, full);
            tc.verifyEqual(out.p_coherent(5, 5, 5), 0, 'AbsTol', 1e-12);       % destructive
            tc.verifyEqual(out.p_incoherent(5, 5, 5), sqrt(2), 'AbsTol', 1e-12);
            tc.verifyEqual(out.p_max(5, 5, 5), 1, 'AbsTol', 1e-12);
            tc.verifyEqual(out.p_coherent(2, 2, 2), 1, 'AbsTol', 1e-12);       % T1 only
            tc.verifyEqual(out.coverage(5, 5, 5), uint8(2));
            tc.verifyEqual(out.coverage(1, 1, 1), uint8(0));
            tc.verifyEqual(out.freq_hz, 3e5);

            flipped = superpose_complex_fields({f1, f2}, Weights=[1, -1]);
            tc.verifyEqual(flipped.p_coherent(5, 5, 5), 2, 'AbsTol', 1e-12);   % constructive
            scaled = superpose_complex_fields({f1, f2}, Scale=[2, 1]);
            tc.verifyEqual(scaled.p_coherent(5, 5, 5), 1, 'AbsTol', 1e-12);
        end

        function test_superpose_rejects_bad_inputs(tc)
            d = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            a = write_pc(fullfile(d, 'a.mat'), ones(3, 3, 3), [], [], 3e5);
            b = write_pc(fullfile(d, 'b.mat'), ones(4, 4, 4), [], [], 3e5);
            c = write_pc(fullfile(d, 'c.mat'), ones(3, 3, 3), [], [], 5e5);
            p_complex = complex(ones(3, 3, 3)); %#ok<NASGU>
            legacy = fullfile(d, 'legacy.mat');
            save(legacy, 'p_complex');
            tc.verifyError(@() superpose_complex_fields({a, b}), 'superpose_complex_fields:gridMismatch');
            tc.verifyError(@() superpose_complex_fields({a, c}), 'superpose_complex_fields:frequencyMismatch');
            tc.verifyError(@() superpose_complex_fields({a, legacy}), ...
                'superpose_complex_fields:mixedPhaseReference');
            tc.verifyWarning(@() superpose_complex_fields({legacy, legacy}), ...
                'superpose_complex_fields:legacyPhaseReference');
        end

        function test_phase_reference_t0(tc)
            % The correction acoustic_wrapper applies: the phase extractAmpPhase
            % returns for a window starting at t0, minus 2*pi*f*t0, is the phase of
            % the signal at t = 0 — whatever t0 and dt are.
            f = 3e5; phi = 0.7;
            for dt = [1.9e-7, 2.3e-7]
                n_per = round(1 / (f * dt));
                dt = 1 / (f * n_per); %#ok<FXSET> integer samples per period
                t0 = 37 * n_per * dt + 5 * dt;
                t  = t0 + (0:3 * n_per - 1) * dt;
                x  = cos(2 * pi * f * t + phi);
                [~, ph] = extractAmpPhase(x, 1 / dt, f, 'Dim', 2);
                tc.verifyEqual(angle(exp(1i * (ph - 2 * pi * f * t0))), phi, 'AbsTol', 0.02);
            end
        end
    end
end

% -------------------------------------------------------------------------
function p = fov_params(n_tx)
p.grid.dims = [60 60 80];
p.grid.resolution_mm = 1;
p.grid.acoustic_fov_diameter_mm = 20;
p.grid.acoustic_fov_length_mm   = 40;
p.transducer = struct('trans_pos', [30 30 5], 'focus_pos', [30 30 40], 'freq_hz', 3e5);
for ti = 2:n_tx
    p.transducer(ti) = p.transducer(1);
end
end

function [pf, km, mm, t_fov, f_fov, off] = run_fov(p)
km.sound_speed = 1500 * ones(p.grid.dims);
masks = zeros(p.grid.dims, 'uint8');
[pf, km, mm, t_fov, f_fov, off] = acoustic_grid_fov(p, km, masks, ...
    p.transducer(1).trans_pos, p.transducer(1).focus_pos);
end

function f = write_cache(f, p, offset, full)
% Minimal acoustic cache as combine_async_intensity reads it (no planimg: no affine).
sensor_data.p_max_all = p;
sensor_data.p_min_all = -p;
kgrid = kWaveGrid(size(p, 1), 1e-3, size(p, 2), 1e-3, size(p, 3), 1e-3);   % the simulated (cropped) grid
kgrid.setTime(10, 1e-7);
kwave_medium = struct('sound_speed', 1500); %#ok<NASGU>
acoustic_provenance = struct();
if ~isempty(offset)
    acoustic_provenance.fov_offset_ac = offset;
    acoustic_provenance.full_ac_dims  = full;
end
acoustic_info.parameters.transducer = struct('target_isppa_wcm2', NaN); %#ok<STRNU>
save(f, 'sensor_data', 'kgrid', 'kwave_medium', 'acoustic_provenance', 'acoustic_info');
end

function f = write_pc(f, p, offset, full, freq)
p_complex = complex(p); %#ok<NASGU>
p_complex_info = struct('phase_reference', 'simulation_t0', 'freq_hz', freq, ...
    'record_t0_s', 0, 'dt_s', 1e-7, 'drive_level', 'simulated', 'grid_dims', size(p), ...
    'fov_offset_ac', offset, 'full_ac_dims', full); %#ok<NASGU>
save(f, 'p_complex', 'p_complex_info');
end
