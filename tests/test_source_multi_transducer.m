classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'water'}) test_source_multi_transducer < matlab.unittest.TestCase
% TEST_SOURCE_MULTI_TRANSDUCER  SOURCE_CREATE keeps each transducer's own drive.
%
%   Regression test: transducers 2..N were built from transducer 1's settings
%   (drive, phases, geometry), with only their positions their own. Builds a
%   two-transducer kWaveArray source on the small water test grid (~1 min).
%
%   Run with:   results = runtests('tests/test_source_multi_transducer.m');

    methods (Test, TestTags = {'source'})

        function test_each_transducer_keeps_its_drive(tc)
            p = make_minimal_parameters(fullfile(tempdir, 'prestus_test_source_multi'));
            p.grid.dims = p.grid.default_dims;
            p.io.dir_cache = fullfile(tempdir, 'prestus_test_source_multi');
            if ~isfolder(p.io.dir_cache), mkdir(p.io.dir_cache); end
            p.io.overwrite_files = 'always';
            p.io.save_source_matrices = 0;
            t2 = p.transducer;
            half = round(p.grid.default_dims / 2);
            t2.trans_pos = [half(1), half(2), p.grid.default_dims(3) - 4];
            t2.focus_pos = [half(1), half(2), half(3) + 6];
            t2.annular.elem_p_amp = 60000;
            t2.annular.elem_phase_rad = 1.0;
            p.transducer(2) = t2;
            p.transducer(1).annular.elem_p_amp = 0;   % only transducer 2 drives

            [~, source, ~, ~, w] = source_sensor_setup(p, 1500, p.transducer(1).trans_pos, ...
                p.transducer(1).focus_pos, [], 1500);

            tc.verifyEqual(w.drive(:), [0; 60000 * exp(1i * 1.0)], 'AbsTol', 1e-6);
            P = source.p(:, 1:20);
            tc.verifyEqual(nnz(any(P ~= 0, 2)), nnz(w.W(:, 2)), ...
                'Only transducer 2''s source points should carry a signal.');
        end

    end
end
