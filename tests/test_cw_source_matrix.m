classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_cw_source_matrix < matlab.unittest.TestCase
% TEST_CW_SOURCE_MATRIX  Unit tests for functions/source/CWSourceMatrix.m
%
% Checks that the lazy rank-2 CW source behaves like the equivalent dense
% num_src x Nt matrix for every operation k-Wave applies to source.p
% (indexing in the time loop and input checks, and the scaling in
% kspaceFirstOrder_scaleSourceTerms).
%
%   Run with:   results = runtests('tests/test_cw_source_matrix.m');

    properties
        Lazy    % CWSourceMatrix
        Dense   % equivalent dense matrix
    end

    methods (TestMethodSetup)
        function make_source(tc)
            rng(1);
            n_src = 7; Nt = 50; omega = 2*pi*250e3; dt = 1e-7;
            C = randn(n_src, 1) + 1i * randn(n_src, 1);
            tc.Lazy  = CWSourceMatrix(real(C), imag(C), omega, dt, Nt);
            t = (0:Nt-1) * dt;
            tc.Dense = single(real(C)) * single(cos(omega*t)) ...
                     - single(imag(C)) * single(sin(omega*t));
        end
    end

    % ------------------------------------------------------------------ %
    %% Indexing
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'cw_source'})

        function test_size_and_rows(tc)
            tc.verifyEqual(size(tc.Lazy), size(tc.Dense));
            tc.verifyEqual(length(tc.Lazy(:, 1)), size(tc.Dense, 1));
            tc.verifyEqual(length(tc.Lazy(1, :)), size(tc.Dense, 2));
        end

        function test_time_column_matches_dense(tc)
            verify_close(tc, tc.Lazy(:, 17), tc.Dense(:, 17), 'AbsTol', 1e-5);
            verify_close(tc, tc.Lazy([2 5], 3), tc.Dense([2 5], 3), 'AbsTol', 1e-5);
        end

        function test_factors_reproduce_dense(tc)
            % jwave_wrapper writes the source as coeffs * basis
            [coeffs, basis] = factors(tc.Lazy);
            tc.verifyEqual(size(coeffs), [size(tc.Dense, 1), 2]);
            tc.verifyEqual(size(basis),  [2, size(tc.Dense, 2)]);
            verify_close(tc, single(coeffs * basis), tc.Dense, 'AbsTol', 1e-5);
        end

        function test_full_rows_match_dense(tc)
            verify_close(tc, tc.Lazy(3, :), tc.Dense(3, :), 'AbsTol', 1e-5);
            verify_close(tc, full(tc.Lazy), tc.Dense, 'AbsTol', 1e-5);
        end

    end

    % ------------------------------------------------------------------ %
    %% Scaling (kspaceFirstOrder_scaleSourceTerms)
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'cw_source'})

        function test_scalar_times_and_divide(tc)
            p = tc.Lazy .* 0.3;
            tc.verifyClass(p, 'CWSourceMatrix');
            verify_close(tc, full(p), tc.Dense * 0.3, 'AbsTol', 1e-5);
            p = tc.Lazy ./ 4;
            verify_close(tc, full(p), tc.Dense / 4, 'AbsTol', 1e-5);
            p = 2 .* tc.Lazy;
            verify_close(tc, full(p), tc.Dense * 2, 'AbsTol', 1e-5);
        end

        function test_column_factor_scales_rows(tc)
            k = (1:7)' / 3;
            verify_close(tc, full(tc.Lazy .* k), tc.Dense .* single(k), 'AbsTol', 1e-5);
        end

        function test_per_row_loop_as_in_kwave(tc)
            % Mirrors the heterogeneous-medium branch of
            % kspaceFirstOrder_scaleSourceTerms.
            c0 = linspace(1400, 2800, 7);
            p  = tc.Lazy;
            d  = tc.Dense;
            for p_index = 1:length(p(:, 1))
                s = 2e-7 / (3 * c0(p_index) * 5e-4);
                p(p_index, :) = p(p_index, :) .* s;
                d(p_index, :) = d(p_index, :) .* s;
            end
            tc.verifyClass(p, 'CWSourceMatrix');
            rel_err = norm(double(full(p)) - double(d), 'fro') / norm(double(d), 'fro');
            tc.verifyLessThan(rel_err, 1e-5);
        end

        function test_row_58_is_not_colon(tc)
            % ':' has character code 58; row index 58 must not be read or
            % assigned as "all rows".
            rng(2);
            C = randn(64, 1) + 1i * randn(64, 1);
            p = CWSourceMatrix(real(C), imag(C), 2*pi*250e3, 1e-7, 30);
            d = full(p);
            tc.verifyEqual(size(p(58, :)), [1 30]);
            verify_close(tc, p(58, :), d(58, :), 'AbsTol', 1e-5);
            for p_index = 1:64
                p(p_index, :) = p(p_index, :) .* 0.5;
            end
            verify_close(tc, full(p), d * 0.5, 'AbsTol', 1e-5);
        end

        function test_non_scaling_assignment_errors(tc)
            p = tc.Lazy;
            function assign_other()
                p(2, :) = fliplr(p(2, :));
            end
            tc.verifyError(@assign_other, 'CWSourceMatrix:notAScaling');
        end

        function test_unsupported_operand_errors(tc)
            tc.verifyError(@() tc.Lazy .* ones(3, 3), 'CWSourceMatrix:unsupportedOperation');
            tc.verifyError(@() 2 ./ tc.Lazy, 'CWSourceMatrix:unsupportedOperation');
        end

    end

end

% ====================================================================== %
%% Local helpers
% ====================================================================== %

function verify_close(tc, actual, expected, tol_type, tol)
% Compare in double: verifyEqual ignores numeric tolerances for singles.
    tc.verifyEqual(double(actual), double(expected), tol_type, tol);
end
