# PRESTUS test suite

The levels, what each guarantees, required data and run times are documented in
[doc/doc_testing.md](../doc/doc_testing.md) (published under *Reference → Testing*).
This file only lists the commands.

```matlab
cd tests
run_all_tests                  % level 1: unit tests (default; no data needed)
run_all_tests('physics')       % + free-water physics benchmarks (2D, ~5 min CPU)
run_all_tests('itrusst')       % ITRUSST BM1 vs. reference (PRESTUS_ITRUSST_DATA)
run_all_tests('physics_3d')    % 3D O'Neil benchmark (GPU node)
run_all_tests('water')         % + 3D water pipeline and multi-ISPPA
run_all_tests('head')          % + subject smoke tests (PRESTUS_TEST_DATA)
run_all_tests('acoustic')
run_all_tests('all')
run_all_tests('demo_inputs')   % demo subject levels (PRESTUS_DEMO_DATA) ...
run_all_tests('demo_acoustic')
run_all_tests('demo_thermal')
run_all_tests('demo_golden')   % key metrics vs. stored reference
```

From a shell (CI): `matlab -batch "cd tests; run_all_tests('unit')"`.

Tests are selected by level tag (`prestus_test_levels.m`): give a new test class
`SharedTestFixtures = {prestus_paths}` and one level tag, e.g. `TestTags = {'unit'}`;
`test_suite_wiring` fails for a test without exactly one level tag.
