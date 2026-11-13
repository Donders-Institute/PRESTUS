function levels = prestus_test_levels()
% PRESTUS_TEST_LEVELS  Test levels of RUN_ALL_TESTS and the level tags each one runs.
%
% Every test carries exactly one level tag (class- or method-level TestTags);
% TEST_SUITE_WIRING checks this. A level runs every test whose level tag is
% in its list. 'standalone' marks tests no level runs (run them with runtests).
%
% Use as:
%   levels = prestus_test_levels()
%
% Output:
%   levels - containers.Map: level name -> cell array of level tags
%
% See also: RUN_ALL_TESTS

    main  = {'unit', 'water', 'head', 'acoustic', 'thermal'};
    demo0 = {'demo_inputs', 'demo_config', 'demo_localite', 'demo_pct', 'demo_head'};

    levels = containers.Map();
    levels('unit')          = main(1);
    levels('water')         = main(1:2);
    levels('head')          = main(1:3);
    levels('acoustic')      = main(1:4);
    levels('all')           = main;
    levels('physics')       = {'unit', 'physics_axisym'};
    levels('physics_3d')    = {'physics_3d'};
    levels('itrusst')       = {'itrusst'};
    levels('slurm_live')    = {'slurm_live'};
    levels('demo_inputs')   = demo0(1);
    levels('demo_config')   = demo0(1:2);
    levels('demo_localite') = demo0(1:3);
    levels('demo_pct')      = {'demo_inputs', 'demo_pct'};
    levels('demo_head')     = demo0;
    levels('demo_acoustic') = [demo0, {'demo_acoustic'}];
    levels('demo_thermal')  = {'demo_inputs', 'demo_config', 'demo_acoustic', 'demo_thermal'};
    levels('demo_golden')   = {'demo_inputs', 'demo_golden'};
end
