classdef prestus_paths < matlab.unittest.fixtures.Fixture
% PRESTUS_PATHS  Shared test fixture: PRESTUS functions, config, externals and test fixtures on the path.
%
% Test classes declare it once and need no addpath of their own:
%
%   classdef (SharedTestFixtures = {prestus_paths}) test_x < matlab.unittest.TestCase
%
% The path is set up once per test run (shared across classes) and restored
% afterwards, so a test run leaves the MATLAB path as it found it.
%
% See also: RUN_ALL_TESTS

    methods
        function setup(fixture)
            here = fileparts(mfilename('fullpath'));
            repo = fileparts(here);
            old_path = path;
            fixture.addTeardown(@() path(old_path));
            addpath(genpath(fullfile(repo, 'functions')));
            addpath(fullfile(repo, 'config'));
            addpath(genpath(fullfile(repo, 'external')));
            addpath(fullfile(here, 'fixtures'));
            fixture.SetupDescription = 'Add PRESTUS functions, config, external and test fixtures to the path';
            fixture.TeardownDescription = 'Restore the MATLAB path';
        end
    end
end
