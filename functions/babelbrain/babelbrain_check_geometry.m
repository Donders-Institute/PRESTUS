function babelbrain_check_geometry(parameters, trans_pos, focus_pos)

% BABELBRAIN_CHECK_GEOMETRY  Error for geometries the BabelBrain backend cannot simulate
%
% BabelBrain simulates one transducer per run, with its beam along the z axis
% of the domain: in its own workflow it resamples the head into a grid aligned
% to the trajectory. PRESTUS hands it the PRESTUS simulation grid unrotated
% (WRITE_BABELBRAIN_INPUT only resamples and flips the depth axis), and
% BabelBrain then places the transducer on that grid's z axis through the
% focus marker. A beam that is not along z would therefore be simulated in the
% wrong direction, aimed at the right focus. Two checks follow from this:
%
%   1. One transducer only. A shared grid for several transducers puts at
%      least one beam off z (grid.mode 'focal_axis' tilts all of them), and
%      BabelBrain's thermal step runs on its own single-transducer acoustic
%      output, so it cannot take a combined heat source either.
%      Multi-transducer runs (coherent or async) need acoustic_backend
%      'kwave' or 'jwave' with thermal_backend 'kwave'.
%   2. Beam along the grid z axis (within 2 degrees), checked when trans_pos
%      and focus_pos are given: excludes e.g. grid.mode 'ras_plus'.
%
% Use as:
%   babelbrain_check_geometry(parameters)                       % check 1
%   babelbrain_check_geometry(parameters, trans_pos, focus_pos) % checks 1 and 2
%
% Input:
%   parameters - PRESTUS config (uses .transducer)
%   trans_pos  - (optional) [1x3] transducer position, grid voxels
%   focus_pos  - (optional) [1x3] focus position, grid voxels
%
% See also: BABELBRAIN_WRAPPER, WRITE_BABELBRAIN_INPUT, LOAD_PARAMETERS

    arguments
        parameters (1,1) struct
        trans_pos  = []
        focus_pos  = []
    end

    max_tilt_deg = 2;

    n_tx = 0;
    if isfield(parameters, 'transducer'), n_tx = numel(parameters.transducer); end
    if n_tx > 1
        error('PRESTUS:babelbrain:multiTransducer', ...
            ['The BabelBrain backend supports one transducer per run (config has %d). ' ...
             'BabelBrain simulates a beam only along the z axis of its domain, and a grid ' ...
             'shared by several transducers puts at least one beam off that axis. Use ' ...
             'acoustic_backend "kwave" or "jwave" and thermal_backend "kwave" for ' ...
             'multi-transducer runs. See doc_babelbrain.md, "Geometry limits".'], n_tx);
    end

    if numel(trans_pos) == 3 && numel(focus_pos) == 3
        u = double(focus_pos(:)' - trans_pos(:)');
        if norm(u) == 0
            error('PRESTUS:babelbrain:zeroBeam', 'trans_pos and focus_pos coincide.');
        end
        tilt = acosd(min(1, abs(u(3)) / norm(u)));
        if tilt > max_tilt_deg
            mode = 'focal_axis';
            if isfield(parameters, 'grid') && isfield(parameters.grid, 'mode'), mode = parameters.grid.mode; end
            error('PRESTUS:babelbrain:obliqueBeam', ...
                ['The beam is %.1f deg off the grid z axis (grid.mode "%s"). BabelBrain would ' ...
                 'simulate it along z instead. Use grid.mode "focal_axis" or ' ...
                 '"transducer_axis" (one transducer: beam along z). See doc_babelbrain.md, ' ...
                 '"Geometry limits".'], tilt, mode);
        end
    end
end
