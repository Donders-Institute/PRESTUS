function parameters = babelbrain_match_grid(parameters)
% BABELBRAIN_MATCH_GRID  Build the PRESTUS grid at BabelBrain's own step (babelbrain.match_grid)
%
% By default PRESTUS builds the head model on grid.resolution_mm, and
% WRITE_BABELBRAIN_INPUT resamples the labels and pseudo-CT onto BabelBrain's
% step (nearest neighbour), which can thicken thin layers such as the skull.
% With babelbrain.match_grid true, grid.resolution_mm is set to BabelBrain's
% step for this run (BABELBRAIN_GRID_STEP_MM at the transducer frequency and
% babelbrain.ppw), so the medium is built once on BabelBrain's grid and not
% resampled. For comparisons with native BabelBrain runs. Called at the start of
% PRESTUS_PIPELINE, before any stage reads grid.resolution_mm.
%
% Use as:
%   parameters = babelbrain_match_grid(parameters)
%
% Input:
%   parameters - PRESTUS parameters (babelbrain.match_grid, babelbrain.ppw,
%                acoustic_backend, transducer(1).freq_hz or its equipment entry)
%
% Output:
%   parameters - with grid.resolution_mm set to BabelBrain's step when
%                babelbrain.match_grid is true; unchanged otherwise
%
% See also: BABELBRAIN_GRID_STEP_MM, WRITE_BABELBRAIN_INPUT

arguments
    parameters (1,1) struct
end

    if ~logical(get_field(parameters, {'babelbrain', 'match_grid'}, false))
        return
    end
    if ~strcmpi(char(string(get_field(parameters, {'acoustic_backend'}, ''))), 'babelbrain')
        error('babelbrain_match_grid:backend', ...
            'babelbrain.match_grid needs acoustic_backend ''babelbrain''.');
    end
    freq_hz = [];
    if isfield(parameters, 'transducer') && isfield(parameters.transducer, 'freq_hz')
        freq_hz = parameters.transducer(1).freq_hz;
    end
    if isempty(freq_hz)
        % frequency from the equipment entry; the resolved transducer is not kept here
        resolved = load_transducer_parameters(parameters);
        freq_hz = resolved.transducer(1).freq_hz;
    end
    if isempty(freq_hz) || ~isnumeric(freq_hz) || freq_hz(1) <= 0
        error('babelbrain_match_grid:noFrequency', ...
            'babelbrain.match_grid needs the transducer frequency (transducer.freq_hz or its equipment entry).');
    end
    ppw = get_field(parameters, {'babelbrain', 'ppw'}, 12);
    [step_mm, sos] = babelbrain_grid_step_mm(parameters, freq_hz(1), ppw);
    old = parameters.grid.resolution_mm;
    parameters.grid.resolution_mm = step_mm;
    fprintf(['[babelbrain_match_grid] grid.resolution_mm %.4g -> %.6f mm: BabelBrain''s step ' ...
        '(%.3f m/s / %.0f Hz / %g PPW), so the medium is not resampled between grids.\n'], ...
        old, step_mm, sos, freq_hz(1), ppw);
end
