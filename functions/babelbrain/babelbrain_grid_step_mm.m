function [step_mm, smallest_sos_ms] = babelbrain_grid_step_mm(parameters, freq_hz, ppw)
% BABELBRAIN_GRID_STEP_MM  BabelBrain's own spatial step for a run
%
% BabelBrain simulates at SpatialStep = SmallestSOS / freq / PPW. SmallestSOS is
% the smallest positive sound speed in its material table, floored by its
% Pichardo shear speed at 1000 kg/m^3 (1000 * 0.422 + 680.515 = 1102.515 m/s;
% GetSmallestSOS / DensityToSSoSPichardo in babel_integration_helpers.py). With
% PRESTUS's fluid-bone model that floor is the binding term, so the medium's
% sound speeds only matter when one is slower.
%
% Use as:
%   [step_mm, smallest_sos_ms] = babelbrain_grid_step_mm(parameters, freq_hz, ppw)
%   [~, smallest_sos_ms] = babelbrain_grid_step_mm(parameters)
%
% Input:
%   parameters - PRESTUS parameters (medium_properties.*.sound_speed)
%   freq_hz    - drive frequency [Hz]; [] returns only smallest_sos_ms (step_mm = NaN)
%   ppw        - points per wavelength; [] as for freq_hz
%
% Output:
%   step_mm         - BabelBrain's spatial step [mm]
%   smallest_sos_ms - the SmallestSOS it uses [m/s]
%
% See also: BABELBRAIN_MATCH_GRID, WRITE_BABELBRAIN_INPUT, BABELBRAIN_WRAPPER

arguments
    parameters (1,1) struct
    freq_hz = []
    ppw = []
end

    smallest_sos_ms = 1000.0 * 0.422 + 680.515;
    if isfield(parameters, 'medium_properties') && isstruct(parameters.medium_properties)
        names = fieldnames(parameters.medium_properties);
        for i = 1:numel(names)
            mp = parameters.medium_properties.(names{i});
            if isstruct(mp) && isfield(mp, 'sound_speed') && isnumeric(mp.sound_speed) && ...
                    isscalar(mp.sound_speed) && mp.sound_speed > 0
                smallest_sos_ms = min(smallest_sos_ms, mp.sound_speed);
            end
        end
    end
    step_mm = NaN;
    if ~isempty(freq_hz) && ~isempty(ppw)
        step_mm = smallest_sos_ms / freq_hz / ppw * 1e3;
    end
end
