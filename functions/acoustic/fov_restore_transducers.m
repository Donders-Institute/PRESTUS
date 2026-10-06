function parameters = fov_restore_transducers(parameters, fov_state)
% FOV_RESTORE_TRANSDUCERS  Put every transducer back at its full-grid position after an FOV crop.
%
% ACOUSTIC_GRID_FOV moves every transducer's trans_pos/focus_pos into FOV
% voxels. Stages that return to the full grid (NIfTI export, thermal) call
% this to undo that for all transducers; transducer(1) uses
% fov_state.trans_pos_full/focus_pos_full as before.
%
% Use as:
%   parameters = fov_restore_transducers(parameters, fov_state)
%
% Input:
%   parameters - PRESTUS config whose transducer positions are in FOV voxels
%   fov_state  - struct from ACOUSTIC_FOV_CROP (fields trans_pos_full,
%                focus_pos_full, and optionally transducer_pos_full)
%
% Output:
%   parameters - the same config with full-grid transducer positions
%
% See also: ACOUSTIC_FOV_CROP, ACOUSTIC_GRID_FOV

parameters.transducer(1).trans_pos = fov_state.trans_pos_full;
if isfield(parameters.transducer(1), 'focus_pos')
    parameters.transducer(1).focus_pos = fov_state.focus_pos_full;
end
if ~isfield(fov_state, 'transducer_pos_full') || isempty(fov_state.transducer_pos_full)
    return
end
for ti = 2:min(numel(parameters.transducer), numel(fov_state.transducer_pos_full))
    parameters.transducer(ti).trans_pos = fov_state.transducer_pos_full(ti).trans_pos;
    parameters.transducer(ti).focus_pos = fov_state.transducer_pos_full(ti).focus_pos;
end
end
