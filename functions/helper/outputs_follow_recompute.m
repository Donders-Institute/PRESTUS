function parameters = outputs_follow_recompute(parameters, stage)
% OUTPUTS_FOLLOW_RECOMPUTE  Let a recomputed stage replace its existing output files
%
% Under io.overwrite_files 'never' (or a declined 'ask'), existing NIfTI
% outputs are kept, while the CSV table of the same run is always rewritten.
% When the stage's results were recomputed in this run rather than loaded
% from cache, keeping the old NIfTIs would leave them describing an earlier
% run next to a CSV that describes this one. This switches the output
% writers of that stage to 'always' for this call, so all outputs of a run
% are written together.
%
% The acoustic and thermal stages record whether they recomputed in
% parameters.state.recomputed.<stage> (ACOUSTIC_STAGE, THERMAL_STAGE).
%
% Use as:
%   parameters = outputs_follow_recompute(parameters, 'acoustic')
%
% Input:
%   parameters - PRESTUS config
%   stage      - 'acoustic' or 'thermal'
%
% Output:
%   parameters - with io.overwrite_files = 'always' when the stage recomputed
%
% See also: CONFIRM_OVERWRITING, NIFTI_ACOUSTIC, NIFTI_THERMAL

    if ~isfield(parameters, 'state') || ~isfield(parameters.state, 'recomputed') || ...
            ~isfield(parameters.state.recomputed, stage) || ~parameters.state.recomputed.(stage)
        return
    end
    if any(strcmp(parameters.io.overwrite_files, {'always', 'force'}))
        return
    end
    fprintf(['[outputs] io.overwrite_files is ''%s'', but the %s results were recomputed in this run: ' ...
             'replacing existing %s outputs so NIfTIs and CSV describe the same run.\n'], ...
            parameters.io.overwrite_files, stage, stage);
    parameters.io.overwrite_files = 'always';
end
