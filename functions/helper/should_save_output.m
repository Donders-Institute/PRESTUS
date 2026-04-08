function tf = should_save_output(io, specific_field)
% SHOULD_SAVE_OUTPUT  Decide whether a particular output file should be saved.
%
% Precedence (highest to lowest):
%   1. io.<specific_field>   — per-step flag (e.g. save_source_matrices)
%   2. io.save_matrices      — global flag covering all simulation outputs
%   3. true                  — default: save everything
%
% Usage:
%   tf = should_save_output(parameters.io, 'save_source_matrices')
%   tf = should_save_output(parameters.io, 'save_acoustic_matrices')
%   tf = should_save_output(parameters.io, 'save_thermal_matrices')
%   tf = should_save_output(parameters.io, 'save_grid_cache')
%
% Inputs:
%   io             - parameters.io struct
%   specific_field - name of the per-step flag field (char)
%
% Output:
%   tf  - logical scalar; true = save, false = skip

    if isfield(io, specific_field)
        tf = logical(io.(specific_field));
    elseif isfield(io, 'save_matrices')
        tf = logical(io.save_matrices);
    else
        tf = true;
    end
end
