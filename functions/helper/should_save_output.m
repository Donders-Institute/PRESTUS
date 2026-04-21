function tf = should_save_output(io, specific_field)
% SHOULD_SAVE_OUTPUT  Decide whether a particular pipeline output should be saved
%
% Precedence (highest to lowest):
%   1. io.<specific_field>  — per-step flag (e.g., save_source_matrices)
%   2. io.save_matrices     — global flag covering all simulation outputs
%   3. true                 — default: save everything
%
% Use as:
%   tf = should_save_output(parameters.io, 'save_source_matrices')
%   tf = should_save_output(parameters.io, 'save_acoustic_matrices')
%   tf = should_save_output(parameters.io, 'save_thermal_matrices')
%   tf = should_save_output(parameters.io, 'save_grid_cache')
%
% Input:
%   io             - parameters.io struct
%   specific_field - name of the per-step flag field
%
% Output:
%   tf - logical scalar; true = save, false = skip
%
% See also: SIMULATION_NIFTI, PRESTUS_PIPELINE

arguments
    io             (1,1) struct
    specific_field (1,:) char
end

    if isfield(io, specific_field)
        tf = logical(io.(specific_field));
    elseif isfield(io, 'save_matrices')
        tf = logical(io.save_matrices);
    else
        tf = true;
    end
end
