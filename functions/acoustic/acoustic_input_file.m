function file = acoustic_input_file(parameters)
% ACOUSTIC_INPUT_FILE  Pressure map a run starts from instead of simulating ('' if none)
%
% path.acoustic_input, or the deprecated io.external_acoustic_nifti that
% LOAD_PARAMETERS migrates to it (also read here for parameter structs built
% without LOAD_PARAMETERS).
%
% Use as:
%   file = acoustic_input_file(parameters)
%
% Input:
%   parameters - PRESTUS parameter struct
%
% Output:
%   file - char path, or '' when the acoustic stage simulates
%
% See also: LOAD_ACOUSTIC_INPUT, ACOUSTIC_STAGE

    file = '';
    if isfield(parameters, 'path') && isfield(parameters.path, 'acoustic_input') && ...
            ~isempty(parameters.path.acoustic_input)
        file = char(parameters.path.acoustic_input);
    elseif isfield(parameters, 'io') && isfield(parameters.io, 'external_acoustic_nifti') && ...
            ~isempty(parameters.io.external_acoustic_nifti)
        file = char(parameters.io.external_acoustic_nifti);
    end
end
