function [sidecar_file, json_file] = find_sidecar(data_file)
% FIND_SIDECAR  JSON sidecar of a data file: same name with .json (BIDS convention)
%
% For 'x/sub-001_pressure.nii.gz' the sidecar is 'x/sub-001_pressure.json'.
%
% Use as:
%   [sidecar_file, json_file] = find_sidecar(data_file)
%
% Input:
%   data_file - path to the data file (.nii, .nii.gz, .mat, .h5)
%
% Output:
%   sidecar_file - the sidecar when it exists, else ''
%   json_file    - where the sidecar of DATA_FILE belongs (for writing)
%
% See also: WRITE_SIDECAR, READ_SIDECAR

arguments
    data_file {mustBeTextScalar}
end

    stem = regexprep(char(data_file), '(\.nii\.gz|\.nii|\.mat|\.h5)$', '', 'ignorecase');
    json_file = [stem '.json'];
    sidecar_file = '';
    if isfile(json_file)
        sidecar_file = json_file;
    end
end
