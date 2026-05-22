function result = loadFile(filePath, options)
% LOADFILE Read YAML file.
%   DATA = YAML.LOADFILE(FILE) reads a YAML file and converts it to
%   appropriate data types DATA.
%
%   DATA = YAML.LOADFILE(STR, "ConvertToArray", true) additionally converts
%   sequences to 1D or 2D non-cell arrays if possible.
%
%   The YAML types are convert to MATLAB types as follows:
%
%       YAML type                  | MATLAB type
%       ---------------------------|------------
%       Sequence                   | cell or array if possible and
%                                  | "ConvertToArray" is enabled
%       Mapping                    | struct
%       Floating-point number      | double
%       Integer                    | double
%       Boolean                    | logical
%       String                     | string
%       Date (yyyy-mm-ddTHH:MM:SS) | datetime
%       Date (yyyy-mm-dd)          | datetime
%       null                       | yaml.Null
%
%   Example:
%       >> DATA.a = 1
%       >> DATA.b = {"text", false}
%       >> FILE = ".\test.yaml"
%       >> yaml.dumpFile(FILE, DATA)
%       >> yaml.loadFile("test.yaml")
%
%         struct with fields:
%
%           a: 1
%           b: {["text"]  [0]}
%
%   See also YAML.LOAD, YAML.DUMP, YAML.DUMPFILE, YAML.ISNULL

arguments
    filePath (1, 1) string
    options.ConvertToArray (1, 1) logical = false
end

content = string(fileread(filePath));
result = yaml.load(content, "ConvertToArray", options.ConvertToArray);

end
