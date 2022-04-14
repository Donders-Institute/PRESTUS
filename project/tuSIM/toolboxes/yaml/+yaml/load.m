function result = load(s, options)
%YAML.LOAD Parse YAML string
%   DATA = YAML.LOAD(STR) parses a YAML string STR and converts it to
%   appropriate data types DATA.
%
%   DATA = YAML.LOAD(STR, "ConvertToArray", true) additionally converts
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
%       >> STR = "{a: 1, b: [text, false]}";
%       >> DATA = yaml.load(STR)
%
%         struct with fields:
%           a: 1
%           b: {["text"]  [0]}
%
%   See also YAML.LOADFILE, YAML.DUMP, YAML.DUMPFILE, YAML.ISNULL

arguments
    s (1, 1) string
    options.ConvertToArray (1, 1) logical = false
end

initSnakeYaml
import org.yaml.snakeyaml.*;
try
    rootNode = Yaml().load(s);
catch cause
    MException("yaml:load:Failed", "Failed to load YAML string.").addCause(cause).throw
end

try
    result = convert(rootNode);
catch exc
    if startsWith(exc.identifier, "yaml:load:")
        error(exc.identifier, exc.message);
    end
    exc.rethrow;
end

    function result = convert(node)
        switch class(node)
            case "double"
                if ~isempty(node)
                    result = node;
                else
                    result = yaml.Null;
                end
            case "char"
                result = string(node);
            case "logical"
                result = logical(node);
            case "java.util.LinkedHashMap"
                result = convertMap(node);
            case "java.util.ArrayList"
                result = convertList(node);
            case "java.util.Date"
                long = node.getTime;
                result = datetime(long, "ConvertFrom", "epochtime", "TicksPerSecond", 1000, "TimeZone", "UTC", "Format", "dd-MMM-uuuu HH:mm:ss.SSS z");
            otherwise
                error("yaml:load:TypeNotSupported", "Data type '%s' is not supported.", class(node))
        end
    end

    function result = convertMap(map)
        result = struct();

        keys = string(map.keySet().toArray())';
        fieldNames = matlab.lang.makeValidName(keys);
        fieldNames = matlab.lang.makeUniqueStrings(fieldNames);

        for i = 1:map.size()
            value = map.get(java.lang.String(keys(i)));
            result.(fieldNames(i)) = convert(value);
        end
    end

    function result = convertList(list)

        % Convert Java list to cell array.
        result = cell(1, list.size());
        for i = 1:list.size()
            result{i} = convert(list.get(i - 1));
        end

        if ~options.ConvertToArray
            return; end

        % Convert to non-cell array if possible
        if isempty(result)
            result = [];
            return
        elseif ~elementsHaveConsistentType(result)
            return
        elseif elementsAreScalar(result)
            result = horzcat(result{:});
        elseif elementsAreRowOrEmpty(result) && elementsHaveConsistentLength(result)
            result = vertcat(result{:});
        end
    end
end

function initSnakeYaml
snakeYamlFile = fullfile(fileparts(mfilename('fullpath')), 'snakeyaml', 'snakeyaml-1.30.jar');
if ~ismember(snakeYamlFile, javaclasspath('-dynamic'))
    javaaddpath(snakeYamlFile);
end
end

function result = elementsHaveConsistentType(c)
result = all(cellfun(@(x) strcmp(class(x), class(c{1})), c));
end

function result = elementsAreScalar(c)
result = all(cellfun(@isscalar, c));
end

function result = elementsAreRowOrEmpty(c)
result = all(cellfun(@(x) isrow(x) || isempty(x), c));
end

function result = elementsHaveConsistentLength(c)
result = all(cellfun(@(x) length(x) == length(c{1}), c));
end
