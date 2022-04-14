function result = dump(data, style)
%DUMP Convert data to YAML string
%   STR = YAML.DUMP(DATA) converts DATA to a YAML string STR.
%
%   STR = YAML.DUMP(DATA, STYLE) uses a specific output style.
%   STYLE can be "auto" (default), "block" or "flow".
%
%   The following types are supported for DATA:
%       MATLAB type             | YAML type
%       ------------------------|----------------------
%       1D cell array           | Sequence
%       1D non-scalar array     | Sequence
%       2D/3D cell array        | Nested sequences
%       2D/3D non-scalar array  | Nested sequences
%       struct                  | Mapping
%       scalar single/double    | Floating-point number
%       scalar int8/../int64    | Integer
%       scalar logical          | Boolean
%       scalar string           | String
%       char vector             | String
%       scalar yaml.Null        | null
%
%   Array conversion can be ambiguous. To ensure consistent conversion
%   behaviour, consider manually converting array data to nested 1D cells
%   before converting it to YAML.
%
%   Example:
%       >> DATA.a = 1
%       >> DATA.b = {"text", false}
%       >> STR = yaml.dump(DATA)
%
%         "a: 1.0
%         b: [text, false]
%         "
%
%   See also YAML.DUMPFILE, YAML.LOAD, YAML.LOADFILE, YAML.ISNULL

arguments
    data
    style {mustBeMember(style, ["flow", "block", "auto"])} = "auto"
end

NULL_PLACEHOLDER = "$%&?"; % Should have 4 characters for correct line breaks.

initSnakeYaml
import org.yaml.snakeyaml.*;

try
    javaData = convert(data);
catch exception
    if string(exception.identifier).startsWith("yaml:dump")
        error(exception.identifier, exception.message);
    end
    exception.rethrow;
end
dumperOptions = getDumperOptions(style);
result = Yaml(dumperOptions).dump(javaData);
result = string(result).replace(NULL_PLACEHOLDER, "null");

    function result = convert(data)
        if iscell(data)
            result = convertCell(data);
        elseif ischar(data) && isvector(data)
            result = convertString(data);
        elseif ~isscalar(data)
            result = convertArray(data);
        elseif isstruct(data)
            result = convertStruct(data);
        elseif isfloat(data)
            result = java.lang.Double(data);
        elseif isinteger(data)
            result = java.lang.Integer(data);
        elseif islogical(data)
            result = java.lang.Boolean(data);
        elseif isstring(data)
            result = convertString(data);
        elseif yaml.isNull(data)
            result = java.lang.String(NULL_PLACEHOLDER);
        else
            error("yaml:dump:TypeNotSupported", "Data type '%s' is not supported.", class(data))
        end
    end

    function result = convertString(data)
        if contains(data, NULL_PLACEHOLDER)
            error("yaml:dump:NullPlaceholderNotAllowed", "Strings must not contain '%s' since it is used as a placeholder for null values.", NULL_PLACEHOLDER)
        end
        result = java.lang.String(data);
    end

    function result = convertStruct(data)
        result = java.util.LinkedHashMap();
        for key = string(fieldnames(data))'
            value = convert(data.(key));
            result.put(key, value);
        end
    end

    function result = convertCell(data)
        data = nest(data);
        result = java.util.ArrayList();
        for i = 1:length(data)
            result.add(convert(data{i}));
        end
    end

    function result = convertArray(data)
        result = convertCell(num2cell(data));
    end

    function result = nest(data)
        if isvector(data) || isempty(data)
            result = data;
            return
        end
        n = size(data, 1);
        nDimensions = length(size(data));
        result = cell(1, n);
        if nDimensions == 2
            for i = 1:n
                result{i} = data(i, :);
            end
        elseif nDimensions == 3
            for i = 1:n
                result{i} = squeeze(data(i, :, :));
            end
        else
            error("yaml:dump:HigherDimensionsNotSupported", "Arrays with more than three dimensions are not supported. Use nested cells instead.")
        end
    end

    function initSnakeYaml
        snakeYamlFile = fullfile(fileparts(mfilename('fullpath')), 'snakeyaml', 'snakeyaml-1.30.jar');
        if ~ismember(snakeYamlFile, javaclasspath('-dynamic'))
            javaaddpath(snakeYamlFile);
        end
    end

    function opts = getDumperOptions(style)
        import org.yaml.snakeyaml.*;
        opts = DumperOptions();
        classes = opts.getClass.getClasses;
        styleFields = classes(4).getDeclaredFields();
        styleIndex = find(style == ["flow", "block", "auto"]);
        opts.setDefaultFlowStyle(styleFields(styleIndex).get([]));
    end

end
