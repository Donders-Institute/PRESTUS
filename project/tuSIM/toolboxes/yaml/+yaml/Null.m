classdef Null
    % YAML.NULL Represents a YAML null value

    methods
        function obj = Null(n, m)
            arguments
                n (1, 1) double {mustBeInteger, mustBePositive} = 1
                m (1, 1) double {mustBeInteger, mustBePositive} = n
            end

            obj(n, m) = obj;
        end
    end
end
