function result = isNull(x)
%YAML.ISNULL Check if value is yaml.Null
%   TF = YAML.ISNULL(X) returns true if X is a scalar yaml.Null object.

result = isequal(x, yaml.Null);

end
