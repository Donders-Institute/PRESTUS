function s = cast_struct(s, cast_type)
% CAST_STRUCT  Recursively cast all numeric fields in a struct to a given type
%
%   Traverses every field of a (possibly nested) scalar struct and casts
%   numeric arrays to the requested MATLAB type. GPU arrays are gathered
%   to the CPU before casting.
%
% Use as:
%   s = cast_struct(s, cast_type)
%
% Input:
%   s         - [1x1] struct, may contain nested structs and numeric arrays
%   cast_type - target MATLAB type name (e.g., 'single', 'double')
%
% Output:
%   s         - [1x1] struct with all numeric fields cast to cast_type
%
% See also: CAST, GATHER

arguments
    s         (1,1) struct
    cast_type (1,:) char
end

    fn = fieldnames(s);
    for i = 1:length(fn)
        if isstruct(s.(fn{i}))
            s.(fn{i}) = cast_struct(s.(fn{i}), cast_type);
        elseif isnumeric(s.(fn{i}))
            s.(fn{i}) = gather(cast(s.(fn{i}), cast_type));
        end
    end
end