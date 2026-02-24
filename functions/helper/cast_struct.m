function s = cast_struct(s, cast_type)
    % Recursively cast all numeric fields in struct to given type
    fn = fieldnames(s);
    for i = 1:length(fn)
        if isstruct(s.(fn{i}))
            s.(fn{i}) = cast_struct(s.(fn{i}), cast_type);
        elseif isnumeric(s.(fn{i}))
            s.(fn{i}) = gather(cast(s.(fn{i}), cast_type));
        end
    end
end