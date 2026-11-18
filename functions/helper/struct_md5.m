function hex = struct_md5(s)
% STRUCT_MD5  MD5 of a value's content, independent of struct field order
%
% Struct fields are sorted (recursively, also inside struct arrays and cells)
% and the value is encoded as JSON before hashing, so two structs with the
% same fields and values give the same hash whatever order their fields were
% set in. NaN and Inf encode as null.
%
% Use as:
%   hex = struct_md5(s)
%
% Input:
%   s   - struct, cell, numeric, logical, or char value (no objects)
%
% Output:
%   hex - 32-character lowercase hexadecimal MD5
%
% See also: ACOUSTIC_STORE_KEY, SEGMENTATION_FINGERPRINT

    txt = jsonencode(canonical(s));
    md  = java.security.MessageDigest.getInstance('MD5');
    md.update(uint8(unicode2native(txt, 'UTF-8')));
    bytes = typecast(md.digest(), 'uint8');
    hex = lower(reshape(dec2hex(bytes, 2)', 1, []));
end

function v = canonical(v)
    if isstruct(v)
        if ~isempty(fieldnames(v))
            v = orderfields(v);
        end
        f = fieldnames(v);
        for i = 1:numel(v)
            for k = 1:numel(f)
                v(i).(f{k}) = canonical(v(i).(f{k}));
            end
        end
    elseif iscell(v)
        for i = 1:numel(v)
            v{i} = canonical(v{i});
        end
    elseif isstring(v)
        v = cellstr(v);
    end
end
