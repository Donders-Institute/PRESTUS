function fp = segmentation_fingerprint(parameters)
% SEGMENTATION_FINGERPRINT  Content identity of the SimNIBS segmentation a stage reads
%
% MD5 checksum and size of <path.seg>/m2m_sub-XXX/final_tissues.nii.gz.
% Stored in the signatures of caches and stage outputs built from the
% segmentation (CACHE_SIGNATURE, PLACEMENT_SIGNATURE, PCT_SIGNATURE), so a
% re-run of charm (another SimNIBS version or other settings) is detected on
% reuse. The checksum follows the file's content, not its date, so copying or
% moving the data keeps it. Checksums are kept per file, size and date for
% the MATLAB session, so the file is read once.
%
% Use as:
%   fp = segmentation_fingerprint(parameters)
%   fp = segmentation_fingerprint(file)
%
% Input:
%   parameters - PRESTUS parameter struct (subject_id, path.seg)
%   file       - or the final_tissues.nii.gz file itself
%
% Output:
%   fp - struct with fields md5 (hex) and bytes; [] when there is no segmentation
%
% See also: SEGMENTATION_FILE, CACHE_SIGNATURE, CHECK_CACHE_SIGNATURE

    fp = [];
    if ischar(parameters) || isstring(parameters)
        f = char(parameters);
    else
        f = segmentation_file(parameters);
    end
    if isempty(f)
        return
    end
    d = dir(f);
    if isempty(d)
        return
    end

    persistent known
    if isempty(known)
        known = containers.Map('KeyType', 'char', 'ValueType', 'char');
    end
    key = sprintf('%s|%d|%.10f', f, d(1).bytes, d(1).datenum);
    if isKey(known, key)
        md5 = known(key);
    else
        md5 = file_md5(f);
        known(key) = md5;
    end
    fp = struct('md5', md5, 'bytes', d(1).bytes);
end


function hex = file_md5(file)
    md = java.security.MessageDigest.getInstance('MD5');
    fid = fopen(file, 'r');
    if fid < 0
        error('segmentation_fingerprint:read', 'Cannot read %s', file);
    end
    closer = onCleanup(@() fclose(fid));
    while true
        chunk = fread(fid, 2^24, '*uint8');
        if isempty(chunk)
            break
        end
        md.update(typecast(chunk, 'int8'));
    end
    hex = sprintf('%02x', typecast(md.digest(), 'uint8'));
end
