function assemble_pfov_fields(files_in, file_out, targets_wcm2)
% ASSEMBLE_PFOV_FIELDS  Deprecated alias of ASSEMBLE_LIMITED_FOV_FIELDS.
%
% The pFOV ("partial field of view") naming was replaced by "limited FOV". This
% alias forwards all arguments unchanged so existing scripts keep working; it
% will be removed in a future release.
%
% Use as:
%   assemble_pfov_fields(files_in, file_out)
%   assemble_pfov_fields(files_in, file_out, targets_wcm2)
%
% See also: ASSEMBLE_LIMITED_FOV_FIELDS

    warning('prestus:deprecated:assemble_pfov_fields', ...
        'assemble_pfov_fields is deprecated; use assemble_limited_fov_fields (same arguments).');
    if nargin < 3
        assemble_limited_fov_fields(files_in, file_out);
    else
        assemble_limited_fov_fields(files_in, file_out, targets_wcm2);
    end
end
