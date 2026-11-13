function frac = focal_reference_to_frac(focal_reference)
% FOCAL_REFERENCE_TO_FRAC  Map a focal-depth convention name to an intensity threshold fraction
%
% Different manufacturers label transducer "focal depth" / steering-table
% settings using different axial-profile conventions. This maps the
% human-readable convention name (as set in calibration.focal_reference)
% to the intensity fraction of the profile peak used to find the FLHM
% (focal length at half maximum) span in COMPUTE_FLHM_CENTER / COMPUTE_FOCAL_METRICS.
%
% Use as:
%   frac = focal_reference_to_frac(focal_reference)
%
% Input:
%   focal_reference - 'flhm_3db' (intensity >= peak/2; ITRUSST/Martin et al.
%                       2025 default)
%                     | 'flhm_6db' (intensity >= peak/4, i.e. pressure >= peak/2;
%                       e.g. Imasonic/IGT steering tables)
%                     | 'peak' has no threshold and raises an error; callers
%                       handle it first (see COMPUTE_FLHM_CENTER,
%                       COMPUTE_FOCAL_METRICS)
%
% Output:
%   frac - intensity threshold as a fraction of peak, e.g. 0.5 or 0.25
%
% See also: COMPUTE_FLHM_CENTER, COMPUTE_FOCAL_METRICS

arguments
    focal_reference (1,:) char
end

    switch focal_reference
        case 'flhm_3db'
            frac = 0.5;
        case 'flhm_6db'
            frac = 0.25;
        case 'peak'
            error('focal_reference_to_frac:peakHasNoFrac', ...
                ['focal_reference_to_frac: ''peak'' has no associated threshold fraction -- ' ...
                 'callers must branch on focal_reference == ''peak'' before calling this function.']);
        otherwise
            error('focal_reference_to_frac:unknownReference', ...
                'focal_reference_to_frac: unknown focal_reference ''%s''; expected ''peak'', ''flhm_3db'', or ''flhm_6db''.', ...
                focal_reference);
    end

end
