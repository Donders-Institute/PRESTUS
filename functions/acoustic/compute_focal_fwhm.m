function fw = compute_focal_fwhm(intensity, peak_pos, beam_vec, res_mm)
% COMPUTE_FOCAL_FWHM  Axial and lateral focal widths at half intensity (−3 dB)
%
% Measures the full width at half maximum of the focal intensity through the
% focal peak: along the transducer axis (axial) and along directions in the
% plane perpendicular to it (lateral). Profiles are sampled by linear
% interpolation every quarter voxel, and each −3 dB crossing (I = I_peak/2) is
% placed by linear interpolation between the two samples that bracket it, so
% the widths resolve below the grid spacing. The lateral width is measured
% every 5° around the axis; the narrowest and widest are returned, which
% coincide for a round focus.
%
% Use as:
%   fw = compute_focal_fwhm(intensity, peak_pos, beam_vec, res_mm)
%
% Input:
%   intensity - (:,:) or (:,:,:) intensity map [W/cm²] (any map proportional
%               to it gives the same widths, e.g. Itp or Ipa)
%   peak_pos  - [1×nd] voxel position of the focal peak (e.g. the
%               brain-masked peak of a layered run)
%   beam_vec  - [1×nd] transducer axis in voxels (focus_pos - trans_pos)
%   res_mm    - scalar grid resolution [mm/voxel] (isotropic)
%
% Output:
%   fw  struct with fields [mm]:
%     axial_mm        FWHM along the transducer axis
%     lateral_min_mm  narrowest FWHM perpendicular to the axis
%     lateral_max_mm  widest FWHM perpendicular to the axis
%   A width is NaN when the profile does not fall below half the peak on
%   both sides inside the grid (focal lobe truncated; enlarge the FOV).
%   On a 2D grid the single in-plane perpendicular gives both lateral fields.
%
% See also: COMPUTE_FOCAL_METRICS, ACOUSTIC_ANALYSIS

arguments
    intensity  {mustBeNumeric}
    peak_pos   (1,:) double
    beam_vec   (1,:) double
    res_mm     (1,1) double {mustBePositive}
end

    fw = struct('axial_mm', NaN, 'lateral_min_mm', NaN, 'lateral_max_mm', NaN);
    nd = ndims(intensity);
    peak_pos = peak_pos(1:nd);
    beam_vec = beam_vec(1:nd);
    if norm(beam_vec) == 0 || any(~isfinite(peak_pos))
        return
    end
    u = beam_vec(:)' / norm(beam_vec);
    intensity = double(intensity);
    thr = interp_at(intensity, peak_pos) / 2;
    if ~(thr > 0)
        return
    end

    fw.axial_mm = profile_fwhm(intensity, peak_pos, u, thr) * res_mm;

    if nd == 3
        [~, k] = min(abs(u));                  % grid axis least aligned with u
        e = zeros(1, 3); e(k) = 1;
        v1 = cross(u, e); v1 = v1 / norm(v1);
        v2 = cross(u, v1);
        angles = (0:5:175) * pi / 180;
        widths = NaN(size(angles));
        for a = 1:numel(angles)
            d = cos(angles(a)) * v1 + sin(angles(a)) * v2;
            widths(a) = profile_fwhm(intensity, peak_pos, d, thr);
        end
    else
        widths = profile_fwhm(intensity, peak_pos, [-u(2), u(1)], thr);
    end
    if all(isfinite(widths))
        fw.lateral_min_mm = min(widths) * res_mm;
        fw.lateral_max_mm = max(widths) * res_mm;
    end
end

function w = profile_fwhm(intensity, p0, d, thr)
% Width [voxels] between the first crossings of thr on either side of p0 along d.
    step = 0.25;
    n_max = ceil(norm(size(intensity)) / step);
    edges = NaN(1, 2);
    for side = 1:2
        sgn = 3 - 2 * side;                    % +1, then -1
        prev = interp_at(intensity, p0);
        for i = 1:n_max
            cur = interp_at(intensity, p0 + sgn * i * step * d);
            if isnan(cur)
                break                          % left the grid: truncated
            end
            if cur < thr
                frac = (prev - thr) / (prev - cur);
                edges(side) = (i - 1 + frac) * step;
                break
            end
            prev = cur;
        end
    end
    w = sum(edges);                            % NaN if either side is truncated
end

function v = interp_at(intensity, p)
    if ndims(intensity) == 3
        v = interpn(intensity, p(1), p(2), p(3), 'linear', NaN);
    else
        v = interpn(intensity, p(1), p(2), 'linear', NaN);
    end
end
