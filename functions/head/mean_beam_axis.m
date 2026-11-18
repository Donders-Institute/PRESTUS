function [axis_vec, axis_norm] = mean_beam_axis(beam_units, min_norm)
% MEAN_BEAM_AXIS  Mean axis of several transducer beams, taken as axes
%
% Averages the beams' unit vectors as axes (lines), not directions: a beam
% pointing against the first one is flipped before averaging. Beams that
% converge from the same side average to their bisector; opposed
% (bilateral) beams average to the bisector of the first beam's axis and
% the other beam's reversed axis, instead of cancelling. The result points
% along the first beam.
%
% Use as:
%   [axis_vec, axis_norm] = mean_beam_axis(beam_units)
%   [axis_vec, axis_norm] = mean_beam_axis(beam_units, min_norm)
%
% Input:
%   beam_units - [N x 3] beam unit vectors (transducer -> focus); row 1 is
%                the reference beam (transducer 1)
%   min_norm   - (optional, default 0.5) if the norm of the mean axis
%                vector is below this, warn and return beam 1's axis. For
%                two beams the norm is always >= cos(45 deg) = 0.71.
%
% Output:
%   axis_vec   - [1 x 3] unit vector of the mean axis
%   axis_norm  - norm of the mean axis vector before normalisation
%                (1 = all beams parallel or antiparallel)
%
% See also: PREPROC_ALIGN_TO_MEAN_FOCAL_AXIS

    arguments
        beam_units (:,3) double
        min_norm   (1,1) double = 0.5
    end

    u = beam_units ./ vecnorm(beam_units, 2, 2);
    flip = (u * u(1,:)') < 0;
    u(flip,:) = -u(flip,:);

    u_mean    = mean(u, 1);
    axis_norm = norm(u_mean);
    if axis_norm < min_norm
        warning('PRESTUS:gridAlignment:ambiguousFocalAxis', ...
            ['Mean beam axis norm %.3f is below threshold %.3f (beam axes heavily spread). ' ...
             'Falling back to the focal axis of transducer(1).'], axis_norm, min_norm);
        axis_vec = u(1,:);
    else
        axis_vec = u_mean / axis_norm;
    end
end
