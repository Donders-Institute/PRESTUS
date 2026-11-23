function [P_ref, info] = remap_complex_field(P, M, out_size, opts)
% REMAP_COMPLEX_FIELD  Resample a complex pressure field onto another grid without losing its phase
%
% A complex steady-state field oscillates once per wavelength along the beam,
% so interpolating it directly loses amplitude between the samples: linear
% interpolation scales the amplitude by cos(k*h/2) at the midpoints (0.87 at
% 6 points per wavelength) and distorts the interference of fields that are
% summed afterwards. The carrier is removed first:
%
%   P(v) = E(v) * exp(-1i * k * s(v)),   s(v) = (v - Origin) . Axis
%
% with k the wavenumber per voxel and Axis the beam direction in the source
% grid's voxels. The envelope E varies slowly and is resampled under M
% (AFFINE_RESAMPLE_3D, real and imaginary part); the carrier is then put back
% at the source position of every output voxel, so it is exact for a plane
% wave along Axis. Where the field departs from that (skull, reflections,
% standing waves, other directions), the envelope varies faster and the
% interpolation error grows; check against a run on one shared grid.
%
% Use as:
%   P_ref = remap_complex_field(P, M, out_size, KPerVoxel=k, Axis=n, Origin=o)
%   P_ref = remap_complex_field(P, M, out_size, Demodulate=false)   % plain interpolation
%
% Input:
%   P          - complex array [D1 x D2 x D3] on the source grid
%   M          - 4x4 forward affine, source voxel -> output voxel, row-vector
%                convention ([out 1] = [in 1] * M), as COMBINE_ASYNC_INTENSITY
%                composes from the planimg sidecars
%   out_size   - [1x3] output grid size
%   KPerVoxel  - wavenumber times the source voxel size [rad per voxel]
%                (2*pi*f*dx/c, c the sound speed the carrier is removed with)
%   Axis       - [1x3] beam direction in source voxels (normalised here)
%   Origin     - [1x3] a point on the beam axis in source voxels (e.g. trans_pos)
%   Method     - interpolation, 'linear' (default) or 'cubic'
%   Demodulate - false: resample real and imaginary parts directly
%
% Output:
%   P_ref - complex array [out_size]; 0 outside the source grid
%   info  - struct with fields demodulated, k_per_voxel, axis, method
%
% See also: COMBINE_ASYNC_INTENSITY, AFFINE_RESAMPLE_3D

arguments
    P        {mustBeNumeric}
    M        (4,4) double
    out_size (1,3) double {mustBeInteger, mustBePositive}
    opts.KPerVoxel  (1,1) double = NaN
    opts.Axis       (1,3) double = [0 0 1]
    opts.Origin     (1,3) double = [1 1 1]
    opts.Method     (1,:) char {mustBeMember(opts.Method, {'linear', 'cubic'})} = 'linear'
    opts.Demodulate (1,1) logical = true
end

info = struct('demodulated', opts.Demodulate, 'k_per_voxel', opts.KPerVoxel, ...
              'axis', opts.Axis, 'method', opts.Method);
P = double(P);
if ~opts.Demodulate
    P_ref = complex(affine_resample_3d(real(P), M, out_size, opts.Method, 0), ...
                    affine_resample_3d(imag(P), M, out_size, opts.Method, 0));
    return
end
if ~isfinite(opts.KPerVoxel) || opts.KPerVoxel <= 0
    error('remap_complex_field:wavenumber', 'KPerVoxel must be a positive wavenumber per voxel.');
end
n = opts.Axis / norm(opts.Axis);
k = opts.KPerVoxel;
info.axis = n;

% Envelope on the source grid
sz = size(P, 1:3);
[x, y, z] = ndgrid(1:sz(1), 1:sz(2), 1:sz(3));
s_src = (x - opts.Origin(1)) * n(1) + (y - opts.Origin(2)) * n(2) + (z - opts.Origin(3)) * n(3);
clear x y z
E = P .* exp(1i * k * s_src);
clear s_src P
E_ref = complex(affine_resample_3d(real(E), M, out_size, opts.Method, 0), ...
                affine_resample_3d(imag(E), M, out_size, opts.Method, 0));
clear E

% Carrier at the source position of every output voxel: v_src = [v_ref 1] * inv(M),
% so s is linear in v_ref
Minv = affine_inverse(M);
g = Minv(1:3, 1:3) * n(:);
c = (Minv(4, 1:3) - opts.Origin) * n(:);
[x, y, z] = ndgrid(1:out_size(1), 1:out_size(2), 1:out_size(3));
s_ref = x * g(1) + y * g(2) + z * g(3) + c;
clear x y z
P_ref = E_ref .* exp(-1i * k * s_ref);
end
