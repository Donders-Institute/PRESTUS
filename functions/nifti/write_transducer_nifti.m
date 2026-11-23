function write_transducer_nifti(parameters, t1, out_file, back_mm)
% WRITE_TRANSDUCER_NIFTI  Write a transducer model on the exact voxel grid of a T1
%
% Writes the transducer bowl(s) as a float32 NIfTI with the dimensions, voxel
% size and affine of the T1. Voxel values:
%   1    emitting face: the bowl surface across the active aperture, from the
%        smallest elem_id_mm to the largest elem_od_mm
%   0.5  the display thickness behind the face (back_mm)
%   0    elsewhere, including the central opening inside the smallest elem_id_mm
%
% The bowl surface comes from makeBowl with the bowl radius TRANSDUCER_SETUP
% uses. Each surface voxel is kept when the distance of its centre from the
% beam axis lies within the aperture, using the exact mm limits of the profile,
% so the outer edge and the central opening are not rounded to whole voxels.
%
% makeBowl needs isotropic voxels, so the bowl is built on a padded isotropic
% grid at the smallest T1 voxel size that shares the T1 axes and voxel [1 1 1].
% Every model voxel is then mapped onto the T1 voxel that contains it; where
% values meet, the larger wins. Parts outside the T1 field of view are dropped,
% and the share dropped is printed.
%
% back_mm thickens each bowl backwards along its beam axis (away from the
% focus) for display. The emitting face stays where PRESTUS places it.
% back_mm = 0 gives the face alone.
%
% Use as:
%   write_transducer_nifti(parameters, t1, out_file)
%   write_transducer_nifti(parameters, t1, out_file, back_mm)
%
% Input:
%   parameters - PRESTUS config struct; parameters.transducer(i) must be annular,
%                with trans_pos (bowl rear centre) and focus_pos as voxel
%                indices of t1 (1-based, as PRESTUS uses them)
%   t1         - planning T1: a NIfTI file path, or a niftiinfo header struct
%                (e.g. planimg.t1_header)
%   out_file   - output path, ending in '.nii.gz'
%   back_mm    - display thickness behind the emitting face [mm] (default 3)
%
% See also: NIFTI_TRANSDUCER, TRANSDUCER_SETUP, NIFTI_TO_T1W

if nargin < 4 || isempty(back_mm), back_mm = 3; end

if ischar(t1) || isstring(t1)
    hdr = niftiinfo(char(t1));
else
    hdr = t1;
end
sz_t1  = double(hdr.ImageSize(1:3));
pixdim = double(hdr.PixelDimensions(1:3));
vox    = min(pixdim);
scale  = pixdim / vox;                       % one T1 voxel step in isotropic voxel steps

% T1 voxel indices -> isotropic grid indices (both grids share axes and voxel [1 1 1])
n    = numel(parameters.transducer);
tpos = zeros(n, 3);
fpos = zeros(n, 3);
for ti = 1:n
    tpos(ti,:) = round(1 + (reshape(double(parameters.transducer(ti).trans_pos), 1, 3) - 1) .* scale);
    fpos(ti,:) = round(1 + (reshape(double(parameters.transducer(ti).focus_pos), 1, 3) - 1) .* scale);
end

% Padded isotropic working grid so the whole bowl fits: r reaches from trans_pos
% to the farthest bowl point (rim of the makeBowl diameter, at its sagitta) plus back_mm
reach = 0;
for ti = 1:n
    reach = max(reach, bowl_reach_mm(parameters.transducer(ti), vox));
end
sz       = ceil(sz_t1 .* scale);
r        = ceil((reach + back_mm) / vox) + 1;
pad_pre  = max(0, r + 1 - min(tpos, [], 1));
pad_post = max(0, max(tpos, [], 1) + r - sz);
sz_pad   = sz + pad_pre + pad_post;

model = zeros(sz_pad, 'single');
for ti = 1:n
    face  = aperture_face(parameters.transducer(ti), tpos(ti,:) + pad_pre, ...
        fpos(ti,:) + pad_pre, sz_pad, vox);
    layer = 0.5 * single(thicken_back(face, tpos(ti,:), fpos(ti,:), back_mm / vox));
    layer(face) = 1;
    model = max(model, layer);
end

% Map every model voxel onto the T1 voxel that contains it; the larger value wins
idx = find(model);
[ii, jj, kk] = ind2sub(sz_pad, idx);
P  = round(([ii jj kk] - pad_pre - 1) ./ scale + 1);
ok = all(P >= 1 & P <= sz_t1, 2);
if any(~ok)
    fprintf('[write_transducer_nifti] %.0f%% of transducer voxels lie outside the T1 field of view and are dropped\n', ...
        100 * mean(~ok));
end
lin = sub2ind(sz_t1, P(ok,1), P(ok,2), P(ok,3));
out = single(reshape(accumarray(lin, double(model(idx(ok))), [prod(sz_t1) 1], @max), sz_t1));

% Write with the T1 header: same dimensions, voxel size and affine
hdr.ImageSize             = sz_t1;
hdr.PixelDimensions       = pixdim;
hdr.Datatype              = 'single';
hdr.BitsPerPixel          = 32;
hdr.MultiplicativeScaling = 1;
hdr.AdditiveOffset        = 0;
niftiwrite_gz(out, regexprep(char(out_file), '\.nii\.gz$', ''), hdr);
end

function reach = bowl_reach_mm(tr, vox)
% Distance [mm] from the bowl's rear centre to the rim of the bowl that
% APERTURE_FACE builds (same makeBowl diameter), at the rim's sagitta.
R     = double(tr.annular.curv_radius_mm);
r_out = max(double(tr.annular.elem_od_mm)) / 2;
rho   = min((2 * floor(2 * r_out / vox / 2) + 3) * vox / 2, R);
reach = hypot(rho, R - sqrt(R^2 - rho^2));
end

function face = aperture_face(tr, tpos, fpos, sz, vox)
% Bowl surface voxels whose centre lies within the active aperture, from the
% smallest elem_id_mm to the largest elem_od_mm (distance from the beam axis).
if ~strcmp(tr.type, 'annular')
    error('write_transducer_nifti:notAnnular', 'Only annular transducers are supported.');
end
r_in   = min(double(tr.annular.elem_id_mm)) / 2;     % central opening radius [mm]
r_out  = max(double(tr.annular.elem_od_mm)) / 2;     % aperture radius [mm]
radius = round(tr.annular.curv_radius_mm / vox);     % bowl radius in voxels, as TRANSDUCER_SETUP
diam   = 2 * floor(2 * r_out / vox / 2) + 3;         % odd; extends past the aperture edge
surf   = makeBowl(sz, tpos, radius, diam, fpos) > 0;
idx = find(surf);
[ii, jj, kk] = ind2sub(sz, idx);
X   = [ii jj kk] - tpos;
u   = (fpos - tpos) / norm(fpos - tpos);
rho = vecnorm(X - (X * u') * u, 2, 2) * vox;         % distance from the beam axis [mm]
face = false(sz);
face(idx(rho >= r_in & rho <= r_out)) = true;
end

function mask = thicken_back(mask, tpos, fpos, back_vox)
% Sweep a surface mask backwards along the beam axis, away from the focus.
u = (fpos - tpos) / norm(fpos - tpos);
[ii, jj, kk] = ind2sub(size(mask), find(mask));
P = [ii jj kk];
for d = 0.5:0.5:back_vox
    Q  = round(P - d * u);
    ok = all(Q >= 1 & Q <= size(mask), 2);
    mask(sub2ind(size(mask), Q(ok,1), Q(ok,2), Q(ok,3))) = true;
end
end

function niftiwrite_gz(data, filepath, hdr)
% Write NIfTI then compress with system gzip, as NIFTI_TO_T1W does: MATLAB's
% Java-based compression fails under macOS Full Disk Access restrictions, and
% the full .nii name stops niftiwrite reading a dot in the affix as the extension.
nii_path = [filepath '.nii'];
gz_path  = [filepath '.nii.gz'];
niftiwrite(data, nii_path, hdr);
if exist(gz_path, 'file'), delete(gz_path); end
[status, msg] = run_command({'gzip', '-f', nii_path});
if status ~= 0
    error('write_transducer_nifti:gzipFailed', 'gzip failed for %s: %s', nii_path, msg);
end
end
