function [out_path, ct_out_path] = write_babelbrain_input(medium_masks, focus_pos, parameters, planimg, out_dir, target_name, ppw_override, pseudoCT, trans_pos)
% WRITE_BABELBRAIN_INPUT  Write BabelViscoInput NIfTI from PRESTUS medium masks
%
% Converts PRESTUS medium_masks (layer index volume in grid space) to a
% single NIfTI with BabelBrain tissue labels and a focal voxel marker,
% using the required BabelBrain filename convention.
%
% BabelBrain tissue label scheme (integer, uint8):
%   0 = background / water
%   1 = skin
%   2 = cortical bone  (skull_cortical, or skull when no cortical/trabecular split)
%   3 = trabecular bone
%   4 = brain
%   5 = focal point    (single voxel — the only accepted way to pass the target)
%
% Required BabelBrain filename:
%   <out_dir>/<ID>/<target_name>_<freq_kHz>kHz_<PPW>PPW_BabelViscoInput.nii.gz
%
% The NIfTI is written in the PRESTUS simulation grid space, using the
% affine derived from planimg (same machinery as nifti_to_t1w). BabelBrain
% reads this directly; no T1-space back-projection is needed for the mask.
%
% Use as:
%   out_path = write_babelbrain_input(medium_masks, focus_pos, parameters, planimg, out_dir, target_name)
%   [out_path, ct_out_path] = write_babelbrain_input(medium_masks, focus_pos, parameters, planimg, ...
%       out_dir, target_name, ppw_override, pseudoCT, trans_pos)
%
% Input:
%   medium_masks  - [Nx x Ny x Nz] integer label volume from preproc_medium_mask
%                   Values are 1-based indices into parameters.medium_properties
%   focus_pos     - [1x3] focus voxel position in grid coordinates (1-based)
%   parameters    - PRESTUS config struct; uses:
%                     .layers           (layer-name → index mapping)
%                     .transducer(1).freq_hz
%                     .grid.resolution_mm
%                     .subject_id  (used in BabelBrain ID field)
%   planimg       - struct from grid_tissue_setup; uses origin_ras_mm and
%                   t1_header for affine reconstruction (may be empty for
%                   non-layered phantom grids)
%   out_dir       - output base directory; file is written under
%                   <out_dir>/sub-<ID>/
%   target_name   - string used as the BabelBrain "target" identifier
%                   (appears in the filename, e.g. 'M1_left')
%   ppw_override  - (optional) points per wavelength resolved by the caller
%                   (babelbrain_wrapper); 0 (default) = derive from
%                   grid.resolution_mm
%   pseudoCT      - (optional) [Nx x Ny x Nz] Hounsfield-equivalent volume,
%                   same grid as medium_masks. Written as CT.nii.gz with the
%                   same resampling and header as the tissue mask, as raw HU;
%                   run_babelbrain_case.py quantizes it and writes CT-cal.npz.
%   trans_pos     - (optional) [1x3] transducer voxel in the same grid as
%                   focus_pos; sets the beam direction for the depth flip
%                   (default: beam +z assumed)
%
% Output:
%   out_path      - full path to the written BabelViscoInput NIfTI file
%   ct_out_path   - full path to the written CT NIfTI file, or '' if pseudoCT
%                   was not supplied
%
% See also: BABELBRAIN_WRAPPER, BABELBRAIN_TRANSDUCER_MAP, EXPORT_BABELBRAIN_WEIGHTS

arguments
    medium_masks  (:,:,:) {mustBeNumeric}
    focus_pos     (1,3)   {mustBeNumeric}
    parameters    (1,1)   struct
    planimg       (1,1)   struct
    out_dir       (1,:)   char
    target_name   (1,:)   char
    ppw_override  (1,1)   double = 0   % pass resolved PPW from wrapper; 0 = auto-compute
    pseudoCT               {mustBeNumeric} = []
    trans_pos              {mustBeNumeric} = []   % [1x3] transducer voxel in the same grid as focus_pos; sets beam direction for the depth flip
end

ct_out_path = '';
if ~isempty(pseudoCT) && ~isequal(size(pseudoCT), size(medium_masks))
    error('write_babelbrain_input: pseudoCT size [%s] must match medium_masks size [%s].', ...
        num2str(size(pseudoCT)), num2str(size(medium_masks)));
end

%% Derive PPW and build filename
freq_hz  = parameters.transducer(1).freq_hz;
freq_kHz = round(freq_hz / 1e3);
res_mm   = parameters.grid.resolution_mm;

% SmallestSOS BabelBrain will actually use for its own spatial step (see below)
% -- computed up front since the PPW fallback (when the caller doesn't already
% resolve one, e.g. ppw_override=0) should target the same reference, not a
% separate/inconsistent water- or brain-based guess.
pichardo_shear_floor_ms = 1000.0 * 0.422 + 680.515;
smallest_sos_ms = pichardo_shear_floor_ms;
if isfield(parameters, 'medium_properties')
    mp_names = fieldnames(parameters.medium_properties);
    for mi = 1:numel(mp_names)
        mp = parameters.medium_properties.(mp_names{mi});
        if isfield(mp, 'sound_speed') && mp.sound_speed > 0
            smallest_sos_ms = min(smallest_sos_ms, mp.sound_speed);
        end
    end
end

if ppw_override > 0
    ppw = ppw_override;
else
    % Chosen to make BabelBrain's own resolution land as close as possible to
    % parameters.grid.resolution_mm (see bb_native_res_mm note below) --
    % callers that already resolve PPW themselves (babelbrain_wrapper.m) pass
    % it in via ppw_override instead of relying on this fallback.
    ppw = round(smallest_sos_ms / (freq_hz * res_mm * 1e-3));
end

%% BabelBrain's own native spatial step (voxel size it will actually simulate at)
% BabelBrain uses the input mask's voxel count directly as its FDTD grid size
% and assumes its own step for each voxel, ignoring the NIfTI PixelDimensions:
%   SpatialStep = SmallestSOS / freq / PPW
% SmallestSOS is the smallest positive sound speed in its material table,
% floored by the Pichardo shear speed at 1000 kg/m3 (1000*0.422 + 680.515 =
% 1102.515 m/s; GetSmallestSOS/DensityToSSoSPichardo in BabelIntegrationBASE.py).
% With the fluid-bone model (_patch_noncat_skull_attenuation in
% run_babelbrain_case.py) that floor is the binding term. A mask written at any
% other step would be read as a smaller or larger physical domain, so the mask
% is resampled to this step below.
bb_native_res_mm = (smallest_sos_ms / freq_hz) / ppw * 1e3;
resample_factor  = res_mm / bb_native_res_mm;   % > 1 => BabelBrain's native grid is finer

fprintf(['[write_babelbrain_input] BabelBrain native resolution: %.4f mm (PRESTUS grid: %.4f mm) -- ' ...
    'resampling mask by factor %.4f so BabelBrain''s voxel-count-based domain sizing preserves physical size.\n'], ...
    bb_native_res_mm, res_mm, resample_factor);

subj_id  = sprintf('sub-%03d', parameters.subject_id);
sub_dir  = fullfile(out_dir, subj_id);
if ~exist(sub_dir, 'dir')
    mkdir(sub_dir);
end

fname    = sprintf('%s_%dkHz_%dPPW_BabelViscoInput.nii.gz', target_name, freq_kHz, ppw);
out_path = fullfile(sub_dir, fname);

%% Map PRESTUS medium layer indices → BabelBrain labels
% medium_masks stores field positions from fieldnames(parameters.medium_properties),
% as assigned by preproc_medium_mask (e.g. find(strcmp(medium_labels,'skull_cortical'))).
% parameters.layers stores SimNIBS tissue indices — a different numbering system.
% We must use field-position lookup (consistent with tissuemask_binary) here.

label_map = zeros(size(medium_masks), 'uint8');  % default 0 = background/water

% Helper: get medium_masks value for a given tissue name (field position in medium_properties)
medium_labels = fieldnames(parameters.medium_properties);
layer_idx = @(name) find(strcmp(medium_labels, name));

idx_brain = layer_idx('brain');
if ~isempty(idx_brain)
    label_map(ismember(medium_masks, idx_brain)) = 4;
end

idx_skin = layer_idx('skin');
if ~isempty(idx_skin)
    label_map(ismember(medium_masks, idx_skin)) = 1;
end

% Prefer split cortical/trabecular if available; fall back to generic skull.
idx_cortical    = layer_idx('skull_cortical');
idx_trabecular  = layer_idx('skull_trabecular');
idx_skull_gen   = layer_idx('skull');

if ~isempty(idx_cortical)
    label_map(ismember(medium_masks, idx_cortical))   = 2;
end
if ~isempty(idx_trabecular)
    label_map(ismember(medium_masks, idx_trabecular)) = 3;
end
if ~isempty(idx_skull_gen)
    unassigned_skull = ismember(medium_masks, idx_skull_gen) & label_map == 0;
    label_map(unassigned_skull) = 2;  % treat as cortical when undifferentiated
end

% Resample the tissue label map (marker not yet placed) from the PRESTUS grid
% to BabelBrain's native step with nearest-neighbour interpolation (integer
% labels). Voxel (1,1,1) maps to itself in both grids. Two mutually inverse
% matrices are needed, and must not be swapped:
%   M_old_to_new - passed to affine_resample_3d, which takes the forward
%                  transform [out 1] = [in 1] * M (OLD PRESTUS grid -> NEW
%                  BabelBrain grid) and inverts it internally:
%                  diag(k) with translation (1-k).
%   M_new_to_old - used for the header affine (NEW voxel -> OLD voxel):
%                  diag(1/k) with translation (1-1/k).
old_dims = size(label_map);
new_dims = max(1, round(old_dims * resample_factor));
if abs(resample_factor - 1) > 1e-6
    M_old_to_new = diag([resample_factor, resample_factor, resample_factor, 1]);
    M_old_to_new(4, 1:3) = 1 - resample_factor;
    label_map = uint8(round(affine_resample_3d(single(label_map), M_old_to_new, new_dims, 'nearest', 0)));

    M_new_to_old = diag([1/resample_factor, 1/resample_factor, 1/resample_factor, 1]);
    M_new_to_old(4, 1:3) = 1 - 1/resample_factor;
end

% Mark focal voxel as label 5 — BabelBrain requires exactly one such voxel.
% Rescale focus_pos into the (possibly resampled) grid using the same
% old<->new voxel relationship as M_new_to_old, anchored at voxel (1,1,1).
fp = round((focus_pos - 1) * resample_factor) + 1;
fp = max(1, min(fp, new_dims));  % clamp to grid
label_map(fp(1), fp(2), fp(3)) = 5;

% Orient the depth axis so BabelBrain's acoustic source lands on the TRANSDUCER
% side. After BabelBrain's own np.flip(axis=2) on load, its z=0 face must be the
% transducer side, with z increasing toward the focus. The beam direction along
% dim 3 depends on the head pose:
%   focus_pos(3) >= trans_pos(3)  (beam +z): pre-flip here, so BabelBrain's load
%                                  flip cancels it
%   focus_pos(3) <  trans_pos(3)  (beam -z): no pre-flip; BabelBrain's load flip
%                                  alone gives the required orientation
% Without trans_pos the direction is unknown and beam +z is assumed.
do_flip = true;
if ~isempty(trans_pos) && numel(trans_pos) >= 3
    % >= so the tie (focus_pos(3)==trans_pos(3)) matches read_babelbrain_output's
    % beam_sign default (+1), keeping input flip and output resample in sync.
    do_flip = focus_pos(3) >= trans_pos(3);
    if do_flip
        beam_str = 'beam +z (pre-flip depth)';
    else
        beam_str = 'beam -z (do NOT pre-flip)';
    end
    fprintf(['[write_babelbrain_input] Beam direction along dim 3: focus_pos(3)=%.1f vs trans_pos(3)=%.1f ' ...
        '-> %s.\n'], focus_pos(3), trans_pos(3), beam_str);
else
    fprintf(['[write_babelbrain_input] trans_pos not provided; assuming beam +z (pre-flipping depth). ' ...
        'If the transducer ends up on the wrong side, pass trans_pos.\n']);
end
if do_flip
    label_map = flip(label_map, 3);
end

fprintf('[write_babelbrain_input] Label counts:\n');
for lbl = [0 1 2 3 4 5]
    n = sum(label_map(:) == lbl);
    names = {'background','skin','cortical','trabecular','brain','focal'};
    if n > 0
        fprintf('  %d (%s): %d voxels\n', lbl, names{lbl+1}, n);
    end
end

%% Build NIfTI header with PRESTUS grid affine
% Write a temporary NIfTI to obtain a valid base header struct (niftiwrite
% requires specific fields like Qfactor that cannot be hand-assembled).
tmp_file = fullfile(tempdir, sprintf('prestus_bb_hdr_%d', randi(1e6)));
niftiwrite(uint8(label_map), tmp_file);
hdr = niftiinfo([tmp_file '.nii']);
delete([tmp_file '.nii']);

hdr.Datatype        = 'uint8';
hdr.BitsPerPixel    = 8;
hdr.PixelDimensions = repmat(bb_native_res_mm, 1, 3);
hdr.SpaceUnits      = 'Millimeter';
hdr.Description     = sprintf('PRESTUS BabelViscoInput: %s  freq=%dkHz  ppw=%d', ...
    subj_id, freq_kHz, ppw);

% Reconstruct affine. For layered simulations, compose the acoustic grid's
% actual head-alignment ROTATION (planimg.inv_transf; the grid is generally
% NOT axis-aligned with the T1 — grid.mode defaults to a beam/focal-axis
% alignment that rotates the head) with the T1 NIfTI's own voxel->RAS affine,
% and any acoustic-FOV-crop voxel offset recorded on planimg.fov_offset_ac —
% via grid_to_ras_affine. For non-layered (phantom) grids, which genuinely
% have no rotation relative to RAS, the diagonal affine built from
% origin_ras_mm remains correct (same as nifti_to_t1w).
%
% Both branches also compose M_new_to_old (defined above, when resampling to
% BabelBrain's native resolution was applied) so the header maps THIS
% (resampled) array's voxels to RAS, not the original PRESTUS-grid voxels.
if ~isempty(planimg) && isfield(planimg, 't1_header') && ~isempty(planimg.t1_header) ...
        && isfield(planimg, 'inv_transf') && ~isempty(planimg.inv_transf)
    if isfield(planimg, 'fov_offset_ac') && ~isempty(planimg.fov_offset_ac)
        fov_offset = planimg.fov_offset_ac;
    else
        fov_offset = [1 1 1];
    end
    T_grid_to_ras = grid_to_ras_affine(planimg, fov_offset);
    if abs(resample_factor - 1) > 1e-6
        T_grid_to_ras = M_new_to_old * T_grid_to_ras;
    end
    % When label_map was pre-flipped along dim 3 (see do_flip above) to match
    % BabelBrain's own load-time flip, compose that same flip into the affine so
    % the header describes THIS (flipped) array. When it was NOT pre-flipped
    % (beam −z), the array is in original grid order and no flip is composed.
    if do_flip
        Nz = new_dims(3);
        M_flip3 = eye(4);
        M_flip3(3,3) = -1;
        M_flip3(4,3) = Nz + 1;
        hdr.Transform = affine3d(nifti_affine_from_1based(M_flip3 * T_grid_to_ras));
    else
        hdr.Transform = affine3d(nifti_affine_from_1based(T_grid_to_ras));
    end
elseif ~isempty(planimg) && isfield(planimg, 't1_header') && ~isempty(planimg.t1_header)
    % Fallback (inv_transf unavailable): origin-only, no rotation — better
    % than nothing but will misalign whenever the grid is rotated relative
    % to the T1 (e.g. grid.mode ~= 'ras_plus' with a non-trivial head pose).
    warning('write_babelbrain_input:noInvTransf', ...
        ['planimg.inv_transf not available; BabelViscoInput.nii.gz header will use T1 origin only ' ...
         '(no rotation/FOV-offset correction). Overlay alignment against the T1/medium_masks may be off.']);
    t1hdr = planimg.t1_header;
    orig_origin = t1hdr.Transform.T(4, 1:3);    % [x y z] mm; affine3d uses row-vector convention: translation is at T(4,1:3)
    T = diag([bb_native_res_mm bb_native_res_mm bb_native_res_mm 1]);
    T(4, 1:3) = orig_origin;
    hdr.Transform = affine3d(T);
elseif ~isempty(planimg) && isfield(planimg, 'origin_ras_mm') && ~isempty(planimg.origin_ras_mm)
    origin = planimg.origin_ras_mm(:)';
    T = diag([bb_native_res_mm bb_native_res_mm bb_native_res_mm 1]);
    T(4, 1:3) = origin;
    hdr.Transform = affine3d(T);
% else: leave MATLAB default identity affine; warn user
else
    warning('write_babelbrain_input: planimg has no affine info. NIfTI will have identity transform.');
end

%% Write NIfTI
niftiwrite(label_map, strrep(out_path, '.nii.gz', ''), hdr, 'Compressed', true);
fprintf('[write_babelbrain_input] Written: %s\n', out_path);

%% Write companion CT NIfTI (raw HU, pre-quantization), same grid as the mask
if ~isempty(pseudoCT)
    ct_vol = single(pseudoCT);
    if abs(resample_factor - 1) > 1e-6
        ct_vol = affine_resample_3d(ct_vol, M_old_to_new, new_dims, 'linear', 0);
    end
    ct_vol = flip(ct_vol, 3);   % same load-time-flip pre-compensation as label_map

    ct_fname    = strrep(fname, 'BabelViscoInput', 'CT');
    ct_out_path = fullfile(sub_dir, ct_fname);

    ct_tmp_file = fullfile(tempdir, sprintf('prestus_bb_ct_hdr_%d', randi(1e6)));
    niftiwrite(ct_vol, ct_tmp_file);
    ct_hdr = niftiinfo([ct_tmp_file '.nii']);
    delete([ct_tmp_file '.nii']);

    ct_hdr.PixelDimensions = repmat(bb_native_res_mm, 1, 3);
    ct_hdr.SpaceUnits      = 'Millimeter';
    ct_hdr.Description     = sprintf('PRESTUS BabelBrain CT (raw HU, pre-quantization): %s  freq=%dkHz  ppw=%d', ...
        subj_id, freq_kHz, ppw);
    ct_hdr.Transform       = hdr.Transform;   % identical grid/geometry as the tissue mask above

    niftiwrite(ct_vol, strrep(ct_out_path, '.nii.gz', ''), ct_hdr, 'Compressed', true);
    fprintf(['[write_babelbrain_input] CT NIfTI written (raw HU, NOT yet the quantized index-map format ' ...
        'BabelBrain''s CT path requires -- quantization happens in run_babelbrain_case.py): %s\n'], ct_out_path);
end

end
