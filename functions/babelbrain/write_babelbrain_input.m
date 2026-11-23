function [out_path, ct_out_path, materials_path] = write_babelbrain_input(medium_masks, focus_pos, parameters, planimg, out_dir, target_name, ppw_override, pseudoCT, trans_pos, medium_table)
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
% With medium_table (babelbrain.medium 'prestus') the labels are material
% classes instead: 0 = water, 5 = focal point, 5 + k = class k of
% BABELBRAIN_MEDIUM_TABLE (uint16). BabelBrain reads labels >= 6 as its
% brain-segmentation materials and uses label - 1 as the material index, so
% <stem>_PRESTUSMaterials.json lists the materials in that order: water, three
% unused rows (water), the focal voxel's class, then classes 1..K.
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
%   medium_table  - (optional) struct with cls, table and water from
%                   BABELBRAIN_MEDIUM_TABLE (same grid as medium_masks); [] =
%                   tissue labels
%
% Output:
%   out_path      - full path to the written BabelViscoInput NIfTI file
%   ct_out_path   - full path to the written CT NIfTI file, or '' if pseudoCT
%                   was not supplied
%   materials_path - full path to the written <stem>_PRESTUSMaterials.json,
%                   or '' without medium_table
%
% See also: BABELBRAIN_WRAPPER, BABELBRAIN_TX_GEOMETRY, EXPORT_BABELBRAIN_WEIGHTS

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
    medium_table           = []                   % struct from babelbrain_medium_table; [] = tissue labels
end

materials_path = '';
use_table = ~isempty(medium_table);
if use_table
    if ~isequal(size(medium_table.cls), size(medium_masks))
        error('write_babelbrain_input:tableSize', 'medium_table.cls size [%s] must match medium_masks size [%s].', ...
            num2str(size(medium_table.cls)), num2str(size(medium_masks)));
    end
    if ~isempty(pseudoCT)
        error('write_babelbrain_input:tableWithCT', ...
            'medium_table carries the pseudo-CT properties itself; pass no pseudoCT with it.');
    end
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
[~, smallest_sos_ms] = babelbrain_grid_step_mm(parameters);   % BabelBrain's SmallestSOS (Pichardo shear floor or slower medium)

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
% 1102.515 m/s; GetSmallestSOS/DensityToSSoSPichardo in babel_integration_helpers.py).
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

% babelbrain.medium 'prestus': the material classes replace the tissue labels
% (same geometry: the classes come from the medium built on these masks)
label_class = 'uint8';
if use_table
    label_class = 'uint16';
    label_map = zeros(size(medium_masks), label_class);
    tissue_vox = medium_table.cls > 0;
    label_map(tissue_vox) = uint16(medium_table.cls(tissue_vox)) + 5;
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
    label_map = cast(round(affine_resample_3d(single(label_map), M_old_to_new, new_dims, 'nearest', 0)), label_class);

    M_new_to_old = diag([1/resample_factor, 1/resample_factor, 1/resample_factor, 1]);
    M_new_to_old(4, 1:3) = 1 - 1/resample_factor;
end

% Mark focal voxel as label 5 — BabelBrain requires exactly one such voxel.
% Rescale focus_pos into the (possibly resampled) grid using the same
% old<->new voxel relationship as M_new_to_old, anchored at voxel (1,1,1).
fp = round((focus_pos - 1) * resample_factor) + 1;
fp = max(1, min(fp, new_dims));  % clamp to grid
focus_class = 0;
if use_table && label_map(fp(1), fp(2), fp(3)) > 5
    focus_class = double(label_map(fp(1), fp(2), fp(3))) - 5;
end
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

% Keep the whole bowl behind BabelBrain's source plane. BabelBrain injects the
% Rayleigh field on the FIRST domain slice on the transducer side (the high-k
% end of this array, after its load flip = its z=0), and requires every source
% element to lie behind that plane: if the bowl rim protrudes into the domain it
% silently snaps the transducer back in whole voxels ("Need to apply correction
% to reposition Tx"), which moves the source off trans_pos (observed: -5.5 mm
% for a pFOV crop that starts at the apex). The slab between the transducer-side
% face and the rim is water inside the bowl concavity, so trim it: the source
% then sits exactly where PRESTUS placed it, and the trimmed slab is simply
% zero-filled on the way back (read_babelbrain_output). The high-k end is
% trimmed, so voxel indices 1..Nz-n (marker, affine) are unchanged.
Nz_full = size(label_map, 3);
n_trim  = 0;
bowl_depth_mm = [];
if isfield(parameters, 'transducer') && isfield(parameters.transducer(1), 'annular') && ...
        all(isfield(parameters.transducer(1).annular, {'curv_radius_mm', 'elem_od_mm'}))
    an = parameters.transducer(1).annular;
    R_mm = an.curv_radius_mm;
    a_mm = max(an.elem_od_mm(:)) / 2;
    if a_mm < R_mm
        bowl_depth_mm = R_mm - sqrt(R_mm^2 - a_mm^2);
    end
end
if ~isempty(trans_pos) && numel(trans_pos) >= 3 && ~isempty(bowl_depth_mm)
    kt = (trans_pos(3) - 1) * resample_factor + 1;           % apex index, native grid, pre-flip
    if do_flip
        d_apex_mm = (kt - 1) * bb_native_res_mm;             % apex -> transducer-side face
    else
        d_apex_mm = (Nz_full - kt) * bb_native_res_mm;
    end
    % The rim lies bowl_depth on the TARGET side of the apex, so the bowl spans
    % d_apex .. d_apex + bowl_depth into the domain (d_apex < 0: apex outside the
    % face). Trim past the rim; +1 voxel margin since BabelBrain tests
    % max(element z) >= source-plane z.
    n_trim = max(0, ceil((d_apex_mm + bowl_depth_mm) / bb_native_res_mm) + 1);
    if n_trim > 0
        slab = label_map(:, :, Nz_full-n_trim+1:Nz_full);
        if any(slab(:) ~= 0)
            warning('write_babelbrain_input:tissueInBowl', ...
                ['%d non-water voxel(s) in the %.1f mm slab between the transducer apex and the bowl rim; ' ...
                 'they are dropped from the BabelBrain domain (BabelBrain cannot model tissue inside the bowl).'], ...
                nnz(slab), n_trim * bb_native_res_mm);
        end
        if do_flip; k_marker = Nz_full + 1 - fp(3); else; k_marker = fp(3); end
        if k_marker > Nz_full - n_trim
            error('write_babelbrain_input:trimTooDeep', ...
                'Bowl-rim trim (%d slices) would remove the focal marker; check trans_pos/focus_pos.', n_trim);
        end
        label_map = label_map(:, :, 1:Nz_full-n_trim);
        fprintf(['[write_babelbrain_input] Trimmed %d transducer-side slice(s) (%.2f mm) so the bowl rim ' ...
            '(depth %.2f mm, apex %.2f mm from the domain face) lies behind BabelBrain''s source plane.\n'], ...
            n_trim, n_trim * bb_native_res_mm, bowl_depth_mm, d_apex_mm);
    end
end

fprintf('[write_babelbrain_input] Label counts:\n');
if use_table
    fprintf('  0 (water): %d voxels, 5 (focal): 1 voxel, classes 6..%d: %d voxels in %d classes\n', ...
        nnz(label_map == 0), 5 + size(medium_table.table, 1), nnz(label_map > 5), numel(unique(label_map(label_map > 5))));
else
    for lbl = [0 1 2 3 4 5]
        n = sum(label_map(:) == lbl);
        names = {'background','skin','cortical','trabecular','brain','focal'};
        if n > 0
            fprintf('  %d (%s): %d voxels\n', lbl, names{lbl+1}, n);
        end
    end
end

%% Build NIfTI header with PRESTUS grid affine
% Write a temporary NIfTI to obtain a valid base header struct (niftiwrite
% requires specific fields like Qfactor that cannot be hand-assembled).
tmp_file = fullfile(tempdir, sprintf('prestus_bb_hdr_%d', randi(1e6)));
niftiwrite(uint8(label_map), tmp_file);
hdr = niftiinfo([tmp_file '.nii']);
delete([tmp_file '.nii']);

hdr.Datatype        = label_class;
hdr.BitsPerPixel    = 8 * (1 + use_table);
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

%% Write the material table (babelbrain.medium 'prestus')
% BabelBrain material index = label - 1 for labels >= 5 (CreateMaterialMaps,
% brain-segmentation branch): 0 water, 1-3 unused, 4 the focal voxel (label 5),
% 4 + k class k (label 5 + k). Columns: density, long. SoS, shear SoS, long.
% attenuation [Np/m], shear attenuation (BabelBrain's AddMaterial order).
if use_table
    w = medium_table.water;
    if focus_class > 0; f_row = medium_table.table(focus_class, :); else; f_row = w; end
    rows = [repmat(w, 4, 1); f_row; medium_table.table];
    mats = [rows(:,1), rows(:,2), zeros(size(rows,1),1), rows(:,3), zeros(size(rows,1),1)];
    materials_path = strrep(out_path, 'BabelViscoInput.nii.gz', 'PRESTUSMaterials.json');
    fid = fopen(materials_path, 'w');
    if fid < 0
        error('write_babelbrain_input:materialsWrite', 'Cannot open %s for writing.', materials_path);
    end
    fwrite(fid, jsonencode(struct('columns', {{'density', 'long_sos', 'shear_sos', 'long_att_np_m', 'shear_att_np_m'}}, ...
        'materials', mats, 'focal_class', focus_class)), 'char');
    fclose(fid);
    fprintf('[write_babelbrain_input] Material table: %d materials (%d classes) -> %s\n', ...
        size(mats, 1), size(medium_table.table, 1), materials_path);
end

%% Write companion CT NIfTI (raw HU, pre-quantization), same grid as the mask
if ~isempty(pseudoCT)
    ct_vol = single(pseudoCT);
    if abs(resample_factor - 1) > 1e-6
        ct_vol = affine_resample_3d(ct_vol, M_old_to_new, new_dims, 'linear', 0);
    end
    if do_flip
        ct_vol = flip(ct_vol, 3);   % same load-time-flip pre-compensation as label_map
    end
    ct_vol = ct_vol(:, :, 1:Nz_full-n_trim);   % same bowl-rim trim as label_map

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
