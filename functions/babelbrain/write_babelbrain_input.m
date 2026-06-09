function out_path = write_babelbrain_input(medium_masks, focus_pos, parameters, planimg, out_dir, target_name)
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
%
% Output:
%   out_path  - full path to the written NIfTI file
%
% See also: BABELBRAIN_WRAPPER, BABELBRAIN_TRANSDUCER_MAP, WRITE_BABELBRAIN_WEIGHTS

arguments
    medium_masks  (:,:,:) {mustBeNumeric}
    focus_pos     (1,3)   {mustBeNumeric}
    parameters    (1,1)   struct
    planimg       (1,1)   struct
    out_dir       (1,:)   char
    target_name   (1,:)   char
end

%% Derive PPW and build filename
freq_hz  = parameters.transducer(1).freq_hz;
freq_kHz = round(freq_hz / 1e3);
res_mm   = parameters.grid.resolution_mm;

% BabelBrain convention: basePPW = c_min / (freq * resolution_mm)
% Use minimum soft-tissue speed (brain, ~1540 m/s) as reference.
c_min = 1500;  % conservative lower bound matching BabelBrain default
if isfield(parameters, 'medium_properties') && isfield(parameters.medium_properties, 'brain')
    c_min = parameters.medium_properties.brain.sound_speed;
end
ppw = round(c_min / (freq_hz * res_mm * 1e-3));

subj_id  = sprintf('sub-%03d', parameters.subject_id);
sub_dir  = fullfile(out_dir, subj_id);
if ~exist(sub_dir, 'dir')
    mkdir(sub_dir);
end

fname    = sprintf('%s_%dkHz_%dPPW_BabelViscoInput.nii.gz', target_name, freq_kHz, ppw);
out_path = fullfile(sub_dir, fname);

%% Map PRESTUS medium layer indices → BabelBrain labels
% parameters.layers maps human-readable names to PRESTUS layer indices.
% We must convert those indices to BabelBrain's fixed label scheme.

label_map = zeros(size(medium_masks), 'uint8');  % default 0 = background/water

% Helper: get all medium indices for a given layer-name pattern
layer_idx = @(name) getidx(parameters.layers, name);

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

% Mark focal voxel as label 5 — BabelBrain requires exactly one such voxel.
fp = round(focus_pos);
fp = max(1, min(fp, size(medium_masks)));  % clamp to grid
label_map(fp(1), fp(2), fp(3)) = 5;

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
hdr.PixelDimensions = repmat(res_mm, 1, 3);
hdr.SpaceUnits      = 'Millimeter';
hdr.Description     = sprintf('PRESTUS BabelViscoInput: %s  freq=%dkHz  ppw=%d', ...
    subj_id, freq_kHz, ppw);

% Reconstruct affine. For layered simulations planimg.t1_header carries the
% original T1 affine which defines the RAS origin; for non-layered grids the
% diagonal RAS affine is built from origin_ras_mm (same as nifti_to_t1w).
if ~isempty(planimg) && isfield(planimg, 't1_header') && ~isempty(planimg.t1_header)
    % Layered simulation: use T1 affine origin directly. The grid is already
    % in the same orientation but may have different voxel size — rebuild
    % the sform with the correct pixel dimensions and origin.
    t1hdr = planimg.t1_header;
    orig_origin = t1hdr.Transform.T(4, 1:3);  % [x y z] mm, row-major MATLAB affine3d
    T = diag([res_mm res_mm res_mm 1]);
    T(4, 1:3) = orig_origin;
    hdr.Transform = affine3d(T');
elseif ~isempty(planimg) && isfield(planimg, 'origin_ras_mm') && ~isempty(planimg.origin_ras_mm)
    origin = planimg.origin_ras_mm(:)';
    T = diag([res_mm res_mm res_mm 1]);
    T(4, 1:3) = origin;
    hdr.Transform = affine3d(T');
% else: leave MATLAB default identity affine; warn user
else
    warning('write_babelbrain_input: planimg has no affine info. NIfTI will have identity transform.');
end

%% Write NIfTI
niftiwrite(label_map, strrep(out_path, '.nii.gz', ''), hdr, 'Compressed', true);
fprintf('[write_babelbrain_input] Written: %s\n', out_path);

end
