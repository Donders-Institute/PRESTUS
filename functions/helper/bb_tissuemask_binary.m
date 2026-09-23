function mask = bb_tissuemask_binary(material_map)
% BB_TISSUEMASK_BINARY  Build per-tissue binary masks from a BabelBrain material map
%
% BabelBrain assigns tissue labels with a fixed convention that differs from
% PRESTUS's medium_masks (which indexes parameters.medium_properties field
% order). This helper builds the same mask struct as TISSUEMASK_BINARY but from
% BabelBrain's own material map, so a BabelBrain field can be masked with the
% segmentation it was actually computed on. After resampling, the map matches
% medium_masks except for single-voxel differences at tissue boundaries; using
% it attributes boundary voxels (e.g. the skull-brain reflection hotspot) the
% way the solver treated them.
%
% BabelBrain label convention (BabelIntegrationBASE.py / material map):
%   0 = water, 1 = skin, 2 = cortical bone, 3 = trabecular bone, 4 = brain
%
% Use as:
%   mask = bb_tissuemask_binary(material_map)
%
% Input:
%   material_map - BabelBrain material index volume, resampled onto the PRESTUS
%                  acoustic grid (see read_babelbrain_output.m)
%
% Output:
%   mask - struct with binary fields: skull, skull_cortical, skull_trabecular,
%          brain, skin, water  (same fields as TISSUEMASK_BINARY)
%
% See also: TISSUEMASK_BINARY, READ_BABELBRAIN_OUTPUT, ACOUSTIC_ANALYSIS

arguments
    material_map {mustBeNumericOrLogical}
end

    m = round(material_map);
    mask.water            = (m == 0);
    mask.skin             = (m == 1);
    mask.skull_cortical   = (m == 2);
    mask.skull_trabecular = (m == 3);
    mask.skull            = (m == 2) | (m == 3);
    mask.brain            = (m == 4);
end
